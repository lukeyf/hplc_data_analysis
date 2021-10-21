# %%
import pathlib
from datetime import datetime
import fnmatch
import os
import pandas as pd
from watcher import Watcher
import numpy as np
# from configuration  import deck
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from aghplctools.data.sample import HPLCSample, HPLCSampleInfo, DADSignalInfo

DEAD_VOLUME_TIME = 0.6


# def extract_data(filename: str = 'DAD1 A 210 nm blank.csv'):
#     """
#     parse the existing csv file data
#     :param filename: the name of the csv file
#     :return: the time points and intensities of the data
#     """
#     data = pd.read_csv(filename).to_numpy()
#     print(data)
#     retention_time = data[1:, 0]
#     intensity_data = data[1:, 1]
#     return retention_time, intensity_data

def extract_data(filename: str, wavelength: int = 210):
    """
    get the data of time and intensity from the home folder
    :param filename: the name of the home
    :return: the time points and intensities of the data
    """
    folder = pathlib.Path(filename)
    data = HPLCSample.create_from_D_file(folder)
    signal = data.signals[0]
    for s in data.signals:
        if int(s.wavelength.real) == wavelength:
            signal = s
    retention_time = signal.retention_times
    intensity_data = signal.mean_unreferenced_intensities
    return retention_time, intensity_data


def extract_time(filename: str = 'Automatically_Generated_Report00.CSV'):
    """
    extract the injection time from the generated report
    :param filename: the name of the report csv file
    :return: a datetime object of the injection time
    """
    f = pd.read_csv(filename, encoding='utf-16').to_numpy()
    date = f[8, 1]
    return datetime.strptime(date, '%d-%b-%y, %H:%M:%S')


def integration(x, y):
    """
    Integrate y=f(x) over interval x using a degree five spline interpolation.
    :param x: x-axis data
    :param y: y-axis data
    :return: the value integrating y over x
    """

    # interpolation
    print(x,y)
    fn = UnivariateSpline(x, y, k=5)

    results = quad(fn, x[0], x[-1])[0]
    return results


def peak_properties(blank_data, reaction_data, internal_standard_retention_time: float = 1.77,
                    peak_width: float = 0.04, detection_limit: float = 5, plot: bool = False):
    """
    The peak properties from a given blank data, peak raw data
    :param blank_data: The time and intensity data of the blank (List of lists)
    :param reaction_data: The time and intensity data of an injection (List of lists)
    :param internal_standard_retention_time: the internal standard retention time
    :param peak_width: the shift tolerance of identifying a peak
    :param detection_limit: The minimum intensity to
    :param plot: True means show and save the current plot of peaks, False otherwise
    :return: a list of peak maximum intensities, areas, retention time (min)
    """

    retention_time, intensity_blank_data = blank_data[0], blank_data[1]
    INTERNAL_STANDARD_RETENTION_TIME = internal_standard_retention_time  # min
    PEAK_WIDTH = peak_width  # min
    DETECTION_LIMIT = detection_limit
    intensity_raw_data = reaction_data[1][:len(retention_time)]

    # Calculate baseline
    processed_data = intensity_raw_data - intensity_blank_data
    base_data = []
    for i in range(len(processed_data)):
        if processed_data[i] - np.average(processed_data) <= 3:
            base_data.append(processed_data[i])
    base_intensity = np.average(base_data)

    ## discard dead volume
    time_index = 0
    for i in range(len(retention_time)):
        time_index = i
        if retention_time[i] >= DEAD_VOLUME_TIME:
            break

    retention_time = retention_time[time_index:]
    processed_data = processed_data[time_index:]

    # %%
    ## Calculate the differences between intensity data points (dy)
    diff = np.diff(processed_data)
    scale = max(processed_data) / max(diff)
    diff = diff * scale


    # calculate base derivatives
    base_diff_data = []
    for i in range(len(diff)):
        if abs(diff[i]) <= 1:
            base_diff_data.append(diff[i])
    base_diff = np.average(base_diff_data)


    # picking out peaks who has at least one point above the detection limit
    data_process_index = 0
    maxima = []
    maximum_retention_time = []
    peak_area_min = []

    for i in range(len(processed_data)):
        already_processed = data_process_index >= i
        if already_processed:
            pass
        else:
            if processed_data[i] - base_intensity >= DETECTION_LIMIT:
                # find a peak who has intensity greater than DETECTION LIMIT
                search_index_left = i
                search_index_right = i
                #search along the left slide until it goes below the baseline/is flat
                while abs(diff[search_index_left]) >= base_diff + 0.3 and processed_data[
                    search_index_left] >= base_intensity:
                    search_index_left -= 1
                # search along the right slide until it goes below the baseline/is flat
                while (abs(diff[search_index_right]) >= base_diff + 0.3 or processed_data[
                    search_index_right] >= DETECTION_LIMIT) \
                        and processed_data[search_index_right] >= base_intensity:
                    search_index_right += 1
                data_process_index = search_index_right

                # current peak analysis
                current_peak_time = retention_time[search_index_left:search_index_right + 1]
                current_peak_intensity = processed_data[search_index_left:search_index_right + 1]


                if plot:
                    fn = UnivariateSpline(current_peak_time, current_peak_intensity, k=5)

                    current_peak_intensity_spline = fn(current_peak_time)
                    plt.plot(current_peak_time, current_peak_intensity_spline)
                # find the peak max and time

                max_int_point = 0
                max_int_time = 0

                for i in range(len(current_peak_intensity)):
                    if current_peak_intensity[i] >= max_int_point:
                        max_int_time = current_peak_time[i]
                        max_int_point = current_peak_intensity[i]

                maximum_retention_time.append(max_int_time)
                maxima.append(max_int_point)

                #
                # #integrate and append the peak information into the current
                current_peak_intensity -= base_intensity

                peak_area_min.append(integration(current_peak_time, current_peak_intensity))


    # convert peak per min to peak per second
    peak_area_sec = []
    for area in peak_area_min:
        peak_area_sec.append(area * 60)

    # plot the peak line
    if plot:
        plt.plot(retention_time, processed_data, '--')
        plt.xlim(1, 2.5)
        plt.show()

    # %%

    # Recognize the internal standard peak
    internal_standard_peak_area = 1
    for i in range(len(peak_area_min)):
        if abs(maximum_retention_time[i] - INTERNAL_STANDARD_RETENTION_TIME) <= PEAK_WIDTH:
            internal_standard_peak_area = peak_area_min[i]
        else:
            pass

    # ratio calculation based on the internal standard
    peak_ratio = []
    for area in peak_area_min:
        peak_ratio.append(area / internal_standard_peak_area)

    return maximum_retention_time, peak_ratio, peak_area_min


# %%


def experimentally_monitored_data(folder: str,
                                  peak_width: float = 0.04,
                                  max_data_point_amount: int = 200):
    """
    From the top level directory, search the blank and all data contained in the folder
    :param folder: the folder contains data. Usually Names LJL + Datatime
    :param peak_width: the tolerance of peak shifting
    :param max_data_point_amount: maximum analysis data point
    :return: peak retention times, time point , peak ratios for plotting
    """

    blank_data_210 = []  # in list form. First is retention time, and second is intensity, third datetime

    reaction_data_210 = []  ## in list of list form
    is_there_blank = False
    reaction_number = '0'  ## implemented but not in use right now
    peak_width = peak_width
    folder = folder
    max_data_point_amount = max_data_point_amount

    # search for blank data
    for file in os.listdir(folder):

        if fnmatch.fnmatch(file, '*blank*'):
            is_there_blank = True

            blank_csv_file = os.path.join(folder, file)
            time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
            retention_time, blank_intensity = extract_data(blank_csv_file)
            blank_data_210.append(retention_time)
            blank_data_210.append(blank_intensity)
            blank_data_210.append(extract_time(time_csv_file))

    # search for injection data by sequence
    for i in range(max_data_point_amount):
        for file in os.listdir(folder):
            if fnmatch.fnmatch(file, "*-" + str(i) + "-*"):
                if is_there_blank:
                    reaction_number = file[file.find(' ') + 1:file.find(' ') + 4]

                reaction_csv_file = os.path.join(folder, file)
                time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
                try:
                    retention_time, intensity = extract_data(reaction_csv_file)
                    time = extract_time(time_csv_file)
                except FileNotFoundError:
                    retention_time, intensity, time = [], [], []
                reaction_data_210.append([retention_time, intensity, time])

    data_amount = len(reaction_data_210)

    # %%

    # calculate the timeline of each data point
    time_point = [0]
    for i in range(len(reaction_data_210) - 1):
        time_point.append((reaction_data_210[i][2] - blank_data_210[2]).seconds / 3600)

    # %%

    distinct_peak_retention_time = [] # Overall, which peaks appeared?
    peak_ratio = [] # What are their ratios with respect to the internal standard?

    # Picks all all appeared peaks
    for i in range(data_amount):
        d = [reaction_data_210[i][0], reaction_data_210[i][1]]
        times, areas, _ = peak_properties(blank_data_210, d)
        for time in times:
            is_in = False
            for t in distinct_peak_retention_time:
                if abs(t - time) <= peak_width:
                    is_in = True

            if not is_in:
                distinct_peak_retention_time.append(time)

    # create 5 empty arrays
    for p in distinct_peak_retention_time:
        peak_ratio.append([])

    # %%
    # fill the spot of each data according to each injection
    # If the current injection does not have this peak, insert 0 in that spot
    for z in range(data_amount):
        times, areas, _ = peak_properties(blank_data_210, reaction_data_210[z])
        for i in range(len(times)):
            for j in range(len(distinct_peak_retention_time)):
                if abs(distinct_peak_retention_time[j] - times[i]) < peak_width:
                    peak_ratio[j].append(areas[i])
        for k in range(len(peak_ratio)):
            if len(peak_ratio[k]) == z:
                peak_ratio[k].append(0)

    return distinct_peak_retention_time, time_point, peak_ratio


def get_last_experimental_data(folder: str = r"/Users/luke/Desktop/LJL 2021-01-31 01-51-50 3/",
                               max_data_point_amount=200):
    """
    Find the last injection's data point
    :param folder: The top level directory
    :param max_data_point_amount: the maximum of the data points
    :return: the retention times and peak ratios of all detected peaks
    """
    blank_data_210 = []  # in list form. First is retention time, and second is intensity, third datetime

    reaction_data_210 = []  ## in list of list form
    is_there_blank = False

    for file in os.listdir(folder):

        if fnmatch.fnmatch(file, '*blank*'):
            is_there_blank = True

            blank_csv_file = os.path.join(folder, file, 'DAD1 A 210 nm (4 nm).csv')
            time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
            retention_time, blank_intensity = extract_data(blank_csv_file)
            blank_data_210.append(retention_time)
            blank_data_210.append(blank_intensity)
            blank_data_210.append(extract_time(time_csv_file))

    for i in range(max_data_point_amount):
        for file in os.listdir(folder):
            if fnmatch.fnmatch(file, "*-" + str(i) + "-hydra*"):
                if is_there_blank:
                    reaction_number = file[file.find(' ') + 1:file.find(' ') + 4]

                reaction_csv_file = os.path.join(folder, file, 'DAD1 A 210 nm (4 nm).csv')
                time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
                try:
                    retention_time, intensity = extract_data(reaction_csv_file)
                    time = extract_time(time_csv_file)
                except FileNotFoundError:
                    retention_time, intensity, time = [], [], []
                reaction_data_210.append([retention_time, intensity, time])

    i = len(reaction_data_210) - 1
    while reaction_data_210[len(reaction_data_210) - 1][0] == []:
        i -= 1

    retention_time, peak_ratio, _ = peak_properties(blank_data_210, reaction_data_210[i])
    return retention_time, peak_ratio


def is_there_peak_in_collection(retention_times: dict, peak_retention_time, shift_tolerance=0.04):
    """
    If the peak is already labelled
    :param retention_times: The dictionary containing all labels and peaks. Normally inside of the deck
    :param peak_retention_time: The target peak retention time
    :param shift_tolerance: The error of detecting that peak
    :return: True if the peak is already identified
    """
    for peaks in retention_times:
        if abs(retention_times[peaks] - peak_retention_time) <= shift_tolerance:
            return True
    return False

def is_there_peak(peak,shift_tolerance=0.04):
    """
    If the last data contains a give peak
    REQUIRES: deck.reaction_folder is updated
    :param peak: target retention time of a peak
    :param shift_tolerance: room for peak shift
    :return: True if the last experimental data contains the give peak
    """
    retention_times,_ = get_last_experimental_data(deck.reaction_folder)
    for time in retention_times:
        if abs(peak-time)<=shift_tolerance:
            return True
    return False

def identify_appeared_peak(peaks: dict, label: str, shift_tolerance=0.04):
    """
    Name one appearing peak with a given label,
    Use only exactly one new peak appeared
    :param peaks: The dictionary containing all labels and peaks. Normally inside of the deck
    :param label: The desired label
    :param shift_tolerance: room for peak shift
    :return: None
    """
    retention_time, _ = get_last_experimental_data(deck.reaction_folder)
    for time in retention_time:
        if not is_there_peak_in_collection(peaks, time, shift_tolerance=shift_tolerance):
            peaks[label] = time


def label_one_peak(peaks: dict, label: str, peak_retention_time, shift_tolerance=0.04):
    """
    label one peak with given peak retention time
    :param peaks: The dictionary containing all labels and peaks. Normally inside of the deck
    :param label: The desired label
    :param peak_retention_time: The target peak retention time
    :param shift_tolerance: room for peak shift
    :return: None
    """
    retention_time, _ = get_last_experimental_data(deck.reaction_folder)

    for time in retention_time:
        if abs(peak_retention_time - time) <= shift_tolerance:
            peak_retention_time = time

    peaks[label] = peak_retention_time

def label_appeared_increasing_peak(peaks: dict,label: str,data_away_from_last :int=5):
    """
    Label a peak if it is increasing and have not been labelled before
    :param peaks: The dictionary containing all labels and peaks. Normally inside of the deck
    :param label: The desired label
    :param data_away_from_last: distance of data point from the last data to determine the changing behaviour
    :return: None
    """
    distinct_peak_retention_time, time_point, peak_ratio= experimentally_monitored_data(deck.reaction_folder)
    for i in range(len(distinct_peak_retention_time)):
        current_peak = peak_ratio[i]
        if not is_there_peak_in_collection(peaks, distinct_peak_retention_time[i]):
            print(current_peak)
            if current_peak[len(current_peak)-data_away_from_last-1] - current_peak[len(current_peak)-1] <= 0:
                peaks[label] = distinct_peak_retention_time[i]
                return


def label_appeared_decreasing_peak(peaks: dict, label: str, data_away_from_last: int = 5):
    """
        Label a peak if it is increasing and have not been labelled before
        :param peaks: The dictionary containing all labels and peaks. Normally inside of the deck
        :param label: The desired label
        :param data_away_from_last: distance of data point from the last data to determine the changing behaviour
        :return: None
        """
    distinct_peak_retention_time, time_point, peak_ratio = experimentally_monitored_data(deck.reaction_folder)
    for i in range(len(distinct_peak_retention_time)):
        current_peak = peak_ratio[i]
        if not is_there_peak_in_collection(peaks, distinct_peak_retention_time[i]):
            if current_peak[len(current_peak) - data_away_from_last - 1] - current_peak[len(current_peak) - 1] >= 0:

                peaks[label] = distinct_peak_retention_time[i]
                return

def update_for_hplc_folder():
    """
    watch for new appearing files so that the program knows the current
    hplc folder
    :return: None
    """
    w = Watcher()
    w.run()