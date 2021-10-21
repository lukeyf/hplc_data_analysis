# %%
import datetime
import fnmatch
import os
import pathlib
from statistics import mean

import pandas as pd
from aghplctools.data.sample import HPLCSample

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline

DEAD_VOLUME_TIME = 0.6


def extract_data(filename: str, wavelength_nm=210):
    # f_blank = pd.read_csv(filename).to_numpy()
    # retention_time = f_blank[1:, 0]
    # intensity_blank_data = f_blank[1:, 1]
    """
        get the data of time and intensity from the home folder
        :param filename: the name of the home
        :return: the time points and intensities of the data
        """
    folder = pathlib.Path(filename)
    data = HPLCSample.create_from_D_file(folder)
    signal = None
    for s in data.signals:
        if int(s.wavelength) == wavelength_nm:
            signal = s
    if signal is not None:
        retention_time = signal.retention_times
        intensity_data = signal.mean_unreferenced_intensities
        return retention_time, intensity_data
    else:
        raise ('Wavelength not found.')


def extract_time(filename: str = 'Automatically_Generated_Report00.CSV'):
    f = pd.read_csv(filename, encoding='utf-16').to_numpy()

    date = f[8][1]
    return datetime.datetime.strptime(date, '%d-%b-%y, %H:%M:%S')


def integration(x, y):
    # interpolation
    fn = UnivariateSpline(x, y, k=5)

    results = quad(fn, x[0], x[-1])[0]
    return results


def peak_properties(blank_data: list, reaction_data: list, internal_standard_retention_time: float = 1.81,
                    peak_width: float = 0.04, detection_limit: float = 20, plot: bool = True, save_plot=False,
                    fig_path='', labels={}, plot_range: list = None):
    '''
    The peak properties from a given blank data, peak raw data
    :return:
    a list of peak maximum intensities, areas, retention time (min)
    '''
    retention_time, intensity_blank_data = blank_data[0], blank_data[1]
    INTERNAL_STANDARD_RETENTION_TIME = internal_standard_retention_time  # min
    PEAK_WIDTH = peak_width  # min
    DETECTION_LIMIT = detection_limit
    intensity_raw_data = reaction_data[1][:len(retention_time)]
    legends = []

    # Calculate baseline
    processed_data = intensity_raw_data - intensity_blank_data
    base_data = []
    for i in range(len(processed_data)):
        if abs(processed_data[i] - np.average(processed_data)) <= 3:
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

    diff = np.diff(processed_data)
    scale = max(processed_data) / max(diff)
    diff = diff * scale

    # calculate base derivatives

    base_diff_data = []
    for i in range(len(diff)):
        if abs(diff[i]) <= 1:
            base_diff_data.append(diff[i])

    base_diff = np.average(base_diff_data)

    # %%

    # picking out peaks

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
                # find a peak who has intensity greater than 5
                search_index_left = i
                search_index_right = i
                while abs(diff[search_index_left]) >= base_diff + 0.3 and processed_data[
                    search_index_left] >= base_intensity:
                    search_index_left -= 1
                while (abs(diff[search_index_right]) >= base_diff + 0.3 or processed_data[
                    search_index_right] >= DETECTION_LIMIT) \
                        and processed_data[search_index_right] >= base_intensity:
                    search_index_right += 1
                data_process_index = search_index_right

                # current peak analysis
                current_peak_time = retention_time[search_index_left:search_index_right + 1]
                current_peak_intensity = processed_data[search_index_left:search_index_right + 1]

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
                # #integrate
                try:
                    current_peak_intensity -= base_intensity

                    peak_area_min.append(integration(current_peak_time, current_peak_intensity))
                    fn = UnivariateSpline(current_peak_time, current_peak_intensity)
                    current_peak_intensity_spline = fn(current_peak_time)
                    if plot or save_plot:
                        plt.plot(current_peak_time, current_peak_intensity_spline)
                        for label in labels:
                            if abs(max_int_time - labels[label]) <= peak_width:
                                legends.append(label)
                                break
                except Exception:
                    maximum_retention_time = maximum_retention_time[:-1]
    peak_area_sec = []
    for area in peak_area_min:
        peak_area_sec.append(area * 60)
    if plot or save_plot:
        plt.plot(retention_time, processed_data, '--')
        if plot_range is not None:
            plt.xlim(plot_range[0], plot_range[1])
        plt.legend(legends)
        if plot:
            plt.show()
        if save_plot:
            plt.savefig(os.path.join(fig_path))
            plt.close()

    # %%

    # Recognize peaks
    internal_standard_peak_area = 1
    for i in range(len(peak_area_min)):
        if abs(maximum_retention_time[i] - INTERNAL_STANDARD_RETENTION_TIME) <= PEAK_WIDTH:
            internal_standard_peak_area = peak_area_min[i]
        else:
            pass

    # ratio calculation
    peak_ratio = []
    for area in peak_area_min:
        peak_ratio.append(area / internal_standard_peak_area)

    return maximum_retention_time, peak_ratio, peak_area_min


# %%


def experimentally_monitored_data(folder: str,
                                  peak_width: float = 0.04,
                                  max_data_point_amount: int = 200,
                                  plot: bool = False, ):
    """
    :param plot:
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
    for file in os.listdir(folder):
        if fnmatch.fnmatch(file, '*blank*'):
            is_there_blank = True

            blank_csv_file = os.path.join(folder, file)
            time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
            retention_time, blank_intensity = extract_data(blank_csv_file)
            blank_data_210.append(retention_time)
            blank_data_210.append(blank_intensity)
            blank_data_210.append(extract_time(time_csv_file))

    for i in range(max_data_point_amount):
        print(i)
        for file in os.listdir(folder):
            if fnmatch.fnmatch((file[:3]), '{0:03}'.format(i + 2)):
                if is_there_blank:
                    reaction_number = file[file.find(' ') + 1:file.find(' ') + 4]

                reaction_csv_file = os.path.join(folder, file)
                time_csv_file = os.path.join(folder, file, '../sample data/Automatically_Generated_Report00.CSV')
                try:
                    retention_time, intensity = extract_data(reaction_csv_file)
                    time = extract_time(time_csv_file)

                    reaction_data_210.append([retention_time, intensity, time])
                except OSError:
                    pass
                    # retention_time, intensity, time = [], [], []

    data_amount = len(reaction_data_210)
    # %%

    # %%

    time_point = []
    for i in range(len(reaction_data_210)):
        time_point.append((reaction_data_210[i][2] - blank_data_210[2]).seconds / 3600)

    # %%

    distinct_peak_retention_time = []
    peak_ratio = []
    peak_concentration = []
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

    for p in distinct_peak_retention_time:
        peak_ratio.append([])
        peak_concentration.append([])

    # %%

    for z in range(data_amount):
        times, areas, concentration = peak_properties(blank_data_210, reaction_data_210[z], plot=plot)
        for i in range(len(times)):
            for j in range(len(distinct_peak_retention_time)):
                if abs(distinct_peak_retention_time[j] - times[i]) < peak_width:
                    peak_ratio[j].append(areas[i])
        for k in range(len(peak_ratio)):
            if len(peak_ratio[k]) == z:
                peak_ratio[k].append(0)

        for i in range(len(times)):
            for j in range(len(distinct_peak_retention_time)):
                if abs(distinct_peak_retention_time[j] - times[i]) < peak_width:
                    peak_concentration[j].append(concentration[i])
        for k in range(len(peak_concentration)):
            if len(peak_concentration[k]) == z:
                peak_concentration[k].append(0)
    return distinct_peak_retention_time, time_point, peak_ratio, peak_concentration


def get_last_experimental_data(folder: str = r"/Users/luke/Desktop/LJL 2021-01-31 01-51-50 3/",
                               max_data_point_amount=200):
    blank_data_210 = []  # in list form. First is retention time, and second is intensity, third datetime

    reaction_data_210 = []  ## in list of list form
    is_there_blank = False

    for file in os.listdir(folder):

        if fnmatch.fnmatch(file, '*-NV-*'):
            is_there_blank = True

            blank_csv_file = os.path.join(folder, file)
            time_csv_file = os.path.join(folder, file, '../sample data/Automatically_Generated_Report00.CSV')
            retention_time, blank_intensity = extract_data(blank_csv_file)
            blank_data_210.append(retention_time)
            blank_data_210.append(blank_intensity)
            blank_data_210.append(extract_time(time_csv_file))
    max = 1
    for i in range(max_data_point_amount):
        for file in os.listdir(folder):
            if fnmatch.fnmatch((file[:3]), '{0:03}'.format(i + 2)):
                if is_there_blank:
                    reaction_number = file[file.find(' ') + 1:file.find(' ') + 4]

                reaction_csv_file = os.path.join(folder, file)
                time_csv_file = os.path.join(folder, file, '../sample data/Automatically_Generated_Report00.CSV')
                try:
                    retention_time, intensity = extract_data(reaction_csv_file)
                    time = extract_time(time_csv_file)
                    reaction_data_210.append([retention_time, intensity, time])
                except:
                    pass
                    # retention_time, intensity, time = [], [], []

                max = i

    retention_time, peak_ratio, _ = peak_properties(blank_data_210, reaction_data_210[max - 2])
    return retention_time, peak_ratio


def get_the_experimental_data(folder: str = r"/Users/luke/Desktop/LJL 2021-01-31 01-51-50 3/",
                              inj_number=1,
                              save_fig=False,
                              labels={}):
    blank_data_210 = []  # in list form. First is retention time, and second is intensity, third datetime

    reaction_data_210 = []  ## in list of list form

    for file in os.listdir(folder):

        if fnmatch.fnmatch(file, '*-NV-*'):
            blank_csv_file = os.path.join(folder, file)
            time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
            retention_time, blank_intensity = extract_data(blank_csv_file)
            blank_data_210.append(retention_time)
            blank_data_210.append(blank_intensity)
            blank_data_210.append(extract_time(time_csv_file))

    for file in os.listdir(folder):
        if fnmatch.fnmatch((file[:3]), '{0:03}'.format(inj_number + 1)):

            reaction_csv_file = os.path.join(folder, file)
            time_csv_file = os.path.join(folder, file, 'Automatically_Generated_Report00.CSV')
            try:
                retention_time, intensity = extract_data(reaction_csv_file)
                time = extract_time(time_csv_file)
                reaction_data_210.append([retention_time, intensity, time])
            except:
                pass

    temp_folder = os.path.join(folder, 'temp')
    fig_name = os.path.join(temp_folder, 'chromatogram_' + str(inj_number) + '.png')
    retention_time, peak_ratio, _ = peak_properties(blank_data_210, reaction_data_210[0], save_plot=save_fig,
                                                    fig_path=fig_name, labels=labels)
    return retention_time, peak_ratio


def range_integration(folderpath: str, range_of_interest: list, wavelength_nm: int = 310):
    score = 0
    x, y = extract_data(folderpath, wavelength_nm=wavelength_nm)
    time, _, integration = peak_properties([x, [mean(y)] * len(y)], [x, y])

    for i in range(len(time)):
        if range_of_interest[0] < time[i] < range_of_interest[1]:
            score += integration[i]

    return score


def evaluate_performance(parent_folder: str,
                         keyword: str,
                         range_of_interest: list,
                         wavelength_nm: int = 310):
    """

    :param parent_folder: The HPLC data folder path
    :param keyword: the data keyword (set in HPLC software) for finding the desired data
    :param range_of_interest: The integration range in min
                              (e.g. [8.5,9] means integrating peaks from 8.5 min to 9 min)
    :param wavelength_nm: The wavelength of data that the integration takes place in
                              Default looking at 310 nm for Au13 integration
    :return: a set of numbers that indicates the integrations of experiments in a given range
    """

    results = []
    contents = os.listdir(parent_folder)
    for i in range(len(contents)):
        for file in contents:
            if fnmatch.fnmatch(file[:3], '{0:03}'.format(i + 1)) and fnmatch.fnmatch(file, '*' + keyword + '*'):
                child_folder = os.path.join(parent_folder, file)
                result = range_integration(child_folder, range_of_interest, wavelength_nm)
                results.append(result)
                break

    return results
