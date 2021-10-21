import os
import pickle
import plotting_test as deck
from hplc_data_anal.backend.analysis_method import get_last_experimental_data, experimentally_monitored_data
from watcher import Watcher
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def is_there_peak_in_collection(retention_times: dict, peak_retention_time, shift_tolerance=0.04):
    for peaks in retention_times:
        if abs(retention_times[peaks] - peak_retention_time) <= shift_tolerance:
            return True
    return False


def is_there_peak(peak, shift_tolerance=0.04):
    retention_times, _ = get_last_experimental_data(deck.reaction_folder)
    for time in retention_times:
        if abs(peak - time) <= shift_tolerance:
            return True
    return False


def identify_appeared_peak(peaks: dict, label: str, shift_tolerance=0.04):
    """
    Name one appearing peak with a given label,
    Use only exactly one new peak appeared
    :param peaks:
    :param label:
    :param shift_tolerance:
    :return:
    """
    retention_time, _ = get_last_experimental_data(deck.reaction_folder)
    for time in retention_time:
        if not is_there_peak_in_collection(peaks, time, shift_tolerance=shift_tolerance):
            peaks[label] = time


def label_one_peak(peaks: dict, label: str, peak_retention_time, shift_tolerance=0.04):
    retention_time, _ = get_last_experimental_data(deck.reaction_folder)

    for time in retention_time:
        if abs(peak_retention_time - time) <= shift_tolerance:
            peak_retention_time = time

    peaks[label] = peak_retention_time


def label_appeared_increasing_peak(peaks: dict, label: str, data_away_from_last: int = 5):
    distinct_peak_retention_time, time_point, peak_ratio, _ = experimentally_monitored_data(deck.reaction_folder)
    for i in range(len(distinct_peak_retention_time)):
        current_peak = peak_ratio[i]
        if not is_there_peak_in_collection(peaks, distinct_peak_retention_time[i]):
            if current_peak[len(current_peak) - data_away_from_last - 1] - current_peak[len(current_peak) - 1] <= 0:
                peaks[label] = distinct_peak_retention_time[i]
                return


def label_appeared_decreasing_peak(peaks: dict, label: str, data_away_from_last: int = 5):
    distinct_peak_retention_time, time_point, peak_ratio, _ = experimentally_monitored_data(deck.reaction_folder)
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
    :return:
    """
    w = Watcher()
    w.run()


def save_json_file(json, file):
    with open(file, 'wb') as fp:
        pickle.dump(json, fp, protocol=pickle.HIGHEST_PROTOCOL)


def load_json_file(file):
    with open(file, 'rb') as fp:
        data = pickle.load(fp)
    return data


def plot_reaction_data_ratio(reaction_folder, labels, save_plot: bool = True, save_data: bool = True, show_plot=False):
    times, time, ratio, _ = experimentally_monitored_data(reaction_folder, peak_width=0.1)
    for i in range(len(times)):
        plt.scatter(time, ratio[i])
    legends = []
    for time in times:
        for label in labels:
            if abs(labels[label] - time) <= 0.06:
                legends.append(label)
    plt.legend(legends)
    if save_plot:
        plt.savefig(os.path.join(reaction_folder, 'reaction_progress_ratio.pdf'))
    if save_data:
        inverse_ratio = np.array(ratio).T.tolist()
        df1 = pd.DataFrame(inverse_ratio,
                           columns=legends)
        df1.to_excel(os.path.join(reaction_folder, 'progress_data_ratio.xlsx'))
    if show_plot:
        plt.show()


def plot_reaction_data_conc(reaction_folder, labels, save_plot: bool = True, save_data: bool = True, show_plot=False):
    times, time, _, conc = experimentally_monitored_data(reaction_folder, peak_width=0.1)
    for i in range(len(times)):
        plt.scatter(time, conc[i])
    legends = []
    for time in times:
        for label in labels:
            if abs(labels[label] - time) <= 0.06:
                legends.append(label)
    plt.legend(legends)
    if save_plot:
        plt.savefig(os.path.join(reaction_folder, 'reaction_progress_conc_min.pdf'))
    if save_data:
        inverse_conc = np.array(conc).T.tolist()
        df1 = pd.DataFrame(inverse_conc,
                           columns=legends)
        df1.to_excel(os.path.join(reaction_folder, 'progress_data_conc_min.xlsx'))
    if show_plot:
        plt.show()
