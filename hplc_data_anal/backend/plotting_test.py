import fnmatch
import os

from hplc_data_anal.backend.analysis_method import experimentally_monitored_data,get_the_experimental_data
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import imageio

reaction_folder = r'D:\Chemstation\1\Data\2021-02-23\LJL 2021-02-23 17-48-47'

labels = {'Internal Standard': 1.81, 'Acid': 1.3179166666666666, 'Intermediate': 1.5554166666666667,
          'Amine': 1.204166666666667, 'Amide': 2.00375}

def plot_referenced_chromatograph(folder,inj_num,labels):
    get_the_experimental_data(folder,inj_num,save_fig=False,labels = labels)


def plot_reaction_data_ratio(reaction_folder, labels = labels, save_plot: bool = True, save_data: bool = True, show_plot=False,
                             plot_until=200):
    times, time, ratio, _ = experimentally_monitored_data(reaction_folder, peak_width=0.05,
                                                          max_data_point_amount=plot_until)
    for i in range(len(times)):
        plt.scatter(time, ratio[i])
    legends = []
    for time in times:
        for label in labels:
            if abs(labels[label] - time) <= 0.06:
                legends.append(label)
    plt.legend(legends)
    if save_plot:
        if plot_until < 200:
            temp_folder = os.path.join(os.path.join(reaction_folder, 'temp'))
            try:
                os.mkdir(temp_folder)
            except:
                pass
            plt.savefig(os.path.join(os.path.join(reaction_folder, 'temp'), 'progress_' + str(plot_until) + '.png'))
        else:
            plt.savefig(os.path.join(reaction_folder, 'reaction_progress_conc_min.pdf'))
    if save_data:
        inverse_ratio = np.array(ratio).T.tolist()
        df1 = pd.DataFrame(inverse_ratio,
                           columns=legends)
        df1.to_excel(os.path.join(reaction_folder, 'progress_data_ratio.xlsx'))
    if show_plot:
        plt.show()
    else:
        plt.close()


def plot_reaction_data_conc(reaction_folder, labels, save_plot: bool = True, save_data: bool = True, show_plot=False,
                            plot_until=200):
    times, time, _, conc = experimentally_monitored_data(reaction_folder, peak_width=0.1,
                                                         max_data_point_amount=plot_until)
    for i in range(len(times)):
        plt.scatter(time, conc[i])
    legends = []
    for time in times:
        for label in labels:
            if abs(labels[label] - time) <= 0.06:
                legends.append(label)
    plt.legend(legends)
    if save_plot:
        if plot_until < 200:
            temp_folder = os.path.join(os.path.join(reaction_folder, 'temp'))
            try:
                os.mkdir(temp_folder)
            except:
                pass
            plt.savefig(os.path.join(os.path.join(reaction_folder, 'temp'), 'progress_' + str(plot_until) + '.png'))
        else:
            plt.savefig(os.path.join(reaction_folder, 'reaction_progress_conc_min.pdf'))

    if save_data:
        inverse_conc = np.array(conc).T.tolist()
        df1 = pd.DataFrame(inverse_conc,
                           columns=legends)
        df1.to_excel(os.path.join(reaction_folder, 'progress_data_conc_min.xlsx'))
    if show_plot:
        plt.show()
    else:
        plt.close()


def get_reaction_progress_gif(reaction_folder, labels, max_data_amount=200, type='ratio',seconds_per_frame=0.5,loop =1):
    data_amount = get_data_amount(max_data_amount, reaction_folder)
    for i in range(data_amount):
        if type == 'ratio':
            try:
                plot_reaction_data_ratio(reaction_folder, labels, save_plot=True, save_data=False, plot_until=i + 1)
            except:
                break
        if type == 'conc':
            try:
                plot_reaction_data_conc(reaction_folder, labels, save_plot=True, save_data=False, plot_until=i + 1)
            except:
                break

    # GIF Generation
    temp_folder = os.path.join(reaction_folder,'temp')
    images = []
    filenames = []
    for i in range(data_amount):
        try:
            filenames.append(os.path.join(temp_folder,"progress_" + str(i+1) + ".png"))
        except:
            pass
    for filename in filenames:
        try:
            images.append(imageio.imread(filename))
        except:
            pass
    imageio.mimsave(os.path.join(reaction_folder,'progress_'+type+'.gif'), images,duration = seconds_per_frame,loop=loop)


def get_data_amount(max_data_amount, reaction_folder):
    data_amount = 0
    for i in range(max_data_amount):
        does_data_point_exist = False
        for file in os.listdir(reaction_folder):
            if fnmatch.fnmatch((file[:3]), '{0:03}'.format(i + 2)):
                does_data_point_exist = True
        if does_data_point_exist:
            data_amount = i + 2
        else:
            break
        # %%
    return data_amount


def get_chromatogram_gif(reaction_folder, labels, max_data_amount=200,seconds_per_frame=0.5,loop =1):
    data_amount = get_data_amount(max_data_amount, reaction_folder)
    for i in range(data_amount):
        try:
            plot_referenced_chromatograph(reaction_folder,i,labels)
        except:
            pass
    temp_folder = os.path.join(reaction_folder, 'temp')
    images = []
    filenames = []
    for i in range(data_amount):
        filenames.append(os.path.join(temp_folder, "chromatogram_" + str(i + 1) + ".png"))
    for filename in filenames:
        try:
            images.append(imageio.imread(filename))
        except:
            pass

    imageio.mimsave(os.path.join(reaction_folder, 'chromatogram' '.gif'), images, duration=seconds_per_frame,
                    loop=loop)

# if __name__ == '__main__':

    # get_reaction_progress_gif(r'D:\Chemstation\1\Data\2021-03-19\LJL 2021-03-19 11-03-50',
    #                      {'Internal Standard': 1.8129166666666667, 'Acid': 1.31625,
    #                       'Intermediate': 1.5704166666666666,
    #                       'Amide': 2.012,
    #                       'Amine': 1.229})
    # get_chromatogram_gif(r'D:\Chemstation\1\Data\2021-01-11\LJL 2021-01-11 11-21-35',
    #                      {'Internal Standard': 1.770166666666667, 'Acid': 1.585,
    #                       'Intermediate': 2.025,
    #                       'Amide': 2.012,
    #                       'Amine': 1.229})
    # plot_referenced_chromatograph(r'D:\Chemstation\1\Data\2021-03-19\LJL 2021-03-19 11-03-50',19,
    #                               {'Internal Standard': 1.8129166666666667, 'Acid': 1.31625, 'Intermediate': 1.5704166666666666,
    #                                'Amide' :2.012,
    #                                'Amine' :1.229})\