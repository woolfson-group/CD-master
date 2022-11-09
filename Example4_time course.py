import math

import CD_master
import datafiles_opener as dfo
import seaborn as sns
import colorcet as cc
# Set up the working directory
input_dir = './input4'
files = dfo.open_files(input_dir, 'txt')

# create an example of the class
exp = CD_master.cdExp()
exp.process_experiment(files, mode='kinet_scans')

# plot the data
labels = [str(p) + " min" for p in exp.get_samples().keys()]
exp.set_labels(labels)
every_nth = 5
exp.spec_plotter(ylim=((-50, 30, 20),(200,1200, 200)), xlim=(190,260,10),fig_title='Test figure\n',
                 palette=sns.color_palette('crest', math.ceil(exp.get_n_samples()/every_nth)), every_nth=every_nth)
# exp.print_off_dict(exp.get_samples())