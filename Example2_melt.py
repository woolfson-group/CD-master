import CD_master
import datafiles_opener as dfo
import seaborn as sns
import colorcet as cc
# Set up the working directory
input_dir = './input12'
files = dfo.open_files(input_dir, 'txt')

# create an example of the class
exp = CD_master.cdExp()
exp.process_experiment(files, mode='var_temp')

# plot the data
labs = ['Cool', 'Melt']

# give it some settings how to plot the data
exp.t_plotter(ylim=((-10, 0, 2),(400, 600, 100)), xlim=(5,85,10),fig_title='25uM AH213\n',
                 labels=labs, palette=sns.color_palette(cc.glasbey, n_colors=24))
