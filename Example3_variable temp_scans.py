import CD_master
import datafiles_opener as dfo
import seaborn as sns
import colorcet as cc
# Set up the working directory
input_dir = './input3'

files = dfo.open_files(input_dir, 'txt')

# create an example of the class
exp = CD_master.cdExp()

# use one of the default methods for processing conc dependencies
exp.process_experiment(files, mode='var_temp_scans')

labels = [str(p) + u"\N{DEGREE SIGN}C" for p in exp.get_samples().keys()]
exp.set_labels(labels)
# exp.print_off_dict(exp.get_samples())
exp.spec_plotter(ylim=((-20,10,5),(200,1200,200)), palette=sns.color_palette("crest", exp.get_n_samples()),
fig_title='',dpi=100, legend_bool=False, lw=2)


                 
