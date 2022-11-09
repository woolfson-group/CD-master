import CD_master
import datafiles_opener as dfo

# Set up the working directory
input_dir = './input1'
files = dfo.open_files(input_dir, 'txt')

# create an example of the class
exp = CD_master.cdExp()

# use one of the default methods for processing conc dependencies
exp.process_experiment(files=files, mode='single_scans')

# set up the order of samples how they'll be plotted
# ordered_samples = [1,4,5,3]
# exp.set_plot_order(ordered_samples)

# set up labels for the curves
conc = [val['Conc_uM'] for val in exp.get_samples().values()]
labels = [str(c) + r' $\mu$M' for c in conc[:]]  # make labels with uM units
exp.set_labels(labels)

# plot the data
# give it some settings how to plot the data
exp.spec_plotter(ylim=((-20, 5, 5),(200,1200, 200)), fig_title='Test figure\n', dpi=100, lw=2)

# export the final datasets [WL, MRE, HT]
# exp.export_all_MRE()