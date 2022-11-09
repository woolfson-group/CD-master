import CD_master
import datafiles_opener as dfo

# Set up the working directory
# input_dir = './input1'
files = dfo.open_files(input_dir, 'txt')

# create an example of the class
exp = CD_master.cdExp()

# use one of the default methods for plotting datasets of the format [WL, MRE, HT]
exp.read_MRE_datafiles(files=files)

# set up the order of samples how they'll be plotted
# ordered_samples = [1,2,3]
# exp.set_plot_order(ordered_samples)

# set up labels for the curves
conc = [val['Conc_uM'] for val in exp.get_samples().values()]
labels = [str(c) + r' $\mu$M' for c in conc[:]]

exp.set_labels(labels)

# plot the data
exp.spec_plotter(ylim=((-20, 5, 5),(200,1000, 200)), xlim=(200, 260, 10), fig_title='Test figure\n', dpi=100, lw=2)
