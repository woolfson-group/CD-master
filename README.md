# CD-master
A class for processing CD data of various formats

Dependenices:
- numpy
- matplotlib
- seaborn
- colorcet

Tested with:
Python-3.9, numpy-1.21.2, matplotlib-3.4.3, seaborn-0-11.2, colorcet-3.0.0
Windows 10 64bit


# How to use

## Input data

Your input files have to be named correctly in order to be processed with the class.
The naming rules are:
<date>_<sample name>_<concentration in uM, mM, or M>_<pathlength in mm>_<temeperature with units>_<pH>_<WL>_<time point>.txt
- File name is parsed expecting the underscore (_) as a separator, so it'll fail if any other seperator is used
- Parameters shouldn NOT have underscores within themselves
- Data and sample name can have any format
- Concentration should be with units, e.g. 250uM
- Pathlength should be with units in mm, e.g. 1mm
- A single temperature (e.g. for spectra) should be with units, either K or C, e.g. 5C
- Temperature range (e.g. for thermal denaturation experiments) should be with units, either K or C, e.g. 5-90C
- pH should be specified like this: pH7.5

Depending on the experiment type, the required parameters are different.
- For spectra: date, name, concentration, pathlength, temperature, pH
- For variable temperature experiments, like thermal denaturation (NOT spectra): date, name, concentration, pathlength, pH, temperature  range, wavelength
- For multiple spectra recorded at different temepratures as a part of a variable temperature experiment: date, name, concentration, pathlength, pH, temperature  range
- For multiple spectra recorded over time (kinetics) at a fixed temperature: date, name, concentration, pathlength, pH, temperature, time range

You need to provide the input data as a list and specify which type of the experiment you are processing by passing an argument mode like this:
exp.process_experiment(files=files, mode='single_scans')
The modes availbale are:
- 'single_scans' for a single or multiple spectra extracted from INDIVIDUAL files
- 'var_temp' for variable temperature experiments (NOT multiple spectra recorded at different temepratures)
- 'var_temp_scans' for multiple spectra recorded at different temepratures as a part of a variable temperature experiment
- 'kinet_scans' for multiple spectra recorded over time
In the examples, you may find the module datafiles_opener that is used to find the text files in the specified directory, therefore, you need to make sure you create a separate directory only for the files that will be processed. Also, the class takes only the datasets in a simple txt format, so the data needs to be exported into txt files using the instrument software.

## Extracted data

All the data extracted from a single or multiple files will be stored as a dictionary that you can access by calling
exp.get_samples())  # where exp is an example of the class that you've created
or print off by doing:
exp.print_off_dict(exp.get_samples())  # where exp is an example of the class that you've created

## Plotting

You can plot the data either manually or using one of the methods that already present in the class.
exp.spec_plotter()  # is made for plotting a single or multiple spectra
You can pass it parameters that you would normally set up for matplotlib or/and seaborn
exp.t_plotter()  # is made for plotting for thermal denaturation datasets

## Labels

By setting up labels for your datasets, you can control how they will be shown in the legend.
For example, for concentration dependencies you may find it convenient to use concentrations in the legend:
conc = [val['Conc_uM'] for val in exp.get_samples().values()]  # make a list of concenetrations by extracting the values
labels = [str(c) + r' $\mu$M' for c in conc[:]]  # make labels with uM units
exp.set_labels(labels)  # set up the labels

