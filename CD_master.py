#-------------------------------------------------------------------------------
# Name:        CDtools
# Purpose:
#
# Author:      Andrey Romanyuk
#
# Created:     02/11/2021
# Copyright:   (c) Andrey Romanyuk 2021
# Licence:     <your licence>
#-------------------------------------------------------------------------------
import numpy as np
import re
import matplotlib.pyplot as plt
import seaborn as sns
import colorcet as cc
import itertools


class cdExp(object):
    """

    """

    def __init__(self):
        # self.mode =
        self.samples = {}  # let's store all the data for the samples of interest as a dict
        self.blanks = {}  # let's store blanks separately
        self.plot_order = []
        self.title = ''
        self.mode = 'single_scans'

    def get_samples(self):
        return self.samples

    def sample_param(self, sample, par_key, par_val):
        self.get_samples()[sample][par_key] = par_val

    def all_sample_param(self, par_key, par_vals):
        for s, v in zip(self.get_samples().values(), par_vals):
            s[par_key] = v

    def get_blanks(self):
        return self.blanks

    def get_mode(self):
        return self.mode

    def set_mode(self, mode):
        self.mode = mode

    def get_title(self):
        return self.title

    def set_title(self, title):
        self.title = title

    def set_peptide_bonds(self, peptide_bonds=[]):
        # Depending on the mode, the number of peptide bonds can either be entered manually or
        # calculated for the mixtures based on the ratios of the components
        if not peptide_bonds:
            peptide_bonds = [int(i) for i in input('Specify the number of peptide bonds for each sample:\t').split(',')]
        if len(peptide_bonds)  == 1 and len(peptide_bonds) < len(self.get_samples()):
            peptide_bonds = peptide_bonds * len(self.get_samples())
        for v, p in zip(self.get_samples().values(), peptide_bonds):
            v['Pep_bonds'] = p

    def get_labels(self):
        for val in self.get_samples().values():
            if 'Label' not in val.keys():
                self.set_labels(
                    [f"{val['Title']}_{val['Conc_uM']} uM_{val['Pathlength_mm']} mm_{val['Temp_C']}C" for val in
                     self.get_samples().values()])
        labs = [val['Label'] for val in self.get_samples().values()]
        return labs

    def set_labels(self, labels):
        # if self.get_labels():
        for val, lab in zip(self.get_samples().values(), labels):
            val['Label'] = lab

    def get_plot_order(self):
        return self.plot_order

    def set_plot_order(self, _list=[]):
        if not _list:
            _list = list(self.get_samples().keys())
        self.plot_order = _list

    def set_precise_conc(self, conc, sample_number):
        self.samples[sample_number]['Conc_uM'] = conc
        self.calculate_MRE()

    def get_n_samples(self):
        return len(self.get_samples())

    def calc_pep_bonds(self, fractions=[], individual_pb=[]):
        assert len(individual_pb) == len(fractions), "The number of components and provided fractions don't match up"
        fractions = np.asarray(fractions)
        individual_pb = np.asarray(individual_pb)
        prod = [pb * fr for pb, fr in zip(individual_pb, fractions)]
        peptide_bonds = np.sum(prod, axis=0)
        self.set_peptide_bonds(peptide_bonds)

    def calc_helical_fraction(self, MRE222, temp, n_pep):
        return 100*(MRE222 - 640+45*temp) / (-42500 * (1 - 3/n_pep) - 640+45*temp)

    def calc_all_helical_fractions(self):
        for val in self.get_samples().values():
            val['Hel_fr'] = self.calc_helical_fraction(MRE222=val['MRE'][np.where(val['WL'] == 222)][0], temp=val['Temp_C'], n_pep=val['Pep_bonds'])

    def print_off_dict(self, _dict):
        for k, v in _dict.items():

            if isinstance(v, dict):
                print(k)
                self.print_off_dict(v)
            else:
                print(f'{k}\t{v}')

    def input_datasets(self, files, bl_cor=1):
        count_blanks = 1
        count_samples = 1
        for f in files:
            # if bl_cor==1:
            if 'blank'.casefold() in f.lower():
                self.blanks[count_blanks] = {'Name': f}
                count_blanks += 1
            else:
                self.samples[count_samples] = {'Name': f}
                count_samples += 1

        print('\nSamples:\n')
        for k, v in self.get_samples().items():
            print(k, v['Name'])
        if bl_cor == 1:
            print('\nBlanks:\n')
            for k, v in self.get_blanks().items():
                print(k, v['Name'])
            if len(self.get_blanks()) == 1:
                blanks_order = [1] * len(self.get_samples())
            else:
                blanks_order = [int(i) for i in input(
                    'Please, list the numbers of blanks in the order corresponding to the samples they should be used with:\t').split(',')]

            for v, b in zip(self.get_samples().values(), blanks_order):
                v['Blank'] = b

    def process_filenames(self, mode='single_scans'):
        if mode:
            self.set_mode(mode)
        for val in self.get_samples().values():
            params = {'Info': []}
            filename = val['Name'].split('_')

            for el in filename:
                if filename.index(el) == 0:
                    params['Date'] = el
                if filename.index(el) == 1:
                    params['Title'] = el
                if filename.index(el) == 2:
                    if el.endswith('mM'):
                        conc_uM = float(el[:-2]) * 1000
                    elif el.endswith('uM'):
                        conc_uM = float(el[:-2])
                    elif el.endswith('M'):
                        conc_uM = float(el[:-2]) * 1e6
                    params['Conc_uM'] = conc_uM
                elif el.endswith('mm'):
                    if re.match('0+\d+', el):
                        pathlength_mm = float(el[:-2]) / 10 ** (len(el[:-2]) - 1)
                    else:
                        pathlength_mm = float(el[:-2])
                    params['Pathlength_mm'] = pathlength_mm
                elif el.endswith('C'):
                    temp_C = np.array([float(te) for te in re.split('-', el[:-1])])
                    if self.get_mode() == 'single_scans' or self.get_mode() == 'kinet_scans':
                        params['Temp_C'] = temp_C[0]
                    elif self.get_mode() == 'var_temp' or self.get_mode() == 'var_temp_scans':
                        params['Temp_C_range'] = temp_C
                elif el.endswith('nm') and (self.get_mode() == 'var_temp' or  self.get_mode() == 'var_temp_scans') :
                    params['WL'] = float(el[:-2])
                elif (el.endswith('min') or el.endswith('s') or el.endswith('h') ) and\
                        (self.get_mode() == 'kin' or self.get_mode() == 'kinet_scans'):
                    if el.endswith('min'):
                        params['Time_min'] = float(el[:-3])
                    elif el.endswith('s'):
                        params['Time_min'] = float(el[:-1])/60
                    elif el.endswith('h'):
                        params['Time_min'] = float(el[:-1]) / 3600

                elif el.startswith('pH'):
                    pH = float(el[2:])
                    if pH > 14:
                        pH /= 10
                    params['pH'] = pH
                else:
                    params['Info'].append(el)

            params['Info'] = '_'.join(params['Info'])

            if self.get_mode() == 'single_scans':
                par_checklist = ['Date', 'Title', 'Conc_uM', 'Pathlength_mm', 'Temp_C', 'pH']
            elif self.get_mode() == 'var_temp':
                par_checklist = ['Date', 'Title', 'Conc_uM', 'Pathlength_mm', 'Temp_C_range', 'pH', 'WL']
            elif self.get_mode() == 'var_temp_scans':
                par_checklist = ['Date', 'Title', 'Conc_uM', 'Pathlength_mm', 'Temp_C_range', 'pH']
            elif self.get_mode() == 'kinet_scans':
                par_checklist = ['Date', 'Title', 'Conc_uM', 'Pathlength_mm', 'Temp_C', 'pH', 'Time_min']
            for par in par_checklist:
                assert par in params.keys(), f'Error: some parameters are missing in the file name:\t {par}'
            for k, v in params.items():
                val[k] = v

    def read_jasco_datafiles(self, mode='single_scans'):
        if mode:
            self.set_mode(mode)
        read_set = { 'single_scans':
                        {'main_var': '',
                         'samp_col_name': ['WL', "mdeg", "HT"],
                         'blank_col_name': ['WL', "mdeg", "HT"],
                         'n_head_lines': 19,
                         'n_foot_lines': 0 },
                     'var_temp':
                         {'main_var': '',
                          'samp_col_name': ['t', "mdeg", "HT"],
                          'blank_col_name': ['WL', "mdeg", "HT"],
                          'n_head_lines': 19,
                          'n_foot_lines': 0},
                     'var_temp_scans':
                         {'main_var': 'Time_s',
                          'samp_col_name': ['t', "mdeg", "HT"],
                          'blank_col_name': ['WL', "mdeg", "HT"],
                          'n_head_lines': 20,
                          'n_foot_lines': 0},
                     'kinet_scans':
                         {'main_var': 'Temp_C',
                          'samp_col_name': ['t', "mdeg", "HT"],
                          'blank_col_name': ['WL', "mdeg", "HT"],
                          'n_head_lines': 20,
                          'n_foot_lines': 0}
                     }
        blank_n_head_lines = 19

        for kk, val in self.get_samples().items():
            with open(val['Name'] + '.txt', 'r') as f:
                val['Metadata'] = f.readlines()[:20]
            for line in val['Metadata']:
                line = line.split('\t')
                if line[0] == 'DELTAX':
                    val['Step_x'] = abs(float(line[1]))
                elif line[0] == 'FIRSTX':
                    val['x0'] = float(line[1])
                elif line[0] == 'LASTX':
                    val['xe'] = float(line[1])
                elif line[0] == 'NPOINTS':
                    val['n_points'] = int(line[1])
                elif line[0] == 'NPOINTST':  # t is either temperature or time
                    val['tn'] = int(line[1])
                elif line[0] == 'TUNITS':  # t is either temperature or time
                    val['t_units'] = (line[1])
                    # if 'sec' in line[1]:
                    #     self.set_mode('kinet')
                elif line[0] == 'FIRSTT':
                    val['t0'] = int(line[1])
                elif line[0] == 'LASTT':
                    val['te'] = int(line[1])

            if self.get_mode() == 'var_temp_scans' or self.get_mode() == 'kinet_scans':
                val['t'] = (np.linspace(val['t0'], val['te'], val['tn'], endpoint=1, dtype=int))
                if 'sec' in line[1]:
                    print('True')
                    val['t'] = val['t'] / 60
                read_set[self.get_mode()]['n_foot_lines'] = val['n_points'] + 1
                read_set[self.get_mode()]['samp_col_name'] = ["WL"] + [str(t) for t in val['t']]
                sub_samples = {}
                for key in read_set[self.get_mode()]['samp_col_name'][1:]:
                    sub_samples[key] = {}

            # generate numpy arrays from data
            array = np.genfromtxt(val['Name'] + '.txt', skip_header=read_set[self.get_mode()]['n_head_lines'],
                                  skip_footer=read_set[self.get_mode()]['n_foot_lines'],
                                        names=read_set[self.get_mode()]['samp_col_name'], delimiter=''
                                  )

            if self.get_mode() == 'var_temp_scans' or self.get_mode() == 'kinet_scans':
                array_HT = np.genfromtxt(val['Name'] + '.txt',
                                       skip_header=read_set[self.get_mode()]['n_head_lines']+val['n_points']+1,
                                       names=read_set[self.get_mode()]['samp_col_name'], delimiter='')

            for count, key in enumerate(read_set[self.get_mode()]['samp_col_name'][:]):
                # sub_samples[key] = {}
                if self.get_mode() == 'single_scans' or self.get_mode() == 'var_temp':

                    val[key] = array[key]
                elif self.get_mode() == 'var_temp_scans' or self.get_mode() == 'kinet_scans':
                    if count != 0:
                        sub_samples[key]['WL'] = array['WL']
                        sub_samples[key]['mdeg'] = array[key]
                        sub_samples[key]['HT'] = array_HT[key]
                        for k in val:
                            sub_samples[key][k] = val[k]

        if self.get_mode() == 'var_temp_scans' or self.get_mode() == 'kinet_scans':
            delete_keys = list(self.get_samples().keys())[:]
            for key, val in sub_samples.items():
                self.samples[key] = val
            if kk in delete_keys:
                self.samples.pop(kk)

        for val in self.get_blanks().values():
            # generate numpy arrays from data
            array = np.genfromtxt(val['Name'] + '.txt', skip_header=blank_n_head_lines,
                                        names=read_set[self.get_mode()]['blank_col_name'], delimiter='')
            # val['Raw_data'] = array
            for key in read_set[self.get_mode()]['blank_col_name']:
                val[key] = array[key]
                # self.print_off_dict(self.get_blanks())

    def calculate_MRE(self):
        for val in self.get_samples().values():
            Divide_by = val['Pathlength_mm'] * val['Conc_uM'] * val['Pep_bonds'] * 10 ** -6
            if self.get_mode() == 'single_scans' or self.get_mode() == 'var_temp_scans' or self.get_mode() == 'kinet_scans':
                # subtract blank array from data array
                scan_blanked = np.subtract(val["mdeg"], self.get_blanks()[val['Blank']]["mdeg"])
                # convert mdeg to MRE
                MRE = scan_blanked[:] / Divide_by
                val['MRE'] = MRE

            if self.get_mode() == 'var_temp':
                condition = self.get_blanks()[val['Blank']]["WL"] == val['WL']
                blank_mdeg_at_WL = np.extract(condition, self.get_blanks()[val['Blank']]["mdeg"])
                melt_blanked = np.subtract(val["mdeg"], blank_mdeg_at_WL)
                MRE_melt = melt_blanked / Divide_by
                val['MRE'] = MRE_melt

    def process_experiment(self, files, mode):
        self.input_datasets(files)
        self.set_mode(mode=mode)
        self.process_filenames(mode=mode)
        self.set_peptide_bonds()
        self.read_jasco_datafiles(mode=mode)
        self.calculate_MRE()
        self.set_plot_order()

    # TODO: add graying out parts of the graph where HT exceeds 700

    def spec_plotter(self, fig_title='', legend_pos='best', xlim=(190, 260, 10), every_nth = 1,
                         ylim=((-50, 90, 20), (150, 1050, 200)), palette="Set2", title_size=16, leg_size=16, axlab_size=16, ticklab_size=14, lw=1,
                     legend_bool=True, figsize=(9,9), dpi=300, spines=2):
        if fig_title:
            self.set_title(fig_title)
        self.get_labels()  # makes sure Labels are in place
        # plot the data
        sns.set_theme()
        sns.set_style("ticks")
        # my_palette = sns.color_palette(cc.glasbey, n_colors=24)
        my_palette = sns.color_palette(palette)

        # ax1 = plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=3)
        fig, axs = plt.subplots(2,1,
                                figsize=figsize,
                                dpi=dpi,
                                sharex=True,
                                gridspec_kw=dict(height_ratios=[2, 0.5]))

        palette = itertools.cycle(sns.color_palette(my_palette))
        coef = 1e-3  #
        
        for count, sample_n in enumerate(self.get_plot_order()):
            if count % every_nth == 0:
                data = self.get_samples()[sample_n]

                color=next(palette)
                # if self.get_mode() == 'var_temp_scans' or self.get_mode() == 'kinet_scans':

                sns.lineplot(x=data['WL'], y=data['MRE'] * coef, linewidth=lw,
                             legend=legend_bool, ax=axs[0], label=data['Label'], color=color)
                sns.lineplot(x=data['WL'], y=data['HT'], linewidth=lw, label=data['Label'],
                             legend=False, ax=axs[1], color=color)
        zero_line_x = np.linspace(min(data['WL']), max(data['WL']), len(data['WL']))
        sns.lineplot(x=zero_line_x, y=[0] * len(zero_line_x), color='k', legend=False, ax=axs[0], linewidth=lw)
        axs[0].lines[-1].set_linestyle("--")
        axs[0].set_ylabel('$\mathrm{MRE*10^{3}}$\n  $\mathrm{(deg\ cm^{2}\ dmol^{-1}\ res^{-1}}$)', size=axlab_size)
        # axs[0].legend(loc=legend_pos, prop={'size': leg_size})
        if fig_title:
            axs[0].set_title(self.get_title(), size=title_size)
        axs[1].set_ylabel('HT (V)', size=axlab_size)
        axs[1].set_xlabel('Wavelength ($\mathrm{nm}$)\n\n', size=axlab_size)

        for i, ax in enumerate(axs):
            ax.set_ylim(ylim[i][:2])
            ax.set_xlim(xlim[:2])
            major_yticks = np.arange(ylim[i][0], ylim[i][1]+0.001, ylim[i][2])
            major_xticks = np.arange(xlim[0], xlim[1]+0.001, xlim[2])
            ax.set_xticks(major_xticks)
            ax.set_yticks(major_yticks)
            ax.tick_params(axis='x', labelsize=ticklab_size)
            ax.tick_params(axis='y', labelsize=ticklab_size)
        # change all spines
        for axis in ['top', 'bottom', 'left', 'right']:
            for ax in [axs[0], axs[1]]:
                ax.spines[axis].set_linewidth(spines)

        box = axs[0].get_position()
        axs[0].set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        box = axs[1].get_position()
        axs[1].set_position([box.x0, box.y0 + box.height * 0.5, box.width, box.height * 0.9])

        # fig.legend(bbox_to_anchor=(0.05,0.0), loc="lower left", fancybox=True, shadow=True, prop={'size': 14},
        #         bbox_transform=fig.transFigure, ncol=1)
        plt.show()

    def t_plotter(self, labels=[], fig_title='', legend_pos='best', xlim=(5, 95, 10),
                     ylim=((-50, 90, 20), (150, 1050, 200)), palette="Set2", title_size=16, leg_size=16, axlab_size=16,
                     ticklab_size=14, marker_size=20):
        if fig_title:
            self.set_title(fig_title)
        if labels:
            self.set_labels(labels)
        # plot the data
        sns.set_theme()
        sns.set_style("ticks")
        # my_palette = sns.color_palette(cc.glasbey, n_colors=24)
        my_palette = sns.color_palette(palette)

        fig, axs = plt.subplots(2, 1,
                                figsize=(9, 9),
                                sharex=True,
                                gridspec_kw=dict(height_ratios=[2, 0.5]))

        palette = itertools.cycle(sns.color_palette(my_palette))
        coef = 1e-3  #
        for sample_n in self.get_plot_order():
            data = self.get_samples()[sample_n]
            color = next(palette)
            sns.scatterplot(x=data['t'], y=data['MRE'] * coef, color=color, s=marker_size,
                         legend=True, ax=axs[0], label=data['Label'])
            sns.lineplot(x=data['t'], y=data['HT'], color=color, linewidth=3, label=data['Label'],
                         legend=False, ax=axs[1])
        axs[0].set_ylabel('$\mathrm{MRE*10^{3}}$\n  $\mathrm{(deg\ cm^{2}\ dmol^{-1}\ res^{-1}}$)', size=axlab_size)
        axs[0].legend(loc=legend_pos, prop={'size': leg_size})
        axs[0].set_title(self.get_title(), size=title_size)
        axs[1].set_ylabel('HT (V)', size=axlab_size)
        axs[1].set_xlabel('Temperature\n', size=axlab_size)

        for i, ax in enumerate(axs):
            ax.set_ylim(ylim[i][:2])
            ax.set_xlim(xlim[:2])
            major_yticks = np.arange(ylim[i][0], ylim[i][1] + 0.001, ylim[i][2])
            major_xticks = np.arange(xlim[0], xlim[1] + 0.001, xlim[2])
            ax.set_xticks(major_xticks)
            ax.set_yticks(major_yticks)
            ax.tick_params(axis='x', labelsize=ticklab_size)
            ax.tick_params(axis='y', labelsize=ticklab_size)

        box = axs[0].get_position()
        axs[0].set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        box = axs[1].get_position()
        axs[1].set_position([box.x0, box.y0 + box.height * 0.5, box.width, box.height * 0.9])
        plt.show()

    # def flexi_plotter(self, X, Y, fig_title='', legend_pos='best', xlim=(190, 260, 10), every_nth=1,
    #                  ylim=((-50, 90, 20), (150, 1050, 200)), palette="Set2", title_size=16, leg_size=16, axlab_size=16,
    #                  ticklab_size=14):

    def export(self, datasets, filename='output.txt'):
        with open(filename, 'w+') as f:
            lines = ''

            for i in range(len(datasets[0])):
                for ds in datasets:
                    lines += str(ds[i]) + '\t'
                lines += '\n'
            f.write(lines)

    def export_all_MRE(self):
        for count, val in self.get_samples().items():
            self.export([val['WL'], val['MRE'], val['HT']], filename=f"{val['Name']}_MRE-{count}.txt")

    def process_experiment(self, files, mode):
        self.input_datasets(files)
        self.set_mode(mode=mode)
        self.process_filenames(mode=mode)
        self.set_peptide_bonds()
        self.read_jasco_datafiles(mode=mode)
        self.calculate_MRE()
        self.set_plot_order()

    def read_MRE_datafiles(self, files, mode='single_scans'):
        self.input_datasets(files, bl_cor=0)
        self.set_mode(mode=mode)
        self.process_filenames(mode=mode)

        for kk, val in self.get_samples().items():
            # with open(val['Name'] + '.txt', 'r') as f:
            array = np.genfromtxt(val['Name'] + '.txt', names=['WL', 'MRE', 'HT'], delimiter='')
            for k in ['WL', 'MRE', 'HT']:
                val[k] = array[k]


        self.set_plot_order()


def main():
    pass

if __name__ == '__main__':
    main()
