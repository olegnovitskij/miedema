from Enthalpy_calculation import H_miedema, Morse_potential, L_J_potential, calculate_enthalpy
import copy
from ase.neighborlist import NeighborList, neighbor_list
from Enthalpy_calculation import convert_to_float
from pyxtal import pyxtal
from ase.spacegroup import crystal
import numpy as np
from kivy.app import App
from kivy.lang import Builder
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.uix.textinput import TextInput
import weakref
import pandas as pd
from kivy.garden.matplotlib.backend_kivyagg import FigureCanvasKivyAgg
from ase.io import write
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


class VisualizationWindow(Screen):
    need_add = 1

    def __init__(self, **kwargs):
        super(VisualizationWindow, self).__init__(**kwargs)
        self.ids.extended_grouped_dist_arr = np.array([])
        self.ids.dist_arr = np.array([])

    def on_press_visu(self, *args):
        if len(self.ids) == 0:
            Clock.schedule_once(self.plot_graph)
        else:
            self.plot_graph()

    def plot_graph(self):
        box = self.ids.fl_visu
        if self.need_add == 1:
            pic = FigureCanvasKivyAgg(plt.gcf())
            self.ids['plot_visu'] = weakref.ref(pic)
            box.add_widget(pic)
            self.need_add = 0
        elif self.need_add == 0:
            box.remove_widget(self.ids.plot_visu)
            self.need_add = 1

    def add_element_button(self):

        self.ids.n += 1
        el = BeutyTI()
        x = BeutyTI(text="-")
        y = BeutyTI(text="-")
        z = BeutyTI(text="-")
        wp = BeutyTI()
        sym = BeutyTI(text="-")
        occ = BeutyTI(text="1")

        self.ids.el_table.ids[f'el_{self.ids.n}'] = weakref.ref(el)
        self.ids.el_table.ids[f'x_{self.ids.n}'] = weakref.ref(x)
        self.ids.el_table.ids[f'y_{self.ids.n}'] = weakref.ref(y)
        self.ids.el_table.ids[f'z_{self.ids.n}'] = weakref.ref(z)
        self.ids.el_table.ids[f'wp_{self.ids.n}'] = weakref.ref(wp)
        self.ids.el_table.ids[f'sym_{self.ids.n}'] = weakref.ref(sym)
        self.ids.el_table.ids[f'occ_{self.ids.n}'] = weakref.ref(occ)

        self.ids.el_table.add_widget(el)
        self.ids.el_table.add_widget(x)
        self.ids.el_table.add_widget(y)
        self.ids.el_table.add_widget(z)
        self.ids.el_table.add_widget(wp)
        self.ids.el_table.add_widget(sym)
        self.ids.el_table.add_widget(occ)

    def remove_element_button(self):
        if self.ids.n != 0:
            for label in ['el', 'x', 'y', 'z', 'wp', 'sym', 'occ']:
                exec(f'self.ids.el_table.remove_widget(self.ids.el_table.ids.{label}_{self.ids.n})')

            self.ids.n -= 1

    def press_restart_button(self):
        self.h.text = f'Doing...'
        self.ids.extended_grouped_dist_arr = copy.deepcopy(self.manager.ids.main.ids.grouped_dist_arr)

        type_compound = int(self.manager.ids.main.ids.compound_type.text)
        num_atoms = int(self.manager.ids.main.ids.num_atoms)
        el_list = []
        for i in range(1, self.ids.n + 1):
            exec(
                f'el_list.append([self.ids.el_table.ids.el_{i}.text, self.ids.el_table.ids.wp_{i}.text, float(self.ids.el_table.ids.occ_{i}.text)])')

        for el, wp, occ in el_list:
            cond1 = ((self.ids.extended_grouped_dist_arr[:, 1] == wp) & (
                        self.ids.extended_grouped_dist_arr[:, 0] == el)).any().any()

            cond2 = ((self.ids.extended_grouped_dist_arr[:, 3] == wp) & (
                        self.ids.extended_grouped_dist_arr[:, 2] == el)).any().any()

            if cond1 or cond2:
                for i, _ in enumerate(self.ids.extended_grouped_dist_arr):
                    if (self.ids.extended_grouped_dist_arr[i, 1] == wp) & (
                            self.ids.extended_grouped_dist_arr[i, 0] == el):
                        self.ids.extended_grouped_dist_arr[i, 5] = occ

                    if (self.ids.extended_grouped_dist_arr[i, 3] == wp) & (
                            self.ids.extended_grouped_dist_arr[i, 2] == el):
                        self.ids.extended_grouped_dist_arr[i, 6] = occ

        dist_arr = np.array([[[A], A_site, [B], B_site, Distance, [occ_A], [occ_B], cnt] for
                             (A, A_site, B, B_site, Distance, occ_A, occ_B, cnt) in self.ids.extended_grouped_dist_arr])

        for el, wp, occ in el_list:
            cond1 = ((self.ids.extended_grouped_dist_arr[:, 1] == wp) & (
                        self.ids.extended_grouped_dist_arr[:, 0] == el)).any().any()

            cond2 = ((self.ids.extended_grouped_dist_arr[:, 3] == wp) & (
                        self.ids.extended_grouped_dist_arr[:, 2] == el)).any().any()
            if not (cond1 or cond2):
                for row in dist_arr[(self.ids.extended_grouped_dist_arr[:, 1] == wp)]:
                    row[0].append(el)
                    row[5].append(occ)

                for row in dist_arr[(self.ids.extended_grouped_dist_arr[:, 3] == wp)]:
                    row[2].append(el)
                    row[6].append(occ)

        print('transformed crystal: \n', dist_arr)

        calculate_enthalpy(dist_arr, n_jobs=1)

        self.h.text = f'Done!'


class BeutyTI(TextInput):
    pass


class MainWindow(Screen):
    def __init__(self, n=4, **kwargs):
        super(MainWindow, self).__init__(**kwargs)
        self.ids.n = n
        self.ids.num_atoms = 1
        self.ids.grouped_dist_arr = np.array([])

    i = 0

    def press_button(self):
        
        columns_el = ['el', 'x', 'y', 'z', 'occ']
        columns_crystal = ['a', 'b', 'c',
                           'alpha', 'beta', 'gamma',
                           'spacegroup', 'compound_type']
        el_list = []
        for i in range(1, self.ids.n + 1):
            exec(
                f'el_list.append([self.ids.el_{i}.text, self.ids.x_{i}.text, self.ids.y_{i}.text, self.ids.z_{i}.text, self.ids.occ_{i}.text])')

        print(f'el_list: \n {el_list}')
        crystal_list = [[self.ids.a.text,
                         self.ids.b.text,
                         self.ids.c.text,
                         self.ids.alpha.text,
                         self.ids.beta.text,
                         self.ids.gamma.text,
                         self.ids.spacegroup.text,
                         self.ids.compound_type.text]]

        print(f'crystal_list: \n {crystal_list}')

        self.crystal_info = pd.DataFrame(crystal_list,
                                         columns=columns_crystal).replace(r'^\s*$', np.nan, regex=True).dropna().astype(
            'float64').astype({'compound_type': 'int64', 'spacegroup': 'int64'})

        self.element_info = pd.DataFrame(el_list, columns=columns_el).replace(r'^\s*$', np.nan,
                                                                              regex=True).dropna().astype({'el': str, })

        self.element_info[['x', 'y', 'z', 'occ']] = self.element_info[['x', 'y', 'z', 'occ']].applymap(
            lambda x: convert_to_float(x))

        type_compound = int(self.ids.compound_type.text)

        self.params = {'spacegroup': int(self.crystal_info.loc[:, 'spacegroup'].values[0]),
                       'basis': [tuple(x) for x in
                                 self.element_info.loc[:, ['x', 'y', 'z']].values],

                       'symbols': self.element_info.loc[:, 'el'].values,
                       'cellpar': self.crystal_info.loc[:, ['a',
                                                            'b',
                                                            'c',
                                                            'alpha',
                                                            'beta',
                                                            'gamma']].values[0],
                       'occupancies': self.element_info.loc[:, 'occ'].values,

                       }

        compound_type = self.crystal_info.loc[:, 'compound_type'].values[0]
        print(f'self.params: \n{self.params}')
        crystal_ase = crystal(**self.params, pbc=True)
        print(f'crystal_ase: \n{crystal_ase}')
        crystal_pyxtal = pyxtal()

        crystal_pyxtal.from_seed(crystal_ase)
        print(f'crystal_pyxtal: \n{crystal_pyxtal}')
        self.ids.num_atoms = crystal_ase.get_global_number_of_atoms()

        sites = []

        for el_name in crystal_pyxtal.get_site_labels():
            for el_site in crystal_pyxtal.get_site_labels()[el_name]:
                sites.append((el_name, el_site))

        ordered_site = [sites[i] for i in crystal_ase.arrays['spacegroup_kinds']]

        if self.ids.cutoff.text == '-1':
            cutoff = ((crystal_ase.get_volume()) ** (1 / 3)) / 2
            print('CUTOFF: ', cutoff)
        else:
            cutoff = float(self.ids.cutoff.text)

        uniq_elem = list(crystal_pyxtal.get_site_labels().keys())

        num_Ions = crystal_pyxtal.numIons

        el_to_c = dict(zip(uniq_elem, num_Ions / sum(num_Ions)))

        neigh_list = neighbor_list('ijd', crystal_ase, cutoff)

        i, j, d = neigh_list

        i_decoded = [ordered_site[num] for num in i]

        j_decoded = [ordered_site[num] for num in j]

        dist = []

        for n in range(len(i_decoded)): 
            dist.append([*i_decoded[n], *j_decoded[n]] + [round(d[n], 3), 1, 1])

        try:
            self.h.text = f'Doing...'
            grouped_dist_df = \
            pd.DataFrame(dist, columns=['A', 'A_site', 'B', 'B_site', 'Distance', 'occ_A', 'occ_B']).groupby(
                ['A', 'A_site', 'B', 'B_site', 'Distance', 'occ_A', 'occ_B'])['Distance'].count().to_frame(
                'cnt').reset_index()

            print('dist', grouped_dist_df)
            grouped_dist_df['dupl_set'] = grouped_dist_df.apply(lambda x: frozenset([x['A'],
                                                                                     x['A_site'],
                                                                                     x['B'],
                                                                                     x['B_site'],
                                                                                     x['Distance']]),
                                                                axis=1)

            grouped_dist_df = grouped_dist_df.drop_duplicates(subset='dupl_set')

            grouped_dist_df = grouped_dist_df.drop(columns='dupl_set')

            self.ids.grouped_dist_arr = grouped_dist_df.values

            write('image.png', crystal_ase)
            img = mpimg.imread('image.png')
            plt.imshow(img)
            print('initial crystal: \n', grouped_dist_df)
            self.manager.ids.visu.ids.el_table.clear_widgets()
            output = []
            for num, s in enumerate(crystal_pyxtal.atom_sites):
                spl = str(s).split()
                el = BeutyTI(text=spl[0])
                x = BeutyTI(text=spl[3])
                y = BeutyTI(text=spl[4])
                z = BeutyTI(text=spl[5][:-2:])
                wp = BeutyTI(text=spl[7][1:-1:])
                sym = BeutyTI(text=spl[9][1:-1:])
                occ = BeutyTI(text="1")

                self.manager.ids.visu.ids.el_table.ids[f'el_{num + 1}'] = weakref.ref(el)
                self.manager.ids.visu.ids.el_table.ids[f'x_{num + 1}'] = weakref.ref(x)
                self.manager.ids.visu.ids.el_table.ids[f'y_{num + 1}'] = weakref.ref(y)
                self.manager.ids.visu.ids.el_table.ids[f'z_{num + 1}'] = weakref.ref(z)
                self.manager.ids.visu.ids.el_table.ids[f'wp_{num + 1}'] = weakref.ref(wp)
                self.manager.ids.visu.ids.el_table.ids[f'sym_{num + 1}'] = weakref.ref(sym)
                self.manager.ids.visu.ids.el_table.ids[f'occ_{num + 1}'] = weakref.ref(occ)

                self.manager.ids.visu.ids.el_table.add_widget(el)
                self.manager.ids.visu.ids.el_table.add_widget(x)
                self.manager.ids.visu.ids.el_table.add_widget(y)
                self.manager.ids.visu.ids.el_table.add_widget(z)
                self.manager.ids.visu.ids.el_table.add_widget(wp)
                self.manager.ids.visu.ids.el_table.add_widget(sym)
                self.manager.ids.visu.ids.el_table.add_widget(occ)

            self.manager.ids.visu.ids.n = len(crystal_pyxtal.atom_sites)

            dist_arr = np.array([[[A], A_site, [B], B_site, Distance, [occ_A], [occ_B], cnt] for
                                 (A, A_site, B, B_site, Distance, occ_A, occ_B, cnt) in self.ids.grouped_dist_arr])

            calculate_enthalpy(dist_arr, n_jobs=1)
            self.h.text = f'Done!'
        except BaseException as e:
            print(e)
            self.h.text = f'error! please, try to change cutoff or structure'

    def add_element_button(self):
        self.ids.n += 1
        el_i = BeutyTI()
        x_i = BeutyTI()
        y_i = BeutyTI()
        z_i = BeutyTI()
        occ_i = BeutyTI(text='1')

        self.ids[f'el_{self.ids.n}'] = weakref.ref(el_i)
        self.ids[f'x_{self.ids.n}'] = weakref.ref(x_i)
        self.ids[f'y_{self.ids.n}'] = weakref.ref(y_i)
        self.ids[f'z_{self.ids.n}'] = weakref.ref(z_i)
        self.ids[f'occ_{self.ids.n}'] = weakref.ref(occ_i)

        self.el_table.add_widget(el_i)
        self.el_table.add_widget(x_i)
        self.el_table.add_widget(y_i)
        self.el_table.add_widget(z_i)
        self.el_table.add_widget(occ_i)

    def remove_element_button(self):
        # add element row
        if self.ids.n != 0:
            for label in ['el', 'x', 'y', 'z', 'occ']:
                exec(f'self.el_table.remove_widget(self.ids.{label}_{self.ids.n})')

            self.ids.n -= 1


class WindowManager(ScreenManager):
    pass


kv = Builder.load_file("my.kv")


class MyMainApp(App):
    def build(self):
        return kv


if __name__ == "__main__":
    MyMainApp().run()