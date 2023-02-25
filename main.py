from Enthalpy_calculation import H_miedema

# import numpy as np, pandas as pd
from Enthalpy_calculation import H_miedema
import copy
import numpy as np, pandas as pd

from pyxtal.symmetry import Group
from pyxtal.molecule import pyxtal_molecule
from pymatgen.core import Molecule
from pyxtal import pyxtal
import ase
from ase.spacegroup import crystal
from ase.visualize import view
from ase.neighborlist import get_distance_matrix
from ase.geometry import get_distances
from ase.neighborlist import NeighborList, neighbor_list
from ase.geometry.analysis import Analysis
from pyxtal.lattice import Lattice
import ase.io as io
from ase.build import cut
from ase.spacegroup import crystal
from pyxtal.viz import display_crystals
# import pandas as pd
from ase.spacegroup import Spacegroup
from Enthalpy_calculation import H
from pyxtal.symmetry import Group
from pyxtal.molecule import pyxtal_molecule
from pymatgen.core import Molecule
from pyxtal import pyxtal
import ase
from ase.spacegroup import crystal
from ase.visualize import view
from ase.neighborlist import get_distance_matrix
from ase.geometry import get_distances
from ase.atoms import Atoms
from ase.geometry.analysis import Analysis
from pyxtal.lattice import Lattice
import ase.io as io
from ase.build import cut
from ase.spacegroup import crystal
from pyxtal.viz import display_crystals
# import pandas as pd
from ase.spacegroup import Spacegroup
import numpy as np
from kivy.app import App
from kivy.lang import Builder
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.properties import ObjectProperty, NumericProperty
from kivy.uix.popup import Popup
from kivy.uix.label import Label
from database import DataBase
from kivy.uix.textinput import TextInput
from kivy.uix.floatlayout import FloatLayout
import weakref
import pandas as pd
from kivy.garden.matplotlib.backend_kivyagg import FigureCanvasKivyAgg
import matplotlib.pyplot as plt
from ase.io import write
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# total = 0
# for A, B, d, cnt in neigh_df.values:
#     total += H(A, B, 1, 1, 1, 'miedema_coefficients.xlsx', 0, 0, 1, 1, d, cnt)
class VisualizationWindow(Screen):
    # pass

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
            
        # box = self.ids.fl_visu
        # box.add_widget(FigureCanvasKivyAgg(plt.gcf()))        
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
            #сделать один суперкласс и от него отнаследовать методы для таблицы
    def add_element_button(self):
        #add element row
        self.ids.n += 1
        print(self.ids)
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
    #add element row
        if self.ids.n != 0:
            for label in ['el', 'x', 'y', 'z', 'wp', 'sym', 'occ']:
                exec(f'self.ids.el_table.remove_widget(self.ids.el_table.ids.{label}_{self.ids.n})')

            self.ids.n -= 1
            print(self.ids.n)
            
            print(self.manager.ids.main.ids)
            
            
    def press_restart_button(self):       
        # считать строчки новые
        # добавить элементы и оккупации в новые позиции скопированного массива
        # отнормализовать оккупации
        # пересчитать энтальпию
        # надо считать со строки под номером main ids n по self ids n 
        # скопировать dist df из main window
        self.ids.extended_grouped_dist_arr = copy.deepcopy(self.manager.ids.main.ids.grouped_dist_arr)
        # self.id.dist_arr = np.array([[[A], A_site, [B], B_site, Distance, [occ_A], [occ_B], cnt] for (A, A_site, B, B_site, Distance, occ_A, occ_B, cnt) in grouped_dist_arr])
        
        print(self.ids.extended_grouped_dist_arr)
        print('visu: \n', self.ids.n)
        print('main: \n', self.manager.ids.main.ids.n)
        print('df : \n', self.ids.extended_grouped_dist_arr)
        el_list = []
        for i in range(1, self.ids.n + 1):
            exec(f'el_list.append([self.ids.el_table.ids.el_{i}.text, self.ids.el_table.ids.wp_{i}.text, float(self.ids.el_table.ids.occ_{i}.text)])')
        print(el_list)
        for el, wp, occ in el_list:
            #проверка, что occ для существуещего элемента
            cond1 = ((self.ids.extended_grouped_dist_arr[:, 1] == wp) & (self.ids.extended_grouped_dist_arr[:, 0] == el)).any().any()
            
            cond2 = ((self.ids.extended_grouped_dist_arr[:, 3] == wp) & (self.ids.extended_grouped_dist_arr[:, 2] == el)).any().any()
            
            print('cond1: ', cond1)
            
            print('cond2: ', cond2)
            
            if cond1 or cond2:
                for i,_ in enumerate(self.ids.extended_grouped_dist_arr):
                    if (self.ids.extended_grouped_dist_arr[i, 1] == wp) & (self.ids.extended_grouped_dist_arr[i, 0] == el):
                        self.ids.extended_grouped_dist_arr[i, 5] = occ

                    if (self.ids.extended_grouped_dist_arr[i, 3] == wp) & (self.ids.extended_grouped_dist_arr[i, 2] == el):
                        self.ids.extended_grouped_dist_arr[i, 6] = occ
        
        dist_arr = np.array([[[A], A_site, [B], B_site, Distance, [occ_A], [occ_B], cnt] for (A, A_site, B, B_site, Distance, occ_A, occ_B, cnt) in self.ids.extended_grouped_dist_arr])
        
        for el, wp, occ in el_list:
            #проверка, что occ для существуещего элемента
            cond1 = ((self.ids.extended_grouped_dist_arr[:, 1] == wp) & (self.ids.extended_grouped_dist_arr[:, 0] == el)).any().any()
            
            cond2 = ((self.ids.extended_grouped_dist_arr[:, 3] == wp) & (self.ids.extended_grouped_dist_arr[:, 2] == el)).any().any()
            
            if not (cond1 or cond2):
                for row in dist_arr[(self.ids.extended_grouped_dist_arr[:, 1] == wp)]:
                    row[0].append(el)
                    row[5].append(occ)
                
                for row in dist_arr[(self.ids.extended_grouped_dist_arr[:, 3] == wp)]:
                    row[2].append(el)
                    row[6].append(occ)
                
        for i, row in enumerate(dist_arr):
            norm_sum = sum(dist_arr[i][6])
            dist_arr[i][6] = [occ_/norm_sum for occ_ in dist_arr[i][6]]
            print(dist_arr[i][6])
            norm_sum = sum(dist_arr[i][5])
            dist_arr[i][5] = [occ_/norm_sum for occ_ in dist_arr[i][5]]
            
            
        
        
        print(dist_arr) 
        
        self.total = 0
        for A_ar, _, B_ar, _, dist, occ_A_ar, occ_B_ar, cnt in dist_arr:
            for i, A in enumerate(A_ar):
                for j, B in enumerate(B_ar):
                    self.total += H(A, B, 1, 1, 1, 'miedema_coefficients.xlsx', 0, 0, occ_A_ar[i], occ_B_ar[j], dist, cnt)

        print('H= ', self.total)
        self.h.text = f'H = {self.total}'
# calculate H by dist_arr

                
              
        
        
class BeutyTI(TextInput):
    pass

class MainWindow(Screen):
    def __init__(self, n=4, **kwargs):
        super(MainWindow, self).__init__(**kwargs)
        self.ids.n = n
        self.ids.grouped_dist_arr = np.array([])
    # self.manager.ids.visu.ids.n = len(crystal_pyxtal.atom_sites) 
    i = 0
    #подсвечивать и выделять пустые значения
    # объединить каждую строку в виджет который будет являться лэяутом и добавлять их динамически
    def press_button(self):
        #return enthalpy


        columns_el = ['el', 'x', 'y', 'z', 'occ']
        columns_crystal = ['a', 'b', 'c', 
                           'alpha', 'beta', 'gamma', 
                           'spacegroup', 'cutoff', 'compound_type']
        el_list = []
        for i in range(1, self.ids.n + 1):
            exec(f'el_list.append([self.ids.el_{i}.text, self.ids.x_{i}.text, self.ids.y_{i}.text, self.ids.z_{i}.text, self.ids.occ_{i}.text])')
        
        crystal_list = [[self.ids.a.text,
         self.ids.b.text,
         self.ids.c.text, 
         self.ids.alpha.text, 
         self.ids.beta.text,
         self.ids.gamma.text,
         self.ids.spacegroup.text,
         self.ids.cutoff.text,
         self.ids.compound_type.text]]

        self.crystal_info = pd.DataFrame(crystal_list, 
                                         columns=columns_crystal).replace(r'^\s*$', np.nan,                                               regex=True).dropna().astype('float64').astype({'compound_type': 'int64', 'spacegroup': 'int64'}) 
        self.element_info = pd.DataFrame(el_list, columns=columns_el).replace(r'^\s*$', np.nan,                                                    regex=True).dropna().astype({'el': str, 
                                                                              'x': 'float64',
                                                                              'y': 'float64',
                                                                              'z': 'float64',
                                                                              'occ': 'float64'}) 
        
        self.params = {'spacegroup': int(self.crystal_info.loc[:, 'spacegroup'].values[0]),
                  'basis': [tuple(x) for x in 
                              self.element_info.loc[:, ['x','y','z']].values],
                                                                                                                   
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
        
        crystal_ase = crystal(**self.params, pbc=True)
        print(crystal_ase)
        crystal_pyxtal = pyxtal()

        crystal_pyxtal.from_seed(crystal_ase)
        
        sites = []
        
        for el_name in crystal_pyxtal.get_site_labels():
            for el_site in crystal_pyxtal.get_site_labels()[el_name]:
                sites.append((el_name, el_site))
                
        ordered_site = [sites[i] for i in crystal_ase.arrays['spacegroup_kinds']]
        
        cutoff = ((crystal_ase.get_volume())**(1/3))/2
        
        neigh_list = neighbor_list('ijd', crystal_ase, cutoff)
        
        i, j, d = neigh_list
        
        i_decoded = [ordered_site[num] for num in i]
        
        j_decoded = [ordered_site[num] for num in j]
        
        dist = []
        
        for n in range(len(i_decoded)): #сюда добавить окупации
            dist.append([*i_decoded[n], *j_decoded[n]] + [round(d[n], 3), 1, 1])
        
        
        
        grouped_dist_df = pd.DataFrame(dist, columns=['A', 'A_site', 'B', 'B_site', 'Distance', 'occ_A', 'occ_B']).groupby(['A', 'A_site', 'B', 'B_site', 'Distance', 'occ_A', 'occ_B'])['Distance'].count().to_frame('cnt').reset_index()
        
        self.ids.grouped_dist_arr = grouped_dist_df.values
        

        
        write('image.png', crystal_ase)
        img = mpimg.imread('image.png')
        plt.imshow(img)
        print(self.ids.grouped_dist_arr)
        # Что делать с атомами в одинаковых позициях но где != 0 по миедеме
        # я хочу составить табличку попарных расстояний между парами частица-пст = частица-пст
        # после этого я смогу добавить occ, и множить базовые строки, добавляя эффективно атом
        #здесь я уже могу загнать энтальпию в формулу
        # добавить эту инфо на вторую стр и дать возможность добавлять в позиции атомы. добавляя или уменьшая dist df   
        #тут надо пересоздавать grid layout
        self.manager.ids.visu.ids.el_table.clear_widgets()
        output = []
        #здесь просто вызвать метод класс окна визуализаций
        for num, s in enumerate(crystal_pyxtal.atom_sites):
            spl = str(s).split()
            el = BeutyTI(text=spl[0])
            x = BeutyTI(text=spl[3])
            y = BeutyTI(text=spl[4])
            z = BeutyTI(text=spl[5][:-2:])
            wp = BeutyTI(text=spl[7][1:-1:])
            sym = BeutyTI(text=spl[9][1:-1:])
            occ = BeutyTI(text="1")
            
            self.manager.ids.visu.ids.el_table.ids[f'el_{num+1}'] = weakref.ref(el)
            self.manager.ids.visu.ids.el_table.ids[f'x_{num+1}'] = weakref.ref(x)
            self.manager.ids.visu.ids.el_table.ids[f'y_{num+1}'] = weakref.ref(y)
            self.manager.ids.visu.ids.el_table.ids[f'z_{num+1}'] = weakref.ref(z)
            self.manager.ids.visu.ids.el_table.ids[f'wp_{num+1}'] = weakref.ref(wp)
            self.manager.ids.visu.ids.el_table.ids[f'sym_{num+1}'] = weakref.ref(sym)
            self.manager.ids.visu.ids.el_table.ids[f'occ_{num+1}'] = weakref.ref(occ)
            
            self.manager.ids.visu.ids.el_table.add_widget(el)
            self.manager.ids.visu.ids.el_table.add_widget(x)
            self.manager.ids.visu.ids.el_table.add_widget(y)
            self.manager.ids.visu.ids.el_table.add_widget(z)
            self.manager.ids.visu.ids.el_table.add_widget(wp)
            self.manager.ids.visu.ids.el_table.add_widget(sym)
            self.manager.ids.visu.ids.el_table.add_widget(occ)
            
        self.manager.ids.visu.ids.n = len(crystal_pyxtal.atom_sites)    
       
        print(self.element_info, self.crystal_info)
        
        dist_arr = np.array([[[A], A_site, [B], B_site, Distance, [occ_A], [occ_B], cnt] for (A, A_site, B, B_site, Distance, occ_A, occ_B, cnt) in self.ids.grouped_dist_arr])
        
        self.total = 0
        for A_ar, _, B_ar, _, dist, occ_A_ar, occ_B_ar, cnt in dist_arr:
            for i, A in enumerate(A_ar):
                for j, B in enumerate(B_ar):
                    self.total += H(A, B, 1, 1, 1, 'miedema_coefficients.xlsx', 0, 0, occ_A_ar[i], occ_B_ar[j], dist, cnt)

        print('H= ', self.total)
        self.h.text = f'H = {self.total}'
            
    def add_element_button(self):
        #add element row
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
    #add element row
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