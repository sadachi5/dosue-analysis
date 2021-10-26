#!/bin/env python
import argparse;

import numpy as np;
np.set_printoptions(threshold=10);
from matplotlib import pyplot as plt;

from read_data import read_data

g_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:olive','tab:cyan','tab:gray','red','royalblue','turquoise','darkolivegr    een', 'magenta', 'blue', 'green']*5;


class Yfactor:
    def __init__(self, filepath_77K, filepath_300K, filepath_plate):
        self.data_77K   = read_data(filepath_77K)
        self.data_300K  = read_data(filepath_300K)
        self.data_plate = read_data(filepath_plate)
        pass

    def 
    def calculate(self):



