from matplotlib import animation
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np
import h5py

def get_smilei_number(file_name):
    file = h5py.File(file_name, 'r')
    l = list(file.keys())
    return len(l)