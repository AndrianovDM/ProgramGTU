import numpy as np
import scipy as sp
import itertools
from math import *
from fuel import *
from scipy.interpolate import interp2d
from scipy.optimize import fsolve
from scipy.interpolate import make_interp_spline
from scipy.optimize import fsolve

#                                                MUSATKIN

                                                # C_a = 0.5
alpha_ro_01_05 = [55.0, 56.25, 57.188, 58.281, 59.375, 60.625, 62.344, 63.75, 65.0, 66.719, 68.75, 70.781, 72.188, 73.594, 75.0, 76.406, 78.125,
79.219, 80.769, 82.0, 83.538, 85.077, 86.769, 88.308, 90.781, 93.125, 94.688, 96.406, 98.281, 99.844, 101.385, 102.615, 104.308, 105.692, 107.385, 108.769, 109.846]
Y_ro_01_05 = [0.429, 0.433, 0.434, 0.437, 0.439, 0.442, 0.446, 0.448, 0.451, 0.454, 0.46, 0.465, 0.467, 0.471, 0.473, 0.477, 0.481, 0.484, 0.487, 0.491, 0.495, 0.499, 0.503, 0.506, 
0.513, 0.519, 0.523, 0.526, 0.531, 0.536, 0.541, 0.545, 0.549, 0.554, 0.559, 0.564, 0.566] 
ro_01_05 = [0.1 for i in range(len(alpha_ro_01_05))]

alpha_ro_02_05 = [55.0, 57.031, 58.281, 59.844, 61.563, 62.969, 64.688, 66.875, 68.906, 71.094, 72.813, 75.0, 76.875, 78.438, 80.462, 82.154, 83.846, 85.385, 86.615, 88.615, 
90.781, 92.188, 93.594, 95.781, 97.813, 98.906, 98.906, 100.923, 102.615, 104.154, 105.692, 107.385, 109.231, 109.846]
Y_ro_02_05 = [0.451, 0.456, 0.458, 0.461, 0.465, 0.468, 0.472, 0.479, 0.482, 0.489, 0.494, 0.499, 0.503, 0.508, 0.513, 0.519, 0.523, 0.526, 0.53, 0.536, 0.544, 0.548, 
0.553, 0.559, 0.565, 0.569, 0.569, 0.576, 0.58, 0.586, 0.591, 0.599, 0.606, 0.609] 
ro_02_05 = [0.2 for i in range(len(alpha_ro_02_05))]

alpha_ro_03_05 = [55.156, 56.719, 58.281, 59.688, 61.406, 62.969, 64.688, 67.031, 69.531, 72.188, 74.688, 77.188, 79.531, 82.154, 84.154, 86.308, 88.462, 90.781, 92.813, 94.375, 97.188, 
98.906, 100.923, 102.769, 104.769, 106.923, 108.462, 109.538]  
Y_ro_03_05 = [0.475, 0.477, 0.482, 0.485, 0.489, 0.492, 0.496, 0.504, 0.51, 0.519, 0.526, 0.534, 0.541, 0.55, 0.558, 0.564, 0.571, 0.579,
 0.586, 0.591, 0.601, 0.609, 0.616, 0.624, 0.631, 0.641, 0.646, 0.654]
ro_03_05 = [0.3 for i in range(len(alpha_ro_03_05))]

alpha_ro_04_05 = [55.0, 57.656, 60.625, 63.125, 65.469, 68.125, 69.688, 72.344, 74.531, 76.719, 78.906, 81.538, 
84.462, 86.923, 89.692, 92.188, 95.156, 97.5, 99.844, 101.692, 103.692, 104.154]
Y_ro_04_05 = [0.497, 0.505, 0.513, 0.52, 0.529, 0.538, 0.543, 0.553, 0.564, 0.571, 
0.58, 0.591, 0.604, 0.615, 0.629, 0.639, 0.653, 0.665, 0.676, 0.686, 0.695, 0.699]
ro_04_05 = [0.4 for i in range(len(alpha_ro_04_05))]

alpha_ro_05_05 = [55.156, 57.188, 59.688, 61.875, 63.906, 65.938, 67.813, 69.531, 71.25, 72.969, 
75.469, 77.344, 79.375, 81.077, 83.077, 84.615, 86.615, 88.308, 89.538, 91.094, 92.188]
Y_ro_05_05 = [0.526, 0.534, 0.541, 0.549, 0.559, 0.566, 0.574, 0.58, 0.589, 0.596, 0.61,
 0.619, 0.629, 0.638, 0.646, 0.656, 0.669, 0.676, 0.684, 0.693, 0.699]
ro_05_05 = [0.5 for i in range(len(alpha_ro_05_05))]

alpha_ro_06_05 = [55.0, 57.5, 59.375, 61.406, 63.906, 65.469, 67.656, 
69.375, 72.031, 74.063, 75.938, 77.969, 79.688, 81.385]
Y_ro_06_05 = [0.566, 0.576, 0.584, 0.591, 0.603, 0.611, 0.623, 0.63, 0.645, 0.655, 0.668, 0.679, 0.688, 0.699]
ro_06_05 = [0.6 for i in range(len(alpha_ro_06_05))]

alpha_ro_07_05 = [55.0, 56.406, 58.438, 59.375, 61.406, 62.969, 64.375, 65.938, 67.813, 69.531, 70.0]
Y_ro_07_05 = [0.601, 0.608, 0.618, 0.623, 0.635, 0.645, 0.655, 0.668, 0.681, 0.695, 0.699]
ro_07_05 = [0.7 for i in range(len(alpha_ro_07_05))]

data_ro_05 = np.array([[55.0, 0.429, 0.1], [56.25, 0.433, 0.1], [57.188, 0.434, 0.1], [58.281, 0.437, 0.1], [59.375, 0.439, 0.1], [60.625, 0.442, 0.1], [62.344, 0.446, 0.1], [63.75, 0.448, 0.1], [65.0, 0.451, 0.1], [66.719, 0.454, 0.1], [68.75, 0.46, 0.1], [70.781, 0.465, 0.1], [72.188, 0.467, 0.1], [73.594, 0.471, 0.1], [75.0, 0.473, 0.1], [76.406, 0.477, 0.1], [78.125, 0.481, 0.1], [79.219, 0.484, 0.1], [80.769, 0.487, 0.1], [82.0, 0.491, 0.1], [83.538, 0.495, 0.1], [85.077, 0.499, 0.1], [86.769, 0.503, 0.1], [88.308, 0.506, 0.1], [90.781, 0.513, 0.1], [93.125, 0.519, 0.1], [94.688, 0.523, 0.1], [96.406, 0.526, 0.1], [98.281, 0.531, 0.1], [99.844, 0.536, 0.1], [101.385, 0.541, 0.1], [102.615, 0.545, 0.1], [104.308, 0.549, 0.1], [105.692, 0.554, 0.1], [107.385, 0.559, 0.1], [108.769, 0.564, 0.1], [109.846, 0.566, 0.1], [55.0, 0.451, 0.2], [57.031, 0.456, 0.2], [58.281, 0.458, 0.2], [59.844, 0.461, 0.2], [61.563, 0.465, 0.2], [62.969, 0.468, 0.2], [64.688, 0.472, 0.2], [66.875, 0.479, 0.2], [68.906, 0.482, 0.2], [71.094, 0.489, 0.2], [72.813, 0.494, 0.2], [75.0, 0.499, 0.2], [76.875, 0.503, 0.2], [78.438, 0.508, 0.2], [80.462, 0.513, 0.2], [82.154, 0.519, 0.2], [83.846, 0.523, 0.2], [85.385, 0.526, 0.2], [86.615, 0.53, 0.2], [88.615, 0.536, 0.2], [90.781, 0.544, 0.2], [92.188, 0.548, 0.2], [93.594, 0.553, 0.2], [95.781, 0.559, 0.2], [97.813, 0.565, 0.2], [98.906, 0.569, 0.2], [98.906, 0.569, 0.2], [100.923, 0.576, 0.2], [102.615, 0.58, 0.2], [104.154, 0.586, 0.2], [105.692, 0.591, 0.2], [107.385, 0.599, 0.2], [109.231, 0.606, 0.2], [109.846, 0.609, 0.2], [55.156, 0.475, 0.3], [56.719, 0.477, 0.3], [58.281, 0.482, 0.3], [59.688, 0.485, 0.3], [61.406, 0.489, 0.3], [62.969, 0.492, 0.3], [64.688, 0.496, 0.3], [67.031, 0.504, 0.3], [69.531, 0.51, 0.3], [72.188, 0.519, 0.3], [74.688, 0.526, 0.3], [77.188, 0.534, 0.3], [79.531, 0.541, 0.3], [82.154, 0.55, 0.3], [84.154, 0.558, 0.3], [86.308, 0.564, 0.3], [88.462, 0.571, 0.3], [90.781, 0.579, 0.3], [92.813, 0.586, 0.3], [94.375, 0.591, 0.3], [97.188, 0.601, 0.3], [98.906, 0.609, 0.3], [100.923, 0.616, 0.3], [102.769, 0.624, 0.3], [104.769, 0.631, 0.3], [106.923, 0.641, 0.3], [108.462, 0.646, 0.3], [109.538, 0.654, 0.3], [55.0, 0.497, 0.4], [57.656, 0.505, 0.4], [60.625, 0.513, 0.4], [63.125, 0.52, 0.4], [65.469, 0.529, 0.4], [68.125, 0.538, 0.4], [69.688, 0.543, 0.4], [72.344, 0.553, 0.4], [74.531, 0.564, 0.4], [76.719, 0.571, 0.4], [78.906, 0.58, 0.4], [81.538, 0.591, 0.4], [84.462, 0.604, 0.4], [86.923, 0.615, 0.4], [89.692, 0.629, 0.4], [92.188, 0.639, 0.4], [95.156, 0.653, 0.4], [97.5, 0.665, 0.4], [99.844, 0.676, 0.4], [101.692, 0.686, 0.4], [103.692, 0.695, 0.4], [104.154, 0.699, 0.4], [55.156, 0.526, 0.5], [57.188, 0.534, 0.5], [59.688, 0.541, 0.5], [61.875, 0.549, 0.5], [63.906, 0.559, 0.5], [65.938, 0.566, 0.5], [67.813, 0.574, 0.5], [69.531, 0.58, 0.5], [71.25, 0.589, 0.5], [72.969, 0.596, 0.5], [75.469, 0.61, 0.5], [77.344, 0.619, 0.5], [79.375, 0.629, 0.5], [81.077, 0.638, 0.5], [83.077, 0.646, 0.5], [84.615, 0.656, 0.5], [86.615, 0.669, 0.5], [88.308, 0.676, 0.5], [89.538, 0.684, 0.5], [91.094, 0.693, 0.5], [92.188, 0.699, 0.5], [55.0, 0.566, 0.6], [57.5, 0.576, 0.6], [59.375, 0.584, 0.6], [61.406, 0.591, 0.6], [63.906, 0.603, 0.6], [65.469, 0.611, 0.6], [67.656, 0.623, 0.6], [69.375, 0.63, 0.6], [72.031, 0.645, 0.6], [74.063, 0.655, 0.6], [75.938, 0.668, 0.6], [77.969, 0.679, 0.6], [79.688, 0.688, 0.6], [81.385, 0.699, 0.6], [55.0, 0.601, 0.7], [56.406, 0.608, 0.7], [58.438, 0.618, 0.7], [59.375, 0.623, 0.7], [61.406, 0.635, 0.7], [62.969, 0.645, 0.7], [64.375, 0.655, 0.7], [65.938, 0.668, 0.7], [67.813, 0.681, 0.7], [69.531, 0.695, 0.7], [70.0, 0.699, 0.7]])


                                            # C_a = 0.8
alpha_ro_01_08 = [55.077, 56.308, 58, 59.538, 61.563, 63.281, 64.688, 66.406, 67.969, 69.688, 71.719, 
72.969, 74.844, 76.875, 78.438, 80.462, 82, 83.538, 85.077, 86.308, 87.692, 89.077, 
90.625, 92.031, 93.438, 95.156, 96.719, 98.281, 100.462, 102, 103.231, 104.769, 106.308, 107.692, 109.231, 110]
Y_ro_01_08 = [0.403,0.405,0.409,0.411,0.416,0.421,0.425,0.43,0.434,0.439,0.445,0.45,0.456,0.463,0.469,0.476,0.481,0.488
,0.493,0.498,0.506,0.509,0.515,0.522,0.528,0.534,0.542,0.549,0.56,0.566,0.573,0.58,0.589,0.595,0.604,0.608]
ro_01_08 = [0.1 for i in range(len(alpha_ro_01_08))]

alpha_ro_02_08 = [55.077,56.769,58.769,60.625,62.344,64.063,65.781,67.5,69.063,70.938,72.656,74.375,76.094,77.656,79.219,80.769,
82.154,83.538,84.923,86.308,87.538,88.923,90.781,91.875,93.125,94.844,96.563,98.438,99.688,100.923,102.769,104.462,
106,107.231,108.462,109.538]
Y_ro_02_08 = [0.42,0.425,0.429,0.434,0.44,0.445,0.451,0.456,0.463,0.47,0.476,0.484,0.49,0.498,0.504,0.511,0.518,0.524,
0.53,0.537,0.543,0.551,0.558,0.565,0.571,0.579,0.587,0.596,0.603,0.61,0.62,0.628,0.636,0.643,0.651,0.656]
ro_02_08 = [0.2 for i in range(len(alpha_ro_02_08))]

alpha_ro_03_08 = [55.077, 56.769, 58.462, 59.692, 61.094, 62.344, 63.594, 65.0, 66.719, 67.969, 69.531, 71.406, 73.281, 75.156, 77.344, 79.219, 80.923, 82.615, 
84.615, 86.154, 87.692, 89.077, 91.563, 93.438, 95.0, 96.875, 98.438, 100.308, 102.0, 103.385, 104.615, 105.846]
Y_ro_03_08 = [0.439, 0.441, 0.446, 0.449, 0.454, 0.459, 0.464, 0.469, 0.476, 0.48, 0.488, 0.496, 0.506, 0.515, 0.525, 0.534, 0.544, 
0.554, 0.565, 0.572, 0.58, 0.589, 0.605, 0.615, 0.625, 0.636, 0.646, 0.659, 0.671, 0.68, 0.689, 0.699]
ro_03_08 = [0.3 for i in range(len(alpha_ro_03_08))]

alpha_ro_04_08 = [55.077, 56.308, 57.692, 59.231, 60.938, 62.5, 63.594, 65.0, 66.875, 68.281, 69.375, 70.469, 
71.563, 72.656, 74.063, 75.156, 76.25, 77.656, 78.75, 80.308, 81.846, 82.923, 84.308, 85.385, 86.462, 87.538, 
88.923, 89.846, 90.625, 91.719, 93.281, 94.375, 95.469, 96.563, 97.656]
Y_ro_04_08 = [0.456, 0.46, 0.465, 0.469, 0.476, 0.483, 0.488, 0.494, 0.503, 0.509, 0.515, 0.522, 0.528, 
0.534, 0.542, 0.548, 0.554, 0.562, 0.57, 0.579, 0.589, 0.595, 0.605, 0.611, 0.619, 0.626, 0.635, 
0.641, 0.648, 0.655, 0.666, 0.674, 0.681, 0.691, 0.699]
ro_04_08 = [0.4 for i in range(len(alpha_ro_04_08))]

alpha_ro_05_08 = [54.923, 56.154, 57.538, 58.923, 60.313, 61.719, 62.969, 64.063, 65.156, 66.25, 67.5, 69.063, 70.781, 72.031, 
73.281, 74.688, 75.781, 76.875, 78.594, 79.688, 80.923, 81.846, 82.923, 84.308, 85.692, 87.077, 88.462, 89.538, 90.781, 91.563]
Y_ro_05_08 = [0.476, 0.48, 0.485, 0.49, 0.496, 0.503, 0.509, 0.515, 0.522, 0.527, 0.534, 0.544, 0.554, 0.562, 0.57, 0.579, 0.585,
 0.592, 0.604, 0.611, 0.62, 0.626, 0.635, 0.644, 0.654, 0.665, 0.675, 0.681, 0.691, 0.696]
ro_05_08 = [0.5 for i in range(len(alpha_ro_05_08))]

alpha_ro_06_08 = [55.231, 56.769, 58.308, 59.538, 61.406, 63.125, 64.531, 65.781, 67.188, 68.906, 70.938, 72.969, 74.688, 
76.094, 77.969, 79.375, 80.923, 82.462, 83.538, 83.846]
Y_ro_06_08 = [0.499, 0.505, 0.51, 0.516, 0.528, 0.538, 0.547, 0.554, 0.566, 0.576, 0.591, 0.608, 0.62, 
0.633, 0.646, 0.659, 0.674, 0.688, 0.696, 0.699]
ro_06_08 = [0.6 for i in range(len(alpha_ro_06_08))]

alpha_ro_07_08 = [55.077, 57.231, 58.154, 59.077, 60.156, 60.938, 62.344, 63.281, 64.219, 65.313, 66.25, 67.344, 68.125, 
69.375, 70.469, 71.25, 72.031, 72.813, 73.594, 74.375, 75.156, 75.781, 76.406]
Y_ro_07_08 = [0.534, 0.543, 0.547, 0.551, 0.557, 0.561, 0.568, 0.575, 0.581, 0.589, 
0.596, 0.605, 0.611, 0.623, 0.631, 0.64, 0.648, 0.656, 0.665, 0.675, 0.685, 0.693, 0.699]
ro_07_08 = [0.7 for i in range(len(alpha_ro_07_08))]

data_ro_08 = np.array([[55.077, 0.403, 0.1], [56.308, 0.405, 0.1], [58, 0.409, 0.1], [59.538, 0.411, 0.1], [61.563, 0.416, 0.1], [63.281, 0.421, 0.1], [64.688, 0.425, 0.1], [66.406, 0.43, 0.1], [67.969, 0.434, 0.1], [69.688, 0.439, 0.1], [71.719, 0.445, 0.1], [72.969, 0.45, 0.1], [74.844, 0.456, 0.1], [76.875, 0.463, 0.1], [78.438, 0.469, 0.1], [80.462, 0.476, 0.1], [82, 0.481, 0.1], [83.538, 0.488, 0.1], [85.077, 0.493, 0.1], [86.308, 0.498, 0.1], [87.692, 0.506, 0.1], [89.077, 0.509, 0.1], [90.625, 0.515, 0.1], [92.031, 0.522, 0.1], [93.438, 0.528, 0.1], [95.156, 0.534, 0.1], [96.719, 0.542, 0.1], [98.281, 0.549, 0.1], [100.462, 0.56, 0.1], [102, 0.566, 0.1], [103.231, 0.573, 0.1], [104.769, 0.58, 0.1], [106.308, 0.589, 0.1], [107.692, 0.595, 0.1], [109.231, 0.604, 0.1], [110, 0.608, 0.1], [55.077, 0.42, 0.2], [56.769, 0.425, 0.2], [58.769, 0.429, 0.2], [60.625, 0.434, 0.2], [62.344, 0.44, 0.2], [64.063, 0.445, 0.2], [65.781, 0.451, 0.2], [67.5, 0.456, 0.2], [69.063, 0.463, 0.2], [70.938, 0.47, 0.2], [72.656, 0.476, 0.2], [74.375, 0.484, 0.2], [76.094, 0.49, 0.2], [77.656, 0.498, 0.2], [79.219, 0.504, 0.2], [80.769, 0.511, 0.2], [82.154, 0.518, 0.2], [83.538, 0.524, 0.2], [84.923, 0.53, 0.2], [86.308, 0.537, 0.2], [87.538, 0.543, 0.2], [88.923, 0.551, 0.2], [90.781, 0.558, 0.2], [91.875, 0.565, 0.2], [93.125, 0.571, 0.2], [94.844, 0.579, 0.2], [96.563, 0.587, 0.2], [98.438, 0.596, 0.2], [99.688, 0.603, 0.2], [100.923, 0.61, 0.2], [102.769, 0.62, 0.2], [104.462, 0.628, 0.2], [106, 0.636, 0.2], [107.231, 0.643, 0.2], [108.462, 0.651, 0.2], [109.538, 0.656, 0.2], [55.077, 0.439, 0.3], [56.769, 0.441, 0.3], [58.462, 0.446, 0.3], [59.692, 0.449, 0.3], [61.094, 0.454, 0.3], [62.344, 0.459, 0.3], [63.594, 0.464, 0.3], [65.0, 0.469, 0.3], [66.719, 0.476, 0.3], [67.969, 0.48, 0.3], [69.531, 0.488, 0.3], [71.406, 0.496, 0.3], [73.281, 0.506, 0.3], [75.156, 0.515, 0.3], [77.344, 0.525, 0.3], [79.219, 0.534, 0.3], [80.923, 0.544, 0.3], [82.615, 0.554, 0.3], [84.615, 0.565, 0.3], [86.154, 0.572, 0.3], [87.692, 0.58, 0.3], [89.077, 0.589, 0.3], [91.563, 0.605, 0.3], [93.438, 0.615, 0.3], [95.0, 0.625, 0.3], [96.875, 0.636, 0.3], [98.438, 0.646, 0.3], [100.308, 0.659, 0.3], [102.0, 0.671, 0.3], [103.385, 0.68, 0.3], [104.615, 0.689, 0.3], [105.846, 0.699, 0.3], [55.077, 0.456, 0.4], [56.308, 0.46, 0.4], [57.692, 0.465, 0.4], [59.231, 0.469, 0.4], [60.938, 0.476, 0.4], [62.5, 0.483, 0.4], [63.594, 0.488, 0.4], [65.0, 0.494, 0.4], [66.875, 0.503, 0.4], [68.281, 0.509, 0.4], [69.375, 0.515, 0.4], [70.469, 0.522, 0.4], [71.563, 0.528, 0.4], [72.656, 0.534, 0.4], [74.063, 0.542, 0.4], [75.156, 0.548, 0.4], [76.25, 0.554, 0.4], [77.656, 0.562, 0.4], [78.75, 0.57, 0.4], [80.308, 0.579, 0.4], [81.846, 0.589, 0.4], [82.923, 0.595, 0.4], [84.308, 0.605, 0.4], [85.385, 0.611, 0.4], [86.462, 0.619, 0.4], [87.538, 0.626, 0.4], [88.923, 0.635, 0.4], [89.846, 0.641, 0.4], [90.625, 0.648, 0.4], [91.719, 0.655, 0.4], [93.281, 0.666, 0.4], [94.375, 0.674, 0.4], [95.469, 0.681, 0.4], [96.563, 0.691, 0.4], [97.656, 0.699, 0.4], [54.923, 0.476, 0.5], [56.154, 0.48, 0.5], [57.538, 0.485, 0.5], [58.923, 0.49, 0.5], [60.313, 0.496, 0.5], [61.719, 0.503, 0.5], [62.969, 0.509, 0.5], [64.063, 0.515, 0.5], [65.156, 0.522, 0.5], [66.25, 0.527, 0.5], [67.5, 0.534, 0.5], [69.063, 0.544, 0.5], [70.781, 0.554, 0.5], [72.031, 0.562, 0.5], [73.281, 0.57, 0.5], [74.688, 0.579, 0.5], [75.781, 0.585, 0.5], [76.875, 0.592, 0.5], [78.594, 0.604, 0.5], [79.688, 0.611, 0.5], [80.923, 0.62, 0.5], [81.846, 0.626, 0.5], [82.923, 0.635, 0.5], [84.308, 0.644, 0.5], [85.692, 0.654, 0.5], [87.077, 0.665, 0.5], [88.462, 0.675, 0.5], [89.538, 0.681, 0.5], [90.781, 0.691, 0.5], [91.563, 0.696, 0.5], [55.231, 0.499, 0.6], [56.769, 0.505, 0.6], [58.308, 0.51, 0.6], [59.538, 0.516, 0.6], [61.406, 0.528, 0.6], [63.125, 0.538, 0.6], [64.531, 0.547, 0.6], [65.781, 0.554, 0.6], [67.188, 0.566, 0.6], [68.906, 0.576, 0.6], [70.938, 0.591, 0.6], [72.969, 0.608, 0.6], [74.688, 0.62, 0.6], [76.094, 0.633, 0.6], [77.969, 0.646, 0.6], [79.375, 0.659, 0.6], [80.923, 0.674, 0.6], [82.462, 0.688, 0.6], [83.538, 0.696, 0.6], [83.846, 0.699, 0.6], [55.077, 0.534, 0.7], [57.231, 0.543, 0.7], [58.154, 0.547, 0.7], [59.077, 0.551, 0.7], [60.156, 0.557, 0.7], [60.938, 0.561, 0.7], [62.344, 0.568, 0.7], [63.281, 0.575, 0.7], [64.219, 0.581, 0.7], [65.313, 0.589, 0.7], [66.25, 0.596, 0.7], [67.344, 0.605, 0.7], [68.125, 0.611, 0.7], [69.375, 0.623, 0.7], [70.469, 0.631, 0.7], [71.25, 0.64, 0.7], [72.031, 0.648, 0.7], [72.813, 0.656, 0.7], [73.594, 0.665, 0.7], [74.375, 0.675, 0.7], [75.156, 0.685, 0.7], [75.781, 0.693, 0.7], [76.406, 0.699, 0.7]])


                                            # C_a = 1.1
alpha_ro_01_11 = [54.844, 57.188, 59.375, 61.875, 64.219, 66.875, 69.063, 71.875, 74.531, 76.875, 78.906, 81.692, 83.692, 85.538, 87.846, 89.692, 91.563, 93.594, 95.156, 
96.719, 98.438, 100.769, 102.308, 103.538, 105.538, 107.231, 109.231]
Y_ro_01_11 = [0.374, 0.378, 0.382, 0.39, 0.397, 0.406, 0.414, 0.425, 0.438, 0.448, 0.457,
0.471, 0.482, 0.492, 0.505, 0.515, 0.525, 0.535, 0.545, 0.555, 0.565, 0.581, 0.59, 0.598, 0.611, 0.623, 0.639]
ro_01_11 = [0.1 for i in range(len(alpha_ro_01_11))]

alpha_ro_02_11 = [55.0, 56.719, 58.906, 60.625, 62.813, 65.156, 67.656, 69.375, 72.031, 74.063, 76.094, 78.125, 79.688, 81.692, 83.538, 
85.231, 87.231, 89.077, 90.938, 92.656, 94.219, 96.094, 97.813, 98.906, 100.769, 102.0, 103.538, 104.923, 106.154, 106.615]
Y_ro_02_11 = [0.388, 0.391, 0.397, 0.404, 0.411, 0.422, 0.432, 0.441, 0.454, 0.466, 0.477, 0.49, 0.499, 0.513, 0.523,
 0.534, 0.549, 0.561, 0.575, 0.588, 0.598, 0.611, 0.625, 0.635, 0.649, 0.658, 0.672, 0.685, 0.694, 0.697]
ro_02_11 = [0.2 for i in range(len(alpha_ro_02_11))]

alpha_ro_03_11 = [55.313, 57.031, 58.906, 61.094, 63.125, 65.469, 67.5, 69.688, 71.231, 73.385, 75.077, 
77.385, 79.538, 81.563, 84.219, 86.719, 88.438, 90.781, 92.813, 94.688, 96.563, 97.969, 99.688, 101.818]
Y_ro_03_11 = [0.401, 0.406, 0.413, 0.42, 0.433, 0.444, 0.453, 0.465, 0.475, 0.487, 
0.497, 0.511, 0.525, 0.539, 0.556, 0.574, 0.588, 0.603, 0.62, 0.634, 0.648, 0.662, 0.676, 0.696]
ro_03_11 = [0.3 for i in range(len(alpha_ro_03_11))]

alpha_ro_04_11 = [55.0, 56.25, 58.281, 60.781, 62.969, 64.844, 66.094, 67.656, 69.219, 70.769, 72.308, 73.538, 75.077,
 76.923, 78.308, 79.692, 81.25, 82.5, 84.219, 85.938, 87.031, 88.594, 89.844, 91.406, 92.969, 94.531, 95.625, 96.875]
Y_ro_04_11 = [0.418, 0.422, 0.429, 0.441, 0.451, 0.462, 0.468, 0.477, 0.487, 0.497, 0.508, 
0.515, 0.526, 0.538, 0.549, 0.558, 0.57, 0.58, 0.591, 0.604, 0.614, 0.627, 0.637, 0.649, 0.662, 0.676, 0.685, 0.697]
ro_04_11 = [0.4 for i in range(len(alpha_ro_04_11))]

alpha_ro_05_11 = [55.156, 57.031, 58.75, 60.781, 62.656, 64.375, 66.25, 67.969, 69.531, 
71.692, 73.538, 75.692, 77.385, 79.077, 80.313, 81.719, 83.594, 85.469, 87.031, 88.75, 90.469, 91.25]
Y_ro_05_11 = [0.434, 0.441, 0.448, 0.459, 0.471, 0.48, 0.492, 0.505, 0.515, 
0.528, 0.543, 0.559, 0.573, 0.586, 0.596, 0.608, 0.624, 0.641, 0.657, 0.673, 0.689, 0.699]
ro_05_11 = [0.5 for i in range(len(alpha_ro_05_11))]

alpha_ro_06_11 = [55.156, 56.875, 59.063, 61.406, 63.125, 65.313, 66.875, 68.594, 70.462, 
71.692, 72.769, 74.308, 75.231, 76.923, 78.0, 79.385, 80.938, 81.719, 82.656, 83.906, 85.0]
Y_ro_06_11 = [0.452, 0.461, 0.471, 0.486, 0.496, 0.511, 0.521, 0.534, 0.549, 0.56, 0.569, 
0.584, 0.593, 0.609, 0.619, 0.634, 0.651, 0.661, 0.671, 0.685, 0.696]
ro_06_11 = [0.6 for i in range(len(alpha_ro_06_11))]

alpha_ro_07_11 = [55.156, 57.188, 59.063, 61.094, 62.656, 64.531, 66.094, 67.813, 69.375, 71.231, 72.923, 74.462, 76.0, 77.385, 78.769, 79.692, 80.469]
Y_ro_07_11 = [0.476, 0.484, 0.495, 0.508, 0.519, 0.533, 0.545, 0.56, 0.574, 0.59, 0.606, 0.623, 0.641, 0.657, 0.672, 0.685, 0.696]
ro_07_11 = [0.7 for i in range(len(alpha_ro_07_11))]

data_ro_11 = np.array([[54.844, 0.374, 0.1], [57.188, 0.378, 0.1], [59.375, 0.382, 0.1], [61.875, 0.39, 0.1], [64.219, 0.397, 0.1], [66.875, 0.406, 0.1], [69.063, 0.414, 0.1], [71.875, 0.425, 0.1], [74.531, 0.438, 0.1], [76.875, 0.448, 0.1], [78.906, 0.457, 0.1], [81.692, 0.471, 0.1], [83.692, 0.482, 0.1], [85.538, 0.492, 0.1], [87.846, 0.505, 0.1], [89.692, 0.515, 0.1], [91.563, 0.525, 0.1], [93.594, 0.535, 0.1], [95.156, 0.545, 0.1], [96.719, 0.555, 0.1], [98.438, 0.565, 0.1], [100.769, 0.581, 0.1], [102.308, 0.59, 0.1], [103.538, 0.598, 0.1], [105.538, 0.611, 0.1], [107.231, 0.623, 0.1], [109.231, 0.639, 0.1], [55.0, 0.388, 0.2], [56.719, 0.391, 0.2], [58.906, 0.397, 0.2], [60.625, 0.404, 0.2], [62.813, 0.411, 0.2], [65.156, 0.422, 0.2], [67.656, 0.432, 0.2], [69.375, 0.441, 0.2], [72.031, 0.454, 0.2], [74.063, 0.466, 0.2], [76.094, 0.477, 0.2], [78.125, 0.49, 0.2], [79.688, 0.499, 0.2], [81.692, 0.513, 0.2], [83.538, 0.523, 0.2], [85.231, 0.534, 0.2], [87.231, 0.549, 0.2], [89.077, 0.561, 0.2], [90.938, 0.575, 0.2], [92.656, 0.588, 0.2], [94.219, 0.598, 0.2], [96.094, 0.611, 0.2], [97.813, 0.625, 0.2], [98.906, 0.635, 0.2], [100.769, 0.649, 0.2], [102.0, 0.658, 0.2], [103.538, 0.672, 0.2], [104.923, 0.685, 0.2], [106.154, 0.694, 0.2], [106.615, 0.697, 0.2], [55.313, 0.401, 0.3], [57.031, 0.406, 0.3], [58.906, 0.413, 0.3], [61.094, 0.42, 0.3], [63.125, 0.433, 0.3], [65.469, 0.444, 0.3], [67.5, 0.453, 0.3], [69.688, 0.465, 0.3], [71.231, 0.475, 0.3], [73.385, 0.487, 0.3], [75.077, 0.497, 0.3], [77.385, 0.511, 0.3], [79.538, 0.525, 0.3], [81.563, 0.539, 0.3], [84.219, 0.556, 0.3], [86.719, 0.574, 0.3], [88.438, 0.588, 0.3], [90.781, 0.603, 0.3], [92.813, 0.62, 0.3], [94.688, 0.634, 0.3], [96.563, 0.648, 0.3], [97.969, 0.662, 0.3], [99.688, 0.676, 0.3], [101.818, 0.696, 0.3], [55.0, 0.418, 0.4], [56.25, 0.422, 0.4], [58.281, 0.429, 0.4], [60.781, 0.441, 0.4], [62.969, 0.451, 0.4], [64.844, 0.462, 0.4], [66.094, 0.468, 0.4], [67.656, 0.477, 0.4], [69.219, 0.487, 0.4], [70.769, 0.497, 0.4], [72.308, 0.508, 0.4], [73.538, 0.515, 0.4], [75.077, 0.526, 0.4], [76.923, 0.538, 0.4], [78.308, 0.549, 0.4], [79.692, 0.558, 0.4], [81.25, 0.57, 0.4], [82.5, 0.58, 0.4], [84.219, 0.591, 0.4], [85.938, 0.604, 0.4], [87.031, 0.614, 0.4], [88.594, 0.627, 0.4], [89.844, 0.637, 0.4], [91.406, 0.649, 0.4], [92.969, 0.662, 0.4], [94.531, 0.676, 0.4], [95.625, 0.685, 0.4], [96.875, 0.697, 0.4], [55.156, 0.434, 0.5], [57.031, 0.441, 0.5], [58.75, 0.448, 0.5], [60.781, 0.459, 0.5], [62.656, 0.471, 0.5], [64.375, 0.48, 0.5], [66.25, 0.492, 0.5], [67.969, 0.505, 0.5], [69.531, 0.515, 0.5], [71.692, 0.528, 0.5], [73.538, 0.543, 0.5], [75.692, 0.559, 0.5], [77.385, 0.573, 0.5], [79.077, 0.586, 0.5], [80.313, 0.596, 0.5], [81.719, 0.608, 0.5], [83.594, 0.624, 0.5], [85.469, 0.641, 0.5], [87.031, 0.657, 0.5], [88.75, 0.673, 0.5], [90.469, 0.689, 0.5], [91.25, 0.699, 0.5], [55.156, 0.452, 0.6], [56.875, 0.461, 0.6], [59.063, 0.471, 0.6], [61.406, 0.486, 0.6], [63.125, 0.496, 0.6], [65.313, 0.511, 0.6], [66.875, 0.521, 0.6], [68.594, 0.534, 0.6], [70.462, 0.549, 0.6], [71.692, 0.56, 0.6], [72.769, 0.569, 0.6], [74.308, 0.584, 0.6], [75.231, 0.593, 0.6], [76.923, 0.609, 0.6], [78.0, 0.619, 0.6], [79.385, 0.634, 0.6], [80.938, 0.651, 0.6], [81.719, 0.661, 0.6], [82.656, 0.671, 0.6], [83.906, 0.685, 0.6], [85.0, 0.696, 0.6], [55.156, 0.476, 0.7], [57.188, 0.484, 0.7], [59.063, 0.495, 0.7], [61.094, 0.508, 0.7], [62.656, 0.519, 0.7], [64.531, 0.533, 0.7], [66.094, 0.545, 0.7], [67.813, 0.56, 0.7], [69.375, 0.574, 0.7], [71.231, 0.59, 0.7], [72.923, 0.606, 0.7], [74.462, 0.623, 0.7], [76.0, 0.641, 0.7], [77.385, 0.657, 0.7], [78.769, 0.672, 0.7], [79.692, 0.685, 0.7], [80.469, 0.696, 0.7]])

def Y_inter(function, parametr, x):
    parametr, residuals, rank, sv, rcond = np.polyfit(parametr, function, 5, full=True)
    function = sp.poly1d(parametr)
    return function(x)

def Ro_inter(function, parametr, x):
    parametr, residuals, rank, sv, rcond = np.polyfit(parametr, function, 4, full=True)
    function = sp.poly1d(parametr)
    return function(x)

def musatkin_Ro(c_2a, Y, alpha2):
    
    if 0 < c_2a <= 0.65: #При C_a=0.5
        r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        Y01 = Y_inter(Y_ro_01_05, alpha_ro_01_05, alpha2) 
        Y02 = Y_inter(Y_ro_02_05, alpha_ro_02_05, alpha2) 
        Y03 = Y_inter(Y_ro_03_05, alpha_ro_03_05, alpha2) 
        Y04 = Y_inter(Y_ro_04_05, alpha_ro_04_05, alpha2) 
        Y05 = Y_inter(Y_ro_05_05, alpha_ro_05_05, alpha2) 
        Y06 = Y_inter(Y_ro_06_05, alpha_ro_06_05, alpha2) 
        Y07 = Y_inter(Y_ro_07_05, alpha_ro_07_05, alpha2) 
        Y_data = [Y01, Y02, Y03, Y04, Y05, Y06, Y07]
        ro = Ro_inter(r, Y_data, Y)

    elif 0.65 < c_2a <= 0.95:
        r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        Y01 = Y_inter(Y_ro_01_08, alpha_ro_01_08, alpha2) 
        Y02 = Y_inter(Y_ro_02_08, alpha_ro_02_08, alpha2) 
        Y03 = Y_inter(Y_ro_03_08, alpha_ro_03_08, alpha2) 
        Y04 = Y_inter(Y_ro_04_08, alpha_ro_04_08, alpha2) 
        Y05 = Y_inter(Y_ro_05_08, alpha_ro_05_08, alpha2) 
        Y06 = Y_inter(Y_ro_06_08, alpha_ro_06_08, alpha2) 
        Y07 = Y_inter(Y_ro_07_08, alpha_ro_07_08, alpha2) 
        Y_data = [Y01, Y02, Y03, Y04, Y05, Y06, Y07]
        ro = Ro_inter(r, Y_data, Y)
    
    else:
        r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
        Y01 = Y_inter(Y_ro_01_11, alpha_ro_01_11, alpha2) 
        Y02 = Y_inter(Y_ro_02_11, alpha_ro_02_11, alpha2) 
        Y03 = Y_inter(Y_ro_03_11, alpha_ro_03_11, alpha2) 
        Y04 = Y_inter(Y_ro_04_11, alpha_ro_04_11, alpha2) 
        Y05 = Y_inter(Y_ro_05_11, alpha_ro_05_11, alpha2) 
        Y06 = Y_inter(Y_ro_06_11, alpha_ro_06_11, alpha2) 
        Y07 = Y_inter(Y_ro_07_11, alpha_ro_07_11, alpha2) 
        Y_data = [Y01, Y02, Y03, Y04, Y05, Y06, Y07]
        ro = Ro_inter(r, Y_data, Y)
    return ro

def musatkin_Y(c_2a, ro, alpha2):
    
    if 0 < c_2a <= 0.65: #При C_a=0.5
        y = [60,70,80,90,100,110] #Величина alfa_2
        x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] #ρ_cр
        z = [[0.44,0.46,0.485,0.512,0.545,0.588,0.628],
              [0.462,0.485,0.512,0.545,0.585,0.636,0.7],   #Y_ст
              [0.485,0.512,0.545,0.588,0.635,0.692,0.772],
              [0.512,0.542,0.578,0.6,0.6783,0.7427,0.844],
              [0.538,0.575,0.615,0.678,0.723,0.7947,0.916],
              [0.57,0.608,0.656,0.7157,0.7683,0.847,0.988]]
        f = interp2d(x, y, z, kind='cubic')
        z = f(ro, alpha2)

    elif 0.65 < c_2a <= 0.95: #При C_a=0.8
        y = [60,70,80,90,100,110] #Величина alfa_2
        x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] #ρ_cр
        z = [[0.411,0.43,0.45,0.472,0.494,0.518,0.554],
              [0.44,0.465,0.488,0.52,0.548,0.584,0.63],   #Y_ст
              [0.474,0.51,0.54,0.578,0.615,0.665,0.706],
              [0.514,0.555,0.595,0.645,0.7,0.736,0.782],
              [0.558,0.606,0.658,0.698,0.7605,0.8095,0.858],
              [0.612,0.658,0.7031,0.7557,0.829,0.883,0.934]]
        f = interp2d(x, y, z, kind='cubic')
        z = f(ro, alpha2)

    else:    #При C_a=1.1
        y0 = [60,70,80,90,100] #Величина alfa_2
        x0 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7] #ρ_cр
        z0 = [[0.368,0.4,0.42,0.438,0.455,0.478,0.5],
            [0.418,0.444,0.466,0.49,0.52,0.55,0.582],   #Y_ст
            [0.464,0.5,0.532,0.564,0.598,0.644,0.69],
            [0.516,0.568,0.6,0.64,0.684,0.723,0.7807],
            [0.576,0.645,0.68,0.703,0.7555,0.8063,0.8757]]
        f = interp2d(x0, y0, z0, kind='cubic')
        z = f(ro, alpha2)

    return z[0]

#                                                SMITH
mu_85 =  [3.122, 3.093, 3.064, 3.029, 2.986, 2.957, 2.913, 2.862, 2.812, 2.768, 2.71, 2.645, 2.602, 2.544, 2.485, 2.441]
c_85 =  [1.005, 1.029, 1.057, 1.083, 1.117, 1.134, 1.157, 1.183, 1.2, 1.217, 1.232, 1.247, 1.255, 1.264, 1.27, 1.274]
etta_85 = [0.85 for i in range(len(mu_85))]

mu_86 =  [3.0, 2.986, 2.964, 2.942, 2.92, 2.891, 2.87, 2.834, 2.797, 2.747, 2.681, 2.623, 2.565, 2.478, 2.419, 2.36, 2.309, 2.243]
c_86 =  [0.939, 0.965, 0.99, 1.012, 1.038, 1.061, 1.083, 1.104, 1.123, 1.143, 1.166, 1.185, 1.2, 1.215, 1.226, 1.234, 1.24, 1.247]
etta_86 = [0.86 for i in range(len(mu_86))]

mu_87 =  [2.819, 2.812, 2.797, 2.783, 2.761, 2.732, 2.689, 2.631, 2.558, 2.485, 2.405, 1.399, 2.324, 1.457, 1.522, 2.243, 1.623, 2.155, 2.044, 1.718, 1.928, 1.826]
c_87 =  [0.859, 0.885, 0.917, 0.945, 0.984, 1.016, 1.053, 1.083, 1.117, 1.14, 1.162, 1.174, 1.179, 1.181, 1.187, 1.189, 1.196, 1.198, 1.202, 1.202, 1.204, 1.204]
etta_87 = [0.87 for i in range(len(mu_87))]

mu_88 =  [2.616, 2.616, 2.623, 2.631, 2.623, 2.616, 2.602, 2.587, 2.558, 2.522, 2.478, 2.427, 2.36, 2.309, 2.235, 1.283, 1.348, 2.147, 1.42, 2.044, 1.486, 1.957, 1.877, 1.594, 1.79, 1.696]
c_88 =  [0.765, 0.782, 0.803, 0.831, 0.861, 0.885, 0.919, 0.954, 0.982, 1.005, 1.031, 1.055, 1.078, 1.098, 1.111, 1.115, 1.121, 1.123, 1.13, 1.134, 1.136, 1.14, 1.143, 1.145, 1.147, 1.147]
etta_88 = [0.88 for i in range(len(mu_88))]

mu_89 =  [2.36, 2.383, 2.405, 2.419, 2.427, 2.441, 2.441, 2.427, 2.397, 2.375, 2.331, 2.287, 2.228, 2.155, 2.088, 1.189, 2.037, 1.261, 1.964, 1.326, 1.899, 1.812, 1.435, 1.725, 1.645, 1.551]
c_89 =  [0.672, 0.696, 0.722, 0.747, 0.767, 0.801, 0.835, 0.878, 0.911, 0.941, 0.967, 0.992, 1.016, 1.04, 1.053, 1.059, 1.066, 1.07, 1.074, 1.078, 1.081, 1.085, 1.087, 1.089, 1.091, 1.091]
etta_89 = [0.89 for i in range(len(mu_89))]

mu_90 =  [1.986, 2.037, 2.081, 2.14, 2.191, 2.235, 2.258, 2.272, 2.272, 2.265, 2.243, 2.206, 2.155, 2.096, 2.037, 1.942, 1.065, 1.102, 1.834, 1.167, 1.261, 1.703, 1.609, 1.493, 1.391]
c_90 =  [0.543, 0.567, 0.588, 0.623, 0.659, 0.694, 0.732, 0.773, 0.799, 0.831, 0.863, 0.889, 0.919, 0.945, 0.965, 0.99, 0.99, 0.997, 1.005, 1.005, 1.018, 1.02, 1.029, 1.029, 1.029]
etta_90 = [0.90 for i in range(len(mu_90))]

mu_91 =  [1.602, 1.638, 1.689, 1.754, 1.812, 1.87, 1.92, 1.971, 2.008, 2.044, 2.066, 2.088, 2.088, 2.074, 2.044, 2.008, 1.942, 1.877, 1.0, 1.051, 1.797, 1.145, 1.232, 1.681, 1.319, 1.587, 1.406, 1.5]
c_91 =  [0.46, 0.472, 0.489, 0.515, 0.539, 0.562, 0.586, 0.614, 0.638, 0.672, 0.7, 0.734, 0.775, 0.816, 0.842, 0.868, 0.894, 0.915, 0.919, 0.926, 0.937, 0.939, 0.945, 0.952, 0.954, 0.958, 0.958, 0.962]
etta_91 = [0.91 for i in range(len(mu_91))]

mu_92 =  [1.319, 1.384, 1.435, 1.471, 1.522, 1.573, 1.623, 1.674, 1.725, 1.776, 1.826, 1.87, 1.899, 1.913, 1.913, 1.891, 1.855, 1.812, 1.761, 0.949, 1.703, 1.036, 1.652, 1.123, 1.594, 1.5, 1.297, 1.21, 1.406]
c_92 =  [0.421, 0.443, 0.462, 0.477, 0.496, 0.515, 0.534, 0.558, 0.582, 0.605, 0.631, 0.661, 0.694, 0.726, 0.758, 0.782, 0.805, 0.825, 0.842, 0.842, 0.857, 0.857, 0.866, 0.868, 0.872, 0.876, 0.876, 0.876, 0.878]
etta_92 = [0.92 for i in range(len(mu_92))]

mu_93 =  [1.094, 1.131, 1.174, 1.21, 1.254, 1.297, 1.341, 1.384, 1.435, 1.471, 1.522, 1.565, 1.594, 1.623, 1.645, 1.66, 1.674, 1.667, 1.645, 1.616, 0.935, 0.986, 1.58, 1.036, 1.102, 1.529, 1.181, 1.464, 1.37, 1.276]
c_93 =  [0.391, 0.404, 0.419, 0.434, 0.451, 0.474, 0.494, 0.515, 0.539, 0.56, 0.586, 0.61, 0.631, 0.655, 0.676, 0.698, 0.719, 0.743, 0.765, 0.786, 0.79, 0.797, 0.803, 0.803, 0.81, 0.814, 0.818, 0.823, 0.827, 0.827]
etta_93 = [0.93 for i in range(len(mu_93))]

mu_94 =  [0.986, 1.029, 1.065, 1.109, 1.152, 1.196, 1.239, 1.268, 1.283, 1.305, 1.319, 1.326, 1.319, 1.305, 0.877, 0.92, 1.283, 0.971, 1.029, 1.232, 1.16, 1.102]
c_94 =  [0.383, 0.398, 0.413, 0.432, 0.455, 0.479, 0.509, 0.534, 0.556, 0.58, 0.605, 0.629, 0.659, 0.687, 0.7, 0.709, 0.711, 0.717, 0.726, 0.728, 0.732, 0.732]
etta_94 = [0.94 for i in range(len(mu_94))]

data_etta = np.array([[1.005, 3.122, 0.85], [1.029, 3.093, 0.85], [1.057, 3.064, 0.85], [1.083, 3.029, 0.85], [1.117, 2.986, 0.85], [1.134, 2.957, 0.85], [1.157, 2.913, 0.85], [1.183, 2.862, 0.85], [1.2, 2.812, 0.85], [1.217, 2.768, 0.85], [1.232, 2.71, 0.85], [1.247, 2.645, 0.85], [1.255, 2.602, 0.85], [1.264, 2.544, 0.85], [1.27, 2.485, 0.85], [1.274, 2.441, 0.85], [0.939, 3.0, 0.86], [0.965, 2.986, 0.86], [0.99, 2.964, 0.86], [1.012, 2.942, 0.86], [1.038, 2.92, 0.86], [1.061, 2.891, 0.86], [1.083, 2.87, 0.86], [1.104, 2.834, 0.86], [1.123, 2.797, 0.86], [1.143, 2.747, 0.86], [1.166, 2.681, 0.86], [1.185, 2.623, 0.86], [1.2, 2.565, 0.86], [1.215, 2.478, 0.86], [1.226, 2.419, 0.86], [1.234, 2.36, 0.86], [1.24, 2.309, 0.86], [1.247, 2.243, 0.86], [0.859, 2.819, 0.87], [0.885, 2.812, 0.87], [0.917, 2.797, 0.87], [0.945, 2.783, 0.87], [0.984, 2.761, 0.87], [1.016, 2.732, 0.87], [1.053, 2.689, 0.87], [1.083, 2.631, 0.87], [1.117, 2.558, 0.87], [1.14, 2.485, 0.87], [1.162, 2.405, 0.87], [1.174, 1.399, 0.87], [1.179, 2.324, 0.87], [1.181, 1.457, 0.87], [1.187, 1.522, 0.87], [1.189, 2.243, 0.87], [1.196, 1.623, 0.87], [1.198, 2.155, 0.87], [1.202, 2.044, 0.87], [1.202, 1.718, 0.87], [1.204, 1.928, 0.87], [1.204, 1.826, 0.87], [0.765, 2.616, 0.88], [0.782, 2.616, 0.88], [0.803, 2.623, 0.88], [0.831, 2.631, 0.88], [0.861, 2.623, 0.88], [0.885, 2.616, 0.88], [0.919, 2.602, 0.88], [0.954, 2.587, 0.88], [0.982, 2.558, 0.88], [1.005, 2.522, 0.88], [1.031, 2.478, 0.88], [1.055, 2.427, 0.88], [1.078, 2.36, 0.88], [1.098, 2.309, 0.88], [1.111, 2.235, 0.88], [1.115, 1.283, 0.88], [1.121, 1.348, 0.88], [1.123, 2.147, 0.88], [1.13, 1.42, 0.88], [1.134, 2.044, 0.88], [1.136, 1.486, 0.88], [1.14, 1.957, 0.88], [1.143, 1.877, 0.88], [1.145, 1.594, 0.88], [1.147, 1.79, 0.88], [1.147, 1.696, 0.88], [0.672, 2.36, 0.89], [0.696, 2.383, 0.89], [0.722, 2.405, 0.89], [0.747, 2.419, 0.89], [0.767, 2.427, 0.89], [0.801, 2.441, 0.89], [0.835, 2.441, 0.89], [0.878, 2.427, 0.89], [0.911, 2.397, 0.89], [0.941, 2.375, 0.89], [0.967, 2.331, 0.89], [0.992, 2.287, 0.89], [1.016, 2.228, 0.89], [1.04, 2.155, 0.89], [1.053, 2.088, 0.89], [1.059, 1.189, 0.89], [1.066, 2.037, 0.89], [1.07, 1.261, 0.89], [1.074, 1.964, 0.89], [1.078, 1.326, 0.89], [1.081, 1.899, 0.89], [1.085, 1.812, 0.89], [1.087, 1.435, 0.89], [1.089, 1.725, 0.89], [1.091, 1.645, 0.89], [1.091, 1.551, 0.89], [0.543, 1.986, 0.9], [0.567, 2.037, 0.9], [0.588, 2.081, 0.9], [0.623, 2.14, 0.9], [0.659, 2.191, 0.9], [0.694, 2.235, 0.9], [0.732, 2.258, 0.9], [0.773, 2.272, 0.9], [0.799, 2.272, 0.9], [0.831, 2.265, 0.9], [0.863, 2.243, 0.9], [0.889, 2.206, 0.9], [0.919, 2.155, 0.9], [0.945, 2.096, 0.9], [0.965, 2.037, 0.9], [0.99, 1.942, 0.9], [0.99, 1.065, 0.9], [0.997, 1.102, 0.9], [1.005, 1.834, 0.9], [1.005, 1.167, 0.9], [1.018, 1.261, 0.9], [1.02, 1.703, 0.9], [1.029, 1.609, 0.9], [1.029, 1.493, 0.9], [1.029, 1.391, 0.9], [0.46, 1.602, 0.91], [0.472, 1.638, 0.91], [0.489, 1.689, 0.91], [0.515, 1.754, 0.91], [0.539, 1.812, 0.91], [0.562, 1.87, 0.91], [0.586, 1.92, 0.91], [0.614, 1.971, 0.91], [0.638, 2.008, 0.91], [0.672, 2.044, 0.91], [0.7, 2.066, 0.91], [0.734, 2.088, 0.91], [0.775, 2.088, 0.91], [0.816, 2.074, 0.91], [0.842, 2.044, 0.91], [0.868, 2.008, 0.91], [0.894, 1.942, 0.91], [0.915, 1.877, 0.91], [0.919, 1.0, 0.91], [0.926, 1.051, 0.91], [0.937, 1.797, 0.91], [0.939, 1.145, 0.91], [0.945, 1.232, 0.91], [0.952, 1.681, 0.91], [0.954, 1.319, 0.91], [0.958, 1.587, 0.91], [0.958, 1.406, 0.91], [0.962, 1.5, 0.91], [0.421, 1.319, 0.92], [0.443, 1.384, 0.92], [0.462, 1.435, 0.92], [0.477, 1.471, 0.92], [0.496, 1.522, 0.92], [0.515, 1.573, 0.92], [0.534, 1.623, 0.92], [0.558, 1.674, 0.92], [0.582, 1.725, 0.92], [0.605, 1.776, 0.92], [0.631, 1.826, 0.92], [0.661, 1.87, 0.92], [0.694, 1.899, 0.92], [0.726, 1.913, 0.92], [0.758, 1.913, 0.92], [0.782, 1.891, 0.92], [0.805, 1.855, 0.92], [0.825, 1.812, 0.92], [0.842, 1.761, 0.92], [0.842, 0.949, 0.92], [0.857, 1.703, 0.92], [0.857, 1.036, 0.92], [0.866, 1.652, 0.92], [0.868, 1.123, 0.92], [0.872, 1.594, 0.92], [0.876, 1.5, 0.92], [0.876, 1.297, 0.92], [0.876, 1.21, 0.92], [0.878, 1.406, 0.92], [0.391, 1.094, 0.93], [0.404, 1.131, 0.93], [0.419, 1.174, 0.93], [0.434, 1.21, 0.93], [0.451, 1.254, 0.93], [0.474, 1.297, 0.93], [0.494, 1.341, 0.93], [0.515, 1.384, 0.93], [0.539, 1.435, 0.93], [0.56, 1.471, 0.93], [0.586, 1.522, 0.93], [0.61, 1.565, 0.93], [0.631, 1.594, 0.93], [0.655, 1.623, 0.93], [0.676, 1.645, 0.93], [0.698, 1.66, 0.93], [0.719, 1.674, 0.93], [0.743, 1.667, 0.93], [0.765, 1.645, 0.93], [0.786, 1.616, 0.93], [0.79, 0.935, 0.93], [0.797, 0.986, 0.93], [0.803, 1.58, 0.93], [0.803, 1.036, 0.93], [0.81, 1.102, 0.93], [0.814, 1.529, 0.93], [0.818, 1.181, 0.93], [0.823, 1.464, 0.93], [0.827, 1.37, 0.93], [0.827, 1.276, 0.93], [0.383, 0.986, 0.94], [0.398, 1.029, 0.94], [0.413, 1.065, 0.94], [0.432, 1.109, 0.94], [0.455, 1.152, 0.94], [0.479, 1.196, 0.94], [0.509, 1.239, 0.94], [0.534, 1.268, 0.94], [0.556, 1.283, 0.94], [0.58, 1.305, 0.94], [0.605, 1.319, 0.94], [0.629, 1.326, 0.94], [0.659, 1.319, 0.94], [0.687, 1.305, 0.94], [0.7, 0.877, 0.94], [0.709, 0.92, 0.94], [0.711, 1.283, 0.94], [0.717, 0.971, 0.94], [0.726, 1.029, 0.94], [0.728, 1.232, 0.94], [0.732, 1.16, 0.94], [0.732, 1.102, 0.94]])

def polyfit2d(x, y, z, order=4):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def smith_etta(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    p = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        p += a * x**i * y**j
    return p

x1 = np.array((np.concatenate((c_85, c_86, c_87, c_88, c_89, c_90, c_91, c_92, c_93, c_94), axis= 0)))
y1 = np.array((np.concatenate((mu_85, mu_86, mu_87, mu_88, mu_89, mu_90, mu_91, mu_92, mu_93, mu_94), axis= 0)))
z1 = np.array((np.concatenate((etta_85, etta_86, etta_87, etta_88, etta_89, etta_90, etta_91, etta_92, etta_93, etta_94), axis= 0)))
val = polyfit2d(x1, y1, z1)

# print(smith_etta(0.8260511580776847, 2.0672147846222697, val))

def I(T, fuel):
    return fuel_inter(fuel.get('h'), fuel.get('temperature_K'), T)

def T(h, fuel):
    return fuel_inter(fuel.get('temperature_K'), fuel.get('h'), h)

def k(T, fuel):
    return fuel_inter(fuel.get('k'), fuel.get('temperature_K'), T)  

def S(T, fuel):
    return fuel_inter(fuel.get('s_0'), fuel.get('temperature_K'), T)

def Cp(T, fuel):
    return fuel_inter(fuel.get('Cp'), fuel.get('temperature_K'), T)

def lamda(velocity, T, fuel, sch):
    return (velocity / sqrt((2e3 * k(T, fuel) * sch[1]['R_g'] * T) / (k(T, fuel) + 1)))

def PI_lamda(lamda, T, fuel):
    return ((1 - (((k(T, fuel) - 1) / (k(T, fuel) + 1)) * lamda**2))**(k(T, fuel) / (k(T, fuel) - 1)))

def q_lamda(lamda, T, fuel):
    return (lamda * (((k(T, fuel) + 1) / 2)**(1 / (k(T, fuel) - 1))) * ((1 - (((k(T, fuel) - 1) / (k(T, fuel) + 1)) * lamda**2))**(k(T, fuel) / (k(T, fuel) - 1))))

def tau_lamda(lamda, T, fuel):
    return(1 - ((k(T, fuel) - 1) / (k(T, fuel) + 1)) * (lamda)**2)

def m_(T, fuel, sch):
    return (sqrt((k(T, fuel) / sch[1]['R_g'] *  10**(-3)) * (2 / (k(T, fuel) + 1))**((k(T, fuel) + 1) / (k(T, fuel) - 1))))

def Mu(T):
    return((44.3*(10**(-6))) * ((T/1073)**0.678))

def Re(P, R, T, C, b ):
    ro = P / (R * T)
    return((C * b * ro) / Mu(T))