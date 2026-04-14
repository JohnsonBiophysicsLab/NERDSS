import numpy as np
import matplotlib.pyplot as plt

# Define the names of the input files
input_files = ["copy_numbers_time_0.dat", "copy_numbers_time_1.dat", 
               "copy_numbers_time_2.dat", "copy_numbers_time_3.dat"]

# Load the data from the input files into numpy arrays
data_list = []
for input_file in input_files:
    data = np.loadtxt(input_file, delimiter=",", skiprows=1)
    data_list.append(data)

# Concatenate the data from the input files into a single numpy array
data = np.concatenate(data_list)

# Sum the A(a) values for each unique time point
unique_times, indices = np.unique(data[:, 0], return_inverse=True)
a_sum = np.zeros_like(unique_times)
np.add.at(a_sum, indices, data[:, 1])

# Plot the data
plt.plot(unique_times, a_sum)
plt.xlabel("Time (s)")
plt.ylabel("A(a)")
plt.show()
