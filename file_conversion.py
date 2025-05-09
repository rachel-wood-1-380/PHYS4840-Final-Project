import numpy as np
import matplotlib.pyplot as plt
import os


##### Any file paths will need to be changed to reflect where you have saved them on your machine,
#####   and where you want them to be saved on your machine. 

# Period of each pulsar found via search.
# Time-step for each calculated by dividing period by 1024 (number of bins per period)

P_1 = 0.693747670469999970      # Period of J0006+1834 (s)
h_1 = P_1 / 1024                # time-step (s)

P_2 = 0.779478492               # Period of J0108-1431 (s)     
h_2 = P_2 / 1024                # time-step (s)

P_3 = 0.8327722180333           # Period of J0152-1637 (s)
h_3 = P_3 / 1024                # time-step (s)

P_4 = 0.630551455               # Period of J0206-4028 (s)
h_4 = P_4 / 1024                # time-step (s)

P_5 = 0.548944259               # Period of J0452-1759 (s)
h_5 = P_5 / 1024                # time-step (s)

P_6 = 1.097990983               # Period of J0946+0951 (s)
h_6 = P_6 / 1024                # time-step (s)

P_7 = 0.0031                    # Period of J1125-5825 (s)
h_7 = P_7 / 1024                # time-step (s)

P_8 = 0.82778                   # Period of J1807-2715 (s)
h_8 = P_8 / 1024                # time-step (s)

P_9 = 0.533                     # Period of J1832+0029 (s)
h_9 = P_9 / 1024                # time-step (s)

P_10 = 0.0160524236697          # Period of J2145-0750 (s)
h_10 = P_10 / 1024              # times-step (s)


h_values = np.array([h_1, h_2, h_3, h_4, h_5, h_6, h_7, h_8, h_9, h_10])

##### Original files can be found at: https://research.csiro.au/pulseatparkes/data-analysis/database/

def file_loop(directory_path, h_values):
	files = sorted([f for f in os.listdir(directory_path) if f.endswith('old.txt')])

	results = {}

	for filename, h in zip(files, h_values):
		file_path = os.path.join(directory_path, filename)

		column1, column2, index, A = np.loadtxt(file_path, unpack=True, usecols=(0,1,2,3))

		count = np.arange(len(A), dtype=int)

		timesteps = []
		for i in count:
			timesteps.append(i * h)

		array = np.column_stack([count, timesteps, A])

		key = os.path.splitext(filename)[0]
		results[key] = array

	return results

# Save files with new columns, now in a format usable for the program

def save_outputs(array_dict, output_dic):
	os.makedirs("/Pulsar_Program/Data_Files", exist_ok=True)

	for original_name, array in array_dict.items():
		new_name = original_name.replace('_old', '')
		output_path = os.path.join("/Pulsar_Program/Data_Files", f"{new_name}.txt")
		np.savetxt(output_path, array, fmt=['%d', '%.10f', '%.10f'])

	return


arrays = file_loop("/Pulsar_Program/Data_Files", h_values)

save_outputs(arrays, "/Pulsar_Program/Data_Files")

