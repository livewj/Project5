import numpy as np

def read(filename):
	#function that reads the x and y coordinate
	infile = open(filename, 'r')
	Timestep = []
	Time =[]
	Temperature = []
	DiffusionConstant = []
	relevant_lines = infile.readlines()[1:] #skips the irrelevant lines
	for line in relevant_lines: 
		data = line.split()
		Timestep.append(float(data[0]))
		Time.append(float(data[1]))
		Temperature.append(float(data[2]))
		DiffusionConstant.append(float(data[3]))
	infile.close()
	Timestep = np.array(Timestep)
	Time = np.array(Time)
	Temperature = np.array(Temperature)
	DiffusionConstant = np.array(DiffusionConstant)

	return Timestep, Time, Temperature, DiffusionConstant

#Timestep, Time, Temperature, DiffusionConstant = read('statistics.txt')

T_init = []
for i in range(1, 1001, 50):
	T_init.append(i)
T_init = np.array(T_init)

T_ratio = np.array([0.669704,0.741992,0.678682,0.690175,0.693458,0.719683,0.759031,0.779948,0.709018,0.740829,0.765328,0.744585,0.761976,0.755197,0.720936,0.738801,0.766893,0.729957,0.729681,0.751274])

T = T_ratio*T_init #temperature at equilibrium

import matplotlib.pyplot as plt

plt.plot(T, T_ratio, 'go')
plt.rcParams.update({'font.size': 14})
plt.ylabel('$T/T_i [m^2/s]$')
plt.xlabel('$T$ [K]')
plt.show()


