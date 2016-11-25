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

D_mean = [3.08224e-12
,1.67068e-10
,3.05171e-10
,4.52251e-10
,6.13243e-10
,7.64216e-10
,1.04179e-09
,1.22794e-09
,1.25638e-09
,1.71514e-09
,1.81977e-09
,1.96159e-09
,2.13156e-09
,2.44635e-09
,2.81205e-09
,3.05605e-09
,3.54855e-09
,3.59755e-09
,3.56333e-09
,4.87054e-09]

import matplotlib.pyplot as plt

plt.plot(T, D_mean, 'go')
plt.rcParams.update({'font.size': 14})
plt.ylabel('$D [m^2/s]$')
plt.xlabel('$T$ [K]')
plt.show()


