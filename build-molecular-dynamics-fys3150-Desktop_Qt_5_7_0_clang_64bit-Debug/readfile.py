import numpy as np

def read(filename):
	#function that reads the x and y coordinate
	infile = open(filename, 'r')
	Timestep = []
	Time =[]
	Temperature = []
	KineticEnergy = []
	PotentialEnergy = []
	TotalEnergy = []
	DiffusionConstant = []
	relevant_lines = infile.readlines()[1:] #skips the irrelevant lines
	for line in relevant_lines: 
		data = line.split()
		Timestep.append(float(data[0]))
		Time.append(float(data[1]))
		Temperature.append(float(data[2]))
		KineticEnergy.append(float(data[3]))
		PotentialEnergy.append(float(data[4]))
		TotalEnergy.append(float(data[5]))
		DiffusionConstant.append(float(data[6]))
	infile.close()
	Timestep = np.array(Timestep)
	Time = np.array(Time)
	Temperature = np.array(Temperature)
	KineticEnergy = np.array(KineticEnergy)
	PotentialEnergy = np.array(PotentialEnergy)
	TotalEnergy = np.array(TotalEnergy)
	DiffusionConstant = np.array(DiffusionConstant)

	return Timestep, Time, Temperature, KineticEnergy, PotentialEnergy, TotalEnergy, DiffusionConstant

Timestep, Time, Temperature, KineticEnergy, PotentialEnergy, TotalEnergy, DiffusionConstant = read('statistics.txt')


import matplotlib.pyplot as plt

plt.plot(Temperature, DiffusionConstant)
plt.rcParams.update({'font.size': 14})
plt.ylabel('Diffusionconstant $[m^2/s]$')
plt.xlabel('Temperature [K]')
plt.show()


