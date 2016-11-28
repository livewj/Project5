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

T = [47.9141, 79.4194, 96.6324, 122.12, 147.195, 177.679, 213.86, 247.495, 247.632 , 280.041, 313.374, 332.292, 355.614, 370.431, 376.699, 422.769, 407.625, 428.792, 460.189, 503.035, 568.394, 612.934, 691.756, 745.181, 826.782, 859.003, 982.932, 1044.38, 1144.01, 1219.79, 1341.79, 1434.83, 1561.96, 1656.9, 1715.72]
D = [1.38188e-10, 2.21717e-10, 2.796e-10, 3.40718e-10, 4.25555e-10,  5.08163e-10, 6.61828e-10, 7.5041e-10, 7.49222e-10, 1.14175e-09, 1.15539e-09, 1.21354e-09, 1.4735e-09, 2.0452e-09, 2.47333e-09, 3.57932e-09, 3.27279e-09, 3.69056e-09, 4.46663e-09, 4.89587e-09, 5.94322e-09, 6.3852e-09, 7.54357e-09, 9.10429e-09, 1.05681e-08, 1.10233e-08, 1.19294e-08, 1.28984e-08, 1.52922e-08, 1.60626e-08, 1.67848e-08, 1.85233e-08, 2.0357e-08, 2.10049e-08, 2.12097e-08]
import matplotlib.pyplot as plt

fit = np.polyfit(T,D,1)
fit_fn = np.poly1d(fit) 
# fit_fn is now a function which takes in x and returns an estimate for y

plt.plot(T,D, 'ro', T, fit_fn(T), '--k')
plt.xlim(400, 1800)
plt.ylim(-0.5e-8, 2.5e-8)
plt.rcParams.update({'font.size': 14})
plt.ylabel('D $[m^2/s]$')
plt.xlabel('T [K]')
plt.legend(['Diffusion Constant', 'Lin fit'], loc='best')
plt.show()

