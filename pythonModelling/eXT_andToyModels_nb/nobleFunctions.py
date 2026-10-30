### functions taken From Austin's 'nEXO' as a lens calculation (in MSc thesis)
import numpy as np

#for normal xenon
#densities in g/cm^3
#photon energy in wavelength

#nEXO properties
height = 120.0 #nEXO height is 1.2 m
LXe_density = 2.9 #g/cm^3
radius = 338/2 
pmin = 1.0 #pressure in bar
pmax = 1.5
Tlow = 165.5 # gabe ortega 2017 simulations
Thigh = 166.5 # one degree from top to bottom
deltaT = Thigh - Tlow

#density scaling
M_xenon = 131.296
M_other = (M_xenon-0.08857*136)/(1 - 0.08857)
M_enr = (1-0.9)*M_other + 0.9*136 #average molar mass of enriched xenon
scaling = M_enr/M_other 

#gives density for LXe as a function of T (within liquid part of phase space)
#using the data from the MJ Terry paper (first measurement, seems to agree with others)
def density_t(T):
	Tr = 161.38
	Vr = 44.312
	A = 0.973/10
	B = 0.325/1e3
	C = 0.539/1e5
	D = -0.297/1e7
	Vm = Vr + A*(T - Tr) + B*(T - Tr)**2 + C*(T-Tr)**3 + D*(T - Tr)**4
	density = M_xenon/Vm
	return density

#get A value
def Avalue(lamb):
	En = 1239.84/lamb
	P = [71.23, 77.75, 1384.89]
	E = [8.4, 8.81, 13.2]

	A = 0
	for i in range(0,3):
		A = A + P[i]/(En**2 - E[i]**2)

	return -A

#only the first term of expansion, density as argument
def n_index_rho(lamb, rho):
	rho_m = rho/M_xenon
	A = Avalue(lamb)

	Fll = A*rho_m
	n = -np.sqrt(4*(1 - Fll)*(2*Fll + 1))/(2*(Fll - 1))
	return n

#only the first term of expansion, temp as argument
def n_index_temp(lamb, temp):
	rho = density_t(temp)
	rho_m = rho/M_xenon
	A = Avalue(lamb)

	Fll = A*rho_m
	n = -np.sqrt(4*(1 - Fll)*(2*Fll + 1))/(2*(Fll - 1))
	return n

######################### nEXO as a lens stuff #########################
#temperature as a function of Z, assume linear gradient
def temp(z):
	grad = deltaT/height
	return Tlow + grad*z

#ignores density fluctuations throughout detector, but minimal. just want to show that density doesnt care about pressure
def extra_pressure(z):
	force_over_area = (height - z)*density*scaling*9.81*10*1e-5 #in Bars
	return force_over_area

def density_height(z):
	temp = temp(z)
	density = density_t(temp)
	return density