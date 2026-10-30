## Keeping the temperature correction functions separate from MathEquations

# you should import this with a namespace because there are a lot of variables

'''
this is implementing the function from the paper
Photoelectric absorption cross section of silicon near the band gap from room temperature to sub-Kelvin temperature
https://doi.org/10.1063/5.0038392

other papers used for the direct bandgap, giving a big better (smoother) threshold behaviour
(from solar cell paper) https://doi.org/10.1016/0038-1101(79)90128-X
which incorrectly uses a different equation from the original paper https://doi.org/10.1063/1.356496 
we use the equation from the later paper, and the same beta, gamma values for bandgap vs temp

See the equations and some example plots in temperature_nb/pdAbs_v1.ipynb
'''

import numpy as np


## Some silicon constants
rho_si = 2.335 # g/cm3
k = 8.617333262145e-5  # Boltzmann constant in eV/K
hc = 1239.84193 # eV nm

## Constants that describe the temperature dependence of the band gap

# Ai - proportionality constant. Some are temp dependent
# c - electron-phonon coupling constant
# gamma, beta - describe bandgap temperature dependence

Ephonons = np.array([18.27, 57.73])*1e-3 # eV, phonon energy. doing things in eV (phonon energy is 18 meV)

# transverse acoustic, optical phonons
Eph_ta = Ephonons[0]    # TA - 212 kelvin
Cta = 5.5               # coupling constant

Eph_to = Ephonons[1] # TO - 675 K (from Solar Cells paper reference
Cto = 4.0

betaL = 7.021e-4 # eV/K
gammaL = 1108.0 # K

# first bandgap (temp dependance)
Eg_01 = 1.134 # eV
c0 = 325 # eV-2/cm
c1 = 1.7e-3 # 1/K
A1_phonon = lambda T:  c0*np.exp(-c1*T)

# second bandgap, no temp dependance
Eg_02 = 2.5 # eV
A2_phonon = 6.0e3 # eV^-2/cm

## Direct Bandgap
# (from solar cell paper) https://doi.org/10.1016/0038-1101(79)90128-X
# Ad = 1.052e6 # eV^-2/cm
# actually no. using function and parameters from the paper https://doi.org/10.1063/1.356496  (Bucher)
# same beta, gamma parameters as indirect bandgap
Egd_0 = 3.2
Adnew = 2.3e7

# no temp in here. function from https://doi.org/10.1063/1.356496
def DirectTerm(Ephoton, EdirectGap, Ad):
    return Ad * (Ephoton - EdirectGap)**1.5 / Ephoton # bucher, what we want. smoother onset and better value wrt optics literature
    # return Ad * (Ephoton - EdirectGap)**0.5 # solar cell paper

## Calculate how the bandgap changes with temperature
def Fbandgap(Eg0, T, beta=betaL, gamma=gammaL):
    return Eg0 - beta*T**2/(T+gamma)

# The term for absorption or emission of a single phonon
def PhononTerm(Egam, T, EgN, Ephonon, Cph):
    
    absorption = Cph * (Egam - EgN + Ephonon)**2/(np.exp(Ephonon/(k*T)) - 1)
    emission = Cph * (Egam - EgN - Ephonon)**2 /(1 - np.exp(-Ephonon/(k*T)))
    # absorbtion only condition
    Edif = Egam - EgN
    conditions = [Egam < EgN - Ephonon, (Egam > EgN-Ephonon) & (Egam < EgN + Ephonon), Egam >= EgN + Ephonon]
    # insufficent energy, absorption, absorption + emission
    outputs = [0, absorption, absorption + emission]    
    return np.select(conditions, outputs)

# The absorption coefficient for both sets of phonons and indirect bandgaps, plus direct bandgap
# include
# E_gamma must be a list or array
def alpha(E_gamma, T, includeDirect = True, fixDirectBandGap =  False): # units of 1/cm
    
    Eg1 = Fbandgap(Eg0=Eg_01, T=T)
    A1t = A1_phonon(T)
    
    Eg2 = Fbandgap(Eg0=Eg_02, T=T)
    A2t = A2_phonon # no temperature dependence
    
    # bandgap 1, TA phonon
    g1ta = PhononTerm(E_gamma, T, Eg1, Eph_ta, Cta)
    # bandgap 1, TO phonon
    g1to = PhononTerm(E_gamma, T, Eg1, Eph_to, Cto)
    bg1 = A1t * (g1ta + g1to)
    
    # bandgap 2, TA phonon
    g2ta = PhononTerm(E_gamma, T, Eg2, Eph_ta, Cta)
    # bandgap 2, TO phonon
    g2to = PhononTerm(E_gamma, T, Eg2, Eph_to, Cto)
    bg2 = A2t * (g2ta + g2to)
    
    # doing direct bandgap
    Edirect_temp = Fbandgap(Eg0=Egd_0, T=T)
    if fixDirectBandGap:
        Edirect_temp = Egd_0
    mask = E_gamma > Edirect_temp
    direct_cutoff = np.zeros_like(E_gamma)
    direct_cutoff[mask] = float(includeDirect) * DirectTerm(E_gamma[mask], Edirect_temp, Adnew)
    
    return bg1 + bg2 + direct_cutoff, A1t*g1ta, A1t*g1to, A2t*g2ta, A2t*g2to, direct_cutoff, Edirect_temp

## make a function that calculates the relative change in alpha for a given temperature change. input is Tdata, Tnew, wlen
# tdata is the temperature the input optical data is taken at (not the PDE data)
# the PDE data is taken at Tnew
# so yourNewAlphaToUse_takenAtTnew = oldAlpha_takenAtTdata * RelativeAlphaChange(Tdata, Tnew, wlen)
def RelativeAlphaChange(Tdata, Tnew, wlen, includeDirect = 1, fixDirectBandGap = False, numericalCorrection=True):
    ephoton = hc / wlen
    alpha_data, *_, EgapDirect_1 = alpha(E_gamma=ephoton, T=Tdata, includeDirect=includeDirect, fixDirectBandGap=fixDirectBandGap)
    alpha_new, *_, EgapDirect_2 = alpha(E_gamma=ephoton, T=Tnew, includeDirect=includeDirect, fixDirectBandGap=fixDirectBandGap)
    
    # print(EgapDirect_1, EgapDirect_2)
    # relativeChange_fixedDirect = alpha(E_gamma=ephoton, T=Tnew, includeDirect=1, fixDirectBandGap=True)[0]/alpha(E_gamma=ephoton, T=Tdata, includeDirect=1, fixDirectBandGap=True)[0]
    
    relative_change = alpha_new / alpha_data
    
    Egap_1_wlen = hc/EgapDirect_1
    Egap_2_wlen = hc/EgapDirect_2
    if numericalCorrection and min(wlen) < min(Egap_1_wlen, Egap_2_wlen):
        # fixing the low spike at the direct bandgap threshold
        # if numericalCorrection and (fixDirectBandGap == False) and (includeDirect == True):
        rightLimit = max(Egap_1_wlen, Egap_2_wlen) # sorting so any T ordering works
        
        if max(wlen) < max(Egap_1_wlen, Egap_2_wlen): # if the wlen range is all within the direct bandgap, we dont have section of the curve to correct with
            rightRelativeChange = RelativeAlphaChange(Tdata, Tnew, np.array(max(Egap_1_wlen, Egap_2_wlen) + 1), includeDirect = 0, fixDirectBandGap = False, numericalCorrection=False)
        else:
            rightRelativeChange = relative_change[wlen > rightLimit][0]
        relative_change[relative_change < rightRelativeChange] = rightRelativeChange
    return relative_change