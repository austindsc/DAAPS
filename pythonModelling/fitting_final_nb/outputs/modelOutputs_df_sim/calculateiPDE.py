import numpy as np

def Pcar_parametrized(V, Vc, A=1.0):
    return A * (1 - np.exp(-V/Vc))

def QuantumPe(wlen, Pe, cutinWl=300, power=0.75): ## 300, 0.75 works for Harry's QY data
    qyield = np.where(wlen < cutinWl, (cutinWl/wlen)**power, 1)
    return 1 - (1-Pe)**qyield

## get Wp, Wn from provided 1D lookup table
## 300, 0.75 works for Harry's QY data. assume HPK and FBK yields are the same
## arXiv:2410.13033
def InternalPDE(wlen, V, Wp, Wn, Ve, Vh, Ae=1.0, Ah=1.0, quantum=True):
    Pe = Pcar_parametrized(Ve, V, Ae)
    Ph = Pcar_parametrized(Vh, V, Ah)
    if quantum:
        qpe = QuantumPe(wlen, Pe)
    else:
        qpe = Pe
    return Wp * qpe + Wn * Ph

## July 2025, values for devices are:
pAvalanche_hpk = {"Ve": 1.79, "Vh": 16.6, "Ae": 1.0, "Ah": 1.0}
pAvalanche_fbk = {"Ve": 1.83, "Vh": 18.77, "Ae": 1.0, "Ah": 1.0}

##### how to use #####
# load in optics data, remember to multiply by FF
# FF = 0.6 for VUV4
# FF = 0.8 for FBK HD3

# when a photon makes it in, you need to do the internal efficiency calculation
# this code/datasheet gives you the internal PDE value for a given wavelength and voltage
# Wp, Wn are the weights for the electron and hole contributions, respectively
# Wp, Wn are taken from a lookup table, Ve and Vh are the avalanche parameters 
# that we think are fixed for a given device, temperature

# this allows the user to calculate iPDE using any voltage values they want (uncertainty analysis?)

# user can toggle quantum on/off in InternalPDE function call
# if studying physics beyond 300 nm, set quantum=False or even write a new function
# without it, to avoid a if/else call and any other confusion


"""
import calculateiPDE as iPDE

V = 5.28  # example voltage
wlen = 169  # example wavelength

#load in Wp, Wn from input table however you want (using df most likely)
Wp_hpk, Wn_hpk = 0.5, 0.5  # example Wp

example_iPDE_HPK = iPDE.InternalPDE(wlen, V, Wp=Wp_hpk, Wn=Wn_hpk,
                               Ve=pAvalanche_hpk["Ve"], Vh=pAvalanche_hpk["Vh"], quantum=True)

print("Yo! this is the iPDE value for HPK at 169 nm and 5.28 V:", example_iPDE_HPK)                  
"""