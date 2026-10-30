import numpy as np
import scipy.integrate as spi

"""
Okay. We will write these 'normally', meaning
these functions don't know about our custom classes.
"""


######## Absorption Functions ########
## Gets passed pre-established values to be used in pure math equations - makes it easier to spot a potential bud in the math
def Mu_fromK(wl, k):
    mu = 4*np.pi*k/wl # AbsLenght = 1/mu
    return mu

def Mu_fromIndex(wl, index):
    k = index.imag
    mu = 4*np.pi*k/wl
    return mu

def kFromMu(wl, mu):
    k = wl * mu / (4 * np.pi)  # k = mu * wl / (4 * pi)
    return k

# its the thickness of each layer that should be passed here
# so abs_layers gets (dp_star, Xpn-dp_star) for calculating Wp on a P on N junction
def Abs_layers(d_above, d_region, mu, theta=0): # d_above and d_region are passed from PNInfo, theta default is 0
    # living through above x dying in region of interest
    living = Abs_singleLayer(d_above, mu, theta)
    dying = (1-Abs_singleLayer(d_region, mu, theta))
    return living * dying

def Abs_singleLayer(d, mu, theta=0): # this function is how much is surviving after distance d
    d_eff = d/np.cos(theta)
    return np.exp(-d_eff * mu)

#### some diffusion models here
# this calculates the absorption in the dead layer multiplies it
# by the flat Pdiffusion probability
# this whole thing gets multiplied by Pe
def Abs_halfBox_Wp(d1, d2, mu, Pdiff = 0.5, theta=0):
    living = Abs_singleLayer(d=d1, mu=mu, theta=theta)
    Wp_hardEdge = Abs_layers(d_above=d1, d_region=d2, mu=mu, theta=theta)
    return Pdiff * (1 - living) + Wp_hardEdge

##### Testing for triple junction formulation
# follow some other convention here
def Pgradient(x, dwidth, Pstart, Pend):
    slope = (Pstart - Pend) / dwidth
    return Pstart - slope * x

def differential_absorption(x, mu, dwidth, Pstart, Pend):
    val = mu * np.exp(-x * mu) * Pgradient(x, dwidth, Pstart, Pend)
    return val

@np.vectorize
def Abs_Pslope_convolution(dx_p, dx_n, mu, Pstart, Pend, theta=0):
    d_above_eff = dx_p/np.cos(theta)
    fracBefore = np.exp(-d_above_eff * mu)
    
    dwidth = dx_n - dx_p
    dwidth_eff = dwidth/np.cos(theta)
    
    pWeightedAbs = spi.quad(differential_absorption, 0, dwidth_eff, args=(mu, dwidth_eff, Pstart, Pend))[0]
    
    return fracBefore * pWeightedAbs
##### Done testing for triple junction formulation

#### Avalanche Production Functions ####
def Pcar_parametrized(V, Vc, A=1.0):
    return A * (1 - np.exp(-V/Vc))

def QuantumPe(wlen, Pe, cutinWl, power): ## 300, 0.75 works for Harry's QY data
    qyield = np.where(wlen < cutinWl, (cutinWl/wlen)**power, 1)
    return 1 - (1-Pe)**qyield

def InternalPDE(wlen, V, Wp, Wn, Ve, Vh, Ae=1.0, Ah=1.0):
    Pe = Pcar_parametrized(Ve, V, Ae)
    Ph = Pcar_parametrized(Vh, V, Ah)
    return Wp * QuantumPe(wlen, Pe, cutinWl=300, power=0.75) + Wn * Ph

######## Optics Functions ########
def SellmeierIndex(wlen, a0, Auv, Air, wluv, wlir):
    n_squared = a0 + Auv * (wlen**2) / (wlen**2 - wluv**2) + Air * (wlen**2) / (wlen**2 - wlir**2)
    return np.array(n_squared**0.5)  # Return complex index n + ik, with k = 0 for Sellmeier. syke just return n

# this keep imaginary part of refractive index
def Snells_Law(n1, n2, theta1):
        """
        Calculates the angle of refraction using Snell's law for complex refractive indices.

        Parameters:
        - n1: complex refractive index of the first medium
        - n2: complex refractive index of the second medium
        - theta1: angle of incidence in radians

        Returns:
        - theta2: angle of refraction in radians
        """
        ## this formulation is correct (keeping imaginary thetas)
        theta2 = np.arcsin((n1 / n2) * np.sin(theta1))
        return theta2

# returns all four field coefficients
def Fresnel_Coefficients(n1, n2, theta1, theta2):
        """
        Calculates the Fresnel coefficients for reflection and transmission at an interface.

        Parameters:
        - n1: complex refractive index of the first medium
        - n2: complex refractive index of the second medium
        - theta1: angle of incidence in radians
        - theta2: angle of refraction in radians

        Returns:
        - r_s: reflection coefficient for s-polarized light
        - r_p: reflection coefficient for p-polarized light
        - t_s: transmission coefficient for s-polarized light
        - t_p: transmission coefficient for p-polarized light
        """
        r_s = (n1 * np.cos(theta1) - n2 * np.cos(theta2)) / (n1 * np.cos(theta1) + n2 * np.cos(theta2))
        r_p = (n2 * np.cos(theta1) - n1 * np.cos(theta2)) / (n2 * np.cos(theta1) + n1 * np.cos(theta2))
        t_s = (2 * n1 * np.cos(theta1)) / (n1 * np.cos(theta1) + n2 * np.cos(theta2))
        t_p = (2 * n1 * np.cos(theta1)) / (n2 * np.cos(theta1) + n1 * np.cos(theta2))
        return r_s, r_p, t_s, t_p

def Poynting_term(n1, n2, theta1, theta2):
    ## pt = np.real_if_close(n2 * np.cos(theta2) / (n1 * np.cos(theta1))) ## maybe do this instead?
    return np.real(n2 * np.cos(theta2)) / np.real(n1 * np.cos(theta1))

## using Poynting vector as well
def Fresnel_RefAndTrans(n1, n2, theta1, theta2):
        pt = Poynting_term(n1, n2, theta1, theta2)
        r_s, r_p, t_s, t_p = Fresnel_Coefficients(n1, n2, theta1, theta2)
        # t_s = r_s + 1.0
        # t_p = (r_p + 1.0) * n1/n2
        R = (np.abs(r_s)**2 + np.abs(r_p)**2) / 2
        T = pt * (np.abs(t_s)**2 + np.abs(t_p)**2) / 2 
        return R, T
    
# for a photon going from medium 1-->2-->3, where we can have internal reflections between 
# 2 and 3, but not between 1 and 2. This comes from the geometric series
def Ray_TwoLayerTransmission(Ttotal, Refftotal):
    return Ttotal / (1 - Refftotal)

############ Wave Transmission Functions ################
# this comes from Kurtis' thesis, chapter 3
# a photon going from medium 1-->2-->3, where we can have internal reflections
# but we do it with wave equations and field coefficients

## define as 4pi instead of 2pi (different than in some references). factor of 2
## carried throughout exponentials below
def twolayer_delta(wlen, n2, d, theta):
    return 4*np.pi*n2*d*np.cos(theta)/wlen

## for single P or S
def Wave_TwoLayerTransmission_lastTerm(t12, t23, r23, r12, delta):
    tt = (t12*t23) / (1 - r23*r12*np.exp(-1j*delta) )
    return np.conjugate(tt) * tt

def Wave_ttTerm(t12, t23, r23, r12, delta):
    """
    Returns:
    - ttt: total transmission field term
    """
    tt = (t12 * t23) * np.exp(1j * delta / 2) / (1 + r12 * r23 * np.exp(1j * delta))
    return np.conjugate(tt) * tt

def Wave_twoLayer_Teff(wlen, thick, n1, n2, n3, theta1, theta2, theta3):
    # equations from this reference
    # https://nano-optics.ac.nz/planar/articles/fresnel.html
    # 1-->2
    r12_s, r12_p, t12_s, t12_p = Fresnel_Coefficients(n1, n2, theta1, theta2)
    # 2-->3
    r23_s, r23_p, t23_s, t23_p = Fresnel_Coefficients(n2, n3, theta2, theta3)
    delta = twolayer_delta(wlen=wlen, n2=n2, d=thick, theta=theta2)
    
    ## P polarized
    # t_p = Wave_tttTerm(t12_p, t23_p, r23_p, r12_p, delta)
    # t_s = (t12_s*t23_s) * np.exp(1j*delta/2)/ (1 + r12_s*r23_s*np.exp(1j*delta) )
    T_p = Wave_ttTerm(t12_p, t23_p, r23_p, r12_p, delta)
    T_s = Wave_ttTerm(t12_s, t23_s, r23_s, r12_s, delta)
    return (T_p + T_s) / 2

## equation we use
## calls Wave_twoLayer_Teff and Wave_ttTerm, with poynting coupling
def Wave_TwoLayerTransmission_total(wlen, thick, n1, n2, n3, theta1, theta2, theta3):
    coupling = Poynting_term(n1, n3, theta1, theta3)
    Teff = Wave_twoLayer_Teff(wlen, thick, n1, n2, n3, theta1, theta2, theta3)
    return coupling * Teff

## different wrapper for above function
def Wave_TwoLayerTransmission_total_noThetas(wlen, thick, n1, n2, n3, theta1):
    theta2 = Snells_Law(n1, n2, theta1)
    theta3 = Snells_Law(n2, n3, theta2)
    return Wave_TwoLayerTransmission_total(wlen, thick, n1, n2, n3, theta1, theta2, theta3)
  
#### for testing some Poynting p component stuff against reflectivity data
## same as below, but returns both polarizations separately (for testing, not most solid framework)
def Wave_TLR_polarized(wlen, theta1, thick, n1, n2, n3):
    theta2 = Snells_Law(n1, n2, theta1)
    theta3 = Snells_Law(n2, n3, theta2)
    
    # 1-->2
    r12_s, r12_p, t12_s, t12_p = Fresnel_Coefficients(n1, n2, theta1, theta2)
    # 2-->3
    r23_s, r23_p, t23_s, t23_p = Fresnel_Coefficients(n2, n3, theta2, theta3)
    delta = twolayer_delta(wlen=wlen, n2=n2, d=thick, theta=theta2)
    
    ## P polarized
    r_p = (r12_p + r23_p * np.exp(1j*delta)) / (1 + r12_p * r23_p * np.exp(1j*delta))
    r_s = (r12_s + r23_s * np.exp(1j*delta)) / (1 + r12_s * r23_s * np.exp(1j*delta))

    R_p = r_p * np.conjugate(r_p)
    R_s = r_s * np.conjugate(r_s)
    
    return R_p, R_s

def kz(wlen, n, theta):
    """
    Calculate the kz component of the wavevector in a medium with refractive index n3 and angle theta3.
    """
    return 2 * np.pi/ wlen * n * np.cos(theta)

## here we are computing bulk absorption vs Z, different than power at interface (equations above)
# def Wave_TwoLayer_bulkAbs(wlen, thick, n1, n2, n3, theta1, theta2, theta3, z=0):
#     Teff = Wave_twoLayer_Teff(wlen, thick, n1, n2, n3, theta1, theta2, theta3)
    
#     preFac = 4*np.pi * n3.real * n3.imag / (wlen * n1.real * np.cos(theta1))
#     kz3 = kz(wlen, n3, theta3)
#     # exp = np.exp(-2 * kz3.imag * z)
#     exp = 1.0
#     exp = 1/(2 * kz3.imag) ## for testing, when doing integral over Z
#     return preFac * Teff * exp

def Wave_TwoLayer_bulkAbs(wlen, thick, n1, n2, n3, theta1, theta2, theta3):
    Teff = Wave_twoLayer_Teff(wlen, thick, n1, n2, n3, theta1, theta2, theta3)
    kz3 = kz(wlen, n3, theta3)
    preFac = np.pi * 2 * n3.real * n3.imag / (wlen * n1.real * np.cos(theta1))
    
    return preFac * Teff / (kz3.imag)  # this is the bulk absorption, not the power at interface
    
def Wave_TwoLayerReflection_total(wlen, thick, n1, n2, n3, theta1, theta2, theta3):
    # equations from this reference
    # https://nano-optics.ac.nz/planar/articles/fresnel.html
    # 1-->2
    r12_s, r12_p, t12_s, t12_p = Fresnel_Coefficients(n1, n2, theta1, theta2)
    # 2-->3
    r23_s, r23_p, t23_s, t23_p = Fresnel_Coefficients(n2, n3, theta2, theta3)
    # 2-->1
    r21_s, r21_p, t21_s, t21_p = Fresnel_Coefficients(n2, n1, theta2, theta1)
    
    delta = twolayer_delta(wlen=wlen, n2=n2, d=thick, theta=theta2)

    # plus in denominator from webpage
    # - in demoninator from nEXO paper
    r_p = (r12_p + r23_p * np.exp(1j*delta)) / (1 + r12_p * r23_p * np.exp(1j*delta))
    r_s = (r12_s + r23_s * np.exp(1j*delta)) / (1 + r12_s * r23_s * np.exp(1j*delta))

    R_p = r_p * np.conjugate(r_p)
    R_s = r_s * np.conjugate(r_s)
    return (R_p + R_s)/2

# wrong equation. included here to test it. taken from nexo paper but has a minus sign off
def Wave_TwoLayerReflection_total_minus(wlen, thick, n1, n2, n3, theta1, theta2, theta3):
    # equations from this reference
    # https://nano-optics.ac.nz/planar/articles/fresnel.html
    # 1-->2
    r12_s, r12_p, t12_s, t12_p = Fresnel_Coefficients(n1, n2, theta1, theta2)
    # 2-->3
    r23_s, r23_p, t23_s, t23_p = Fresnel_Coefficients(n2, n3, theta2, theta3)
    # 2-->1
    r21_s, r21_p, t21_s, t21_p = Fresnel_Coefficients(n2, n1, theta2, theta1)
    
    delta = twolayer_delta(wlen=wlen, n2=n2, d=thick, theta=theta2)

    # plus in denominator from webpage
    # - in demoninator from nEXO paper
    r_p = (r12_p + r23_p * np.exp(1j*delta)) / (1 - r12_p * r23_p * np.exp(1j*delta))
    r_s = (r12_s + r23_s * np.exp(1j*delta)) / (1 - r12_s * r23_s * np.exp(1j*delta))

    R_p = r_p * np.conjugate(r_p)
    R_s = r_s * np.conjugate(r_s)
    return (R_p + R_s)/2