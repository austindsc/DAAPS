import numpy as np

# from 
import PDEModel.MathEquations as ME
import PDEModel.PhotonInfo as PI # this calls DataLoader when defining init
import PDEModel.PNJunction as PNJ

from scipy.stats import norm
from scipy.stats import triang
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import interp1d

from scipy.ndimage import convolve1d

def quadraticSolution(a, b, c):
    return (-b + np.sqrt(b**2 - 4*a*c))/(2*a), (-b - np.sqrt(b**2 - 4*a*c))/(2*a)

# this returns the relative (fraction) of the microstructure width
# relative to pixel size
def W_from_FF(FF):
    sol = quadraticSolution(1, -1, (1-FF)/4)
    return min(sol)

# this returns the absolute width of the microstructure
# this is the same as above multiplied by pitch
def ff_from_FF(FF, pitch):
    sol = quadraticSolution(1, -1, (1-FF)/4)
    return min(sol) * pitch

## simpler functions with objects
def StandardPDE(FF, T, Wp, Wn, Pe, Ph):
    return FF * T * (Wp*Pe + Wn*Ph)

def SimplePDE_fromObjects(pn, pi):
    Wp, Wn = pn.calcW_both(pi)
    T = pi.Wave_TotalTransmission(pn)
    T = 1 - pi.Wave_TotalReflection(pn)
    return StandardPDE(FF=pn.FF, T=T, Wp=Wp, Wn=Wn, Pe=pn.pe, Ph=pn.ph)
 
######### Functions for full PDE equations #########
# need to get the locals() dict, (and remove whatever your x value array is for recursion to only happen once)
def MiscFunction_smoothing(x, sigma, argDict, function, stepPerSigma = 5, sigmaPastEdges = 2.5):
    stepSize = sigma / stepPerSigma
    newXvals = np.arange(min(x)-sigmaPastEdges*sigma, max(x) + sigmaPastEdges*sigma + stepSize, stepSize)
    newFunctionCall = function(newXvals, **argDict)
    
    smoothedFunction = gaussian_filter1d(newFunctionCall, sigma=stepPerSigma)
    
    originalFunction = interp1d(newXvals, smoothedFunction, kind='linear', fill_value=0)
    return originalFunction(x)

## no interpolation on this one (using strategy from relative toxide)
## triangular only
## width here is FWHM for a 'even' triangle with shape factor 0.5. 
def create_centered_triangular_arrays_vected(xarr, width, num_points):
    stepSize = 2 * width / num_points
    xshift = np.arange(-width+stepSize/2, width+stepSize/2, stepSize)
    # yVals = triang.pdf(xshift, c=0.5, loc=-width, scale=2 * width)
    yVals = np.bartlett(num_points)
    
    # xVals = [np.arange(x - width, x + width + stepSize, stepSize) for x in xarr]
    xVals = xarr[:, np.newaxis] + xshift
    return np.array(xVals), yVals

## halfwidth here is half the triangle base, which is (probably) the FWHM
# this function creates a new x array for each x point, averaging over wlen for each x (wlen) data point. expensive and slow
def MiscFunction_smoothing2_triang(x, argDict, function, halfWidth, NperPoint = 16):
    # NperPoint - number of averages done to compute the new value
    # halfWidth - half width of the smoothing triangle
    allXvals, triWeight = create_centered_triangular_arrays_vected(x, halfWidth, NperPoint)

    outputXes = np.array([np.average(function(xVals, **argDict), weights=triWeight) for xVals in allXvals ])
    
    # vectorized way. need function to take 2d input(s)
    # functionValues = function(allXvals, **argDict)
    # outputXes = np.average(functionValues, axis=1, weights=triWeight)
    
    return outputXes

# we make sure the new x values extend below the original array by halfwidth,
# such that the extrapolation that convolve does at the edges doesn't matter
def FuncSmoothing_trian_interp1d(x, argDict, function, halfWidth, NperPoint = 11):
    if NperPoint % 2 == 0:
        print("NperPoint must be odd, adding 1")
        NperPoint += 1
    stepSize = 2 * halfWidth / NperPoint
    allXvals = np.arange(min(x) - halfWidth, max(x) + halfWidth + stepSize, stepSize)
    functionValues = function(allXvals, **argDict)
    
    kernel = np.bartlett(NperPoint)
    kernel = kernel / np.sum(kernel)
    
    # originalPoints = interp1d(allXvals, functionValues, kind='linear', fill_value=0)
    smoothed = convolve1d(functionValues, kernel, mode='reflect')
    
    return interp1d(allXvals, smoothed, kind='cubic', fill_value=0)(x)

def FuncAveraging_triangle(x, argDict, function, halfWidth, NperPoint = 11):
    stepSize = 2 * halfWidth / NperPoint
    xSpan = np.arange(x-halfWidth, x + halfWidth + stepSize, stepSize)
    newFuncValues = np.array([function(xVal, **argDict) for xVal in xSpan])
    weight = triang.pdf(xSpan, c=0.5, loc=x-halfWidth, scale=2*halfWidth)
    return np.average(newFuncValues, weights=weight, axis=0)

############### Our Standard PDE Calculations ###############

## Cases:
## A: wlen > 160, normal
## B: wlen <= 160, normal
## A, B - normal incidence: Field equations work

## C: wlen > 160, oblique
## field equations 'fail', but Poynting gives R+T = 1.0
## In this case, T = 1-R

## D: wlen <= 160, oblique
## in this case, T+R < 1.0 by field or Poynting, as there is absorption in middle layer
## for T = 1 - R to work, would need to account for absorption in quartz. In which case:
## T = 1 - A' - R where A' is absorption in thin film

def LocalTransmissionWrapper(photon, tOxide=17, window=False, case="A"):
    if window:
        if case == "A" or case == "B": # field calculation because normal incidence
            Ttotal = photon.T_total_windowAndOxide(tOxide, absorption=True)
        elif case == "C": # 'poynting' sum because incorrectly R+T < 1.0
            Ttotal = 1 - photon.R_total_windowAndOxide(tOxide, absorption=True)
        elif case == "D":
            print("option D not complete, must account for A in quartz. Using field equations for now")
            Ttotal = photon.T_total_windowAndOxide(tOxide, absorption=True)
    else:
        if case == "A" or case == "B":
            Ttotal = photon.Wave_TotalTransmission(tOxide)
        elif case == "C": # 'poynting' sum because incorrectly R+T < 1.0
            Ttotal = 1 - photon.Wave_TotalReflection(tOxide)
        elif case == "D":
            print("option D not complete, must account for A in quartz. Using field equations for now")
            Ttotal = photon.Wave_TotalTransmission(tOxide)
    return Ttotal

def LocalReflectivityWrapper(photon, tOxide, window):
    if window:
        Rtotal = photon.R_total_windowAndOxide(tOxide, absorption=True)
    else:
        Rtotal = photon.Wave_TotalReflection(tOxide)
    return Rtotal

# this one is for a 1d sigma smoothing in wavelength (smoothSigma)
# x axis as wavelength
def PDE_allArgs_wlen(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, pe, ph, Temperature = None, smoothSigma=None, quantum = False, dh = None, window=False, case = "A"):
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        
    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window, case)

    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp, Wn = pn.calcW_both(photon)
    
    if quantum:
        pe = ME.QuantumPe(wlen, pe, 300, 0.75)

    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        # return MiscFunction_smoothing(x=wlen, sigma=smoothSigma, argDict=theseArgsDict, function=PDE_allArgs, stepPerSigma=5, sigmaPastEdges=2.5)
        # return MiscFunction_smoothing2_triang(x=wlen, argDict=theseArgsDict, function=PDE_allArgs, halfWidth=smoothSigma, NperPoint=51) # slow but no interp1d?
        return FuncSmoothing_trian_interp1d(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_wlen, halfWidth=smoothSigma, NperPoint=13)
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)
    return PDE

# same as above, but Pe and Ph are parametrized by Ve, Vh instead of letting them be free (floating)
# really smooth sigma should be input as the FWHM of the assumed triangle smoothing function
def PDE_allArgs_wlen_parametrized(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, voltage, Ve, Vh, Ae=1.0, Ah=1.0, Temperature = None, smoothSigma=None, quantum = False, dh = None, window=False, case = "A"):
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        
    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window, case)

    pe = ME.Pcar_parametrized(voltage, Ve, Ae)
    ph = ME.Pcar_parametrized(voltage, Vh, Ah)
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp, Wn = pn.calcW_both(photon)
    
    if quantum:
        pe = ME.QuantumPe(wlen, pe, 300, 0.75)
    
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)
    
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        # return MiscFunction_smoothing(x=wlen, sigma=smoothSigma, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized, stepPerSigma=5, sigmaPastEdges=2.5)   
        # return MiscFunction_smoothing2_triang(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized, halfWidth=smoothSigma, NperPoint=8)
        return FuncSmoothing_trian_interp1d(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized, halfWidth=smoothSigma, NperPoint=11)
            
    return PDE
  
def PDE_allArgs_wlen_parametrized_forMedium(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, voltage, Ve, Vh, Ae=1.0, Ah=1.0, Temperature = None, smoothSigma=None, quantum = False, dh = None, window=False, medium="vacuum"):
  if smoothSigma is not None: # expensive but simple/works/robust
      theseArgsDict = locals()
      del theseArgsDict["wlen"]
      
  photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature, mediumStr=medium)
  Ttotal = LocalTransmissionWrapper(photon, tOxide, window)

  pe = ME.Pcar_parametrized(voltage, Ve, Ae)
  ph = ME.Pcar_parametrized(voltage, Vh, Ah)
  pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
  Wp, Wn = pn.calcW_both(photon)
  
  if quantum:
      pe = ME.QuantumPe(wlen, pe, 300, 0.75)
  
  PDE = FF * Ttotal * (Wp*pe + Wn*ph)
  
  if smoothSigma is not None: # expensive but simple/works/robust
      theseArgsDict["smoothSigma"] = None
      return FuncSmoothing_trian_interp1d(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized_forMedium, halfWidth=smoothSigma, NperPoint=11)
          
  return PDE

# this one is for a 1d sigma smoothing in wavelength, with x axis as angular dependance
def PDE_allArgs_Angle(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, pe, ph, smoothSigma=None, quantum = False, dh = None, case = "C"):
    # this is very different from the wavelength smoothing, we are smoothing in a different dimension
    if smoothSigma is not None:
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        theseArgsDict["smoothSigma"] = None

    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window=False, case=case)
    # Ttotal = photon.Wave_TotalTransmission(tOxide)
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp, Wn = pn.calcW_both(photon)
    
    if quantum:
        cutinWl = 300
        nyield = np.where(wlen < cutinWl, np.sqrt(cutinWl/wlen), 1)
        pe = 1-(1-pe)**nyield
    
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)
    if smoothSigma is not None and smoothSigma > 0:
        return FuncAveraging_triangle(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_Angle, halfWidth=smoothSigma, NperPoint=11) 
    return PDE
  
## made in Nov 2025, hopefully the most complete (other than lacking smoothing in either dimension)
def PDE_allArgs_Angle_parametrized_medium_shadowing(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, voltage, Ve, Vh, Ae=1.0, Ah=1.0,
                                                    smoothSigma=None, quantum = False, dh = None, window = None, Temperature = None,case = "C",
                                                    shadowingH = None, mediumStr="vacuum"):
    if smoothSigma is not None:
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        theseArgsDict["smoothSigma"] = None
    
    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature, mediumStr=mediumStr)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window=window, case=case).real
    pe = ME.Pcar_parametrized(voltage, Ve, Ae)
    ph = ME.Pcar_parametrized(voltage, Vh, Ah)
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp, Wn = pn.calcW_both(photon)
    
    ## doing shadowing
    if shadowingH is not None and shadowingH > 0:
        relative_W = W_from_FF(FF)
        FF_theta = FF - 1 * FF_theta_relative_shadow(theta = theta1, W = relative_W, H = shadowingH)
        # no factor of 2 here because we are NOT subtracting shadowing twice due to reflection (incoming and outgoing)
    else:
        FF_theta = FF
    
    ## doing quantum yield
    if quantum:
        pe = ME.QuantumPe(wlen, pe, 300, 0.75)
    
    PDE = FF_theta * Ttotal * (Wp*pe + Wn*ph)
  
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        return FuncSmoothing_trian_interp1d(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_Angle_parametrized_medium_shadowing, halfWidth=smoothSigma, NperPoint=11)
            
    return PDE.real
    

##### Optics equations #####
# x axis is wavelength
# smoothing convolutes with triangular kernel (along x axis)
def Reflectivity_wlen(wlen, theta1, tOxide, FF = 0.6, scaleDown =0.0, shadowingH = None, window = False, smoothSigma=None, dh=None, mediumStr="vacuum"):
    # to use later
    if smoothSigma is not None:
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        theseArgsDict["smoothSigma"] = None

    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, mediumStr=mediumStr)
    Rtotal = LocalReflectivityWrapper(photon, tOxide, window)

    if shadowingH is not None and shadowingH > 0:
        relative_W = W_from_FF(FF)
        FF_theta = FF - 2 * FF_theta_relative_shadow(theta = theta1, W = relative_W, H = shadowingH)
        ## extra factor of 2 here because we are subtracting shadowing twice due to reflection (incoming and outgoing)
    else:
        FF_theta = FF
    
    if smoothSigma is not None and smoothSigma > 0: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        return FuncSmoothing_trian_interp1d(x=wlen, argDict=theseArgsDict, function=Reflectivity_wlen, halfWidth = smoothSigma, NperPoint=7)

    return np.real(Rtotal * FF_theta * (1-scaleDown) )

# x axis is angle
# smoothing just has to average different curves together (not smooth along x axis)
def Reflectivity_angle(wlen, theta1, tOxide, FF = 0.6, scaleDown = 0.0, shadowingH = None, window = False, smoothSigma=None, dh=None, mediumStr="vacuum"):
    # to use later
    if smoothSigma is not None:
        theseArgsDict = locals()
        del theseArgsDict["wlen"]

    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, mediumStr=mediumStr)
    Rtotal = LocalReflectivityWrapper(photon, tOxide, window)
            
    if shadowingH is not None and shadowingH > 0:
        relative_W = W_from_FF(FF)
        FF_theta = FF - 2 * FF_theta_relative_shadow(theta = theta1, W = relative_W, H = shadowingH)
        # extra factor of 2 here because we are subtracting shadowing twice due to reflection (incoming and outgoing)
    else:
        FF_theta = FF
    
    Ref_total = np.real(Rtotal * FF_theta * (1-scaleDown))

    if smoothSigma is not None and smoothSigma > 0:
        theseArgsDict["smoothSigma"] = None
        return FuncAveraging_triangle(x=wlen, argDict=theseArgsDict, halfWidth = smoothSigma, function=Reflectivity_angle, NperPoint=11)
    return Ref_total

def Reflectivity_either(wlen, theta1, tOxide, FF = 0.6, scaleDown = 0.0, shadowingH = None, window = False, smoothSigma=None, dh=None, mediumStr="vacuum"):
    if hasattr(theta1, "__len__"):
        return Reflectivity_angle(wlen, theta1, tOxide, FF, scaleDown, shadowingH, window, smoothSigma, dh, mediumStr)
    elif hasattr(wlen, "__len__"):
        return Reflectivity_wlen(wlen, theta1, tOxide, FF, scaleDown, shadowingH, window, smoothSigma, dh, mediumStr)
    else:
        print("Invalid input for Reflectivity_either")
        return 0

##-------- Shadowing equations --------##
def l_loss(h, theta): # theta in degrees, h in absolute height
    theta = np.abs(theta * np.pi/180)
    return h * np.tan(theta)

def FF_theta_absolute(p, ff, h, theta):
    l = l_loss(h, theta)
    top1 = p - 2*ff
    top2 = p - 2*ff - l
    return (top1 * top2)/(p**2)

def FF_theta_relative(theta, W, H, offset = 0):
    # W defined as ff/p
    # H defined as h/p
    # offset is vertical scaling up or down (degenerate with W)
    theta = np.abs(theta * np.pi/180)
    
    fixed = 4*(W**2) - 4*W + 1
    var = (1 - 2*W)*H*np.tan(theta)
    
    return offset + fixed - var

# only the subtracting part of the equation above
# (for doing reflectivity where we subtract shadow twice)
def FF_theta_relative_shadow(theta, W, H):
    theta= np.abs(theta * np.pi/180)
    var = (1 - 2*W)*H*np.tan(theta)
    return var
# print("testing FF func:", FF_theta_relative(Wfbk, Hfbk, 0))

##----------- Relative equations -----------##
# set H=0 if you want zero angular dependance in the FF
# this only accounts for the optics
# ---- this is the function used for the relative PDE measurements (because we ignore internal absorption)
# x axis is angle, for a fixed input angle, smoothing is done in wavelength, with options between gaussian, triangular and box.
# @np.vectorize
def RelativeTrans_oxide(wlen, tOxide, theta1, W=0.05, H=0, smoothSigma = None, dataHelper = None, window=False, smoothStyle = "gauss", case = "C"):
    
    if smoothSigma is not None and smoothSigma > 0: # removing args for the single recursion call
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        theseArgsDict["smoothSigma"] = None
    
    if smoothStyle not in ["gauss", "triang", "box"]:
        print("Invalid smoothing style, using Gaussian")
        smoothStyle = "gauss"
    
    # normal incidence values (to normalize)
    photon0 = PI.PhotonInfo(wlen, theta_1=0, dataHelper=dataHelper)
    Tzero = LocalTransmissionWrapper(photon0, tOxide, window, case=case)
    FFZer0 = FF_theta_relative(theta=0.0, W=W, H=H)
    
    # angular dependent values
    photonAngle = PI.PhotonInfo(wlen, theta_1=theta1, dataHelper=dataHelper)
    TAngle = LocalTransmissionWrapper(photonAngle, tOxide, window, case=case)
    FFAngle = FF_theta_relative(theta=theta1, W=W, H=H)
    
    # smoothing
    if smoothSigma is not None and smoothSigma > 0:
        NperSigma = 10
        stepSize = 2*smoothSigma / NperSigma
        
        wlenSpan = np.arange(wlen-smoothSigma, wlen + smoothSigma + stepSize, stepSize)
        pdeVals = np.array([RelativeTrans_oxide(wl, **theseArgsDict) for wl in wlenSpan])
        weight = triang.pdf(wlenSpan, c=0.5, loc=wlen-smoothSigma, scale=2*smoothSigma) # triang by default
        
        if smoothStyle == "gauss":
            NsigmaPast = 2 # move two sigma past the edges of the data (x axis)
            wlenSpan = np.arange(wlen-NsigmaPast*smoothSigma, wlen + NsigmaPast*smoothSigma + stepSize, stepSize)
            weight = norm.pdf(wlenSpan, wlen, smoothSigma)
        elif smoothStyle == "box": # uniform
            weight = 1.0 * [(wlenSpan < wlen + smoothSigma/2 - stepSize/2) & (wlenSpan > wlen - smoothSigma/2)][0]
        return np.average(pdeVals, axis=0, weights=weight)
    
    return TAngle * FFAngle / (Tzero * FFZer0)


############## non-updated functions, or with something slightly wrong ##############

## just changing the signature so we can include the medium
## should be the exact same as PDE_allArgs_wlen_parametrized with one extra arg (cause im lazy and dont want to add arguments to all my fitting dicts)
def PDE_allArgs_wlen_parametrized_forMedium(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, voltage, Ve, Vh, Ae=1.0, Ah=1.0, Temperature = None, smoothSigma=None, quantum = False, dh = None, window=False, medium="vacuum"):
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        
    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature, mediumStr=medium)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window)

    pe = ME.Pcar_parametrized(voltage, Ve, Ae)
    ph = ME.Pcar_parametrized(voltage, Vh, Ah)
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp, Wn = pn.calcW_both(photon)
    
    if quantum:
        pe = ME.QuantumPe(wlen, pe, 300, 0.75)
    
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)
    
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        # return MiscFunction_smoothing(x=wlen, sigma=smoothSigma, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized, stepPerSigma=5, sigmaPastEdges=2.5)   
        # return MiscFunction_smoothing2_triang(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized, halfWidth=smoothSigma, NperPoint=8)
        return FuncSmoothing_trian_interp1d(x=wlen, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized_forMedium, halfWidth=smoothSigma, NperPoint=11)
            
    return PDE

## BadWrapper because I pass the photon directly
## here photon has temperature, and n, k of the active medium, and dh
## no smoothing either
def PDE_allArgs_wlen_parametrized_BADWRAPPER(photon, FF, tOxide, dp_star, Xpn, dw_star, voltage, Ve, Vh, Ae=1.0, Ah=1.0, quantum = False, window=False):
        
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window)
    pe = ME.Pcar_parametrized(voltage, Ve, Ae)
    ph = ME.Pcar_parametrized(voltage, Vh, Ah)
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp, Wn = pn.calcW_both(photon)
    
    if quantum:
        pe = ME.QuantumPe(photon.wlen, pe, 300, 0.75)
    
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)
 
    return PDE

# this one is for a 1d sigma smoothing in wavelength (smoothSigma), with x axis as wavelength
# with the halfBox diffusion 
## need to be updated. wrong, re smoothing or something else
def PDE_allArgs_halfBox(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, pe, ph, Pdiff, Temperature = None, smoothSigma=None, quantum = False, dh = None, window=False):
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window)
        
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    Wp = pn.calcHalfBox_Wp(photonInfo=photon, Pdiff=Pdiff)
    Wn = pn.calcWn(photon)
    
    if quantum:
        pe = ME.QuantumPe(wlen, pe, 300, 0.75)
    
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)

    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        return MiscFunction_smoothing(x=wlen, sigma=smoothSigma, argDict=theseArgsDict, function=PDE_allArgs_halfBox, stepPerSigma=5, sigmaPastEdges=2.5)   
        
    return PDE

## xx Testing, doing halfbox diffusion for Wp
## need to be updated. wrong, re smoothing or something else
def PDE_allArgs_wlen_parametrized_halfBox(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, voltage, Ve, Vh, Ae=1.0, Ah=1.0, Pdiff = 0.5, Temperature = None, smoothSigma=None, quantum = False, dh = None, window=False):

    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        
    photon = PI.PhotonInfo(wlen, theta1, dataHelper=dh, materialTemp=Temperature)
    Ttotal = LocalTransmissionWrapper(photon, tOxide, window)
    pe = ME.Pcar_parametrized(voltage, Ve, Ae)
    ph = ME.Pcar_parametrized(voltage, Vh, Ah)
    pn = PNJ.PNJunction(tOxide = tOxide, dp_star = dp_star, Xpn = Xpn, dw_star=dw_star, pe = pe, ph = ph)
    
    Wp = pn.calcHalfBox_Wp(photonInfo=photon, Pdiff=Pdiff)
    Wn = pn.calcWn(photon)
    
    if quantum:
        pe = ME.QuantumPe(wlen, pe, 300, 0.75)
    
    PDE = FF * Ttotal * (Wp*pe + Wn*ph)
    
    if smoothSigma is not None: # expensive but simple/works/robust
        theseArgsDict["smoothSigma"] = None
        return MiscFunction_smoothing(x=wlen, sigma=smoothSigma, argDict=theseArgsDict, function=PDE_allArgs_wlen_parametrized_halfBox, stepPerSigma=5, sigmaPastEdges=2.5)   
            
    return PDE

## --- relative --- ##
## relative PDE barely changes with angle so this function only was used for validating
## the statement that relative PDE ~ relative Transmission
def RelativePDE_angle(wlen, theta1, FF, tOxide, dp_star, Xpn, dw_star, pe, ph, smoothSigma=None, quantum = False, dh = None, window = False):
    
    if smoothSigma is not None and smoothSigma > 0:
        theseArgsDict = locals()
        del theseArgsDict["wlen"]
        theseArgsDict["smoothSigma"] = None
        
    zeroVal = PDE_allArgs_wlen(wlen, theta1=0, FF=FF, tOxide=tOxide, dp_star=dp_star, Xpn=Xpn, dw_star=dw_star, pe=pe, ph=ph, smoothSigma=None, quantum=quantum, dh=dh, window=window)
    angleVal = PDE_allArgs_wlen(wlen, theta1=theta1, FF=FF, tOxide=tOxide, dp_star=dp_star, Xpn=Xpn, dw_star=dw_star, pe=pe, ph=ph, smoothSigma=None, quantum=quantum, dh=dh, window=window)
    
    # gaussian smoothing
    if smoothSigma is not None and smoothSigma > 0:
        NperSigma = 10
        stepSize = smoothSigma / NperSigma
        NsigmaPast = 2
        
        wlenSpan = np.arange(wlen-NsigmaPast*smoothSigma, wlen + NsigmaPast*smoothSigma + stepSize, stepSize)
        pdeVals = np.array([RelativePDE_angle(wl, **theseArgsDict) for wl in wlenSpan])
        weight = norm.pdf(wlenSpan, wlen, smoothSigma)
        # weight = triang.pdf(wlenSpan, c=0.5, loc=wlen-smoothSigma/2, scale=smoothSigma)
        return np.average(pdeVals, axis=0, weights=weight)
    return angleVal / zeroVal