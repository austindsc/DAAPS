import os
import numpy as np
from scipy.interpolate import interp1d
import MathEquations as ME

import sys
sys.path.append('../')
from DAAPS_TOOLKIT import GenPaths as GP

## get the relative path to the data so we can run from anywhere
script_dir = os.path.dirname(os.path.realpath(__file__))

# from DAAPS_TOOLKIT.GenPaths import *
## alternatively can use "PROJECTPATH" from GenPaths

# print("OpticsDataLoader running from", script_dir)  

"""
This is a singleton class for loading and storing optical data for silicon and silicon oxide.
It uses the singleton pattern to ensure that the data is only loaded once, and that the same data is used throughout the program.

if you build a new instance with different data paths, it WILL make a new object. so this is a quasi-singleton class (ha ha ha)

you can set the spline type you want (cubic, linear, etc) by setting the 'interp_Type' variable in the __new__ method.
"""

"""
usage:
from PDEhelpers.OpticsDataLoader import *
od = OpticsData(optionalPathToSiDataPath, optionalPathTooxideData)
od.SiIndex(500)
od.OxideIndex(800)

DEFAULT DATA IS ASPNES IN MIDDLE wlen region ?? not anymore 
silicon_crystalline_lithoAndSchinke_merged.txt is other


be sure to call 
OpticsData.ResetInstance() if you have loaded in data previously and want to load in new default data. 
"""

"""
Interpolation:
we do linear interpolation, to avoid curves (especially curves going to negative n, k values)
We extrapolate past the bounds using the last value (flat past bounds)
"""
## random functions that don't get called all the time ##
def Load_Ham_data(fname = "HamVUV4_6050CQ-02_fromDatasheet.txt"):
  path = os.path.join(GP.VERA_DATA_PATH, "..", "datasheets", fname)
  absPath = os.path.abspath(path)
  x, y = np.genfromtxt(absPath, skip_header= 1, delimiter=",", unpack=True)
  return x, y

class OpticsData:
  _instance = None

  def __new__(cls, siDataPath = script_dir+"/../opticsData/GoodMerged/silicon_merged_AspnesMiddle.txt",
        oxideDataPath = script_dir+"/../opticsData/GoodMerged/silicon_dioxide-merged-clear.txt", 
        
        ## put other arguments here
        
        verbose = 0):
    if cls._instance is None or siDataPath != cls._instance.mPathSi or oxideDataPath != cls._instance.mPathSiO2:
      cls._instance = super().__new__(cls) 
      ###------- init happens here -------###
      # using our fake name for 'self' haha
      slf = cls._instance
      slf.mPathSi = siDataPath
      slf.mPathSiO2 = oxideDataPath
      slf.interp_Type = "linear" # cubic is dangerous, linear is safe
      slf.verbose = verbose
      
      ## for tracking data wavelengths
      slf.dataWlenSi = None
      slf.dataWlenOxide = None
      
      ## LXe stuff
      slf.loadAllData() # load data and create interpolate objects
      
    return cls._instance
    # don't define __init__ because it will be called every time __new__ is called



  # ## more singleton stuff
  # @classmethod
  # def get_instance(cls):
  #     if not cls._instance:
  #         cls._instance = cls()
  #     return cls._instance
  
  @classmethod
  def ResetInstance(self):
    self._instance = None
    return

  def loadAllData(self):
    self.loadData_Silicon()
    self.loadData_SiliconOxide()
    return

  def loadData_Silicon(self):
    try:
      data_nk = np.genfromtxt(self.mPathSi, delimiter=',', skip_header=1, usecols=(0, 1, 2)) # Load data from files
      x_n, n, k = data_nk[:, 0], data_nk[:, 1], data_nk[:, 2] # Extract x, n, k values

      print("-----------------------------")
      if np.mean(x_n) < 100:
        print("Silicon data is in um, converting to nm")
        x_n *= 1e3
      
      self.dataWlenSi = x_n
      self.nSi_interp = interp1d(x_n, n, kind = self.interp_Type, bounds_error=False, fill_value=(n[0], n[-1])) # Create interpolated functions for n(x) and k(x)
      self.kSi_interp = interp1d(x_n, k, kind = self.interp_Type, bounds_error=False, fill_value=(k[0], k[-1])) # fill_value='extrapolate' is used to allow the function to be called outside the range of the data
      
      print(f"Silicon data from\n   {self.mPathSi}\nloaded and interpolated successfully.")
      if self.verbose:
        print("x_n\t\tn\t\tk")
        [print(f"{x_n[i]}\t{n[i]}\t{k[i]}") for i in range(len(x_n))]
      print("-----------------------------")
      
    except Exception as e:
      print("Error loading silicon data:", str(e))

  def loadData_SiliconOxide(self):
    try:
      data_nk = np.genfromtxt(self.mPathSiO2, delimiter=',', skip_header=1, usecols=(0, 1, 2))
      x_n, n, k = data_nk[:, 0], data_nk[:, 1], data_nk[:, 2]
      
      print("-----------------------------")
      if np.mean(x_n) < 100:
        print("Oxide data is in um, converting to nm")
        x_n *= 1e3
      
      self.dataWlenOxide = x_n
      self.nSiO2_interp = interp1d(x_n, n, kind = self.interp_Type, bounds_error=False, fill_value=(n[0], n[-1]))
      self.kSiO2_interp = interp1d(x_n, k, kind = self.interp_Type, bounds_error=False, fill_value=(k[0], k[-1]))

      print(f"Oxide data from\n   {self.mPathSiO2}\nloaded and interpolated successfully.")
      if self.verbose:
        print("x_n\t\tn\t\tk")
        [print(f"{x_n[i]}\t{n[i]}\t{k[i]}") for i in range(len(x_n))]
      print("-----------------------------")
      
    except Exception as e:
      print("Error loading silicon oxide data:", str(e))
  ###########------- End of functions that get called at Init -------########### 
  
  ## random functions that don't get called all the time ##
  # this is now a wrapper function
  def Load_Ham_data(self, fname = "HamVUV4_6050CQ-02_fromDatasheet.txt"):
    return Load_Ham_data(fname)
    
  def Load_LXe_nData(self, fname="lxe_n_merged.txt"):
    # path = os.path.join(script_dir, "../", "opticsData", fname)
    path = os.path.join(GP.OPTICS_DATA_PATH, fname)
    absPath = os.path.abspath(path)
    wlen, n = np.genfromtxt(absPath, skip_header= 1, delimiter="\t", unpack=True)
    self.nLXe_interp = interp1d(wlen, n, kind = self.interp_Type, bounds_error=False, fill_value=(n[0], n[-1]))
    return wlen, n 

  ###########------- Functions to be used repeatedly -------###########
  def LXeIndex(self, wl): ## no k values in LXe
    return self.nLXe_interp(wl) + 1j * 0

  ## we truncate these to within the range of the resonances
  def LArIndex(self, wl): # argon, taken from grace, munroe paper, 90K
    # return (1.26+0.23/(1-(0.1066/x)**2)+0.0023/(1-(0.9083/x)**2))**.5 + 1j * 0
    wluv = 106.6
    uvBuf = 1
    wlir = 908.3
    irBuf = 20
    lowVal = ME.SellmeierIndex(wluv + uvBuf, a0=1.26, Auv=0.23, Air=0.0023, wluv=wluv, wlir=wlir)
    highVal = ME.SellmeierIndex(wlir - irBuf, a0=1.26, Auv=0.23, Air=0.0023, wluv=wluv, wlir=wlir)
    
    # return np.piecewise(wl, 
    #     [wl < wluv + uvBuf, (wl >= wluv + uvBuf) & (wl <= wlir - irBuf), wl > wlir - irBuf],
    #     [lambda x: lowVal, 
    #      lambda x: ME.SellmeierIndex(x, a0=1.26, Auv=0.23, Air=0.0023, wluv=wluv, wlir=wlir), 
    #     lambda x: highVal])
    return ME.SellmeierIndex(wl, a0=1.26, Auv=0.23, Air=0.0023, wluv=wluv, wlir=wlir)
  
  def LXeIndex_Sellmeier(self, wl): # 178 K
    wluv = 146.9
    uvBuf = 1
    wlir = 827.0
    irBuf = 20
    lowVal = ME.SellmeierIndex(wluv + uvBuf, a0=1.4, Auv=0.4, Air=0.002, wluv=wluv, wlir=wlir)
    highVal = ME.SellmeierIndex(wlir - irBuf, a0=1.4, Auv=0.4, Air=0.002, wluv=wluv, wlir=wlir)
    # return np.piecewise(wl,
    #     [wl < wluv + uvBuf, (wl >= wluv + uvBuf) & (wl <= wlir - irBuf), wl > wlir - irBuf],
    #     [lambda x: lowVal,
    #      lambda x: ME.SellmeierIndex(x, a0=1.4, Auv=0.4, Air=0.002, wluv=wluv, wlir=wlir),
    #      lambda x: highVal])
    return ME.SellmeierIndex(wl, a0=1.4, Auv=0.4, Air=0.002, wluv=wluv, wlir=wlir)
  
  def CF4Index(self, wl): # cf4, taken from grace, munroe paper, 300K
    liq_cf4_photonE = np.array([7, 5, 2, 1]) # in eV
    liq_cf4_n = [1.24, 1.225, 1.21, 1.2039]
    liq_cf4_wl = 1239.8 / liq_cf4_photonE # convert to nm
    self.cf4Interp = interp1d(liq_cf4_wl, liq_cf4_n, kind = self.interp_Type, bounds_error=False, fill_value=(liq_cf4_n[0], liq_cf4_n[-1]))
    return self.cf4Interp(wl) + 1j * 0

  def SiIndex(self, wl):
    # Return the complex refractive index n + ik for a given x value
    return self.nSi_interp(wl) + 1j * self.kSi_interp(wl)
  
  def OxideIndex(self, wl):
    return self.nSiO2_interp(wl) + 1j * self.kSiO2_interp(wl)

  def Si_absCoeff(self, wl):
    return ME.Mu_fromK(wl = wl, k = self.kSi_interp(wl)) #Still unsure if this works given it can't find MathEquations
  
  def Oxide_absCoeff(self, wl):
    return ME.Mu_fromK(wl = wl, k = self.kSiO2_interp(wl))
