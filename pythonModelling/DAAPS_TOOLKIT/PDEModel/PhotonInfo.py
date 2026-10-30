import numpy as np
import pandas as pd
import scipy as sp
import math

#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import ticker, cm

## our custom imports
import MathEquations as ME
import OpticsDataLoader as ODL
import PNJunction as PNJ
import TemperatureCorrections as TC

rad2deg = 180/np.pi


## ----- Photon Information Class ----- ##
class PhotonInfo:
  availableMediums = ["vacuum", "air", "LXe", "lxe", "LAr", "lar", "cf4"]
  
  # takes theta1 in degrees! Then converts all to radians
  def __init__(self, wlen = 420, theta_1 = 0, dataHelper = None, mediumStr = "vacuum", materialTemp = None):
    # print("new photon")
    if dataHelper is None:
      dataHelper = ODL.OpticsData()
    
    # input photon information
    self.wlen = wlen
    self.theta_1 = theta_1 / rad2deg
      
    self.theta_2 = None
    self.theta_3 = None
    
    # Reflection/Transmission at each boundary (no geometric bouncebacks here)
    # _ray is for pure 
    self.Ref_12_ray = None
    self.T_12_ray = None
    self.Ref_23_ray = None
    self.T_23_ray = None
    
    # for reverse going photons and wave optics
    self.Ref_21_ray = None
    
    # ignore window for now, but window should be same as boundary 1<-->2
    # note that Ref for a photon going 1-->2 is same as Transmission going 2-->1
    
    # building material information
    self.MediumStr = mediumStr
    self.nVacuum = None
    
    self.nSi = None
    self.nOxide = None
    
    self.muSi = None
    self.muOxide = None
    
    self.mMaterialTemp = materialTemp
    self.mReferenceTemp = 300
    
    ## sanity/safety checks ##
    if (np.array(self.theta_1) > np.pi/2).any():
      print("Warning: Theta1 is greater than pi/2.")
    
    # photon info, dynamics and kinematics
    self.BuildData(helper = dataHelper, Temp = self.mMaterialTemp)
    self.BuildThetas()
    self.BuildPureRT()
    
    # print("len muSi = ", len(self.muSi))
    return
  
  ## ----- End of __init__ ----- ##
  def BuildData(self, helper = None, Temp = None):
    if helper is None:
      helper = ODL.OpticsData()
    
    if self.MediumStr == "vacuum" or self.MediumStr == "air":
      self.nVacuum = 1.00 + 0j
    elif self.MediumStr == "LXe" or self.MediumStr == "lxe":
      helper.Load_LXe_nData()
      self.nVacuum = helper.LXeIndex(self.wlen)
    elif self.MediumStr == "LAr" or self.MediumStr == "lar":
      self.nVacuum = helper.LArIndex(self.wlen)
    elif self.MediumStr == "cf4":
      self.nVacuum = helper.CF4Index(self.wlen)
    else:
      raise ValueError(f"Medium {self.MediumStr} not recognized.")
    
    self.nSi = helper.SiIndex(self.wlen)
    self.nOxide = helper.OxideIndex(self.wlen)
    
    self.muSi = ME.Mu_fromIndex(wl = self.wlen, index = self.nSi)
    self.muOxide = ME.Mu_fromIndex(wl = self.wlen, index = self.nOxide)
    
    if Temp is not None:
      self.mMaterialTemp = Temp
      self.TempCorrection_absorption(Temp, referenceTemp=300)
    return
  
  def TempCorrection_absorption(self, Temp, referenceTemp = 300):
    alphaChange = TC.RelativeAlphaChange(Tdata = referenceTemp, Tnew = Temp, wlen = self.wlen)
    if min(self.wlen) < 300:
      print("temperature scaling not tested/validated below 300nm, be careful")
    if Temp < 150:
      print("PDE becomes affected by charge carrier freezeout below 150k, photo absorption alone may not be enough.")
      
    self.muSi *= alphaChange
    return
  
  def BuildThetas(self):
    self.theta_2 = ME.Snells_Law(self.nVacuum, self.nOxide, self.theta_1)
    self.theta_3 = ME.Snells_Law(self.nOxide, self.nSi, self.theta_2)
    # self.theta_3 = np.real(self.theta_3) # xx Testing
    # stheta2 = np.sin(self.theta_1) * self.nVacuum / self.nOxide
    # self.theta_2 = np.arcsin(stheta2) if (np.abs(stheta2) <= 1).any() else np.arcsin(1j * np.sqrt(np.abs(1 - stheta2**2)))
    # stheta3 = np.sin(self.theta_2) * self.nOxide / self.nSi
    # self.theta_3 = np.arcsin(stheta3) if (np.abs(stheta2) <= 1).any() else np.arcsin(1j * np.sqrt(np.abs(1 - stheta3**2)))
    return
  
  def BuildPureRT(self):
    self.Ref_12_ray, self.T_12_ray = ME.Fresnel_RefAndTrans(n1=self.nVacuum, n2=self.nOxide, theta1=self.theta_1, theta2=self.theta_2)
    self.Ref_23_ray, self.T_23_ray = ME.Fresnel_RefAndTrans(n1=self.nOxide, n2=self.nSi, theta1=self.theta_2, theta2=self.theta_3)
    
    # reverse photons coming up from Si
    self.Ref_21_ray, _, = ME.Fresnel_RefAndTrans(n1=self.nOxide, n2=self.nVacuum, theta1=self.theta_2, theta2=self.theta_1)
    return
  
  # def BuildFresnelFieldCoeffs(self):
  #     # 1-->2
  #     r12_s, r12_p, t12_s, t12_p = ME.Fresnel_Coefficients(self.nVacuum, self.nOxide, self.theta_1, self.theta_2)
  #     # 2-->3
  #     r23_s, r23_p, t23_s, t23_p = ME.Fresnel_Coefficients(self.nOxide, self.nSi, self.theta_2, self.theta_3)
  
  # redundant but nice that we can call from PNinfo or PhotonClass
  def GetWValues_withPN(self, PNInfo):
    Wp = PNInfo.calcWp(self)
    Wn = PNInfo.calcWn(self)
    return Wp, Wn
  
  def criticalAngle(self, n1, n2):
    # return np.arcsin(n2/n1)
    raise ArithmeticError("Total internal reflection not possible for photons entering Silicon!")
    return
  
  ## -- more intense math that relates more directly to PDE -- ##
  
  def OxideLiving(self, thick_oxide = None):
    # how much does light survive past silicon oxide layer
    # Only to be used for classical calculations
    if self.muOxide.all() == 0 or thick_oxide == 0:
      return 1
    elif thick_oxide is None:
      raise ValueError("Need to know the thickness of the passivation layer.")
    elif isinstance(thick_oxide, PNJ.PNJunction):
      thick_oxide = thick_oxide.tOxide
    oxideAtten = ME.Abs_singleLayer(d=thick_oxide, mu=self.muOxide, theta = self.theta_2)
    return oxideAtten

  def Pure_T13(self, thick_oxide = 0):
    OxideLiving = self.OxideLiving(thick_oxide)
    return self.T_12_ray * self.T_23_ray * OxideLiving
  
  def Ray_T13(self, thick_oxide = 0): # includes geometric reflections
    OxideLiving = self.OxideLiving(thick_oxide)
    R_eff = self.Ref_23_ray * self.Ref_21_ray * OxideLiving**2
    T_eff = self.T_12_ray * self.T_23_ray * OxideLiving
    return ME.Ray_TwoLayerTransmission(Ttotal = T_eff, Refftotal = R_eff)
  
  def Wave_TotalTransmission(self, thick_oxide=0):
    if thick_oxide is None:
      raise ValueError("Need to know the thickness of the passivation layer.")
    elif isinstance(thick_oxide, PNJ.PNJunction):
      thick_oxide = thick_oxide.tOxide
    fullT13 = ME.Wave_TwoLayerTransmission_total(self.wlen, thick_oxide, self.nVacuum, self.nOxide, self.nSi, self.theta_1, self.theta_2, self.theta_3)
    # fullR13 = ME.Wave_TwoLayerReflection_total(self.wlen, thick_oxide, self.nVacuum, self.nOxide, self.nSi, self.theta_1, self.theta_2, self.theta_3)
    return fullT13
  
  def Wave_TwoLayer_bulkAbs(self, thick_oxide=0):
    fullA = ME.Wave_TwoLayer_bulkAbs(self.wlen, thick_oxide, self.nVacuum, self.nOxide, self.nSi, self.theta_1, self.theta_2, self.theta_3)
    return fullA
  
  def Wave_TotalReflection(self, thick_oxide=0):
    if thick_oxide is None:
      raise ValueError("Need to know the thickness of the passivation layer.")
    elif isinstance(thick_oxide, PNJ.PNJunction):
      thick_oxide = thick_oxide.tOxide
    fullR13 = ME.Wave_TwoLayerReflection_total(self.wlen, thick_oxide, self.nVacuum, self.nOxide, self.nSi, self.theta_1, self.theta_2, self.theta_3)
    return fullR13
  
  def Wave_TotalReflection_polarized(self, thick_oxide=0):
    ## using ME.Wave_TLR_polarized(wlen, theta1, thick, n1, n2, n3)
    fullR13_s, fullR13_p = ME.Wave_TLR_polarized(self.wlen, self.theta_1, thick_oxide, self.nVacuum, self.nOxide, self.nSi)
    return fullR13_s, fullR13_p
  
  # this equation is wrong, but it's a good starting point
  def Wave_TotalReflection_minus(self, thick_oxide=0):
    if thick_oxide is None:
      raise ValueError("Need to know the thickness of the passivation layer.")
    elif isinstance(thick_oxide, PNJ.PNJunction):
      thick_oxide = thick_oxide.tOxide
    fullR13 = ME.Wave_TwoLayerReflection_total_minus(self.wlen, thick_oxide, self.nVacuum, self.nOxide, self.nSi, self.theta_1, self.theta_2, self.theta_3)
    return fullR13
  
  ######### Quartz Window stuff #########
  ## absorption in quartz window available via option
  ## this is the geometric ray equations for the quartz window
  ## See derivation in notes (somewhere lol)
  def T_R_quartzWindow(self, absorption=True):
    # here n1, n2 are functions of wavelength
    
    # 1 is air, 2 is quartz
    R12, T12 = ME.Fresnel_RefAndTrans(n1=self.nVacuum, n2=self.nOxide, theta1=self.theta_1, theta2=self.theta_2)
    R21, T21 = ME.Fresnel_RefAndTrans(n1=self.nOxide, n2=self.nVacuum, theta1=self.theta_2, theta2=self.theta_1)
    
    living = 1.0
    if absorption==True:
      living = ME.Abs_singleLayer(d=0.5e6, mu=self.muOxide, theta = self.theta_2) # 5mm is the thickness of the quartz window, in nm
      # print(self.muOxide)
      # print("living fraction:", living)
    A= living
    
    Rquartz = R12 + T12*R21*T21*A/(1-(R21*A)**2) ## correct for ray optics
    Tquartz = T12*T21*A/(1-(R21*A)**2) # correct for ray optics
    return Tquartz, Rquartz
  
  def T_total_windowAndOxide(self, thick_oxide = None, absorption=True):
    if thick_oxide is None:
      raise ValueError("Need to know the thickness of the passivation layer.")
    elif isinstance(thick_oxide, PNJ.PNJunction):
      thick_oxide = thick_oxide.tOxide
    
    # quartz window
    Tq, Rq = self.T_R_quartzWindow(absorption = absorption)
    # Tq *= self.OxideLiving(5e6) # 5mm is the thickness of the quartz window
    # print(self.OxideLiving(5e6)) # 5mm is the thickness of the quartz window
    Tsi = self.Wave_TotalTransmission(thick_oxide)
    Rsi = self.Wave_TotalReflection(thick_oxide)
    Ttotal = Tq * Tsi /(1 - Rq*Rsi)
    return Ttotal
  
  def R_total_windowAndOxide(self, thick_oxide = None, absorption = False):
    if thick_oxide is None:
      raise ValueError("Need to know the thickness of the passivation layer.")
    elif isinstance(thick_oxide, PNJ.PNJunction):
      thick_oxide = thick_oxide.tOxide
    
    # quartz window
    Tq, Rq = self.T_R_quartzWindow(absorption=absorption)
    # Tq *= self.OxideLiving(5e6) # 5mm is the thickness of the quartz window
    Tsi = self.Wave_TotalTransmission(thick_oxide)
    Rsi = self.Wave_TotalReflection(thick_oxide)
    Rtotal = Rq + Tq * Rsi * Tq / (1 - Rq*Rsi)
    return Rtotal
  
  ## ----- Printing Methods ----- ##
  def PrintPhotonInfo(self):
    info = [
        f"Theta1: {self.theta_1 * rad2deg} degrees",
        f"Theta2: {self.theta_2 * rad2deg} degrees",
        f"Theta3: {self.theta_3 * rad2deg} degrees"
        ]
    if self.Ref_12_ray is not None:
      info.extend([
        f"Pure R, T 1<->2: {self.Ref_12_ray}{self.T_12_ray}",
        f"Pure R, T 2<->3: {self.Ref_23_ray}, {self.T_23_ray}",
        ])
    
    print(f"-- Photon Information: {self.wlen} nm --")
    for st in info: print(st)
    print("")
    return info
  
  def PrintMaterialInfo(self):
    info = [
      f"nSi: {self.nSi}",
      f"nOxide: {self.nOxide}",
      f"muSi: {self.muSi}",
      f"muOxide: {self.muOxide}"
      ]
    
    print(f"-- Material Information: {self.wlen} nm --")
    for st in info: print(st)
    # print("")
    return info
## ----- End of PhotonInfo Class ----- ##