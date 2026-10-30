import numpy as np
import pandas as pd
import scipy as sp
import math

import inspect

#%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import ticker, cm

## our custom imports
import MathEquations as ME
import PhotonInfo as PC

rad2deg = 180/np.pi

## ----- PN Junction Information Class ----- ##

class PNJunction:
    def __init__(self, dp_star = 4, Xpn = 2e3, dw_star=5e3, dp = None, dw = None, pe = 1.0, ph = 1.0, tOxide = 110, FF = 0.60):
        # the numbers are all in nm
        # these things are measured wrt top of Silicon (NOT from top of passivation layer)
        self.dp_star = dp_star  # effective dept of P region (diffusion)
        self.Xpn = Xpn          # middle of the junction
        self.dw_star = dw_star  # effective depth of n region (diffusion)
        
        ## don't actually need these for now, not used
        # self.dw = dw            # depth of the n region
        # self.dp = dp            # depth of the p region

        # junction geometry
        self.pe = pe            # probability electron causes avalanche
        self.ph = ph            # probability hole causes avalanche
        
        # passivation layer
        self.tOxide = tOxide      # thickness of passivation layer (SiO2 on surface)
        self.FF = FF          # fill factor for VUV4
        self.fitNames = ['dp_star', 'Xpn', 'dw_star', 'tOxide', 'pe', 'ph']
        
        # sanity check on geometry
        # if self.dp_star > self.dp:
        #     print("Warning: Effective depth of P region is deeper than actual depth.")
        # if self.dw_star < self.dw:
        #     print("Warning: Effective depth of N region is shallower than actual depth.")  
        if self.dp_star > self.dw_star:
            print("Warning: Effective depth of P region is deeper than N region.")
        if self.Xpn > self.dw_star or self.Xpn < self.dp_star:
            print("Warning: Middle of junction is not between top of P or bottom of N")

        # some voltage dependance?
        # self.vov_slope = 0.69 # unused for now
        
        return
    

    def fitParList(self):
        if not all(hasattr(self, name) for name in names):
            print("Not all elements in names are attributes of self.")
            return None
        return self.fitNames
    def fitParDict(self):
        return {name: getattr(self, name) for name in self.fitParList()}
    
    # organizing the thicknesses as properties because they are derived from above values
    @property
    def tdead(self):
        return self.dp_star
    @property
    def tp(self):
        return self.Xpn - self.dp_star
    @property
    def tn(self):
        return self.dw_star - self.Xpn
    
    ##### Reading/Helper Methods #####
    def printAllInfo(self, doPrint = False):
        info = ["Device Information:",
                f"Fill Factor: {self.FF}",
                f"Oxide Thickness: {self.tOxide} nm",
                f"Dp* (start of P): {self.dp_star:0.2f} nm",
                f"Xpn (middle): {self.Xpn/1000:0.2f} um",
                f"Dw* (end of N): {self.dw_star/1000:0.2f} um",
                f"P Thickness: {self.tp/1000:0.1f} um",
                f"N Thickness: {self.tn/1000:0.1f} um"]
        
        if doPrint:
            print(f"-- PN Information: --")
            for st in info: print(st)
            print("")
        return info
    
    def printPNInfo(self, doPrint = False):
        info = [
                fr"$dp*=$ {self.dp_star:0.1f} nm",
                r"$X_{PN}  =$" +  f"{self.Xpn/1000:0.1f} um",
                fr"$dw^* =$ {self.dw_star/1000:0.1f} um"
                ]
        
        if doPrint:
            print(f"-- PN Information: --")
            for st in info: print(st)
            print("")
        return info

    
    ## option 'pure path alone'
    def calcWp(self, photonInfo):
        # print("theta 3 = ", photonInfo.theta_3)
        return ME.Abs_layers(d_above=self.tdead, d_region=self.tp, mu=photonInfo.muSi, theta = photonInfo.theta_3)
    
    def calcWn(self, photonInfo):
        return ME.Abs_layers(d_above = self.tdead+self.tp, d_region = self.tn, mu = photonInfo.muSi, theta = photonInfo.theta_3)

    def calcW_both(self, photonInfo):
        return self.calcWp(photonInfo), self.calcWn(photonInfo)
    
    def calcOxideAttenuation(self, photonInfo):
        return ME.Abs_singleLayer(d=self.tOxide, mu=photonInfo.muOxide, theta = photonInfo.theta_2)
    
    # this calculates the absorption in the dead layer multiplies it
    # by the flat Pdiffusion probability
    # this whole thing gets multiplied by Pe
    def calcHalfBox_Wp(self, photonInfo, Pdiff=0.5):
        return ME.Abs_halfBox_Wp(d1=self.tdead, d2=self.tp, mu=photonInfo.muSi, Pdiff=Pdiff, theta=photonInfo.theta_3)
    
## ----- End of PNInfo Class ----- ##  


        
## -- Stuff for plotting -- ##

def Plot_Ws(someDictOfStuff, xRange = (150, 900), theta=0, ax=None):
    
    func_args = inspect.signature(PNJunction.__init__).parameters.keys()

    # Filter the dictionary
    filtered_args = {k: v for k, v in someDictOfStuff.items() if k in func_args}
    
    newPN = PNJunction(**filtered_args)
    
    xVals = np.linspace(xRange[0], xRange[1], 100)

    Wp = newPN.calcWp(PC.PhotonInfo(wlen = xVals, theta_1=theta))
    Wn = newPN.calcWn(PC.PhotonInfo(wlen = xVals, theta_1=theta))
    
    plt.figure(figsize=(5, 3))
    ax=plt.gca() if ax is None else 0
    ax.plot(xVals, Wp, label = 'P Region')
    ax.plot(xVals, Wn, linestyle = "--", label = 'N Region')
    
    pnText =[
        "Device info (from dict):",
        f"Dp* (start of P): {newPN.dp_star:0.2f}nm",
        f"Xpn (middle): {newPN.Xpn/1000:0.2f}um",
        f"Dw* (end of N): {newPN.dw_star/1000:0.2f}um"
    ]
    ax.text(0.1, 0.2, '\n'.join(pnText), fontsize = 8, bbox=dict(facecolor='red', alpha=0.5), transform=ax.transAxes)
    
    ax.legend()
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel(r'$W_n$ or $W_p$')
    plt.show()
    
    return ax