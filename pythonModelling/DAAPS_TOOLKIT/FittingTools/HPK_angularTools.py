import os
import pandas as pd
import numpy as np

import PDEModel.PhotonInfo as PI

import FittingFunctions as FitF

##### used for loading in many csv's at once (see below) #####
def get_csv_files_with_identifier(directory, str_id):
    csv_files = []
    for file in os.listdir(directory):
        if file.endswith(".csv") and str_id in file:
            csv_files.append(os.path.join(directory, file))
    return csv_files

def get_chars_before(toSearch, findBefore, numBefore = 3):
    index1 = toSearch.find(findBefore)
    if index1 == -1:
        return None
    else:
        return toSearch[index1-numBefore:index1]

################ old HPK data loader ################
## for use with HPK data that isn't packaged into a csv (which becomes a df)
class AngularData:
    
    def __init__(self, path, comment="#"):
        self.df = pd.read_csv(path, header = 0, comment=comment)
        self.names = self.df.columns.values
        self.baseName = os.path.basename(path)
        
        self.wlen = float(get_chars_before(self.baseName, "nm_", numBefore=3))
        self.voltage = get_chars_before(self.baseName, "oV", numBefore=1)
      
        self.xRaw = None
        self.Ipeak = None
        self.IpeakErr = None
        
        # dn subtracted
        # self.Isub = None
        # self.IsubErr = None
        return None
    
    def RemoveMax(self):
        maxAngle = self.df[self.names[0]].max()
        # remove all rows with max angle
        self.df = self.df[self.df[self.names[0]] != maxAngle]
        return None
    
    def LoadRawdata(self, firstDark = 1, secondDark = -1):
        self.xRaw = -1*self.df[self.names[0]].values[firstDark:secondDark] # flip sign of angles
        self.Iraw = self.df[self.names[1]].values[firstDark:secondDark]
        self.IrawErr = self.df[self.names[2]].values[firstDark:secondDark]
        return
    
    def GetDark(self):
        first_row = self.df.iloc[0]
        last_row = self.df.iloc[-1]
        
        self.darkAvg = np.mean([last_row[self.names[1]], first_row[self.names[1]]])
        self.darkDif = np.abs(last_row[self.names[1]] - first_row[self.names[1]])/first_row[self.names[1]]
        self.darkError = np.sqrt(last_row[self.names[2]]**2 + first_row[self.names[2]]**2)/2
        return first_row[self.names[1]], first_row[self.names[2]], last_row[self.names[1]], last_row[self.names[2]]
    
    def SubtractDark(self): # must call GetDark first
        self.Iraw = np.array(self.Iraw) - self.darkAvg
        self.IrawErr = np.sqrt(self.IrawErr**2 + self.darkError**2)
        return None
    
    def Normalize(self):
        if self.Ipeak is None:
            
            self.Ipeak = self.Iraw[np.argmin(np.abs(self.xRaw))]
            self.IpeakErr = self.IrawErr[np.argmin(np.abs(self.xRaw))]
        
        self.Inorm = self.Iraw/self.Ipeak
        # self.InormErr = self.Inorm * (self.IrawErr/self.Iraw)
        self.InormErr = self.Inorm * np.sqrt((self.IrawErr/self.Iraw)**2 + (self.IpeakErr/self.Ipeak)**2)
        return None
    
    def LoadData(self, doDark=True, normalize = False, darkDifThreshold = 0.05, skip=[1, -1]):
        self.LoadRawdata(skip[0], skip[1])
        
        if doDark:
            self.GetDark()
            if self.darkDif > darkDifThreshold:
                print("Dark subtraction may not be valid, dark difference is too large")
                print(f"Dark difference: {self.darkDif}")
                
            self.SubtractDark()
        
        if normalize:
           self.Normalize()
    
        return None
    
    
    def Describe(self, verbose = False):
        print("---- Angular Data ----")
        print("Base Name:", self.baseName) if verbose else None
        print("Voltage:", self.voltage) if verbose else None
        print("Wavelength:", self.wlen)
        print(f"peak Current: {self.Ipeak*1e6:0.3f} +/- {self.IpeakErr*1e6:0.4f}") if self.Ipeak is not None else None
        print("----- \n")
        return None
    

######################### For newer data from large dataframes #########################
## also wraps up previous dataclass, everything is finalized here

# for individual data that has been extracted from a larger dataframe 
class AngularObject:
    
    def __init__(self, voltage, wavelength, fwhm, angles, yRaw, yRawErr, yNorm=None, yNormErr=None, yNorm_range=5):
        # fixed variables
        self.voltage = voltage
        self.wlen = wavelength
        self.fwhm = fwhm
        
        # self.df = DF
        # arrays
        self.yNormalMean = None # raw data mean at normal incidence
        self.yAll = None # fill these somewhere else
        self.yAllErr = None
        
        # dark stuff
        self.d1 = None
        self.d2 = None
        self.darkAvg = None
        self.darkError = None
        self.darkDif = None
        
        # Dark subtracted
        self.angles = angles
        self.yRaw = yRaw
        self.yRawErr = yRawErr
        
        # normalized
        self.yNorm = yNorm
        self.yNormErr = yNormErr
        self.yNorm_range = yNorm_range
        self.yMeanstd = None
        
        # factored out
        self.yNorm_noOptics = None
        self.yNormErr_noOptics = None
        self.yNorm_noPDE = None
        self.yNormErr_noPDE = None
        return
    
    ##
    # incorporate these into the 'build angular object from...' function
    def GetDark(self, firstDark = 1, secondDark = -1):
        self.d1 = self.yAll[firstDark]
        self.d1err = self.yAllErr[firstDark]
        self.d2 = self.yAll[secondDark]
        self.d2err = self.yAllErr[secondDark]
        self.darkAvg = np.mean([self.d1, self.d2])
        self.darkError = np.sqrt(self.yAllErr[firstDark]**2 + self.yAllErr[secondDark]**2)/2
        return self.darkAvg, self.darkError
    
    def SubtractDark(self):
        self.angles = self.angles
        self.yRaw = self.yRaw - self.darkAvg
        self.yRawErr = np.sqrt(self.yRawErr**2 + self.darkError**2)
        return None
    ##
    
    def FactorOutTransmission(self, tOxide, dh):
        # trans0 = PI.PhotonInfo(wlen = self.wlen, theta_1=0., dataHelper = dh).T_total_windowAndOxide(thick_oxide=tOxide)
        # transVals1 = PI.PhotonInfo(wlen = self.wlen, theta_1=self.angles, dataHelper = dh).T_total_windowAndOxide(thick_oxide=tOxide)
        # transForm = transVals1/trans0
        Ref0 = PI.PhotonInfo(wlen = self.wlen, theta_1=0., dataHelper = dh).R_total_windowAndOxide(thick_oxide=tOxide)
        RefVals1 = PI.PhotonInfo(wlen = self.wlen, theta_1=self.angles, dataHelper = dh).R_total_windowAndOxide(thick_oxide=tOxide)
        transForm = (1-RefVals1)/(1-Ref0) # relative transmission
        
        self.yNorm_noOptics = self.yNorm / transForm
        self.yNormErr_noOptics = self.yNormErr / transForm
        return self.yNorm_noOptics, self.yNormErr_noOptics, transForm
    
    def FactorOutPDE(self, pnArgs):
        pdeTransform = FitF.RelativePDE_angle(wlen=self.wlen, theta1=self.angles, **pnArgs)
        self.yNorm_noPDE = self.yNorm / pdeTransform
        self.yNormErr_noPDE = self.yNormErr / pdeTransform
        return self.yNorm_noPDE, self.yNormErr_noPDE, pdeTransform
    
    # normalize when only having passed the raw data
    def Normalize_fromSelfArrays(self, meanRange):
        self.yNorm_range = meanRange
        
        # finding middle (normal incidence) angle
        normal_meanArray = self.yRaw[(self.angles >= -meanRange) & (self.angles <= meanRange)]
        normal_meanErr_array = self.yRawErr[(self.angles >= -meanRange) & (self.angles <= meanRange)]
        
        self.yNormalMean = np.mean(normal_meanArray)
        self.yMeanstd = np.std(normal_meanArray)
        normal_meanErr = np.sqrt(np.sum([i**2 for i in normal_meanErr_array])) / len(normal_meanErr_array)        
        
        self.yNorm = self.yRaw / self.yNormalMean # current normalized
        self.yNormErr = self.yNorm * np.sqrt((self.yRawErr/self.yRaw)**2 + (normal_meanErr/self.yNormalMean)**2)
        
        return self.yNorm, self.yNormErr
    
    def GetNormalizedSTD(self, shoulderRange = 30):
        noShoulder = self.yNorm[(self.angles >= -shoulderRange) & (self.angles <= shoulderRange)]
        return np.std(noShoulder)
        
    
    def DataQuality(self):
        self.sortbyX()
        self.removeDuplicatesinX()
        return None
    
    def sortbyX(self):
        idx = np.argsort(self.angles)
        self.angles = self.angles[idx]
        self.yRaw = self.yRaw[idx]
        self.yRawErr = self.yRawErr[idx]
        self.yNorm = self.yNorm[idx]
        self.yNormErr = self.yNormErr[idx]
        return None
    def removeDuplicatesinX(self):
        self.angles, idx = np.unique(self.angles, return_index=True)
        self.yRaw = self.yRaw[idx]
        self.yRawErr = self.yRawErr[idx]
        self.yNorm = self.yNorm[idx]
        self.yNormErr = self.yNormErr[idx]
        return None

def BuildAngularObject_fromDF(df, strWlen, strFWHM, strV, strAngle, strCurrent, strCurrentErr, angleMean_range = 5, darkInd=[1, -1], doDark = True, normalize = True):
    voltage = df[strV].values[0]
    wlen = df[strWlen].values[0]
    fwhm = df[strFWHM].values[0]
    
    xVals = df[strAngle].values[darkInd[0]:darkInd[1]]
    # angleNormal = df.loc[df[strAngle].abs().idxmin()][strAngle] # get angle for minimum angle
    # currentNormArray = df[(df[strAngle] >= angleNormal - angleMean_range) & (df[strAngle] <= angleNormal + angleMean_range)][strCurrent].values
    # currentErrorNormArray = df[(df[strAngle] >= angleNormal - angleMean_range) & (df[strAngle] <= angleNormal + angleMean_range)][strCurrentErr].values
    
    yRaw = df[strCurrent].values[darkInd[0]:darkInd[1]]
    yRawErr = df[strCurrentErr].values[darkInd[0]:darkInd[1]]
    
    yAll = df[strCurrent].values
    yAllErr = df[strCurrentErr].values
    
    ao = AngularObject(voltage=voltage, wavelength=wlen, yRaw=yRaw, yRawErr=yRawErr, angles=xVals, fwhm=fwhm, yNorm_range=angleMean_range)
    
    ao.yAll = yAll # just to hold them so we can call doDark anywhere
    ao.yAllErr = yAllErr
    
    if doDark:
        ao.GetDark(darkInd[0], darkInd[1])
        ao.SubtractDark()
    if normalize:
        ao.Normalize_fromSelfArrays(angleMean_range) 
    return ao

# must pass Dark current subtracted values to here
def BuildAngularObject_fromAngularData(angData: AngularData, fwhm=None, angleMean_range = None):
    
    ao = AngularObject(angles = angData.xRaw, voltage=angData.voltage, wavelength=angData.wlen, yRaw=angData.Iraw, yRawErr=angData.IrawErr, fwhm = fwhm)
    ao.Normalize_fromSelfArrays(angleMean_range)
    ao.fwhm = fwhm
    return ao


