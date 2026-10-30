import inspect
from typing import Callable
import numpy as np
import iminuit

# from iminuit import Minuit # need to install this!

####### Handling Errors and Outputs (to get errors) #######
# class myParameter():
#     def __init__(self, name, value, error):
#         self.name = name
#         self.value = value
#         self.bounds = None
#         self.error = error

# def ResultDict(fitResult, x0Names):
#     resultDict = {}
#     for i, key in enumerate(x0Names):
#         value = fitResult.x[i]
#         error = np.sqrt(fitResult.hess_inv[i, i])
#         mpar = myParameter(name=key, value = value, error = error)
#         resultDict[key] = mpar


def GetResDict(fitResult, x0Names, rounding = 4):
    return {x0Names[i]: round(fitResult.x[i], rounding) for i in range(len(x0Names))}

def GetResErrorDict(fitResult, x0Names, rounding = 4):
    errors = GetErrorFromInvHess(fitResult.hess_inv)
    return {x0Names[i]: round(errors[i], rounding) for i in range(len(x0Names))}

def GetErrorFromInvHess(inv_hess): # only use with L-BFGS-B
    # inv_hess = resHessProduct.todense()
    # hess = np.linalg.inv(inv_hess)
    return np.sqrt(np.diag(inv_hess))

def GetCorrelationMatrix(inv_hess):
    std_devs = GetErrorFromInvHess(inv_hess)
    corrMatrix = inv_hess / np.outer(std_devs, std_devs)
    return corrMatrix

##-------- Misc. Printing and Formatting --------##
def ShortenDict_toPrint(d):
    soutArrs = []
    for key, value in d.items():
        if hasattr(value, '__len__') and len(value) > 3:
            value = list(value)[:3] + ['...']
        sEntry= f"'{key}': {value}"
        soutArrs.append(sEntry)
    return "{" + ", ".join(soutArrs) + "}"

##-------- List sorting and checking tools --------##
def GetNames(func):
    return inspect.getfullargspec(func).args
GetFuncArgs = GetNames

def ListSubtraction_orderd(list1, list2):
    return [x for x in list1 if x not in list2]

def NewDictFromKeysInList(dict1, keyList):
    return {key: val for key, val in dict1.items() if key in keyList}
def NewDictFromKeysNOTInList(dict1, keyList):
    return {key: val for key, val in dict1.items() if key not in keyList}

def is_list_of_dicts(obj):
    return isinstance(obj, list) and all(isinstance(item, dict) for item in obj)

# returns a copy of dict1 minus dict2 keys
# def dictKeySubtraction(dict1, dict2):
#     return {key: dict1[key] for key in dict1.keys() - dict2.keys()}


##------------------------ Sorting and Building Inputs ------------------------##
# order x0 has to match the ordering in functionToFit
# therefore, this sorts the floating parameters into the correct order
def BuildX0(floatDict, functionToFit):
    x0 = []
    x0Names = []
    for key in GetNames(functionToFit):
        if key in floatDict.keys():
            x0.append(floatDict[key])
            x0Names.append(key)
    return x0, x0Names

# ditto as above, but sorts the bounds
def BuildX0_bounds(floatBoundsDict, functionToFit):
    bounds = []
    for key in GetNames(functionToFit):
        if key in floatBoundsDict.keys():
            bounds.append((floatBoundsDict[key][0], floatBoundsDict[key][1]))
    return bounds

def AddNonesToDict(dAddTo, dRef):
    for key in dRef.keys():
        if key not in dAddTo:
            dAddTo[key] = (None, None)
    return None

##------------------------ Checking Data and Key organization ------------------------##
###########
# make sure all keys in the dataset are consistent (same headers, same dimensionality of data)
def CheckDataSetKey_consistent(dataObj, verbose = False):
    print(f"-- CheckDataSetKey_consistent() --") if verbose else None
    keys1 = set(dataObj[0].keys())
    all_same_keys = all(set(d.keys()) == keys1 for d in dataObj)
    if all_same_keys:
        print("\tDataset(s) have same keys, proceed") if verbose else None
        return all_same_keys
    else:
        print("\tData sets have different keys, cannot proceed") if verbose else None
        print([list(d.keys()) for d in dataObj])
        return all_same_keys

###########
# check keys are complete and consistent (nothing missing, nothing extra)
def CheckInputKeySets_completeness(func, dataObj, dFloat, dFixed, yNames=("yVals", "yErrs"), verbose = False):
    print(f"-- CheckInputKeySets_completeness(func={func.__name__}) --") if verbose else None
    #### First check everything is present ####
    funcKeys = set(GetNames(func))
    dataKeys = set(dataObj[0].keys())
    
    extraData = dataKeys - funcKeys
    if len(extraData) == 1 and yNames[0] in extraData:
        print("\tSuccess: one extra key in data array:", extraData, "\t(no uncertainties included)") if verbose else None
    elif len(extraData) == 2 and set(yNames) == extraData:
        print("\tSuccess: two extra keys in data array:", extraData, "\t(uncertainties included)") if verbose else None
    else:
        print('\tExtra data keys:', extraData, "\n\tBAD, this is wrong")
        print(f"\textra data keys should have len 1 or 2, and include headers {yNames[0]} and/or {yNames[1]}")
        return False
    noYdata_Keys = dataKeys - extraData
    
    ## check for extra keys in float or fixed dictionaries
    set_extra_inputs = set(dFloat.keys()).union(set(dFixed.keys())) - funcKeys - noYdata_Keys
    if len(set_extra_inputs) != 0:
        print("\tMistake! Extra keys passed to dFloat or dFixed that aren't in the function or dataset(func args):  ", set_extra_inputs)
        return False
    
    ## check if float and fixed pars have any of the SAME variables
    sharedKeys = set(dFloat.keys()).intersection(set(dFixed.keys()))
    if len(sharedKeys) != 0:
        print("\tMistake! Float and Fixed dictionaries have the same key:", sharedKeys)
        return False
    
    ## check for keys missing from the function signature
    set_missing_fromFunc = funcKeys - noYdata_Keys - set(dFloat.keys()) - set(dFixed.keys())
    if len(set_missing_fromFunc) != 0:
        print("\tMistake! Keys missing. In function but not included in the dataset or float/fixed:  ", set_missing_fromFunc)
        return False
    
    print("\tFloat/fixed dicts and data are complete, all keys are accounted for.") if verbose else None
    return True

###########
# check that the data is in the right format, that everything maps up 1:1 between float, fixed, x0, data and the function
def CheckInputs_forResiduals(x0, fixedArgs_dict, dataObj, functionToFit, yNames = ("yVals", "yErrs"), verbose = True):
    print(f"-- CheckInputs_forResiduals(function={functionToFit.__name__}) --") if verbose else None
    
    ## AA labels match labelling in function below
    # x0 must be sorted appropriately
    dataHeaders = list(dataObj[0].keys()) # already checked with CheckDataSetKey_consistent
    yName = yNames[0]
    yErrName = yNames[1]

    ## AA
    # check the data format is at least long enough
    if len(dataHeaders) < 2:
        print("\tError. Not enough data headers supplied!\nneed at least 2 (x and y)")
        return False
    
    # sorted, ordered lists of args
    allArgs = GetNames(functionToFit)           # everything in function signature
    nonDataArgs = ListSubtraction_orderd(allArgs, dataHeaders) # without data headers that are in function
    missingFloatArgs = ListSubtraction_orderd(nonDataArgs, fixedArgs_dict.keys()) # this is the ordered keys that x0 correspond to
    # missingFloatArgs = allFunctionArgs - dataHeaders - fixedArgs = x0 parameters
    
    # xx troubleshooting
    # print("\t all func args minus data", nonDataArgs)
    # print("\t function missing args (x0 keys)", missingFloatArgs)
    
    ## BB
    ## length testing, ensuring enough vars are passed to x0 (float) and fixed, compared to function - relevant data Args
    if len(x0) + len(fixedArgs_dict) != len(nonDataArgs):
        print("\tError:missing or extra args passed, len(x0) + len(fixedDict) != len(argsNeeded)")
        print("\t\tlen(x0 + fixed):", len(x0) + len(fixedArgs_dict), "  len(nonDataArgs)", len(nonDataArgs))
        return False
    
    # check that number of missing keys in fixed parameter list equal x0 len (these keys correspond x0)
    # in other words, do we have enough keys missing from the function signature to match the x0 values?
    if len(x0) != len(missingFloatArgs):
        print("\tError: len of args missing from function does not match x0!")
        print("\t\t", len(missingFloatArgs), len(x0))
        return False
    
    # create a new dictionary, add the keys that are missing from the function signature
    # and set the values to the x0 values (match x0 to the float parameter names)
    # imo this is the most dangerous step if the ordering of x0 gets messed up somehow
    floatDict = {}
    for i, key in enumerate(missingFloatArgs): # x0 must be sorted appropriately
        floatDict[key] = x0[i]
    fullArgDict = {**fixedArgs_dict, **floatDict} # artifact. keep it
    
    ## Data loop! Most important part
    yErrsIncluded = []
    for i, oneD_dict in enumerate(dataObj):
        # get the keys from one-dimensional data dict that are also in the function signature (everything except Y, Yerr)
        dataDict_forFunction = NewDictFromKeysInList(dict1=oneD_dict, keyList=allArgs)
        dataDict_forData = NewDictFromKeysNOTInList(dict1=oneD_dict, keyList=allArgs) # just Ys and possibly Yerrs
        dataLen = len(dataDict_forData)
        
        ## CC
        # check if we have at "y" or "yErrs" in the data dict
        if dataLen==1 and yName not in dataDict_forData: # no errors, missing y key
            print(f"\tError: No y values in data dict, need {yName} key/value pair")
            print("\t\tData dict:", dataDict_forData)
            return False
        elif dataLen==1 and yName in dataDict_forData: # no errors, good conditions
            yErrsIncluded.append(False)
        elif dataLen==2 and (yName not in dataDict_forData or yErrName not in dataDict_forData): # errors, missing a key:
            print(f"\t Error: No y values or y-uncertainties in data dict, need *{yName}* and *{yErrName}* key/value pairs")
            print("\t\tData dict:", dataDict_forData)
            return False
        elif dataLen==2 and (yName in dataDict_forData or yErrName in dataDict_forData): # errors and y key included
            yErrsIncluded.append(True)
        else:
            print("\tremaining stuff in Data dict is not the right length, need 1 value (Y) or 2 values (Y, Yerrs)")
            print("\tData dict for data:", dataDict_forData.keys())
            return False
        
        # print outputs
        if verbose:
            print(f"\t---------------- data block {i} --------------------")
            print("\t\tfloat dict (from x0):", floatDict)
            print("\t\tfixed dict:", fixedArgs_dict)
            print(f"\t\tdataDict for function, (keys):", ShortenDict_toPrint(dataDict_forFunction))
            print(f"\t\tdataDict for data (y):", ShortenDict_toPrint(dataDict_forData))
    print("\t----------------") if verbose else None
    # done forloop over 1d data sets
    

    # output status of errors being good or not
    if all(yErrsIncluded):
        print("\t(All data has uncertainties, Chi2 is cost function)") if verbose else None
        print("\tSafe to proceed with GenericResiduals() function") if verbose else None
        print("\t----------------") if verbose else None
        return True
    elif not any(yErrsIncluded):
        print("\t(No datasets have uncertainties, using Residual^2 as cost function)") if verbose else None
        print("\tSafe to proceed with GenericResiduals() function") if verbose else None
        print("\t----------------") if verbose else None
        return True
    else:
        print("\tERROR: Data sets have mixed inclusion of uncertainties, should not proceed")
        return False

#### Wrap conditional checks into single function ####
def CheckAllConditions_forResiduals(func, dataObj, dFloat, dFixed, dBounds=None, yNames = ("yVals", "yErrs"), verbose=False):
    print(f"------- CheckAllConditions_forResiduals(function={func.__name__}) -------") if verbose else None
    x0, x0Names = BuildX0(dFloat, func)
    
    if dBounds is not None:
        x0bounds = BuildX0_bounds(dBounds, func)
        if len(x0) != len(x0bounds):
            print("Error: x0 and x0bounds are not the same length")
            return False
    
    
    ## check data
    test1 = CheckDataSetKey_consistent(dataObj, verbose=verbose)
    if not test1:
        CheckDataSetKey_consistent(dataObj, verbose=True)
        print(f"\tError: ({CheckDataSetKey_consistent.__name__}()) failed, should not proceed")
        return False
    
    ## check input completeness     
    test2 = CheckInputKeySets_completeness(func, dataObj, dFloat, dFixed, yNames = yNames, verbose=verbose)
    if not test2:
        CheckInputKeySets_completeness(func, dataObj, dFloat, dFixed, yNames = yNames, verbose=True)
        print(f" ({CheckInputKeySets_completeness.__name__}) failed, should not proceed")
    
    ## check consistent with residual logic
    test3 = CheckInputs_forResiduals(x0=x0, fixedArgs_dict=dFixed, dataObj=dataObj, functionToFit=func, yNames=yNames, verbose=verbose)
    if not test3:
        CheckInputs_forResiduals(x0=x0, fixedArgs_dict=dFixed, dataObj=dataObj, functionToFit=func, yNames=yNames, verbose=True)
        print(f" ({CheckInputs_forResiduals.__name__}) failed, should not proceed")    
    
    if test1 and test2 and test3:
        print("Passed all tests, ready to proceed with GenericResiduals() function") if verbose else None
        print("------------------------------------") if verbose else None
        return True
    else:
        print("------------------------------------")
        raise ValueError("Failed one or more tests, should not proceed with GenericResiduals() function")

##------------------------ Actual Residual Function to minimize (no errors yet) ------------------------##
# before calling this, you need define your dicts and bounds, and sort:
# x0, x0Names = BuildX0(floatDict, functionToFit)
# x0bounds = BuildX0_bounds(floatBoundsDict, functionToFit)

# then you need check everything with:
# CheckAll(func, dataObj, dFloat, dFixed, yNames = ("yVals, yErrs"))
# this calls:
    # CheckDataSetKey_consistent(dataObj)
    # CheckInputKeySets_completeness(func, dataObj, dFloat, dFixed, yNames = ("yVals, yErrs"))
    # CheckInputs_forResiduals(x0, fixedArgs_dict, dataObj, functionToFit, yNames = ("yVals, yErrs"))

# this adds the unreduced chi2 value, as you properly should for uncertainty calculation
def GenericResiduals(x0, fixedArgs_dict, dataObj, functionToFit, yNames = ("yVals", "yErrs"), verbose = False):
    print(f"------- GenericResiduals(function={functionToFit.__name__}) -------") if verbose else None
    
    # x0 must be sorted appropriately
    dataHeaders = list(dataObj[0].keys()) # already checked with CheckDataSetKey_consistent
    yName = yNames[0]
    yErrName = yNames[1]
    # usingErrors = yErrName in dataHeaders
    ## AA
    
    # sorted, ordered lists of args
    allArgs = GetNames(functionToFit)           # everything in function signature
    nonDataArgs = ListSubtraction_orderd(allArgs, dataHeaders) # without data headers that are in function
    missingFloatArgs = ListSubtraction_orderd(nonDataArgs, fixedArgs_dict.keys())
    # this is the ordered keys that x0 correspond to
    # missingFloatArgs = allFunctionArgs - dataHeaders - fixedArgs = x0 parameters
    if verbose:
        print("\tall input parameters (except data): ", nonDataArgs)
        print("\tfloating args: ", missingFloatArgs)
    
    # BB
    
    # create a new dictionary, add the keys that are missing from the function signature
    # and set the values to the x0 values (match x0 to the float parameter names)
    floatDict = {}
    for i, key in enumerate(missingFloatArgs): # x0 must be sorted appropriately
        floatDict[key] = x0[i]
    fullArgDict = {**fixedArgs_dict, **floatDict}    
    
    ## Data loop! Most important part
    resSum = 0
    for i, oneD_dict in enumerate(dataObj):
        # get the keys from one-dimensional data dict that are also in the function signature (everything except Y, Yerr)
        dataDict_forFunction = NewDictFromKeysInList(dict1=oneD_dict, keyList=allArgs)
        dataDict_forData = NewDictFromKeysNOTInList(dict1=oneD_dict, keyList=allArgs) # just Ys and Yerrs
        dataLen = len(dataDict_forData)
        ## CC
        
        # get the Y data
        y = dataDict_forData[yName]
        yErr = dataDict_forData[yErrName] if (yErrName in dataHeaders or dataLen==2) else np.ones(len(y))
        
        # print outputs
        if verbose:
            print(f"\t---------------- data block {i} --------------------")
            print("\t\tfloat dict (from x0):", floatDict)
            print("\t\tfixed dict:", fixedArgs_dict)
            print(f"\t\tdataDict for function, (keys):", ShortenDict_toPrint(dataDict_forFunction))
            print(f"\t\tdataDict for data (y):", ShortenDict_toPrint(dataDict_forData))

        # calculate the function
        yFit = functionToFit(**fullArgDict, **dataDict_forFunction)
        
        # calculate the residuals (chi2 if errors included)
        resSum += np.sum( ((y - yFit)/yErr)**2)
    
    print("\t----------------") if verbose else None
    # done forloop over 1d data sets
    return resSum

## ------------------------ Doing Nuisance Parameters! ------------------------##



## -- Class for nuisance parameter stuff below -- ##
class OneNuisanceDataSet:
    def __init__(self, oneDData:dict, nuisanceFloat:dict, nuisanceBounds:dict, nuisanceFixed:dict = {}):
        self.oneDData = [oneDData] if type(oneDData) == dict else oneDData # data has to be a list of dictionaries, even if a single dictionary
        self.mDict_nuisFloat = dict(nuisanceFloat) # use dict() because need to copy it
        self.mDict_nuisBounds = dict(nuisanceBounds)
        self.mDict_nuisFixed = dict(nuisanceFixed)
        self.minuitResult = None # output of the fit
        self.resultDict = None
        self.resultErrorsDict = None
        
        # for NDOF calculations
        self.nDOF = None
    
    # def MakeErrorsDict(self):
    #     if self.minuitResult is None:
    #         print("No fit result to make errors from")
    #         return None
    #     else:
    #         return GetResErrorDict(self.minuitResult, GetNames(self.minuitResult))
    
    @staticmethod
    def DescribeThisClass():
        print("This class is used to organize data for nuisance parameter fitting. It is a container for one-dimensional data, and the floating, fixed, and bounds dictionaries for the nuisance parameters.")
        print("Members:", OneNuisanceDataSet.__dict__)
        return None
        
    def CalcNDOF(self, xkey="wlen"):
        self.Nxdata = sum([len(dataDict[xkey]) for dataDict in self.oneDData])
        self.NlocalPars = len(self.mDict_nuisFloat)
        self.nDOF = self.Nxdata - self.NlocalPars
        return self.nDOF


# This function below takes similar input to GenericResiduals, but is designed to handle multiple datasets
# each 'dataset' is held with the class above. The actual data is the same list of dictionaries as before
# but you also pass the specific 'nuisance parameters' for each dataset, as well as any run specificy fixed parameters
# Look at the example notebook for how to package your data to use this function

def GenericResiduals_Nuisance_1(globalFloat:list, globalFixed:dict, dataPackage: list, functionToFit:Callable, yNames:tuple=("yVals", "yErrs") , verbose=False, testing=False, nestedMinimizerKwargs={}, ax=None):
    
    score_allData = 0
    
    for OneD_Class in dataPackage: # This loops on something with a set of nuisance parameters
        print("------------------- new dataset -------------------") if verbose else None
        if not isinstance(OneD_Class, OneNuisanceDataSet):
            raise ValueError("Data package is not a list of OneNuisanceDataSet objects")
        
        ## get what you need from the class. Check data formatting
        dataList_oneD = OneD_Class.oneDData # this is a list of the data wanted to analyze (including some non x/y parameters)
        
        if not is_list_of_dicts(dataList_oneD): # check the actual data is correct as well
            raise ValueError("OneD_Class.oneDData is not a list of dictionaries")
        CheckDataSetKey_consistent(dataObj=dataList_oneD, verbose=testing) # check the data keys
        
        # the local, fixed dict is the combination of the global float and local fixed parameters
        localFloatDict = OneD_Class.mDict_nuisFloat
        localFloatBoundsDict = OneD_Class.mDict_nuisBounds
        localFixedDict = OneD_Class.mDict_nuisFixed
        
        # sort the global float arguments back into a dict of 'fixed' for the local fit
        funcArgs = GetNames(functionToFit)
        remainingKeys = set(funcArgs) - set(globalFixed.keys()) - set(dataList_oneD[0].keys()) - set(localFloatDict.keys()) - set(localFixedDict.keys()) #what keys are missing from function signature
        orderedRemainingKeys = [key for key in funcArgs if key in remainingKeys] # order the keys to match the function signature order (because x0 is also sorted)
        if len(orderedRemainingKeys) != len(globalFloat):
            print("orderedRemainingKeys:", orderedRemainingKeys)
            raise ValueError("Missing required local float parameters do not match the global float parameter length.")
        globalFloat_toLocalFixed_dict = {key: globalFloat[i] for i, key in enumerate(orderedRemainingKeys)} # put the global floating parameters into a dictionary
        
        # this is everything that is 'fixed' for the nuisance parameters
        completeLocalFixedDict = {**globalFloat_toLocalFixed_dict, **globalFixed, **OneD_Class.mDict_nuisFixed}
        
        # check inputs
        localCheck = CheckAllConditions_forResiduals(func=functionToFit, dataObj=dataList_oneD, dFixed=completeLocalFixedDict, dFloat=localFloatDict, dBounds=localFloatBoundsDict, yNames=yNames, verbose=testing)
        
        if testing:
            print("\t---------------------------------")
            print("\tLocal check of inputs:", localCheck)
        
        # organize inputs
        x0local, x0Names_local = BuildX0(localFloatDict, functionToFit)
        x0local_bounds = BuildX0_bounds(localFloatBoundsDict, functionToFit)
        guessRes = GenericResiduals(x0=x0local, fixedArgs_dict=completeLocalFixedDict, dataObj=dataList_oneD, functionToFit=functionToFit, yNames=yNames, verbose=testing)
        
        # don't actually need to do local minimization. this is just a global fit
        if len(localFloatDict)==0:
            score_allData += guessRes
            continue
            
        # do nuisance parameter minimization
        localResult = iminuit.minimize(GenericResiduals, x0=x0local, bounds = x0local_bounds,
                                       args=(completeLocalFixedDict, dataList_oneD, functionToFit), method='migrad', **nestedMinimizerKwargs)
        
        print(localResult) if verbose else None
        
        # check the results?
        if localResult.success == False:
            print("ERROR: Nuisance fit did not converge")
            raise ValueError("Nuisance fit did not converge")
        
        # save the result
        x0Fit_local = localResult.x
        # resErrors = GetErrorFromInvHess(localResult.hess_inv)
        resDict = GetResDict(localResult, x0Names_local)
        resErrorsDict = GetResErrorDict(localResult, x0Names_local)
        
        OneD_Class.minuitResult = localResult
        OneD_Class.resultDict = resDict
        OneD_Class.resultErrorsDict = resErrorsDict
        
        # calculate the residuals, add to total residual
        nuisanceScore = GenericResiduals(x0=x0Fit_local, fixedArgs_dict=completeLocalFixedDict, dataObj=dataList_oneD,
                                            functionToFit=functionToFit, yNames=yNames, verbose=verbose)
        
        if verbose or testing:
            print("\tNuisance score:", nuisanceScore)

        # add scores up
        score_allData += nuisanceScore
    
    print("Total score:", score_allData) if verbose or testing else None
    
    return score_allData