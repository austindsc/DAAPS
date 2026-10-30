import sys
import os

# this is to make everything path agnostic without getting too technical
# this adds the DAAPS_TOOLKIT to the path, and the desired subfolderToAdd
# this also sets the global variables for the paths to the data folderToAdd

# call this from any folder you are working in

#-- utility functions --#
def find_folder_in_path(target_folder):
    current_path = os.path.abspath(__file__)
    base_path = os.path.basename(current_path)
    
    pathComponents = current_path.split(os.sep)
    for i in range(len(pathComponents)):
        if pathComponents[-1] == target_folder:
            return os.sep.join(pathComponents)
        else:
            pathComponents.pop()
    
    print(f"Could not find {target_folder} in the path.")
    sys.exit(1)
    
    # Return None if the target folder is not found
    return None

def printLastNfolderToAdd(path, N=3):
    pathComponents = path.split(os.sep)
    if len(pathComponents) < N:
        return path
    lastThreefolderToAdd = pathComponents[-N:]
    shortPath = os.path.join(*lastThreefolderToAdd)
    return ".../"+shortPath
# --------------------- #

## Control everything from here
BASE_FOLDER = "pythonModelling"
TOOLKIT_NAME = "DAAPS_TOOLKIT"
FOLDERS_TO_ADD = ["FittingTools", "PDEModel"]
veraName = "veraData"

# to add to global variables
DATA_FOLDERS = ["opticsDATA", "literatureDATA"]

## ---- Specific functions ---- ##
def GetBasePath(baseFolder=BASE_FOLDER, toolkitName=TOOLKIT_NAME):
    ##-------- Find the path to the entire project --------##
    basePath = find_folder_in_path(baseFolder)

    ##-------- Make the path to the DAAPS_TOOLKIT --------##
    toolkitPath = os.path.join(basePath, toolkitName)
    # print("toolkit:", toolkitPath)
    if not os.path.exists(toolkitPath):
        print(f"{toolkitPath} does not exist.")
        sys.exit(1)
    return toolkitPath

#---- Set the paths to the toolkit so Python can find them easily ----#
def SetToolkitPaths(verbose=True):
    projectPath = GetBasePath()
    
    #---- adding to python path! ----#
    sys.path.insert(0, projectPath)
    print("Added to path:", projectPath) if verbose else None
    
    projectPath = GetBasePath()
    for folder in FOLDERS_TO_ADD:
        new_path = os.path.join(projectPath, folder)
        if not os.path.exists(new_path):
            print(f"{new_path} does not exist.")
            sys.exit(1)
        sys.path.insert(0, new_path)
        print("Added to path:", os.path.join(projectPath, folder)) if verbose else None
    return

def PrintDAAPSGlobalVariables():
    print("Global variables:")
    print("\tPROJECT_PATH = ", PROJECT_PATH)
    print("\tOPTICS_DATA_PATH = ", OPTICS_DATA_PATH)
    print("\tLITERATURE_DATA_PATH = ", LITERATURE_DATA_PATH)
    return

def SetGlobalPathVariables(verbose=True):
    projectPath = GetBasePath()
    global PROJECT_PATH, OPTICS_DATA_PATH, LITERATURE_DATA_PATH, VERA_DATA_PATH
    #---- Global path variables ----#
    PROJECT_PATH = projectPath
    OPTICS_DATA_PATH = os.path.join(projectPath, "opticsDATA")
    LITERATURE_DATA_PATH = os.path.join(projectPath, "literatureDATA")
    
    ##---- Vera data path ----##
    VERA_DATA_PATH = os.path.join(os.path.dirname(os.path.dirname(PROJECT_PATH)), veraName)

    globalPaths = {
        'PROJECT_PATH': projectPath,
        'OPTICS_DATA_PATH': OPTICS_DATA_PATH,
        'LITERATURE_DATA_PATH': LITERATURE_DATA_PATH,
        'VERA_DATA_PATH': VERA_DATA_PATH
    }

    # Check if the global paths exist
    for path in globalPaths.values():
        if not os.path.exists(path):
            print(f"{path} does not exist.")
            sys.exit(1)

    if verbose:
        print("global path variables:")
        for name, value in globalPaths.items():
            print("\t", name, "=.../"+printLastNfolderToAdd(value, 4))
    return

#-- Verbosity for first time --#
this_path = os.path.abspath(__file__) # Get the full file path of this file
print("calling GenPaths.py from:", printLastNfolderToAdd(this_path, 4))

SetToolkitPaths()
SetGlobalPathVariables()