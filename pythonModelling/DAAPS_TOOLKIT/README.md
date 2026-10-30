## Code Requirements
You will need a few things. Nothing crazy. In general the scientific packages:
- numpy
- scipy
- Iminuit (ROOT minimizer wrapped for python) if you want to use the fitting tools


## Using Python Tools (simple)
To load the data and use the tools, uou can do it manually, or use GenPaths.py. I did NOT write a setup.py file, or __init__.py, or build an environment or export as its own package to keep things simple.

Genpaths.py is a simple script that add the python paths to the system path, and sets up some global variables (strings) to access data in DAAAPS_TOOLKIT folder

```
import sys
sys.path.append('../path to find DAAPS_TOOLKIT or GenPaths.py from wherever you are working')
from DAAPS_TOOLKIT.GenPaths import * # this line sets up the paths to find. All done!
```
now you can import the modules however you want. Look at the notebooks for examples. May need to stay consistent in namespace aliases between modules and your notebooks. Usually it looks like below:

```
import PDEModel.MathEquations as ME
from PDEModel.OpticsDataLoader import *
import PDEModel.PDEStyle as MS
from PDEModel.PhotonInfo import * # this calls DataLoader when defining init
from PNJunction import * 
```

## Using Python Tools (advanced)
someone can make this an enviroment or standalone package. Unecessary overhead for our group.

## Notebook Status
After reorganizing the code, I am not going to modify every notebook to get it to work. ``` fitting_preliminary_nb/angular_v5_HPK_newData.ipynb``` works and that is all i fixed. The rest are not working. You can fix as needed, its just a few paths that must be changed and the imports.

## Tool Description

## Folder Structure in DAAPS_Toolkit
```
pythonModelling
├── DAAPS_TOOLKIT
│   ├── GenPaths.py                 
│   ├── FittingTools
│   │   ├── FillBetween3d.py        # 3d fill between, for visualizing 3d fits
│   │   ├── FittingFunctions.py     # functions to fit to data. This wraps some
│   │                               # PDEModel tools into single functions with
│   │                               # large sets of arguments or parameters
│   │   ├── HPK_angularTools.py     # tools for organizing HPK angular data
│   │   ├── PDE_fittingTools.py     # most useful. Tool for fitting multi dimensional
│   │                               # datasets with flexible fixed and floating input arguments
│   ├── PDEModel
│   │   ├── MathEquations.py        # math functions for PDEModel
│   │   ├── OpticsDataLoader.py     # tools for loading n, k data and interpolating
│   │   ├── PDEStyle.py             # style for plots
│   │   ├── PNJunction.py           # tools for PN junctions, including absorption
│   │   ├── PhotonInfo.py           # Hold photon info, and use math equations
│   │   └── test.py                 # test file for PDEModel
│   ├── README.md                   # this file LOL
│   │
│   └── opticsDATA          # self explanatory
│       ├── GoodMerged      # merged datasets
│       ├── ...
```

## Folder Structure in remainder of pythonModelling
```
pythonModelling
├── main_Angular_PlotTools.ipynb    # plotting angular tools and such
├── main_InputData.ipynb            # showing input n and k data, and how to load it, and merging some
├── main_PlotTools.ipynb            # plotting tools for 1d data (sections of PDE equation)
├── main_PlotTools2d.ipynb          # plotting 2d examples (sections of PDE equation)
│
├── eXT_andToyModels_nb             # modelling external cross-talk using PDE 'backwards'
│   ├── HPK_140nm_Production_And_Transmission_highd.root
│   ├── ext_validations.ipynb
│   └── nobleFunctions.py
│
├── temperature_nb                  # modelling temperature effects (diffusion and absorption)
│   ├── diffusion_v1.ipynb
│   ├── dmpaper_supplements
│   └── peAbs_v1.ipynb
│
│
├── fitting_preliminary_nb         # fitting notebooks that are preliminary. only some working...
│   ├── angular_fitting_v1.ipynb
│   ├── angular_v2_AllData.ipynb
│   ├── .....
│
├── fitting_final_nb                # final fitting notebooks (highest quality etc)
```