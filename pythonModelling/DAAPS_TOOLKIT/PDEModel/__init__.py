# # # Get the directory of the current file
# # current_dir = os.path.dirname(__file__)

# # # Get all .py files in the directory, excluding __init__.py
# # modules = glob.glob(os.path.join(current_dir, "*.py"))
# # __all__ = [os.path.basename(f)[:-3] for f in modules if os.path.isfile(f) and not f.endswith('__init__.py')]


# # # Remove 'test.py' from __all__
# # dontInclude = ['test']
# # __all__ = [module for module in __all__ if module not in dontInclude]

# # # Import all modules
# # for module in __all__:
# #     __import__(f"{__name__}.{module}")

# # import sys
# # import os

# # # Add the path to DAAPS_TOOLKIT
# # sys.path.append(os.path.abspath('/path/to/DAAPS_TOOLKIT'))

# # from ../GenPaths import *

# import numpy as np

# from .PDEStyle import *
# import DAAPS_TOOLKIT.PDEModel.PDEStyle as MS

# from .MathEquations import *
# import DAAPS_TOOLKIT.PDEModel.MathEquations as ME

# from .OpticsDataLoader import *
# import DAAPS_TOOLKIT.PDEModel.OpticsDataLoader as ODL

# from .PhotonInfo import *
# from .PNJunction import *

# # Optionally, you can define what gets imported when using `from my_package import *`
# __all__ = ['MS', 'ME', 'ODL', 'PhotonInfo', 'PNJunction', 'np']


