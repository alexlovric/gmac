# This file defines the gmac.io submodule.

from .gmac import write_vtp, write_stl, write_obj, read_stl, read_obj

# Other python specific functions
def info():
    print("This is the GMAC I/O submodule.")