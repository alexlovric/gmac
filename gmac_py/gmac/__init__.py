# This file makes 'gmac' a package and controls what's in its namespace.

# Import the submodule files
from . import core
from . import morph
from . import io

# Import all functions from core for convenience
from .core import *

# Define what `from gmac import *` does.
__all__ = ["core", "morph", "io"]
