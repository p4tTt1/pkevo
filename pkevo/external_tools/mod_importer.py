import ctypes
import os
import sys

from config.pkevo_config import GeneralConfig

# Make sure that the "out" directory exists
if not os.path.exists("./out"):
    os.makedirs("./out")

# Set the path to the Mød Python bindings
mod_path = GeneralConfig.MOD_PATH
sys.path.insert(0, mod_path)

# Import the mod module
import mod