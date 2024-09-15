# Note Patrick: for reasons unknown to me, RDKit imports can be sometimes troubling (causing errors 
# like "ModuleX.methodY did not match C++ signature:..."), which is most likely related to RDKit using Boost.Python.
# Such errors can occur when the first test module imports for example only "Chem" from rdkit and then in the second
# test module we import not only "Chem", but also "AllChem". Even worse the second test module does not have
# to explicitly import those modules. In order for this error to occur it is enough if that second test module is testing 
# another module which is importing and using "AllChem".
# 
# For that reason, perform all imports from RDKit in the init file to ensure that this gets imported first when executing
# the tests. Now, all test modules will run without errors 

from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem 

