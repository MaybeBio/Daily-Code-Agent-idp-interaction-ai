# Homology modeling by the automodel class
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
import sys

log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments('A', 686)

a = MyModel(env,
              alnfile  = 'aligs.pir',     # alignment filename
              knowns   = ('3w3w',),              # codes of the templates
              sequence = 'ecm29_684_1749',
              assess_methods=assess.normalized_dope)              # code of the target
a.starting_model= 1                 # index of the first model
a.ending_model  = 10                 # index of the last model
                                    # (determines how many models to calculate)
if '--test' in sys.argv: a.ending_model = 1
a.make()                            # do the actual homology modeling

             
