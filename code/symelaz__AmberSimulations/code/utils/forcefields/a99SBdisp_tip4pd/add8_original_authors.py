import parmed
from parmed.tools import addLJType
from parmed.tools import addDihedral
from parmed.amber import AmberParm
import sys

## load prmtop file
PRMTOP=sys.argv[1]
NR = sys.argv[2]
parm = AmberParm(PRMTOP)

## system-specific information
Nresidues =  int(NR)
str_n = NR
N_minus_one = Nresidues - 1
str_n_minus_one = str(N_minus_one)

prot = f"1- {Nresidues}"

prot_2_to_n_minus_one = f":2-{Nresidues -1}"
prot_2_to_n = f":2- {Nresidues}"
prot_3_to_n = f":3- {Nresidues}"
prot_1_to_n_minus_one = f":1-{Nresidues -1}"

f1action1 = addDihedral(parm, prot_2_to_n_minus_one+"@N", prot_2_to_n_minus_one+"@CA", prot_2_to_n_minus_one+"@C", prot_3_to_n+"@N", -0.197, 1, 180.0)
f1action2 = addDihedral(parm, prot_2_to_n_minus_one+"@N", prot_2_to_n_minus_one+"@CA", prot_2_to_n_minus_one+"@C", prot_3_to_n+"@N", -0.255, 1, 90.0)
f1action3 = addDihedral(parm, prot_2_to_n_minus_one+"@N", prot_2_to_n_minus_one+"@CA", prot_2_to_n_minus_one+"@C", prot_3_to_n+"@N",  0.075, 2, 90.0)
f1action4 = addDihedral(parm, prot_2_to_n_minus_one+"@N", prot_2_to_n_minus_one+"@CA", prot_2_to_n_minus_one+"@C", prot_3_to_n+"@N",  0.057, 3, 90.0)
f1action5 = addDihedral(parm, prot_2_to_n_minus_one+"@N", prot_2_to_n_minus_one+"@CA", prot_2_to_n_minus_one+"@C", prot_3_to_n+"@N",  0.025, 4, 90.0)
f1action6 = addDihedral(parm, prot_2_to_n_minus_one+"@N", prot_2_to_n_minus_one+"@CA", prot_2_to_n_minus_one+"@C", prot_3_to_n+"@N", -0.011, 5, 90.0)
f1action8 = addDihedral(parm, prot_1_to_n_minus_one+"@C", prot_2_to_n+"@N", prot_2_to_n+"@CA", prot_2_to_n+"@C",  0.40, 2, 90.0)

actions = [f1action1, f1action2, f1action3, f1action4, f1action5, f1action6, f1action8]

for i in range(len(actions)):
    actions[i].execute()
#    print(str(actions[i]))
parm.save("add8_"+PRMTOP)
