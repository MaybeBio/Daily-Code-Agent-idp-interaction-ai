import argparse
from Protein_RPA.utils.input_parse import input_parse
import pickle
import numpy as np

def read_state(args):
    input_state=args.s
    input_param=args.i

    input_parameters = input_parse(args.i);
    nx = int(input_parameters['nx'])
    ny = int(input_parameters['ny'])
    dx = input_parameters['dx']
    phis=input_parameters['phis']
    kappa=input_parameters['kappa']
    u=input_parameters['temperature_inverse']

    state_file=open(input_state,'rb')
    steps_t,t_t,input_parameters_t,phi_value_t,xi_value_t=pickle.load(state_file) 

    cond1=int(input_parameters_t['nx'])==int(nx)
    cond2=int(input_parameters_t['phis'])==int(phis)
    cond3=int(input_parameters_t['temperature_inverse'])==int(u)
    cond4=input_parameters_t['kappa']==kappa
    rules=[cond1,cond2, cond3, cond4]
    print("rules ",rules)
    state_file.close()

    print(' state file: ')
    print( "nx ny dx ",input_parameters_t['nx']," ",ny," ",dx)
    print( "phis kappa u",input_parameters_t['phis']," ",input_parameters_t['kappa']," ",input_parameters_t['temperature_inverse'])
    print("steps t", steps_t,"  ",t_t)

    print(" average phi ", np.mean(phi_value_t))
    print(" average xi ", np.mean(xi_value_t))


if __name__ == "__main__":
     parser = argparse.ArgumentParser(description='Take output filename to produce Movie')
     parser.add_argument('--i',help="input parameters", required = True);

     parser.add_argument('--s',help="Name of state file", required = True);
     args = parser.parse_args();

     read_state(args);
