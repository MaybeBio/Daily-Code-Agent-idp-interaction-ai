import argparse
from Protein_RPA.utils.input_parse import input_parse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from postprocess import *
import fipy as fp
import os
from os import listdir
from os.path import join, isfile

def key_func(x):
    """
    Input   =   List of state file paths (x)

    Output  =   List of filenames without extensions (to sort)
    """
    return int(x.split('-')[-1].rstrip('.pickle'))

def read_state(args):
    input_state=args.s
    output_dir=args.o
    input_param=args.i

    input_parameters = input_parse(args.i);

    nx = int(input_parameters['nx'])
    ny = int(input_parameters['ny'])
    dx = input_parameters['dx']
    phis=input_parameters['phis']
    kappa=input_parameters['kappa']
    u=input_parameters['temperature_inverse']

    mesh = fp.Grid2D(nx=nx, ny=ny, dx=dx, dy=dx)
    mesh=mesh-float(nx)*dx*0.5


    phi_range=[0,0.2]
    xi_range=[-1.14,-0.91]
    phi=fp.CellVariable(mesh=mesh, name=r'$\phi_{monomer}$')
    xi = fp.CellVariable(mesh=mesh)



    file_names = list((fn for fn in os.listdir(input_state) if fn.endswith('.pickle')))
    file_names = sorted(file_names,key=key_func)

    print("file_names", file_names[0],"  ",file_names[-1])

if __name__ == "__main__":
     parser = argparse.ArgumentParser(description='Take output filename to produce Movie')
     parser.add_argument('--o',help="output of plot density", required = True);
     parser.add_argument('--i',help="input parameters", required = True);
     
     parser.add_argument('--s',help="Directory of state files", required = True);
     args = parser.parse_args();

     read_state(args);
