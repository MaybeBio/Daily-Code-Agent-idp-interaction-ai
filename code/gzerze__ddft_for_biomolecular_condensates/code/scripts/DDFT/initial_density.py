import numpy as np
import matplotlib.pyplot as plt


# Define the initial density profiles and normalize them.
# Inputs:
#   nx : number of grid points in the x-direction
#   ny : number of grid points in the y-direction
# The initial profiles can be modified or extended to match the proposed simulation setup.

def phi_1D(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    mid = int((ny*nx)/2)
    for i in range(nx):
        for j in range(ny):
            if (i<mid) :
                phi_value[i]=0.149949*0.9
            else:
                phi_value[i]=0.00002687

    return phi_value

def phi_2D(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            if (i>10) and (i<25) :
                    if(j>-1) and (j<250):
                        phi_value[i,j]=0.149949*0.5
                        #phi_value[i,j]=0.149949*0.9
                    else:
                        phi_value[i,j]=0.00002687
            else:
                phi_value[i,j]=0.00002687


    for i in range(nx):
        for j in range(ny):
            if (j>10) and (j<25) :
                phi_value[i,j]=0.149949*0.5


    return phi_value

def phi_circle0(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-30)**2)
            if d<15 :
                phi_value[i,j]=0.149949*0.5
            else:
                phi_value[i,j]=0.00002687

    return phi_value



def phi_circle(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-30)**2)
            if d<15 :
                phi_value[i,j]=0.149949*0.5
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949
    return phi_value

def phi_circle2(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-20)**2)
            if d<10 :
                phi_value[i,j]=0.149949*0.5
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949*0.7
    return phi_value

def phi_circle2_2(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-24)**2)
            if d<10 :
                phi_value[i,j]=0.149949*0.9
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<5:
                phi_value[i,j]=0.149949
    return phi_value

def phi_circle2_3(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-24)**2)
            if d<10 :
                phi_value[i,j]=0.149949
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949
    return phi_value

def phi_circle2_4(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-24)**2)
            if d<10 :
                phi_value[i,j]=0.149949*0.9
            else:
                phi_value[i,j]=0.00002687*100


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949*0.9
    return phi_value

def phi_circle2_5(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-24)**2)
            if d<10 :
                phi_value[i,j]=0.149949*0.7
            else:
                phi_value[i,j]=0.00002687*100


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949*0.3
    return phi_value


def phi_circle2_6(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-24)**2)
            if d<10 :
                phi_value[i,j]=0.149949*0.7
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949*0.3
    return phi_value


def phi_circle3(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-30)**2)
            if d<15 :
                phi_value[i,j]=0.149949*0.1
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-40)**2)
            if  d2<10:
                phi_value[i,j]=0.149949*0.7
    return phi_value

def phi_circle4(nx,ny):
    phi_value=np.array( [[0.]*ny]*nx )
    for i in range(nx):
        for j in range(ny):
            d=np.sqrt((i-30)**2+(j-30)**2)
            if d<15 :
                phi_value[i,j]=0.149949*0.5
            else:
                phi_value[i,j]=0.00002687


            d2=np.sqrt((i-30)**2+(j-30)**2)
            if  d2<10:
                phi_value[i,j]=0.149949*0.7
    return phi_value


