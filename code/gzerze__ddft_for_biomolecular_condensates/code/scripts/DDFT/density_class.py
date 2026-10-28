class initial_density:
"""
Define the initial density profiles and normalize them.
Parameters
----------
nx : int
    Number of grid points in the x-direction.
ny : int
    Number of grid points in the y-direction.

The initial profiles can be modified or extended to match
the proposed simulation setup.
"""
    
    def __init__(self,nx,ny,idense,idil,nz=0):
    
        import numpy as np
        import matplotlib.pyplot as plt
        self.nx=nx
        self.ny=ny
        self.nz=nz
        self.idense = idense
        self.idil = idil
        print("dense: " + str(self.idense))
        print("dilute: " + str(self.idil))

    

    def phi_1D(self):
        import numpy as np
        mid = int((self.ny*self.nx)/2)
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                if (i<mid) :
                    phi_value[i]=self.idense*0.9
                else:
                    phi_value[i]=self.idil

        return phi_value

    def phi_2D(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                if (i>10) and (i<25) :
                        if(j>-1) and (j<250):
                            phi_value[i,j]=self.idense*0.5
                            #phi_value[i,j]=0.149949*0.9
                        else:
                            phi_value[i,j]=self.idil
                else:
                    phi_value[i,j]=self.idil


        for i in range(self.nx):
            for j in range(self.ny):
                if (j>10) and (j<25) :
                    phi_value[i,j]=self.idense*0.5


        return phi_value

    def phi_circle0(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-30)**2+(j-30)**2)
                if d<15 :
                    phi_value[i,j]=0.149949*0.9
                else:
                    phi_value[i,j]=self.idil

        return phi_value



    def phi_circle_concentric(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-60)**2)
                if d<=10 :
                    phi_value[i,j]=self.idense*0.9
                elif (d>=15) and (d<=20):
                    phi_value[i,j]=self.idense*0.9
                else:
                    phi_value[i,j]=self.idil

        return phi_value

    def phi_uniform(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                phi_value[i,j]=self.idense*0.1
        
        return phi_value
    
    def phi_uniform2(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                phi_value[i,j]=self.idense*0.1
        
        return phi_value

    def phi_uniform_m(self):
        import numpy as np
        #phi_value=np.array( [[0.]*self.ny]*self.nx )
        phi_value=np.full((self.nx,self.ny),0.01)
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-60)**2)
                if d<=5:
                    phi_value[i,j]=self.idense
        return phi_value

    def phi_uniform_x(self):
        import numpy as np
        #phi_value=np.array( [[0.]*self.ny]*self.nx )
        phi_value=np.full((self.nx,self.ny),self.idil)
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-60)**2)
                if d<=5:
                    phi_value[i,j]=self.idense
        return phi_value

    def phi_uniform_y(self):
        import numpy as np
        #phi_value=np.array( [[0.]*self.ny]*self.nx )
        phi_value=np.full((self.nx,self.ny),self.idil)
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-60)**2)
                if d<=5:
                    phi_value[i,j]=self.idense
        return phi_value

    def phi_uniform_z(self):
        import numpy as np
        #phi_value=np.array( [[0.]*self.ny]*self.nx )
        phi_value=np.full((self.nx,self.ny),self.idil)
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-60)**2)
                if d<=5:
                    phi_value[i,j]=self.idense
        return phi_value

    def phi_circle2_7(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-54)**2)
                if d<10 :
                    phi_value[i,j]=self.idense*0.9
                else:
                    phi_value[i,j]=self.idil


                d2=np.sqrt((i-60)**2+(j-70)**2)
                if  d2<10:
                    phi_value[i,j]=self.idense*0.9
        return phi_value
    
    def phi_sphere2_7(self):
        import numpy as np
        r = self.nx/2
        dr = 10
        dm = 10
        phi_value=np.array( [[[0.]*self.ny]*self.nx]*self.nz )
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    d=np.sqrt((i-r)**2+(j-(r-dr))**2+(k-r)**2)
                    if d<dm :
                        phi_value[i,j,k]=self.idense*0.9
                    else:
                        phi_value[i,j,k]=self.idil


                    d2=np.sqrt((i-r)**2+(j-(r+dr))**2+(k-r)**2)
                    if  d2<dm:
                        phi_value[i,j,k]=self.idense*0.9
        return phi_value

    def phi_circle2_7m(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-54)**2)
                if d<10 :
                    phi_value[i,j]=self.idense*0.1
                else:
                    phi_value[i,j]=self.idil


                d2=np.sqrt((i-60)**2+(j-70)**2)
                if  d2<10:
                    phi_value[i,j]=0.9*self.idense
        return phi_value



    def phi_random(self):
        import numpy as np
        import random
        random.seed(10)
        phi_value=np.array( [[0.]*self.ny]*self.nx )

        phi_value=np.full((self.nx,self.ny),0.01)
        for k in range(20):
            x_rand=random.randint(4,self.nx-4)
            y_rand=random.randint(4,self.ny-4)
            r_rand=random.randint(3,5)
            print("random phi",x_rand,"  ",y_rand," ",r_rand)
            for i in range(self.nx):
                for j in range(self.ny):
                    d=np.sqrt((i-x_rand)**2+(j-y_rand)**2)
                    if d<=r_rand :
                        phi_value[i,j]=self.idense*0.5

        return phi_value
    def phi_random2(self):
        import numpy as np
        import random
        random.seed(10)
        phi_value=np.array( [[0.]*self.ny]*self.nx )

        phi_value=np.full((self.nx,self.ny),self.idil)
        for k in range(20):
            x_rand=random.randint(4,self.nx-4)
            y_rand=random.randint(4,self.ny-4)
            r_rand=random.randint(3,5)
            print("random phi",x_rand,"  ",y_rand," ",r_rand)
            for i in range(self.nx):
                for j in range(self.ny):
                    d=np.sqrt((i-x_rand)**2+(j-y_rand)**2)
                    if d<=r_rand :
                        phi_value[i,j]=self.idense

        return phi_value

    def phi_circle00(self):
        import numpy as np
        phi_value=np.array( [[0.]*self.ny]*self.nx )
        for i in range(self.nx):
            for j in range(self.ny):
                d=np.sqrt((i-60)**2+(j-60)**2)
                if d<15 :
                    phi_value[i,j]=self.idense*0.9
                else:
                    phi_value[i,j]=self.idil

        return phi_value
