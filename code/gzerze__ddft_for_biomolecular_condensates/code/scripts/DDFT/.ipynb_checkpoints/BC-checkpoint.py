class Boundary:
    def __init__(self,System):
        
        import numpy as np
        import fipy as fp
        from initialize import Initialize


        #self.option=option
    
        self.mesh=System.mesh
        self.input_parameters=System.input_parameters

        try:
            self.option=System.input_parameters['boundary_type']
        except:
            self.option="diffusion_proportional"

        try:
            self.h_mul=-System.input_parameters['h_mul']
        except:
            self.h_mul=-1e-3


        self.dx=System.input_parameters['dx']
        self.phi_inf=System.input_parameters['phi_boundary']

        

        #self.dPf = fp.FaceVariable(mesh=self.mesh,value=[0.5*self.dx,0.5*self.dx],rank=1 )
        self.n = self.mesh.faceNormals

        self.dPf = 0.5*self.dx*self.n
        
        self.mask=self.mesh.exteriorFaces

        self.Gamma0=1.0
        #self.a = fp.FaceVariable(mesh=self.mesh, value=[self.h_mul,self.h_mul], rank=1)
        self.a = self.h_mul*self.n
        self.b = fp.FaceVariable(mesh=self.mesh, value= 1, rank=0)
        self.g = fp.FaceVariable(mesh=self.mesh, value= self.h_mul*self.phi_inf, rank=0)
        #self.RobinCoeff_diff = (self.mask * self.Gamma0 * self.n / (-self.dPf.dot(self.a) +self.b))

        
        self.coeff_diff=0.0
        self.coeff_impl=0.0
        self.alpha=0
        self.beta=0
       
        self.xi_now=System.xi_now
        self.xi_boundary=System.xi_boundary 
    def parameters(self):
        import fipy as fp
        import numpy as np
        """"
        Option of Boundary Condition are :
        1. "zero_flux":
        The simplest BC, with no flux of protein going out of/in to exterior boundary : normal.dot.faceGrad(phi)=0
        2. "diffusion_pde":
        This is used when the PDE of ddFT is written as in diffusion equation. This is not going to be used in the future.
        3. "diffusion_proportional"
        Mathematically this is written as : n.dot(nabla(phi_s))=-h_mul(phi_s-phi_inf)
        phi_s: density at cell centers that are facing exterior boundary.
        phi_inf : density of reservoir outside the exterior boundary
        h_mul : parameter that need to be provided
        """
        match self.option: 
        # Boundary Condition 
        # The boundary condition implemented is Robin Boundary Condition
        # There are some modifications because some of subroutine in the manual does not work:
        # dPf is set to rank 1 ( not 0)
        # Use a.dot(n), instead of n.dot(a). The former create a rank 0 term

            


            case "zero_flux":

                self.alpha=0
                self.beta=0

                self.coeff_impl=0.0
                self.coeff_diff=0.0

                self.RobinCoeff=0.0
                self.g=0.0

                self.chem_grad_surf=fp.FaceVariable(mesh=self.mesh,value=(self.xi_now-self.xi_boundary).faceValue/(0.5*self.dx),rank=0)
            case "diffusion_pde":
      
                self.coeff_diff=1.0
                self.coeff_impl=self.a
                self.alpha=1.0
                self.beta=-1.0

            
                self.RobinCoeff = (self.mask * self.Gamma0 * self.n / (-self.dPf.dot(self.a) + self.b))



            case "diffusion_proportional":

                #Below is the initial implementation for RobinCoeff, but is replaced with a better one
                #RobinCoeff_conv=mask*n*(xi_now.faceGrad.dot(n))/(-dPf.dot(a)+b)

                #Below is the newer implementation of RobinCoeff, the improvement is done by first calculating normal gradient  
                #of chemical potential at surface
                #The following is to transform cell variable to face variable

                self.chem_grad_surf=fp.FaceVariable(mesh=self.mesh,value=-(self.xi_now-self.xi_boundary).faceValue/(0.5*self.dx),rank=0)
                #self.RobinCoeff=self.mask*self.n*self.chem_grad_surf/(-self.dPf.dot(self.a)+self.b)
                self.RobinCoeff=self.mask*self.n*self.chem_grad_surf/(-0.5*self.dx*self.h_mul+self.b)
                
                self.c_conv=fp.FaceVariable(mesh=self.mesh,value=0.5*self.dx,rank=0)

                self.coeff_diff=self.c_conv
                self.coeff_impl=self.b
                self.alpha=-1.0
                self.beta=1.0
        
            
        return self.alpha,self.beta,self.coeff_diff,self.coeff_impl,self.g,self.RobinCoeff

    def update_parameters(self,System):
        import fipy as fp
        import numpy as np
        from initialize import Initialize
        #self.xi_now=xi_now
        match self.option:
            case "zero_flux":
                self.chem_grad_surf.setValue( -0*(System.xi_now-System.xi_boundary).faceValue/(0.5*self.dx))
                self.RobinCoeff=self.mask*self.n*self.chem_grad_surf/(-0.5*self.dx*self.h_mul+self.b)
            case "diffusion_pde":
                print("no update")

            case "diffusion_proportional":
                self.chem_grad_surf.setValue( -(System.xi_now-System.xi_boundary).faceValue/(0.5*self.dx))
                self.RobinCoeff=self.mask*self.n*self.chem_grad_surf/(-0.5*self.dx*self.h_mul+self.b)

