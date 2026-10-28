class Initialize:


    """
    This class is to call all needed parameters, and to construct components needed 
    before running DDFT evolution
    """

    def __init__(self,args):

        import os
        import fipy as fp
        import argparse
        import Protein_RPA.utils.free_energy as f_en
        from Protein_RPA import seq_list as sl
        import initial_density as density
        import chempot_table as chempot
        import numpy as np
        from Protein_RPA.utils.input_parse import input_parse
        import pickle
        from density_class import initial_density
        import chempot_table as chempot #very similar inputs as main DDFT function
        self.input_parameters = input_parse(args.i);
        input_parameters = input_parse(args.i);
        seq_name = "Ddx4_N1"  # should be part of input_parameters

        #finally, a list of input parameters for args
        ehs0=input_parameters['enthalpy']
        ehs1=input_parameters['entropy']
        ehs=np.array([ehs0, ehs1])
        ef=False
        sig, N, the_seq = sl.get_the_charge(seq_name,4.95)
        self.HP=f_en.RPAFH(sig, ehs=ehs, epsfun=ef)
        self.FE=f_en.free_energy(self.HP)

        #Salt Concentration
        #Gathering all input parameters
        self.nx = int(input_parameters['nx'])
        self.ny = int(input_parameters['ny'])
        self.dx = input_parameters['dx']
        self.dimension = int(input_parameters['dimension'])
        self.plot_flag = bool(input_parameters['plot_flag'])

        self.n_capture=input_parameters['n_capture']
        self.phis = input_parameters['phis'];


        #Surface tension/ influence parameter 
        self.kappa=input_parameters['kappa'];
        self.dt = input_parameters['dt'];
        self.dt_max = input_parameters['dt_max'];
        self.dt_min = input_parameters['dt_min'];
        self.tolerance = input_parameters['tolerance'];
        self.total_steps = int(input_parameters['total_steps']);
        self.checkpoint = int(input_parameters['checkpoint']);
        if 'text_log' in input_parameters.keys():
            self.text_log = int(input_parameters['text_log'])
        else:
            self.text_log = checkpoint;
        self.duration = input_parameters['duration'];

        self.time_step = fp.Variable(self.dt)

        self.t = fp.Variable(0.0)
        self.dt = input_parameters['dt'];

        self.initial_density=input_parameters['initial_density']
        if self.dimension==1:
            self.mesh=fp.Grid1D(nx=self.nx)
            self.mesh=mesh-float(self.nx)*self.dx*0.5

        if self.dimension==2:
            self.mesh = fp.Grid2D(nx=self.nx, ny=self.ny, dx=self.dx, dy=self.dx)
            self.mesh=self.mesh-float(self.nx)*self.dx*0.5

        elif self.dimension==3:
            self.mesh = fp.Grid3D(nx=self.nx, ny=self.ny,nz=self.nx, dx=self.dx, dy=self.dx,dz=self.dx)
            self.mesh=self.mesh-float(self.nx)*self.dx*0.5


        self.phi_inf=input_parameters['phi_boundary']
    
        #inverse of temperature is u
        self.u=input_parameters['temperature_inverse'];
        self.phi=fp.CellVariable(mesh=self.mesh, name=r'$\phi_{monomer}$', hasOld=True)

        self.phi_value=np.array( [[0.]*self.ny]*self.nx )
        self.xi_value=np.array( [[0.]*self.ny]*self.nx )


        self.output_dir=args.o
        self.state_dir=args.s
        #self.state_file=args.s
        self.elapsed = 0.0
        self.steps0 = 0
        self.t0=0

        #Read initial density
        mydensity=initial_density(self.nx,self.ny)
        func_density=getattr(mydensity,self.initial_density)

        #Read state file of initial density 
        try:
            #self.filename=os.getcwd()+'/'+self.state_dir+'state'+'.pickle'
            self.filename=self.state_dir+'state'+'.pickle'
            self.file = open(self.filename, 'rb')
            self.steps_t,self.t_t,self.input_parameters_t,self.phi_value_t,self.xi_value_t=pickle.load(self.file)
            print('state file is available')
            print(self.filename)

            cond1=int(self.input_parameters_t['nx'])==int(self.nx)
            cond2=int(self.input_parameters_t['phis'])==int(self.phis)
            cond3=int(self.input_parameters_t['temperature_inverse'])==int(self.u)
            cond4=self.input_parameters_t['kappa']==self.kappa
            rules=[cond1,cond2, cond3, cond4]
            print("rules ",rules)
            print(self.input_parameters_t['nx'])
            if ( all(rules) ):
                self.steps0=self.steps_t
                self.t0=self.t_t
                self.total_steps+=self.steps_t
                self.phi_value=self.phi_value_t
                self.initial_state=True
                print(' initial average of density is ', np.mean(self.phi_value))
            self.file.close()
        except:
            print('state file is NOT available')
            print(self.filename)
            
            #self.phi_value=density.phi_circle2_6(self.nx,self.ny)
            self.phi_value=func_density()
            self.initial_state=False

        #Below are variables needed in DDFT equation

        self.mu_spl=chempot.lookup_table(self.FE,self.u,self.phis) 
        self.dmu_spl=chempot.dmu_table(self.FE,self.u,self.phis)

        
        self.mu_value=np.array( [[0.]*self.ny]*self.nx )
        self.dmu_value=np.array( [[0.]*self.ny]*self.nx )

        self.mu_value=self.mu_spl(self.phi_value)
        self.dmu_value=self.dmu_spl(self.phi_value)

        #The chemical potential is transferred into variable xi, which is a variabe in fipy data structure,

        self.phi.setValue(self.phi_value.flatten())
        self.xi = fp.CellVariable(mesh=self.mesh,hasOld=True)
        self.xi.setValue(self.mu_value.flatten())


        self.phi.updateOld()


        #Some temporary variables used in solving DDFT equation

        self.xi_now = fp.CellVariable(mesh=self.mesh)
        self.dxi_now= fp.CellVariable(mesh=self.mesh)
        self.phi_now=fp.CellVariable(mesh=self.mesh)

        #Boundary value
        self.xi_boundary=self.mu_spl(self.phi_inf)
    def calc_chemical(self):
        
        self.mu_value=self.mu_spl(self.phi_value)
        self.dmu_value=self.dmu_spl(self.phi_value)

        self.xi_now.setValue(self.mu_value.flatten())
        self.dxi_now.setValue(self.dmu_value.flatten())
        self.phi_now.setValue(self.phi.value)


    def global_phi(self):
        self.phi_value=self.phi.globalValue.reshape(self.nx,self.ny)        
        self.xi_value=self.xi.globalValue.reshape(self.nx,self.ny)        

    def update_parameters(self,update):
        return 0


