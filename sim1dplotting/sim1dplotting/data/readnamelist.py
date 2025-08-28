class namelist:
    
    def __init__(self, configuration_exp = None, \
                       configuration_rheo = None, 
                       configuration_fig = None,
                       configuration_time = None):
        
        #--- Reading the Model Parameters ---#
        
        #------ Experiment ------#
        self.expno = [int(x) for x in configuration_exp['expno'].split(',')]
        self.savevar = int(configuration_exp['savevar'])
        self.plotfields = int(configuration_exp['plot_fields'])
        self.read_all = int(configuration_exp['read_all'])
        self.plotsingle = int(configuration_exp['plot_single'])
        self.plotresidual = int(configuration_exp['residual_analysis'])
        self.dissipation = int(configuration_exp['energy_dissipation'])
        self.mechanical_energy = int(configuration_exp['mech_energy'])
        self.maxvelocities = int(configuration_exp['maxvelocities'])

        
        #------ Rheology ------#
        self.muphi      = int(configuration_rheo['muphi'])
        self.mu0        = float(configuration_rheo['mu0'])
        self.muinf      = float(configuration_rheo['muinf'])
        self.mub        = [float(x) for x in configuration_rheo['mub'].split(',')]
        self.dx         = float(configuration_rheo['dx'])*1e3
        self.Nx         = float(configuration_rheo['Nx'])
        self.adv        = int(configuration_rheo['advection_scheme'])
        self.imex       = int(configuration_rheo['IMEX'])
        self.solv       = int(configuration_rheo['solver'])
        self.maxcapping = [float(x) for x in configuration_rheo['maxcapping'].split(',')]
        
        #------ Figures ------#
        self.outputdir = str(configuration_fig['outputdir'])
        self.figdir    = str(configuration_fig['figdir'])
        
        #------ Time -------#
        self.dt = float(configuration_time['dt'])
        self.dt_read = float(configuration_time['dt_read'])
        self.startk = int(configuration_time['start_k'])
        self.nstep = int(configuration_time['nstep'])
        
        