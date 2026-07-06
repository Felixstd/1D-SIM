import numpy as np
import sim1dplotting.utils.analysis_utils as utils

from sim1dplotting.data import read_data


def load_experiment_data(Parameters, solver, exp, ind, tstep):
    
    #--------------------------------------------------------------
    # Default variables
    #--------------------------------------------------------------
    
    datadict, datadict_energy = read_data.read_data(
                                exp, 
                                int(Parameters.dt_read), 
                                int(Parameters.dx[ind]/1e3), 
                                solver, 
                                Parameters.imex, 
                                Parameters.adv, 
                                tstep, 
                                Parameters.outputdir, 
                                Dissipation = Parameters.dissipation, 
                                Energy      = Parameters.mechanical_energy, 
                                GC          = Parameters.GC, 
                                strength    = Parameters.strength
        )
    
    data = {"base": datadict, "energy": datadict_energy}

    if Parameters.mechanical_energy:
    #--------------------------------------------------------------
    # Energy related variables
    #--------------------------------------------------------------
        data['energy_avg'] = read_data.read_avg_energy(
                                        exp, 
                                        int(Parameters.dt_read), 
                                        int(Parameters.dx[ind]/1e3), 
                                        Parameters.solv, 
                                        Parameters.imex, 
                                        Parameters.adv, 
                                        Parameters.nstep, 
                                        Parameters.outputdir)
   
        data['num_visc'], data['num_plas'] = read_data.read_visc_plas_file('outputs_VP/num_visc_plas_{}.out'.format(exp))
        
        # number_mixed_viszeta, number_mixed_viseta = read_data.read_visc_plas_file('outputs_VP/num_mixed_visc_plas_{}.out'.format(exp))

    Parameters.dx[ind] = utils.map_dx(Parameters.Nx[ind],Parameters.dx[ind])
    
    x_exp = np.linspace(0, Parameters.Nx[ind]*Parameters.dx[ind], int(Parameters.Nx[ind]))	

    data["x"] = x_exp
    data["A"] = datadict['A_dates']
    data["h"] = datadict['h_dates']
    data["u"] = datadict['u_dates']
    
    return data 
    
    