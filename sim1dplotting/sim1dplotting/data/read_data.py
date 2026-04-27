import numpy as np
import re

def load_data(data_dict, files_info, date, expno, dt, dx, solv, IMEX, adv,outputdir):
    
    
        for prefix, key in files_info:
            try:
                
                filename = f"{outputdir}{prefix}{'_'}{dt:05}{'s_'}{dx:06}{'km_'}{'solv'}{solv}{'_IMEX'}{IMEX}{'_adv'}{adv}{'_BDF20_ts'}{date:09}{'.'}{expno:02d}"

                data_dict[key].append(np.loadtxt(filename, dtype=None))

            
            except FileNotFoundError:
                
                filename = f"{outputdir}{prefix}{'_'}{dt:05}{'s_'}{dx:06}{'km_'}{'solv'}{solv}{'_IMEX'}{IMEX}{'_adv'}{adv}{'_BDF20_ts'}{date:08}{'.'}{expno:02d}"
                data_dict[key].append(np.loadtxt(filename, dtype=None))
        # return data_dict


def read_data(expno, dt, dx, solv, IMEX, adv, dates, outputdir, Dissipation = False, Energy = False, GC = False):
    """
    This function reads the ouput, for the specified dates and experiment number (expno) from the McGIll-SIM and puts
    them in a dictionnary. 
    
    It will read the following variables: 


    Args:
        expno (str): experiment number
        dates (list of str): list of the output dates
        outputdir (str): output directory
        MuPhi (bool, optional): If the expno was computed with MuPhi. Defaults to True.
        
        
    Returns:
        data_dict
    """
        # Initialize lists
    data_dict_base = {
        'divergence_dates': [],
        'h_dates': [],
        'A_dates': [],
        'u_dates':[],
        'sigma_dates':[], 
        'signorm_dates':[], 
        # 'b_dates':[], 
        'eta_dates':[], 
        'zeta_dates':[]

    }
    
    if Dissipation:
        data_dict_base.update({
            'Wsigma_dates':[]
        }
        )
    
    if Energy:
        data_dict_energy = {
            
            'Plat_dates':   [],
            'Pfric_S_dates':[],
            'Pfric_R_dates':[], 
            'Ph_dates':     [],
            'Pw_dates':     [], 
            'Ppot_dates':   []     

        }
        
    if GC: 
        data_dict_base.update({
            'P_dates':[],
            'muI_dates':[]
        }
        )
    
    for k, date in enumerate(dates, start=1):
        # List of file prefixes and associated keys
        print('Reading Base, date: ', date,'expno: ', expno)
            
        files_info = [
            ('div', 'divergence_dates'),
            ('h', 'h_dates'),
            ('A', 'A_dates'),
            ('u', 'u_dates'),
            ('sig', 'sigma_dates'),
            ('sigN', 'signorm_dates'),
            # ('b', 'b_dates'), 
            ('zeta', 'zeta_dates'), 
            ('eta', 'eta_dates')
        ]
        
        load_data(data_dict_base, files_info, date, expno, dt, dx, solv, IMEX, adv, outputdir)
            
        
        if Dissipation:
            for prefix, key in [('Wdissip', 'Wsigma_dates')]:
                filename = f"{outputdir}{prefix}{'_'}{dt:05}{'s_'}{dx:06}{'km_'}{'solv'}{solv}{'_IMEX'}{IMEX}{'_adv'}{adv}{'_BDF20_ts'}{date:09}{'.'}{expno:02d}"
                data_dict_base[key].append(np.loadtxt(filename, dtype=None))
                
        if GC:
            files_info = [('P', 'P_dates'), 
                          ('Mu', 'muI_dates')]
            load_data(data_dict_base, files_info, date, expno, dt, dx, solv, IMEX, adv, outputdir)
        

        if Energy:
            print('Reading Energy: ', date, expno)
            files_energy = [
            ( 'Plat', 'Plat_dates'), 
            ('Pfric_s', 'Pfric_S_dates'),
            ('Pfric_R', 'Pfric_R_dates'), 
            ('Ph', 'Ph_dates'),
            ('Pw', 'Pw_dates'), 
            ('Ppot', 'Ppot_dates')]
            
            
            load_data(data_dict_energy, files_energy, date, expno, dt, dx, solv, IMEX, adv, outputdir)
            
            
    if Energy:
        return data_dict_base, data_dict_energy
    else: 
        return data_dict_base

def read_maxvelocities(namefile):
    """
    This function is used to read the maximum and minimum velocities at each time step by reading the 
    post files. 
    
    In order to obtain the files to read, use the following command 

        grep 'min, max u' postfile > namefile.out
        
        this works also awk '/min, max u ice/ { print  $(NF-7), $(NF-6),  $(NF-5), $(NF-4), $(NF-3), $(NF-2), $(NF-1), $(NF) }'  output_mu_test_64
        
        
    Args:
        namefile (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    number_columns = []
    with open(namefile, "r") as f:
        for line in f:
            parts = line.split("=")[-1].split()
            number_columns.append([float(val) for val in parts])


    number_columns = np.asarray(number_columns)        # shape (N, 3)

    # grab column 1 (max‑values) and a matching x‑axis
    max_velocities = number_columns[:, 1]                           # shape (N,)
    min_velocities = number_columns[:, 0]
    tstep = np.arange(1, len(number_columns)+1)


    return max_velocities, min_velocities, tstep


def read_avg_energy(expno, dt, dx, solv, IMEX, adv, date, outputdir):

    data_dict_energy_avg = {
            
            'Plat_Tavg':   [],
            'Pfric_S_Tavg':[],
            'Pfric_R_Tavg':[], 
            'Pfric_S_visc_Tavg':[],
            'Pfric_R_visc_Tavg':[], 
            'Pfric_S_plas_Tavg':[],
            'Pfric_R_plas_Tavg':[], 
            'P_fric_R_pZ_vE_Tavg':[], 
            'P_fric_R_vZ_pE_Tavg':[],
            'Ph_Tavg':     [],
            'Pw_Tavg':     [], 
            'Ppot_Tavg':   []     

        }
    
    files_energy = [
            ( 'Pl_t', 'Plat_Tavg'), 
            ('PfR_t', 'Pfric_R_Tavg'),
            ('PfS_t', 'Pfric_S_Tavg'), 
            ('PfRv_t', 'Pfric_R_visc_Tavg'),
            ('PfSv_t', 'Pfric_S_visc_Tavg'), 
            ('PfRp_t', 'Pfric_R_plas_Tavg'),
            ('PfSp_t', 'Pfric_S_plas_Tavg'), 
            ('PfRvzpe_t', 'P_fric_R_vZ_pE_Tavg'), 
            ('PfRpzve_t', 'P_fric_R_pZ_vE_Tavg'), 
            ('Ph_t', 'Ph_Tavg'),
            ('Pw_t',   'Pw_Tavg'), 
            ('Pp_t', 'Ppot_Tavg')]


    load_data(data_dict_energy_avg, files_energy, date, expno, dt, dx, solv, IMEX, adv, outputdir)

    return data_dict_energy_avg

def read_visc_plas_file(namefile):
    
    viscous = []
    plastic = []
    
    with open(namefile, "r") as f:
        for line in f:
            # Find all numbers
            numbers = re.findall(r"\d+", line)
            if len(numbers) >= 2:
                plastic.append(int(numbers[0]))
                viscous.append(int(numbers[1]))
    
    return np.asarray(viscous), np.asarray(plastic)
    
    
    
    
    
    