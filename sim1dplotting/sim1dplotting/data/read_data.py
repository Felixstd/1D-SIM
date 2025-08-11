import numpy as np

def read_data(expno, dt, dx, solv, IMEX, adv, dates, outputdir, MuPhi = False, Dissipation = False):
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
    data_dict = {
        'divergence_dates': [],
        'h_dates': [],
        'A_dates': [],
        'u_dates':[],
        # 'b_dates':[], 
        'eta_dates':[], 
        'zeta_dates':[]

    }
    
    if Dissipation:
        data_dict.update({
            'Wsigma_dates':[]
        }
        )
    
    # if MuPhi:
    #     data_dict.update({
    #         'muI_dates': [],
    #         'phi_dates': [],
    #         'I_dates': [],
    #         'e_II_dates':[], 
    #         'Pmax_dates':[], 
    #         'Peq_dates':[]
    #     })
    
    for k, date in enumerate(dates, start=1):
        # List of file prefixes and associated keys
        files_info = [
            ('div', 'divergence_dates'),
            ('h', 'h_dates'),
            ('A', 'A_dates'),
            ('u', 'u_dates'),
            # ('b', 'b_dates'), 
            ('zeta', 'zeta_dates'), 
            ('eta', 'eta_dates')
        

        ]
        
        for prefix, key in files_info:
           
            # filename = f"{outputdir}{prefix}{'_'}{dt:05}{'s_'}{dx:03}{'km_'}{'solv'}{solv}{'_IMEX'}{IMEX}{'_adv'}{adv}{'_BDF20_ts'}{date:08}{'.'}{expno}"
            filename = f"{outputdir}{prefix}{'_'}{dt:05}{'s_'}{dx:06}{'km_'}{'solv'}{solv}{'_IMEX'}{IMEX}{'_adv'}{adv}{'_BDF20_ts'}{date:08}{'.'}{expno:02d}"
            # filename = f"{outputdir}{prefix}{date}{('_k{:04d}'.format(k) + '.' + expno) if 'sig' in prefix or prefix in ['div', 'shear'] else '.' + expno}"
            # print(filename)
            data_dict[key].append(np.loadtxt(filename, dtype=None))
        
        if Dissipation:
            for prefix, key in [('Wdissip', 'Wsigma_dates')]:
                filename = f"{outputdir}{prefix}{'_'}{dt:05}{'s_'}{dx:06}{'km_'}{'solv'}{solv}{'_IMEX'}{IMEX}{'_adv'}{adv}{'_BDF20_ts'}{date:08}{'.'}{expno:02d}"
                data_dict[key].append(np.loadtxt(filename, dtype=None))
        
    
        k+=1
    return data_dict

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



    
    
    