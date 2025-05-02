import numpy as np

def read_data(expno, dt, dx, solv, IMEX, adv, dates, outputdir, MuPhi = True):
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
        'h_dates': [],
        'A_dates': [],
        'u_dates':[]

    }
    
    # if MuPhi:
    #     data_dict.update({
    #         'muI_dates': [],
    #         'phi_dates': [],
    #         'I_dates': [],
    #         'e_II_dates':[], 
    #         'Pmax_dates':[], 
    #         'Peq_dates':[]
    #     })
    
    for k, date in enumerate(dates, start=0):
        print('Reading: ', dates[k])
        # List of file prefixes and associated keys
        files_info = [
            ('h', 'h_dates'),
            ('A', 'A_dates'),
            ('u', 'u_dates'),

        ]
        
        for prefix, key in files_info:
            filename = f"{outputdir}{prefix}{'_'}{dt:05}{'km'}{dx:05}{'_ts'}{date:06}{'.'}{expno}"
            # filename = f"{outputdir}{prefix}{date}{('_k{:04d}'.format(k) + '.' + expno) if 'sig' in prefix or prefix in ['div', 'shear'] else '.' + expno}"
            
            data_dict[key].append(np.loadtxt(filename, dtype=None))
        
        # if MuPhi:
        #     for prefix, key in [('mu_I', 'muI_dates'), ('I', 'I_dates'), ('phi_I', 'phi_dates'), ('e_II', 'e_II_dates'), ('Pmax', 'Pmax_dates'), ('Peq', 'Peq_dates')]:
        #         filename = f"{outputdir}{prefix}{date}.{expno}"
        #         data_dict[key].append(np.loadtxt(filename, dtype=None))
        
    
        k+=1
    return data_dict




    
    
    