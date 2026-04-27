import numpy as np

def init_condition(type_u, type_h, Nx, random = False):
    """
    Function which treats the initial conditions of the 
    simulations. 

    Args:
        type_u (_type_): _description_
        type_h (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    hinit = np.zeros(Nx+1)
    uinit = np.zeros(Nx+1)
    
    x = np.arange(0, Nx+1,1)
    
    if random:
        rng = np.random.default_rng(12345)
        noise = rng.random(size=Nx+1, dtype=np.float64)*1e-2
        
    else:
        noise = np.zeros_like(uinit)
    
    if type_h == 'Step':
        # if Nx == 801:
        #     hinit = np.exp(-((x-Nx/2)**2/10000)**5) + noise
        if (Nx == 401) or (Nx == 801):
            # hinit = np.exp(-((x-Nx/2)**2/2500)**5) + noise
            hinit = np.exp(-((x-Nx/2)**2/1000)**5) + noise
        elif Nx == 101:
            hinit = np.exp(-((x-Nx/2)**2/200)**5) + noise
        elif Nx == 4001:
            hinit = np.exp(-((x-Nx/2)**2/250000)**5) + noise
        elif Nx == 8001:
            hinit = np.exp(-((x-Nx/2)**2/100000)**5) + noise
    
    if type_h == 'Gray':
        if Nx == 401:
            # hinit = np.exp(-((x-Nx/2)**2/2500)**5) + noise
            hinit = (3000-(x-Nx/2)**2)**(1/20)
            hinit[np.isnan(hinit)] = 0
            # np.exp(-((x-Nx/2)**2/5000)**5) + noise
        
    # elif type_h == 'Gray':
    #     hinit = (1 - (x-Nx/2)**2)**(1/20) + noise
    #     hinit[np.isnan(hinit)] = 0
    
    if type_u == 'Step':
        uinit = np.exp(-((x-Nx/2)**2/10000)**5) + noise
    
    elif type_u == 'Gray':
        uinit += noise
        # if Nx == 801:
        uinit[0:int(Nx/4)] = 0
        uinit[int(3*Nx/4)] = 0
            # uinit[int(Nx/4):int(3*Nx/4)] = x[int(Nx/4):int(3*Nx/4)]/1e4
        
        
        uinit[np.where(hinit>0)] = x[np.where(hinit>0)]/2e3
        
    
    elif type_u == 'Gray_smooth':
        idx = np.where(hinit>0.2)[0]
        if Nx == 401:
            uinit = x/1e3*(np.tanh((x - idx[0])/1e1) - np.tanh((x - idx[-1])/1e1))
            
        elif Nx == 101:
            uinit = x/1e3*(np.tanh((x - idx[0])/1e0) - np.tanh((x - idx[-1])/1e0))
            
            
    elif type_u == 'GrayDiv': 
        idx = np.where(hinit>0.2)[0]
        if (Nx == 401):
            # left = -x /1e3 * (np.tanh((x - idx[0])/1e1) - np.tanh((x - idx[int(len(idx)/2)])/1e1))
            # right = x /1e3 * (np.tanh((x - idx[int(len(idx)/2)])/1e1) - np.tanh((x - idx[-1])/1e1))
            # uinit = left+right
            # uinit_1 = np.linspace(-2.5,2.5, Nx+1)
            
            uinit_1 = (x - Nx/2)/20
            # uinit_1[hinit < 0.2] = 0
            
            width = 0.1
            mask_tanh = 0.4 * (1 + np.tanh((hinit - 0.5) / width))
            uinit = uinit_1 * mask_tanh 
        
        
        elif (Nx == 801):
            uinit_1 = (x - Nx/2)/10
            # uinit_1[hinit < 0.2] = 0
            
            width = 0.1
            mask_tanh = 0.4 * (1 + np.tanh((hinit - 0.5) / width))
            uinit = uinit_1 * mask_tanh 
            
        if Nx == 4001:
            uinit_1 = (x - Nx/2)/800
            # uinit_1[hinit < 0.2] = 0
            
            width = 0.1
            mask_tanh = 0.4 * (1 + np.tanh((hinit - 0.5) / width))
            uinit = uinit_1 * mask_tanh 
            # uinit[hinit < 1e-3] = 0
        # uinit[int(2*Nx/5):int(3*Nx/5)] = x[int(2*Nx/5):int(3*Nx/5)]/1e4
    elif type_u == 'zero':
        uinit = np.zeros(Nx+1) + noise
        

        
    uinit[0] = 0
    uinit[-1] = 0
    
    hinit[0] = 0
    hinit[-1] = 0
    

    return uinit, hinit