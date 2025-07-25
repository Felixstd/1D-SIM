import numpy as np

class TimeUtility:
    
    def __init__(self, configuration_time = None, 
                 configuration_fig = None):
        
        self.dt           = int(configuration_time['dt'])
        self.outputdir     = str(configuration_fig['outputdir'])
        
        
    def read_time(self, expno = None): 
        
        tstep = []
        with open(self.outputdir+'time_run.{:02d}'.format(expno)) as f:
            lines = f.readlines()
            for line in lines:
                tstep.append(int(line))
                
        tstep = np.asarray(tstep)
                
        time = tstep*self.dt/(60*60)

        return tstep[:-1], time[:-1]