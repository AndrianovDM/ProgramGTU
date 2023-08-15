import pandas as pd
import numpy as np

def fuel_table(temperature_C, temperature_K, u, h, s_0, Cp, Cv, k):  
     
    data = pd.DataFrame({ 't':pd.Series(np.round((temperature_C),2),index=range(1,len(temperature_C)+1,1)),
                              'T':pd.Series(np.round((temperature_K),2),index=range(1,len(temperature_K)+1,1)),
                              'u':pd.Series(np.round((u),1),index=range(1,len(u)+1,1)),
                              'h':pd.Series(np.round((h),1),index=range(1,len(h)+1,1)),
                              's°':pd.Series(np.round((s_0),3),index=range(1,len(s_0)+1,1)),
                            # #   'lg(pi_0)':pd.Series(self.lg_Pi,index=range(1,len(self.lg_Pi)+1,1)),
                            # #   'pi_0':pd.Series(self._Pi,index=range(1,len(self._Pi)+1,1)),
                              'C_p':pd.Series(np.round((Cp),4),index=range(1,len(Cp)+1,1)),
                              'C_v':pd.Series(np.round((Cv),4),index=range(1,len(Cv)+1,1)),
                              'k':pd.Series(np.round((k),4),index=range(1,len(k)+1,1))},
                            columns=['t','T','u','h','s°','C_p','C_v','k'])
    return data