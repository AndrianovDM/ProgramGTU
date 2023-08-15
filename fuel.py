import streamlit as st
from math import *
from fuelFunction import *
from fuelPlot import *
from fuelTable import *


def properties(temperature_C, coefficient, mu, h_k, s_k):
    
    R_const = 8.31451
    temperature_K = [round((temperature_C[i] + 273.15), 3) for i in range(len(temperature_C))]
    tau = [temperature_K[i] / 1000 for i in range(len(temperature_C))]
    R = R_const / mu
    Cp = [R * (((coefficient[0] * (tau[i])**0) + (coefficient[1] * (tau[i])**1) +
                (coefficient[2] * (tau[i])**2) + (coefficient[3] * (tau[i])**3) +
                (coefficient[4] * (tau[i])**4) + (coefficient[5] * (tau[i])**5) +
                (coefficient[6] * (tau[i])**6)) + ((coefficient[7] * (1 / tau[i])**1) +
                (coefficient[8] * (1 / tau[i])**2) + (coefficient[9] * (1 / tau[i])**3) +
                (coefficient[10] * (1 / tau[i])**4) + (coefficient[11] * (1 / tau[i])**5)+
                (coefficient[12] * (1 / tau[i])**6))) for i in range(len(tau))]

    h = [R * 1000 * ((((coefficient[0] / 1) * tau[i]**1) + ((coefficient[1] / 2) * tau[i]**2) +
                ((coefficient[2] / 3) * tau[i]**3) + ((coefficient[3] / 4) * tau[i]**4) +
                ((coefficient[4] / 5) * tau[i]**5) + ((coefficient[5] / 6) * tau[i]**6) +
                ((coefficient[6] / 7) * tau[i]**7)) + coefficient[7] * log(tau[i]) +
                (((coefficient[8] / (-1)) * (1 / tau[i])**1) + ((coefficient[9] / (-2)) * (1 / tau[i])**2) +
                ((coefficient[10] / (-3)) * (1 / tau[i])**3) + ((coefficient[11] / (-4)) * (1 / tau[i])**4) +
                ((coefficient[12] / (-5)) * (1 / tau[i])**5))) + h_k * R for i in range(len(temperature_K))]

    s_0 = [R * (coefficient[0] * log(tau[i]) + (((coefficient[1] * tau[i]**1) / 1) +
            ((coefficient[2] * tau[i]**2) / 2) + ((coefficient[3] * tau[i]**3) / 3) +
            ((coefficient[4] * tau[i]**4) / 4) + ((coefficient[5] * tau[i]**5) / 5) +
            ((coefficient[6] * tau[i]**6) / 6)) + (((coefficient[7] * (1 / tau[i])**(1)) / (-1)) +
            ((coefficient[8] * (1 / tau[i])**(2)) / (-2)) + ((coefficient[9] * (1 / tau[i])**(3)) / (-3)) + 
            ((coefficient[10] * (1 / tau[i])**(4)) / (-4)) + ((coefficient[11] * (1 / tau[i])**(5)) / (-5)) +
            ((coefficient[12] * (1 / tau[i])**(6)) / (-6)))) + s_k * R for i in range(len(temperature_K))]

    u=[h[i] - R * temperature_K[i] for i in range(len(temperature_K))]
    Cv=[Cp[i] - R for i in range(len(temperature_K))]
    k= [Cp[i] / Cv[i] for i in range(len(temperature_K))]
    return {'temperature_C': temperature_C, 'temperature_K': temperature_K, 'temperature_K': temperature_K, 'u':u, 'h':h, 's_0':s_0, 'Cp':Cp, 'Cv':Cv, 'k':k, 'R':R, 'mu':mu}

@st.cache
def gas(_H_2S, _CO_2, _O_2, _CO, _H_2, _CH_2, _CH_4, _C_2H_4, _C_2H_6, _C_3H_8, _C_4H_10, temperature_C):
        
        _N_2_atm = 100-_H_2S - _CO_2 - _O_2 - _CO - _H_2 - _CH_2 - _CH_4 - _C_2H_4 - _C_2H_6 - _C_3H_8 - _C_4H_10
        temperature_K = [round((temperature_C[i]+273.15),3) for i in range(len(temperature_C))]
        R_const = 8.31451
        excess_air_ratio_in_flue_gases = 1
        d_t = 47.5
        d_air = 8

        heat_of_combustion = 358.2 * _CH_4 + 637.46 * _C_2H_6 + 860.05 * _C_3H_8 + 107.98 * _H_2 + 126.36 * _CO
        stoichiometric_air_flow = ((1 / 21) * (0.5 * _H_2 + 0.5 * _CO + 2 * _CH_2 + 2 * _CH_4 + 3.5 * _C_2H_6 + 5 * _C_3H_8 + 6.5 * _C_4H_10 + 1.5 * _H_2S - _O_2))
        theoretical_volumes_ro_2 = 0.01 * (_CO_2 + _CO + _H_2S + _CH_4 + _CH_2 + 2 * _C_2H_6 + 3 * _C_3H_8 + 4 * _C_4H_10)
        theoretical_volumes_h20 = 0.01 * (_H_2 + 2 * _CH_2 + 2 * _C_2H_4 + _H_2S + 3 * _C_2H_6 + 4 * _C_3H_8 + 5 * _C_4H_10 + 0.124 * (d_t + excess_air_ratio_in_flue_gases * stoichiometric_air_flow * d_air))
        theoretical_volumes_N_2 = 0.79 * stoichiometric_air_flow + 0.01 * _N_2_atm
        full_volumes = theoretical_volumes_ro_2 + theoretical_volumes_N_2 + theoretical_volumes_h20
        r_ro_2 = (theoretical_volumes_ro_2 / full_volumes)
        r_h20 = (theoretical_volumes_h20 / full_volumes)
        r_N_2 = (theoretical_volumes_N_2 / full_volumes)
        
        enthalpy_CO2 = CO_2.get('h')
        enthalpy_H20 = H2O.get('h')
        enthalpy_N2 = N_2_atm.get('h')
        h = [(r_ro_2 * enthalpy_CO2[i] + r_h20 * enthalpy_H20[i] + r_N_2 * enthalpy_N2[i]) for i in range(len(temperature_C))]
       
        heat_capacity_CO2 = CO_2.get('Cp')
        heat_capacity_N2 = N_2_atm.get('Cp')
        heat_capacity_H20 = H2O.get('Cp')
        Cp = [r_ro_2 * heat_capacity_CO2[i] + r_h20 * heat_capacity_H20[i] + r_N_2 * heat_capacity_N2[i] for i in range(len(temperature_C))]
        
        entropy_CO2 = CO_2.get('s_0')
        entropy_N2 = N_2_atm.get('s_0')
        entropy_H20 = H2O.get('s_0')
        s_0 = [r_ro_2 * entropy_CO2[i] + r_h20 * entropy_H20[i] + r_N_2 * entropy_N2[i] for i in range(len(temperature_C))] 

        mu = (mu_CO_2 * r_ro_2 + mu_H2O * r_h20 + mu_N_2_atm * r_N_2)
        R = R_const / mu
        u = [h[i] - R * temperature_K[i] for i in range(len(temperature_K))]
        Cv = [Cp[i] - R for i in range(len(temperature_K))]
        k = [Cp[i] / Cv[i] for i in range(len(temperature_K))]

        return {'temperature_C':temperature_C, 'temperature_K':temperature_K, 'temperature_K':temperature_K, 
                'u':u, 'h':h, 's_0':s_0, 'Cp':Cp, 'Cv':Cv, 'k':k, 'R':R, 'mu':mu,
                'heat_of_combustion':heat_of_combustion,
                'stoichiometric_air_flow':stoichiometric_air_flow}

N_2 = properties(temperature_C = temperature_Cels, coefficient = coefficient_N_2, mu = mu_N_2, h_k = h_N_2, s_k = s_N_2)
O_2 = properties(temperature_C = temperature_Cels, coefficient = coefficient_O_2, mu = mu_O_2, h_k = h_O_2, s_k = s_O_2)
C_O = properties(temperature_C = temperature_Cels, coefficient = coefficient_C_O, mu = mu_C_O, h_k = h_C_O, s_k = s_C_O)
CO_2 = properties(temperature_C = temperature_Cels, coefficient = coefficient_CO_2, mu = mu_CO_2, h_k = h_CO_2, s_k = s_CO_2)
H2O = properties(temperature_C = temperature_Cels, coefficient = coefficient_H2O, mu = mu_H2O, h_k = h_H2O, s_k = s_H2O)
SO2 = properties(temperature_C = temperature_Cels, coefficient = coefficient_SO2, mu = mu_SO2, h_k = h_SO2, s_k = s_SO2)
air = properties(temperature_C = temperature_Cels, coefficient = coefficient_air, mu = mu_air, h_k = h_air, s_k = s_air)
N_2_atm = properties(temperature_C = temperature_Cels, coefficient = coefficient_N_2_atm, mu = mu_N_2_atm, h_k = h_N_2_atm, s_k = s_N_2_atm)
NO = properties(temperature_C = temperature_Cels, coefficient = coefficient_NO, mu = mu_NO, h_k = h_NO, s_k = s_NO)
NO_2 = properties(temperature_C = temperature_Cels, coefficient = coefficient_NO_2, mu = mu_NO_2, h_k = h_NO_2, s_k = s_NO_2)
H_2 = properties(temperature_C = temperature_Cels, coefficient = coefficient_H_2, mu = mu_H_2, h_k = h_H_2, s_k = s_H_2)
Ar = properties(temperature_C = temperature_Cels, coefficient = coefficient_Ar, mu = mu_Ar, h_k = h_Ar, s_k = s_Ar)
Ne = properties(temperature_C = temperature_Cels, coefficient = coefficient_Ne, mu = mu_Ne, h_k = h_Ne, s_k = s_Ne)

# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# fuel_plot(fuel['s_0'], fuel['temperature_K'], 's', 't')

# print(fuel)
# print(fuel_inter(fuel.get('Cp'),fuel.get('temperature_C'), 1060))

# N_2_table = fuel_table(temperature_C = N_2.get('temperature_C'), temperature_K = N_2.get('temperature_K'),
#                   u = N_2.get('u'), h = N_2.get('h'), 
#                   s_0 = N_2.get('s_0'), Cp = N_2.get('Cp'), 
#                   Cv = N_2.get('Cv'), k = N_2.get('k'))

# O_2_table = fuel_table(temperature_C = O_2.get('temperature_C'), temperature_K = O_2.get('temperature_K'),
#                   u = O_2.get('u'), h = O_2.get('h'), 
#                   s_0 = O_2.get('s_0'), Cp = O_2.get('Cp'), 
#                   Cv = O_2.get('Cv'), k = O_2.get('k'))

# C_O_table = fuel_table(temperature_C = C_O.get('temperature_C'), temperature_K = C_O.get('temperature_K'),
#                   u = C_O.get('u'), h = C_O.get('h'), 
#                   s_0 = C_O.get('s_0'), Cp = C_O.get('Cp'), 
#                   Cv = C_O.get('Cv'), k = C_O.get('k'))

# CO_2_table = fuel_table(temperature_C = CO_2.get('temperature_C'), temperature_K = CO_2.get('temperature_K'),
#                   u = CO_2.get('u'), h = CO_2.get('h'), 
#                   s_0 = CO_2.get('s_0'), Cp = CO_2.get('Cp'), 
#                   Cv = CO_2.get('Cv'), k = CO_2.get('k'))

# H2O_table = fuel_table(temperature_C = H2O.get('temperature_C'), temperature_K = H2O.get('temperature_K'),
#                   u = H2O.get('u'), h = H2O.get('h'), 
#                   s_0 = H2O.get('s_0'), Cp = H2O.get('Cp'), 
#                   Cv = H2O.get('Cv'), k = H2O.get('k'))

# SO2_table = fuel_table(temperature_C = SO2.get('temperature_C'), temperature_K = SO2.get('temperature_K'),
#                   u = SO2.get('u'), h = SO2.get('h'), 
#                   s_0 = SO2.get('s_0'), Cp = SO2.get('Cp'), 
#                   Cv = SO2.get('Cv'), k = SO2.get('k'))

# air_table = fuel_table(temperature_C = air.get('temperature_C'), temperature_K = air.get('temperature_K'),
#                   u = air.get('u'), h = air.get('h'), 
#                   s_0 = air.get('s_0'), Cp = air.get('Cp'), 
#                   Cv = air.get('Cv'), k = air.get('k'))

# N_2_atm_table = fuel_table(temperature_C = N_2_atm.get('temperature_C'), temperature_K = N_2_atm.get('temperature_K'),
#                   u = N_2_atm.get('u'), h = N_2_atm.get('h'), 
#                   s_0 = N_2_atm.get('s_0'), Cp = N_2_atm.get('Cp'), 
#                   Cv = N_2_atm.get('Cv'), k = N_2_atm.get('k'))

# NO_table = fuel_table(temperature_C = NO.get('temperature_C'), temperature_K = NO.get('temperature_K'),
#                   u = NO.get('u'), h = NO.get('h'), 
#                   s_0 = NO.get('s_0'), Cp = NO.get('Cp'), 
#                   Cv = NO.get('Cv'), k = NO.get('k'))

# NO_2_table = fuel_table(temperature_C = NO_2.get('temperature_C'), temperature_K = NO_2.get('temperature_K'),
#                   u = NO_2.get('u'), h = NO_2.get('h'), 
#                   s_0 = NO_2.get('s_0'), Cp = NO_2.get('Cp'), 
#                   Cv = NO_2.get('Cv'), k = NO_2.get('k'))

# H_2_table = fuel_table(temperature_C = H_2.get('temperature_C'), temperature_K = H_2.get('temperature_K'),
#                   u = H_2.get('u'), h = H_2.get('h'), 
#                   s_0 = H_2.get('s_0'), Cp = H_2.get('Cp'), 
#                   Cv = H_2.get('Cv'), k = H_2.get('k'))

# Ar_table = fuel_table(temperature_C = Ar.get('temperature_C'), temperature_K = Ar.get('temperature_K'),
#                   u = Ar.get('u'), h = Ar.get('h'), 
#                   s_0 = Ar.get('s_0'), Cp = Ar.get('Cp'), 
#                   Cv = Ar.get('Cv'), k = Ar.get('k'))

# Ne_table = fuel_table(temperature_C = Ne.get('temperature_C'), temperature_K = Ne.get('temperature_K'),
#                   u = Ne.get('u'), h = Ne.get('h'), 
#                   s_0 = Ne.get('s_0'), Cp = Ne.get('Cp'), 
#                   Cv = Ne.get('Cv'), k = Ne.get('k'))

# f_table = fuel_table(temperature_C = fuel.get('temperature_C'), temperature_K = fuel.get('temperature_K'),
#                   u = fuel.get('u'), h = fuel.get('h'), 
#                   s_0 = fuel.get('s_0'), Cp = fuel.get('Cp'), 
#                   Cv = fuel.get('Cv'), k = fuel.get('k'))

# print(N_2_table)
# print(O_2_table)
# print(C_O_table)
# print(CO_2_table)
# print(H2O_table)
# print(SO2_table)
# print(air_table)
# print(N_2_atm_table)
# print(NO_table)
# print(NO_2_table)
# print(H_2_table)
# print(Ar_table)
# print(Ne_table)
# print(f_table)