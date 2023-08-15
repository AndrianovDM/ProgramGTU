from geometry import *
from distributionFunction import *
from distributionPlot import *
from distributionTable import *

@st.cache
def parameters(fuel, sch, geom, alpha_02_i, delta_H, periodicity):
    
    alpha_02_z = 0
    ro_st_z = 0

    n_steps = [(i) / (geom[0]['number_of_steps']) for i in range(geom[0]['number_of_steps'] + 1)]
    delta_C2a_i = [np.round((2.94 * (n_steps[i]**1.753) * np.exp(-1.0693 * n_steps[i])),4) for i in range(geom[0]['number_of_steps'] + 1)]
    delta_C2a_i = np.delete(delta_C2a_i, geom[0]['number_of_steps'])
    delta_C2a_i = np.append(delta_C2a_i, 1)
    delta_C2a_i = np.delete(delta_C2a_i, geom[0]['number_of_steps'] - geom[0]['number_of_steps'])
    C2a_i = [(delta_C2a_i[i] * (geom[0]['axial_speed_outlet'] - geom[0]['axial_speed_input']) + geom[0]['axial_speed_input']) for i in range(geom[0]['number_of_steps'])]
    D_sr_i = [(geom[2]['Dsr_sopl_inl_i'][i] + geom[2]['Dsr_sopl_out_i'][i] + geom[2]['Dsr_rab_inl_i'][i] + geom[2]['Dsr_rab_out_i'][i]) / 4 for i in range(geom[0]['number_of_steps'])]
    U_sr_i = [pi * D_sr_i[i] * periodicity  for i in range(geom[0]['number_of_steps'])]    
    C2a_i_ = [C2a_i[i] / U_sr_i[i] for i in range(geom[0]['number_of_steps'])]
    D_sr = sum([geom[2]['Dsr_sopl_inl_i'][i] + geom[2]['Dsr_sopl_out_i'][i] +
    geom[2]['Dsr_rab_inl_i'][i] + geom[2]['Dsr_rab_out_i'][i] for i in range(geom[0]['number_of_steps'])]) / (4 * geom[0]['number_of_steps'])
    
    Y_t = (pi * D_sr * periodicity) * sqrt(geom[0]['number_of_steps'] * sch[1]['etta_is_t'] / (2000 * sch[1]['H_t']))
    Y_st_z = musatkin_Y(C2a_i_[geom[0]['number_of_steps'] - 1], ro_st_z, alpha_02_z)    
    delta_Y = abs((Y_st_z - Y_t) / Y_st_z * 100)

    if delta_Y >= 1:
        while delta_Y >= 0:
            alpha_02_z += 0.01
            ro_st_z = musatkin_Ro(C2a_i_[geom[0]['number_of_steps'] - 1], Y_st_z, alpha_02_z)
            Y_st_z = musatkin_Y(C2a_i_[geom[0]['number_of_steps'] - 1], ro_st_z, alpha_02_z)
            delta_Y = abs((Y_st_z - Y_t) / Y_st_z * 100)
            if alpha_02_z >= alpha_02_i[geom[0]['number_of_steps']-1]:
                break
    
    alpha_02_i[geom[0]['number_of_steps'] - 1] = alpha_02_z
    alpha_0_i = [90] 
    alpha_0_i.extend(alpha_02_i[0:geom[0]['number_of_steps']]) 
    alpha_0_i = (alpha_0_i[0:geom[0]['number_of_steps']]) 
    
    C0a_i = np.delete(C2a_i, geom[0]['number_of_steps'] - 1)
    C0a_i = np.insert(C0a_i, geom[0]['number_of_steps'] - geom[0]['number_of_steps'], geom[0]['axial_speed_input'])

    mu_st_z = 1 / (2 * Y_st_z**2)
    Ls_st_z_ = U_sr_i[geom[0]['number_of_steps'] - 1]**2 / (2000 * Y_st_z**2)
    etta_st_z = smith_etta(C2a_i_[geom[0]['number_of_steps'] - 1], mu_st_z, val)
    
    Ls_1_ = Ls_st_z_ * delta_H[0] + Ls_st_z_
    Ls_st_i_ = [Ls_1_]
    for i in range(geom[0]['number_of_steps'] - 2):
        Ls_st_i_.append((sch[1]['H_t'] / sch[1]['etta_is_t'] - Ls_st_z_ - Ls_1_) * delta_H[i+1])
    Ls_st_i_.append(Ls_st_z_)
    
    Y_st_i  = [sqrt(U_sr_i[i]**2 / (2000 * Ls_st_i_[i])) for i in range(geom[0]['number_of_steps'])]
    ro_st_i = [musatkin_Ro(C2a_i_[i], Y_st_i[i], alpha_02_i[i]) for i in range(geom[0]['number_of_steps'])]
    mu_st_i = [1 / (2 * Y_st_i[i]**2) for i in range(geom[0]['number_of_steps'])]
    etta_st_i = [float(smith_etta(C2a_i_[i], mu_st_i[i], val)) for i in range(geom[0]['number_of_steps'])]
    L_st_i_ = [Ls_st_i_[i] * etta_st_i[i] for i in range(geom[0]['number_of_steps'])]
    
    I0_i_, T0_i_, P0_i_= [I((sch[1]['t_c'] + 273.15), fuel)], [(sch[1]['t_c'] + 273.15)],[sch[1]['P_c']]
    I2_i_, T2_i_, P2_i_=[],[],[]
    I2s_i_, T2s_i_=[],[]  
    
    for i in range(geom[0]['number_of_steps']):
        if i < geom[0]['number_of_steps'] - 1:
            I0_i_.append(I0_i_[i] - L_st_i_[i])
            T0_i_.append(T(I0_i_[i + 1], fuel))
        I2_i_.append(I0_i_[i] - L_st_i_[i])
        T2_i_.append(T(I2_i_[i], fuel))
        I2s_i_.append(I0_i_[i] - Ls_st_i_[i])
        T2s_i_.append(T(I2s_i_[i], fuel))
        if i< geom[0]['number_of_steps']-1:
            P0_i_.append(P0_i_[i] * ((T2s_i_[i]) / (T0_i_[i]))**(k(T0_i_[i], fuel) / (k(T0_i_[i],fuel) - 1)))
        P2_i_.append(P0_i_[i] * ((T2s_i_[i]) / (T0_i_[i]))**(k(T0_i_[i], fuel) / (k(T0_i_[i], fuel) - 1)))
    
    q_lamda_c2a_i = [(sch[2]['G_t'] * sqrt(T2_i_[i])) / (m_(T2_i_[i], fuel, sch) * P2_i_[i] * geom[2]['F_2_out_i'][i] * sin(radians(alpha_02_i[i]))) for i in range(geom[0]['number_of_steps'])]
    lamda_i=[]
    for n in q_lamda_c2a_i:
        def q(x):
            return x * (1 - (k(T2_i_[i], fuel) - 1) / (k(T2_i_[i], fuel) + 1) * x**2)**(1 / (k(T2_i_[i], fuel) - 1)) - n / ((k(T2_i_[i], fuel) + 1) / 2)**(1 / (k(T2_i_[i], fuel) - 1))
        lamda_i.append(fsolve(q,0.01)) 
    lamda_c2a_i = [i[0] for i in lamda_i]      
    PI_lamda_c2a_i = [PI_lamda(lamda_c2a_i[i], T2_i_[i], fuel) for i in range(geom[0]['number_of_steps'])]
    tau_lamda_c2a_i = [tau_lamda(lamda_c2a_i[i], T2_i_[i], fuel) for i in range(geom[0]['number_of_steps'])]
    
    P2_i = [P2_i_[i] * PI_lamda_c2a_i[i] for i in range(geom[0]['number_of_steps'])]
    T2_i = [T2_i_[i] * tau_lamda_c2a_i[i] for i in range(geom[0]['number_of_steps'])]
    T2s_i = [T0_i_[i] * (P2_i[i] / P0_i_[i])**((k(T0_i_[i], fuel) - 1) / (k(T0_i_[i], fuel))) for i in range(geom[0]['number_of_steps'])]
    Hs_st_i = [I(T0_i_[i], fuel) - I(T2s_i[i], fuel) for i in range(geom[0]['number_of_steps'])]

    return {'alpha_0_i':alpha_0_i, 'alpha_02_i':alpha_02_i, 'ro_st_i':ro_st_i, 'C0a_i':C0a_i ,'C2a_i':C2a_i, 'C2a_i_':C2a_i_,
            'Y_st_i':Y_st_i, 'mu_st_i':mu_st_i, 'etta_st_i':etta_st_i, 'Ls_st_i_':Ls_st_i_, 'L_st_i_':L_st_i_, 
            'Hs_st_i':Hs_st_i,'number_of_steps':geom[0]['number_of_steps'], 'periodicity':periodicity}, {'T0_i_':T0_i_, 'T2_i_':T2_i_, 'T2_i':T2_i, 'T2s_i':T2s_i, 'T2s_i_':T2s_i_, 'P0_i_':P0_i_,
            'P2_i_':P2_i_, 'P2_i':P2_i, 'q_lamda_c2a_i':q_lamda_c2a_i, 'lamda_c2a_i':lamda_c2a_i, 'PI_lamda_c2a_i':PI_lamda_c2a_i,
            'tau_lamda_c2a_i':tau_lamda_c2a_i},(alpha_0_i, alpha_02_i, ro_st_i, C0a_i, C2a_i, C2a_i_, Y_st_i, 
            mu_st_i, etta_st_i, Ls_st_i_, L_st_i_, Hs_st_i),(T0_i_, T2_i_, T2_i, T2s_i, T2s_i_, P0_i_, P2_i_, P2_i, q_lamda_c2a_i, lamda_c2a_i, PI_lamda_c2a_i, tau_lamda_c2a_i)

# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# schemegtu = scheme(fuel = fuel, ele_power = 153*10e5, t_c = 1060,
#                                 t_a = 15, P_a = 100.5*10e2, epsilon = 10.9, coefficient_pressure_loss = 0.97,
#                                 etta_c_c = 0.995, etta_mch = 0.995, etta_e_g = 0.982,
#                                 etta_is_t = 0.89, etta_is_k = 0.87, leak_rate = 0.005, t_w = 850, z = 4, method = 'notcold')

# geometrygtu = geometry(fuel = fuel, sch = schemegtu, number_of_steps = 4, axial_speed_input = 180, axial_speed_outlet = 252,
#              D_sr__h_2_z = 3.5, K_s = 0.065, K_r = 0.06, radial_clearance = 0.0008, method = 'root')

# parametergtu = parameters(fuel = fuel, sch = schemegtu, geom = geometrygtu, alpha_02_i = [50, 85, 85, 90], delta_H = [0.7, 0.5, 0.50], periodicity = 50)

# print(distribution_table(parametergtu[2] , method = 'kinematics'))
# print(distribution_table(parametergtu[3] , method = 'termod'))
# distribution_plot(parametergtu, 'ro')
# distribution_plot(parametergtu, 'alpha')
# distribution_plot(parametergtu, 'velocity')
# distribution_plot(parametergtu, 'temperature')
# distribution_plot(parametergtu, 'pressure')
# distribution_plot(parametergtu, 'etta')
# distribution_plot(parametergtu, 'heatdrop')

# etta_plot_2D([c_85, c_86, c_87, c_88, c_89, c_90, c_91, c_92, c_93, c_94], 
#              [mu_85, mu_86, mu_87, mu_88, mu_89, mu_90, mu_91, mu_92, mu_93, mu_94], parametergtu[0]['C2a_i_'], parametergtu[0]['mu_st_i'], parametergtu[0]['etta_st_i'])



