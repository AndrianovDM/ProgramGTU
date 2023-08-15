from stage import *

@st.cache
def spin_laws_stage(fuel, sch, geom, distr, stg, n, m, method):  
    radius_sopl_i = [ (geom[2]['h_sopl_out_i'][n] / (m - 1)) * i + (geom[2]['Dk_sopl_out_i'][n]) / 2 for i in range(m) ]
    radius_rab_i = [ (geom[2]['h_rab_out_i'][n] / (m - 1)) * i + (geom[2]['Dk_rab_out_i'][n]) / 2 for i in range(m) ]
    
    radius_sopl_i_ = [ radius_sopl_i[i] / ((geom[2]['Dk_sopl_out_i'][n] + geom[2]['h_sopl_out_i'][n]) / 2) for i in range(m) ]  
    radius_rab_i_ = [ radius_rab_i[i] / ((geom[2]['Dk_rab_out_i'][n] + geom[2]['h_rab_out_i'][n]) / 2) for i in range(m) ]
    
    relative_sopl_i =[ (1 / (m - 1))* i for i in range(m)]
    relative_rab_i =[ (1 / (m - 1))* i for i in range(m)]
    fi_i = [stg['fi'][n] for i in range(m)]
    psi_i = [stg['psi'][n] for i in range(m)]

    if method == 'rtgconst':
        alpha_1_i = [ degrees(atan(tan(radians(stg['alpha1'][n])) / radius_sopl_i_[i])) for i in range(m) ]  
        C_1a_i = [ stg['C1a'][n] * (1 + (tan(radians(stg['alpha1'][n])))**2) / (radius_sopl_i_[i]**2 + (tan(radians(stg['alpha1'][n])))**2) for i in range(m) ]
        C_1u_i = [ C_1a_i[i] / (tan(radians(alpha_1_i[i]))) for i in range(m) ]

    if method == 'C1uconst':
        C_1a_i = [ stg['C1a'][n] for i in range(m) ]  
        C_1u_i = [ stg['C1u'][n] / (radius_sopl_i_[i]) for i in range(m) ]  
        alpha_1_i = [degrees(atan(C_1a_i[i] / C_1u_i[i])) for i in range(m)]

    if method == 'alpha1const': 
        alpha_1_i = [ stg['alpha1'][n] for i in range(m)] 
        C_1a_i = [ stg['C1a'][n] / (radius_sopl_i_[i]**((fi_i[i]**2) * (cos(radians(alpha_1_i[i]))**2))) for i in range(m) ]  
        C_1u_i = [ stg['C1u'][n] / (radius_sopl_i_[i]**((fi_i[i]**2) * (cos(radians(alpha_1_i[i]))**2))) for i in range(m) ]

    C_1_i =  [C_1a_i[i] / (sin(radians(alpha_1_i[i]))) for i in range(m) ]
    C_1s_i = [C_1_i[i] / fi_i[i] for i in range(m) ]
    C_2a_i = [stg['C2a'][n] for i in range(m) ]
    U_1_i =  [stg['U1'][n] * radius_sopl_i_[i] for i in range(m) ]
    U_2_i =  [stg['U2'][n] * radius_rab_i_[i]  for i in range(m) ]
    Hs_sa_i = [ (C_1s_i[i]**2) / 2000 for i in range(m)]
    Hs_rab_i = [(stg['Hs_sa'][n] + stg['Hs_rk'][n]) - Hs_sa_i[i] for i in range(m)]        

    # ro_term_i = [1 - ((1 - ro_k)/((radius_sopl_i[i] / radius_sopl_i[0])**(2 * (fi_i[i]**2) * (cos(radians(alpha_1_i[i])))**2))) for i in range(m)]
    ro_term_i = [((stg['Hs_sa'][n] + stg['Hs_rk'][n]) - Hs_sa_i[i]) / (stg['Hs_sa'][n] + stg['Hs_rk'][n]) for i in range(m)]

    betta_1_i = []
    for i in range(m):
        if C_1u_i[i] - U_1_i[i] >= 0:
            betta_1_i_ = degrees(atan((C_1a_i[i]) / (C_1u_i[i] - U_1_i[i])))
            betta_1_i.append(betta_1_i_)    
        else:
            betta_1_i_ = 180 - degrees(atan((C_1a_i[i]) / abs(C_1u_i[i] - U_1_i[i])))
            betta_1_i.append(betta_1_i_)

    W_1_i = [ C_1a_i[i] / (sin(radians(betta_1_i[i]))) for i in range(m) ]
    W_1u_i = [ C_1u_i[i] - U_1_i[i] for i in range(m) ]
    W_1a_i = [ W_1_i[i] * sin(radians(betta_1_i[i])) for i in range(m) ]
    C_2u_i =[ (stg['Lu'][n] * 1e3 / U_2_i[i]) - C_1u_i[i] for i in range(m)]
    W_2u_i = [ C_2u_i[i] + U_2_i[i]  for i in range(m)] 

    alpha_2_i = []
    for i in range(m):
        if W_2u_i[i] >= U_2_i[i]:
            alpha_2_i_ = degrees(atan(C_2a_i[i] / abs((C_2u_i[i]))))
            alpha_2_i.append(alpha_2_i_)  
        else:
            alpha_2_i_ = 180 - degrees(atan(C_2a_i[i] / abs((C_2u_i[i]))))
            alpha_2_i.append(alpha_2_i_)   

    C_2_i = [ C_2a_i[i] / (sin(radians(alpha_2_i[i]))) for i in range(m)]
    C_2a_i = [C_2_i[i] * sin(radians(alpha_2_i[i])) for i in range(m)]
    betta_2_i = [degrees(atan(C_2a_i[i] / (C_2u_i[i] + U_2_i[i]))) for i in range(m)]
    W_2_i = [ C_2a_i[i] / (sin(radians(betta_2_i[i]))) for i in range(m)]
    W_2a_i = [W_2_i[i] * sin(radians(betta_2_i[i])) for i in range(m)]
    ro_k_i = [ 1 - ((C_1u_i[i] - (C_2u_i[i])) / (2 * U_1_i[i])) for i in range(m)]
    
    # Точка 0_
    ################################################################
    T0_i_ = [stg['T0_'][n] for i in range(m)] 
    I0_i_ = [I(T0_i_[i], fuel) for i in range(m)]
    P0_i_ = [stg['P0_'][n] for i in range(m)] 
    V0_i_ = [sch[1]['R_g'] * 1e3 * T0_i_[i] / P0_i_[i] for i in range(m)] 
    S0_i_ = [S(T0_i_[i], fuel) for i in range(m)]
    point0_i_ = {'I0_i_':I0_i_, 'T0_i_':T0_i_, 'S0_i_':S0_i_, 'P0_i_':P0_i_, 'V0_i_':V0_i_}
    
    # Точка 1s
    ################################################################
    I1s_i = [I0_i_[i] - (C_1_i[i]**2 / ((fi_i[i]**2) * 2e3)) for i in range(m)]
    T1s_i = [T(I1s_i[i], fuel) for i in range(m)]
    P1s_i = [P0_i_[i] * (T1s_i[i] / T0_i_[i])**(k(T0_i_[i], fuel) / (k(T0_i_[i], fuel) - 1)) for i in range(m)]
    V1s_i = sch[1]['R_g'] * 1e3 * T1s_i[i] /P1s_i[i]   
    S1s_i = [S(T1s_i[i], fuel) for i in range(m)]
    point1s_i = {'I1s_i':I1s_i, 'T1s_i':T1s_i, 'S1s_i':S1s_i, 'P1s_i':P1s_i, 'V1s_i':V1s_i}
    ################################################################

    # Точка 1
    ################################################################  
    I1_i = [I0_i_[i] - (C_1_i[i]**2 / 2e3) for i in range(m)]
    T1_i  = [T(I1_i[i], fuel) for i in range(m)]
    P1_i = [P0_i_[i] * (T1s_i[i] / T0_i_[i])**(k(T0_i_[i], fuel) / (k(T0_i_[i], fuel) - 1)) for i in range(m)]
    V1_i = [sch[1]['R_g'] * 1e3 * T1_i[i] / P1_i[i] for i in range(m)] 
    S1_i = [S(T1_i[i], fuel) for i in range(m)] 
    point1_i = {'I1_i':I1_i, 'T1_i':T1_i, 'S1_i':S1_i, 'P1_i':P1_i, 'V1_i':V1_i}
    ################################################################   

    # Точка 1w*
    ################################################################     
    I1w_i_ = [I(T1_i[i], fuel) + (W_1_i[i]**2) / 2e3 for i in range(m)]
    T1w_i_ = [T(I1w_i_[i], fuel) for i in range(m)]
    P1w_i_ = [P1_i[i] * (T1w_i_[i] / T1_i[i])**(k(T1_i[i], fuel) / (k(T1_i[i], fuel) - 1)) for i in range(m)]
    V1w_i_ = [sch[1]['R_g'] * 1e3 * T1w_i_[i] / P1w_i_[i] for i in range(m)]
    S1w_i_ = [S(T1w_i_[i], fuel) for i in range(m)]
    point1w_i_ = {'I1w_i_':I1w_i_, 'T1w_i_':T1w_i_, 'S1w_i_':S1w_i_, 'P1w_i_':P1w_i_, 'V1w_i_':V1w_i_}
    ################################################################  

    lambda_W2s_i = [lamda(W_2_i[i] / psi_i[i], T1w_i_[i], fuel, sch) for i in range(m)]
    lambda_W2_i =  [lamda(W_2_i[i], T1w_i_[i], fuel, sch) for i in range(m)]
    PI_lambda_W2s_i = [PI_lamda(lambda_W2s_i[i], T1w_i_[i], fuel) for i in range(m)] 
    PI_lambda_W2_i = [PI_lamda(lambda_W2_i[i], T1w_i_[i], fuel) for i in range(m)]
    sigmma_rk_i = [PI_lambda_W2s_i[i] / PI_lambda_W2_i[i] for i in range(m)]

    # Точка 2
    ################################################################  
    I2_i = [I(T1w_i_[i], fuel) - (W_2_i[i]**2) / 2e3 for i in range(m)]
    T2_i = [T(I2_i[i], fuel) for i in range(m)]
    P2w_i_= [P1w_i_[i] * sigmma_rk_i[i] for i in range(m)]
    P2_i = [P2w_i_[i] * (T2_i[i] / T1w_i_[i])**(k(T1w_i_[i], fuel)/(k(T1w_i_[i], fuel) - 1)) for i in range(m)]
    V2_i = [sch[1]['R_g'] * 1e3 * T2_i[i] / P2_i[i] for i in range(m)]
    S2_i = [S(T2_i[i], fuel) for i in range(m)]
    point2_i = {'I2_i':I2_i, 'T2_i':T2_i, 'S2_i':S2_i, 'P2_i':P2_i, 'V2_i':V2_i}
    ################################################################ 
   
    a1c_i = [sqrt(k(T1_i[i], fuel) * sch[1]['R_g'] * 1e3 * T1_i[i]) for i in range(m)]
    M1c_i = [C_2_i[i] / a1c_i[i] for i in range(m)]
    a1w_i = [sqrt(k(T1_i[i], fuel) * sch[1]['R_g'] * 1e3 * T1_i[i]) for i in range(m)]
    M1w_i = [W_1_i[i] / a1w_i[i] for i in range(m)]

    a2c_i = [sqrt(k(T2_i[i], fuel) * sch[1]['R_g'] * 1e3 * T2_i[i]) for i in range(m)]
    M2c_i = [C_2_i[i] / a2c_i[i] for i in range(m)]
    a2w_i = [sqrt(k(T2_i[i], fuel) * sch[1]['R_g'] * 1e3 * T2_i[i]) for i in range(m)]
    M2w_i = [W_2_i[i] / a2w_i[i] for i in range(m)]

    parametrs = {'radius_sopl_i':np.round((radius_sopl_i),4),'radius_rab_i':np.round((radius_rab_i),4),
                 'radius_sopl_i_':np.round((radius_sopl_i_),4),'radius_rab_i_':np.round((radius_rab_i_),4),
                 'relative_sopl_i':np.round((relative_sopl_i),4),'relative_rab_i':np.round((relative_rab_i),4),
                 'T0_i_':np.round(T0_i_, 2), 'P0_i_':np.round(P0_i_, 2), 
                 'T1_i':np.round(T1_i, 2), 'P1_i':np.round(P1_i, 2), 
                 'T2_i':np.round(T2_i,2), 'P2_i':np.round(P2_i,2),
                 'Hs_sa_i':np.round((Hs_sa_i),2), 'Hs_rab_i':np.round((Hs_rab_i),2),
                 'ro_term_i':np.round((ro_term_i),2), 'ro_k_i':np.round((ro_k_i),2),
                 'C_1_i':np.round((C_1_i),2), 'C_1u_i':np.round((C_1u_i),2),
                 'C_1a_i':np.round((C_1a_i),2), 'W_1_i':np.round((W_1_i),2),
                 'W_1u_i':np.round(W_1u_i,2), 'W_1a_i':np.round(W_1a_i,2),
                 'U_1_i':np.round((U_1_i),2), 'alpha_1_i':np.round((alpha_1_i),2), 
                 'betta_1_i':np.round((betta_1_i), 2), 'C_2_i':np.round((C_2_i),2),
                 'C_2u_i':np.round((C_2u_i),2), 'C_2a_i':np.round((C_2a_i),2),
                 'W_2_i':np.round((W_2_i),2), 'W_2u_i':np.round((W_2u_i),2),
                 'W_2a_i':np.round(W_2a_i,2), 'U_2_i':np.round((U_2_i),2),
                 'alpha_2_i':np.round((alpha_2_i),2), 'betta_2_i':np.round((betta_2_i),2),
                 'M1c_i':np.round(M1c_i, 2), 'M1w_i':np.round(M1w_i, 2),
                 'M2c_i':np.round(M2c_i, 2), 'M2w_i':np.round(M2w_i, 2)}
    param = [radius_sopl_i, radius_rab_i, radius_sopl_i_, radius_rab_i_, relative_sopl_i, relative_rab_i, T0_i_, P0_i_, T1_i, P1_i, T2_i, P2_i, Hs_sa_i, Hs_rab_i, ro_term_i, ro_k_i, C_1_i, C_1u_i, C_1a_i, W_1_i, W_1u_i, W_1a_i, U_1_i, alpha_1_i, betta_1_i, C_2_i, C_2u_i, C_2a_i, W_2_i, W_2u_i, W_2a_i, U_2_i, alpha_2_i, betta_2_i, M1c_i, M1w_i, M2c_i, M2w_i]

    velocity_triangle = [C_1_i, W_1_i, U_1_i, alpha_1_i, betta_1_i,
                         C_2_i, W_2_i, U_2_i, alpha_2_i, betta_2_i]
    
    return parametrs, param, velocity_triangle, n, 

# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# schemegtu = scheme(fuel = fuel, ele_power = 153*1e6, t_c = 1060,
#                                 t_a = 15, P_a = 100.5*1e3, epsilon = 10.9, coefficient_pressure_loss = 0.97,
#                                 etta_c_c = 0.995, etta_mch = 0.995, etta_e_g = 0.982,
#                                 etta_is_t = 0.89, etta_is_k = 0.87, leak_rate = 0.005, t_w = 850, z = 4, method = 'notcold')

# geometrygtu = geometry(fuel = fuel, sch = schemegtu, number_of_steps = 4, axial_speed_input = 180, axial_speed_outlet = 252,
#              D_sr__h_2_z = 3.5, K_s = 0.05, K_r = 0.04, radial_clearance = 0.0008, method = 'root')

# parametergtu = parameters(fuel = fuel, sch = schemegtu, geom = geometrygtu, alpha_02_i = [60, 85, 85, 90], delta_H = [0.5, 0.5, 0.50], periodicity = 50)

# stages = stage(fuel = fuel , sch = schemegtu, geom = geometrygtu, distr = parametergtu, consumption = schemegtu[2]['G_t'],
#     value1_sopl = 2, value1_rab = 2, value2_sopl = 0.79, value2_rab = 1, value3_sopl = 0.001455, value3_rab = 0.00133, 
#     coef_sopl = 0.005, coef_rab = 0.005, B_sopl = 0.25, B_rab = 0.25, SorU_sopl = 0, SorU_rab = 1, ks_sopl = 1e-6, ks_rab = 1e-6, method_losses_sopl = 'SDB', method_losses_rab = 'SDB', method_bandage = False, n = 0)

# sls = spin_laws_stage(fuel = fuel, sch = schemegtu, geom = geometrygtu, distr = parametergtu, stg = stages, m = 5, n = 0, method = 'alpha1const')
# print(stage_section_table(sls[1]))

# v_tr_1 = velocity_triangle_i(C_1_i = sls[2][0], W_1_i = sls[2][1], U_1_i = sls[2][2], alpha_1_i = sls[2][3], betta_1_i = sls[2][4],
#                               C_2_i = sls[2][5], W_2_i = sls[2][6], U_2_i = sls[2][7], alpha_2_i = sls[2][8], betta_2_i = sls[2][9]) 



