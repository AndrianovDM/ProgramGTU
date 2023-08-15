from distribution import *
from stagePlot import *
from stageTable import *
from stageFunction import *
from profiling import *
from Losses import *

@st.cache
def stage(fuel, sch, geom, distr, consumption, value1_sopl = None, value1_rab = None, 
        value2_sopl = None, value2_rab = None, value3_sopl = None, value3_rab = None,
        coef_sopl = None, coef_rab = None, B_sopl = None, B_rab = None,
        SorU_sopl = None, SorU_rab = None, ks_sopl = None, ks_rab = None, 
        method_losses_sopl = None, method_losses_rab = None, method_bandage = None, n = None):
    
    Hs_sa = distr[0]['Hs_st_i'][n] * (1 - distr[0]['ro_st_i'][n])
    Hs_rk = distr[0]['Hs_st_i'][n] - Hs_sa

    C1s = sqrt(2e3 * Hs_sa)
    fi = 1
    C1 = C1s * fi
    U1 = pi * ((geom[2]['Dsr_sopl_inl_i'][n] + geom[2]['Dsr_sopl_out_i'][n]) / 2) * distr[0]['periodicity']
    U2 = pi * ((geom[2]['Dsr_rab_inl_i'][n] + geom[2]['Dsr_rab_out_i'][n]) / 2) * distr[0]['periodicity']
  
    # Точка 0*
    ################################################################
    I0_ = I(distr[1]['T0_i_'][n], fuel)
    T0_ = T(I0_, fuel)
    P0_ = distr[1]['P0_i_'][n]
    V0_ = sch[1]['R_g'] * 1e3 * T0_ / P0_
    S0_ = S(T0_, fuel)
    point0_ = {'I0_':I0_, 'T0_':T0_, 'S0_':S0_, 'P0_':P0_, 'V0_':V0_ }
    ################################################################

    lambda_C1s = lamda(C1s, T0_, fuel, sch)
    lambda_C1 = lamda(C1, T0_, fuel, sch)
    PI_lambda_C1s = PI_lamda(lambda_C1s, T0_, fuel)
    PI_lambda_C1 = PI_lamda(lambda_C1, T0_, fuel)
    q_lambda_C1 = q_lamda(lambda_C1, T0_, fuel)
    q_lambda_C1s = q_lamda(lambda_C1s, T0_, fuel)
    
    # Точка 1s
    ################################################################
    I1s = I(T0_, fuel) - (C1s**2 / 2e3)
    T1s = T(I1s, fuel)
    P1s =  P0_ * (T1s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
    V1s = sch[1]['R_g'] * 1e3 * T1s / P1s   
    S1s = S(T1s, fuel)
    point1s = {'I1s':I1s, 'T1s':T1s, 'S1s':S1s, 'P1s':P1s, 'V1s':V1s}
    ################################################################

    # Точка 1
    ################################################################  
    I1 = I(T0_, fuel) - (C1**2 / 2e3)
    T1  = T(I1, fuel) 
    P1 = P0_ * (T1s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
    V1 = sch[1]['R_g'] * 1e3 * T1 / P1   
    S1 = S(T1, fuel)   
    point1 = {'I1':I1, 'T1':T1, 'S1':S1, 'P1':P1, 'V1':V1}
    ################################################################   

    if n == 0:
        G_gap_sopl_gas = 0
    else:
        D_y_sopl = 0.5 * (((geom[2]['Dk_sopl_inl_i'][n] + geom[2]['Dk_sopl_out_i'][n]) / 2) + ((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2)) + 0.005
        delta_y_sopl = ((geom[2]['h_sopl_inl_i'][n] + geom[2]['h_sopl_out_i'][n]) / 4) * 0.005
        F_gap_sopl = pi * D_y_sopl * delta_y_sopl
        G_gap_sopl_gas = 0.6 * F_gap_sopl * sqrt((P0_**2 - P1**2) / (4 * P0_ * V0_))

    G1 = consumption - G_gap_sopl_gas

    # Точка 1*
    ################################################################  
    P1_ = P1 / PI_lambda_C1
    T1_ = T0_
    I1_ = I(T1_, fuel)
    V1_ = sch[1]['R_g'] * 1e3 * T1_ / P1_ 
    S1_ = S(T1_, fuel)
    point1_ = {'I1_':I1_, 'T1_':T1_, 'S1_':S1_, 'P1_':P1_, 'V1_':V1_}
    ################################################################

    sigmma_sa = PI_lambda_C1s / PI_lambda_C1
    alpha1 = degrees(asin((G1 * sqrt(T0_)) / (m_(T0_, fuel, sch) * P0_ * q_lambda_C1 * geom[2]['F_1_out_i'][n] * sigmma_sa))) 
    C1u = C1 * cos(radians(alpha1))
    C1a = C1 * sin(radians(alpha1))
    W1u = C1u - U1

    if W1u >= 0:
        betta1 = degrees(atan(C1a / W1u))
    else:
        betta1 = 180 - degrees(atan( C1a / abs(W1u))) 
    W1 = C1a / sin(radians(betta1))

    W1a = W1 * sin(radians(betta1))

    # Точка 1w*
    ################################################################     
    I1w_ = I(T1, fuel) + (W1**2) / 2e3
    T1w_ = T(I1w_, fuel)
    P1w_ = P1 * (T1w_ / T1)**(k(T1, fuel) / (k(T1, fuel) - 1))
    V1w_ = sch[1]['R_g'] * 1e3 * T1w_ / P1w_ 
    S1w_ = S(T1w_, fuel)   
    point1w_ = {'I1w_':I1w_, 'T1w_':T1w_, 'S1w_':S1w_, 'P1w_':P1w_, 'V1w_':V1w_}
    ################################################################  

    a1c = sqrt(k(T1, fuel) * sch[1]['R_g'] * 1e3 * T1)
    M1c = C1 / a1c
    a1w = sqrt(k(T1, fuel) * sch[1]['R_g'] * 1e3 * T1)
    M1w = W1 / a1w
    a0c = sqrt(k(T0_, fuel) * sch[1]['R_g'] * 1e3 * T0_)
    M0c = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))) / a1c

    profile_sopl = profiling(D_k = geom[2]['Dk_sopl_out_i'][n] + geom[2]['h_sopl_out_i'][n], width = geom[2]['S_sopl_i'][n], height = geom[2]['h_sopl_out_i'][n], 
    alpha_0 = distr[0]['alpha_0_i'][n], alpha_1 = alpha1, M_c1 = M1c, 
    value_1 = value1_sopl, value_2 = value2_sopl, value_3 = value3_sopl, method_2 = 'sopl')

    fi_sopl = losses(fuel, temperature = T1, pressure = P1, velocity_inlet = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))), 
    velocity_outlet = C1, R = sch[1]['R_g'] * 1e3, M_inlet = M0c, M_outlet = M1c,
    lamda_velocity = lambda_C1, height = geom[2]['h_sopl_out_i'][n], D_k = geom[2]['Dk_sopl_out_i'][n], 
    pitch = profile_sopl['t_sopl'] * 1e-3, width = profile_sopl['B_sopl'] * 1e-3, chord = profile_sopl['b_sopl'] * 1e-3, 
    te_radius = profile_sopl['r_out_sopl'] * 1e-3, C_max = profile_sopl['Cmax_sopl'] * 1e-3, 
    a_inlet = profile_sopl['a_inl_sopl'] * 1e-3, a_outlet = profile_sopl['a_inl_sopl'] * 1e-3, 
    blade_inlet_angle = profile_sopl['alpha0sc_sopl'], inlet_angle = distr[0]['alpha_0_i'][n], outlet_angle = alpha1, 
    design_inc = profile_sopl['alpha0sc_sopl'] - distr[0]['alpha_0_i'][n], 
    coef = coef_sopl, B = B_sopl, SorU = SorU_sopl, nu = Mu(T1s), ks = ks_sopl, method_1 = method_losses_sopl, method_2 = 'sopl')

    while abs(fi - fi_sopl[0]) >= 10e-6:
        fi = fi_sopl[0]
        C1 = C1s * fi
        # Точка 0*
        ################################################################
        I0_ = I(distr[1]['T0_i_'][n], fuel)
        T0_ = T(I0_, fuel)
        P0_ = distr[1]['P0_i_'][n]
        V0_ = sch[1]['R_g'] * 1e3 * T0_ / P0_
        S0_ = S(T0_, fuel)
        point0_ = {'I0_':I0_, 'T0_':T0_, 'S0_':S0_, 'P0_':P0_, 'V0_':V0_ }
        ################################################################

        lambda_C1s = lamda(C1s, T0_, fuel, sch)
        lambda_C1 = lamda(C1, T0_, fuel, sch)
        PI_lambda_C1s = PI_lamda(lambda_C1s, T0_, fuel)
        PI_lambda_C1 = PI_lamda(lambda_C1, T0_, fuel)
        q_lambda_C1 = q_lamda(lambda_C1, T0_, fuel)
        q_lambda_C1s = q_lamda(lambda_C1s, T0_, fuel)
        
        # Точка 1s
        ################################################################
        I1s = I(T0_, fuel) - (C1s**2 / 2e3)
        T1s = T(I1s, fuel)
        P1s =  P0_ * (T1s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
        V1s = sch[1]['R_g'] * 1e3 * T1s / P1s   
        S1s = S(T1s, fuel)
        point1s = {'I1s':I1s, 'T1s':T1s, 'S1s':S1s, 'P1s':P1s, 'V1s':V1s}
        ################################################################
        
        # Точка 1
        ################################################################  
        I1 = I(T0_, fuel) - (C1**2 / 2e3)
        T1  = T(I1, fuel) 
        P1 = P0_ * (T1s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
        V1 = sch[1]['R_g'] * 1e3 * T1 / P1   
        S1 = S(T1, fuel)   
        point1 = {'I1':I1, 'T1':T1, 'S1':S1, 'P1':P1, 'V1':V1}
        ################################################################   

        if n == 0:
            G_gap_sopl_gas = 0
        else:
            D_y_sopl = 0.5 * (((geom[2]['Dk_sopl_inl_i'][n] + geom[2]['Dk_sopl_out_i'][n]) / 2) + ((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2)) + 0.005
            delta_y_sopl = ((geom[2]['h_sopl_inl_i'][n] + geom[2]['h_sopl_out_i'][n]) / 4) * 0.005
            F_gap_sopl = pi * D_y_sopl * delta_y_sopl
            G_gap_sopl_gas = 0.6 * F_gap_sopl * sqrt((P0_**2 - P1**2) / (4 * P0_ * V0_))
        G1 = consumption - G_gap_sopl_gas

        # Точка 1*
        ################################################################  
        P1_ = P1 / PI_lambda_C1
        T1_ = T0_
        I1_ = I(T1_, fuel)
        V1_ = sch[1]['R_g'] * 1e3 * T1_ / P1_ 
        S1_ = S(T1_, fuel)
        point1_ = {'I1_':I1_, 'T1_':T1_, 'S1_':S1_, 'P1_':P1_, 'V1_':V1_}
        ################################################################

        sigmma_sa = PI_lambda_C1s / PI_lambda_C1
        alpha1 = degrees(asin((G1 * sqrt(T0_)) / (m_(T0_, fuel, sch) * P0_ * q_lambda_C1 * geom[2]['F_1_out_i'][n] * sigmma_sa))) 
        C1u = C1 * cos(radians(alpha1))
        C1a = C1 * sin(radians(alpha1))
        W1u = C1u - U1
        
        if W1u >= 0:
            betta1 = degrees(atan(C1a / W1u))
        else:
            betta1 = 180 - degrees(atan( C1a / abs(W1u))) 
        W1 = C1a / sin(radians(betta1))

        W1a = W1 * sin(radians(betta1))

        # Точка 1w*
        ################################################################     
        I1w_ = I(T1, fuel) + (W1**2) / 2e3
        T1w_ = T(I1w_, fuel)
        P1w_ = P1 * (T1w_ / T1)**(k(T1, fuel) / (k(T1, fuel) - 1))
        V1w_ = sch[1]['R_g'] * 1e3 * T1w_ / P1w_ 
        S1w_ = S(T1w_, fuel)   
        point1w_ = {'I1w_':I1w_, 'T1w_':T1w_, 'S1w_':S1w_, 'P1w_':P1w_, 'V1w_':V1w_}
        ################################################################  

        a1c = sqrt(k(T1, fuel) * sch[1]['R_g'] * 1e3 * T1)
        M1c = C1 / a1c
        a1w = sqrt(k(T1, fuel) * sch[1]['R_g'] * 1e3 * T1)
        M1w = W1 / a1w
        a0c = sqrt(k(T0_, fuel) * sch[1]['R_g'] * 1e3 * T0_)
        M0c = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))) / a1c

        profile_sopl = profiling(D_k = geom[2]['Dk_sopl_out_i'][n] + geom[2]['h_sopl_out_i'][n], width = geom[2]['S_sopl_i'][n], height = geom[2]['h_sopl_out_i'][n], 
        alpha_0 = distr[0]['alpha_0_i'][n], alpha_1 = alpha1, M_c1 = M1c, 
        value_1 = value1_sopl, value_2 = value2_sopl, value_3 = value3_sopl, method_2 = 'sopl')
        
        fi_sopl = losses(fuel, temperature = T1, pressure = P1, velocity_inlet = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))), 
        velocity_outlet = C1, R = sch[1]['R_g'] * 1e3, M_inlet = M0c, M_outlet = M1c,
        lamda_velocity = lambda_C1, height = geom[2]['h_sopl_out_i'][n], D_k = geom[2]['Dk_sopl_out_i'][n], 
        pitch = profile_sopl['t_sopl'] * 1e-3, width = profile_sopl['B_sopl'] * 1e-3, chord = profile_sopl['b_sopl'] * 1e-3, 
        te_radius = profile_sopl['r_out_sopl'] * 1e-3, C_max = profile_sopl['Cmax_sopl'] * 1e-3, 
        a_inlet = profile_sopl['a_inl_sopl'] * 1e-3, a_outlet = profile_sopl['a_inl_sopl'] * 1e-3, 
        blade_inlet_angle = profile_sopl['alpha0sc_sopl'], inlet_angle = distr[0]['alpha_0_i'][n], outlet_angle = alpha1, 
        design_inc = profile_sopl['alpha0sc_sopl'] - distr[0]['alpha_0_i'][n], 
        coef = coef_sopl, B = B_sopl, SorU = SorU_sopl, nu = Mu(T1), ks = ks_sopl, method_1 = method_losses_sopl, method_2 = 'sopl')
        
        fi = fi_sopl[0]
    H_sa = I0_ - I1
    psi = 1
    W2s = sqrt(2e3 * Hs_rk + W1**2)
    W2 = W2s * psi
    T2w_ = T1w_

    lambda_W2s = lamda(W2s, T2w_, fuel, sch)
    lambda_W2 =  lamda(W2, T2w_, fuel, sch)
    PI_lambda_W2s = PI_lamda(lambda_W2s, T2w_, fuel)
    PI_lambda_W2 = PI_lamda(lambda_W2, T2w_, fuel)
    q_lambda_W2s = q_lamda(lambda_W2s, T2w_, fuel)
    q_lambda_W2 = q_lamda(lambda_W2, T2w_, fuel)

    sigmma_rk = PI_lambda_W2s / PI_lambda_W2
    P2w_= P1w_ * sigmma_rk

    # Точка 2
    ################################################################  
    I2 = I(T2w_, fuel) - (W2**2) / 2e3 
    T2 = T(I2, fuel)
    P2 = P2w_ * (T2 / T2w_)**(k(T2w_, fuel)/(k(T2w_, fuel) - 1))
    V2 = sch[1]['R_g'] * 1e3 * T2 / P2 
    S2 = S(T2, fuel)   
    point2 = {'I2':I2, 'T2':T2, 'S2':S2, 'P2':P2, 'V2':V2}
    ################################################################  

    betta2 = degrees(asin((G1 * sqrt(T2w_)) / (m_(T2w_, fuel, sch) * P1w_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 
    
    coef_S_rab = 0.5
    if method_bandage == False:
        G_gap_rab_gas = (0.5 * 0.006 * G1) / (coef_S_rab * sin(radians(betta2)))
    else:
        D_y = 0.5 * (((geom[2]['Dp_sopl_inl_i'][n] + geom[2]['Dp_sopl_out_i'][n]) / 2) + ((geom[2]['Dp_rab_inl_i'][n] + geom[2]['Dp_rab_out_i'][n]) / 2)) + 0.005
        delta_y = ((geom[2]['h_rab_inl_i'][n] + geom[2]['h_rab_out_i'][n]) / 2) * 0.005
        F_gap = pi * D_y * delta_y
        G_gap_rab_gas = 0.6 * F_gap * sqrt((P1**2 - P2**2) / (4 * P1 * V1))        
    
    G2 = G1 - G_gap_rab_gas
    betta2 = degrees(asin((G2 * sqrt(T2w_)) / (m_(T2w_, fuel, sch) * P1w_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

    # Точка 2w_
    ################################################################  
    I2w_ = I2 + (W2**2) / 2e3
    T2w_ = T(I2w_, fuel)
    P2w_= P1w_ * sigmma_rk
    V2w_ = sch[1]['R_g'] * 1e3 * T2w_ / P2w_ 
    S2w_ = S(T2w_, fuel)  
    point2w_ = {'I2w_':I2w_, 'T2w_':T2w_, 'S2w_':S2w_, 'P2w_':P2w_, 'V2w_':V2w_}    
    ################################################################

    W2u = W2 * cos(radians(betta2))
    W2a = W2 * sin(radians(betta2))
    C2u = W2u - U2
    C2 = sqrt(C2u**2 + W2a**2)

    # Точка 2s
    ################################################################       
    I2s = I(T0_, fuel) - distr[0]['Hs_st_i'][n]
    T2s = T(I2s, fuel) 
    P2s = P2w_ * (T2 / T2w_)**(k(T2w_, fuel)/(k(T2w_, fuel) - 1))
    V2s = sch[1]['R_g'] * 1e3 * T2s / P2s 
    S2s = S(T2s, fuel)  
    point2s = {'I2s':I2s, 'T2s':T2s, 'S2s':S2s, 'P2s':P2s, 'V2s':V2s}
    ################################################################
 
    # Точка 2*
    ################################################################    
    I2_ = I2 + (C2**2) / 2e3 
    T2_ = T(I2_, fuel) 
    lambda_C2 = lamda(C2, T2_, fuel, sch)
    PI_lambda_C2 = PI_lamda(lambda_C2, T2_, fuel)
    tau_lambda_C2 = tau_lamda(lambda_C2, T2_, fuel) 
    q_lambda_C2 = q_lamda(lambda_C2, T2_, fuel) 
    P2_ = P2 / PI_lambda_C2
    V2_ = sch[1]['R_g'] * 1e3 * T2_ / P2_ 
    S2_ = S(T2_, fuel)  
    point2_ = {'I2_':I2_, 'T2_':T2_, 'S2_':S2_, 'P2_':P2_, 'V2_':V2_}
    ################################################################

    # Точка 2s_
    ################################################################       
    I2s_ = I(T1, fuel) - Hs_rk
    T2s_ = T(I2s_, fuel) 
    P2s_ = P2w_ * (T2 / T2w_)**(k(T2w_, fuel)/(k(T2w_, fuel) - 1))
    V2s_ = sch[1]['R_g'] * 1e3 * T2s_ / P2s_ 
    S2s_ = S(T2s_, fuel)  
    point2s_ = {'I2s_':I2s_, 'T2s_':T2s_, 'S2s_':S2s_, 'P2s_':P2s_, 'V2s_':V2s_}
    ################################################################


    if  W2u >= U2:
        alpha2 = degrees(asin(W2a / C2))
    else:
        alpha2 = 180 - degrees(asin(W2a / C2))
    
    C2u = C2 * cos(radians(alpha2)) 
    C2a = C2 * sin(radians(alpha2)) 

    a2c = sqrt(k(T2, fuel) * sch[1]['R_g'] * 1e3 * T2)
    M2c = C2 / a2c
    a2w = sqrt(k(T2, fuel) * sch[1]['R_g'] * 1e3 * T2) 
    M2w = W2 / a2w

    profile_rab = profiling(D_k = geom[2]['Dk_rab_out_i'][n] + geom[2]['h_rab_out_i'][n], width = geom[2]['S_rab_i'][n], height = geom[2]['h_rab_out_i'][n], 
    alpha_0 = betta1, alpha_1 = betta2, M_c1 = M2w, value_1 = value1_rab, 
    value_2 = value2_rab, value_3 = value3_rab, method_2 = 'rab')

    psi_rab = losses(fuel, temperature = T2, pressure = P2, velocity_inlet = W1, 
    velocity_outlet = W2, R = sch[1]['R_g'] * 1e3, M_inlet = M1w, M_outlet = M2w,
    lamda_velocity = lambda_W2, height = geom[2]['h_rab_out_i'][n], D_k = geom[2]['Dk_rab_out_i'][n], 
    pitch = profile_rab['t_rab'] * 1e-3, width = profile_rab['B_rab'] * 1e-3, chord = profile_rab['b_rab'] * 1e-3, 
    te_radius = profile_rab['r_out_rab'] * 1e-3, C_max = profile_rab['Cmax_rab'] * 1e-3, 
    a_inlet = profile_rab['a_inl_rab'] * 1e-3, a_outlet = profile_rab['a_out_rab'] * 1e-3, 
    blade_inlet_angle = profile_rab['betta0sc_rab'], inlet_angle = betta1, outlet_angle = betta2, 
    design_inc = profile_rab['betta0sc_rab'] - betta1, 
    coef = coef_rab, B = B_rab, SorU = SorU_rab, nu = Mu(T2), ks = ks_rab, method_1 = method_losses_rab, method_2 = 'rab')
        
    while abs(psi - psi_rab[0]) >= 10e-6:
        
        psi = psi_rab[0]     
        W2 = W2s * psi
        T2w_ = T1w_

        lambda_W2s = lamda(W2s, T2w_, fuel, sch)
        lambda_W2 =  lamda(W2, T2w_, fuel, sch)
        PI_lambda_W2s = PI_lamda(lambda_W2s, T2w_, fuel)
        PI_lambda_W2 = PI_lamda(lambda_W2, T2w_, fuel)
        q_lambda_W2s = q_lamda(lambda_W2s, T2w_, fuel)
        q_lambda_W2 = q_lamda(lambda_W2, T2w_, fuel)

        sigmma_rk = PI_lambda_W2s / PI_lambda_W2
        P2w_= P1w_ * sigmma_rk

        # Точка 2
        ################################################################  
        I2 = I(T2w_, fuel) - (W2**2) / 2e3 
        T2 = T(I2, fuel)
        P2 = P2w_ * (T2 / T2w_)**(k(T2w_, fuel)/(k(T2w_, fuel) - 1))
        V2 = sch[1]['R_g'] * 1e3 * T2 / P2 
        S2 = S(T2, fuel)   
        point2 = {'I2':I2, 'T2':T2, 'S2':S2, 'P2':P2, 'V2':V2}
        ################################################################  
        
        betta2 = degrees(asin((G1 * sqrt(T2w_)) / (m_(T2w_, fuel, sch) * P1w_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

        coef_S_rab = (profile_rab['t_rab'] * 1e-3) / (profile_rab['b_rab'] * 1e-3)
        if method_bandage == False:
            G_gap_rab_gas = (0.5 * 0.006 * G1) / (coef_S_rab * sin(radians(betta2)))
        else:
            D_y = 0.5 * (((geom[2]['Dp_sopl_inl_i'][n] + geom[2]['Dp_sopl_out_i'][n]) / 2) + ((geom[2]['Dp_rab_inl_i'][n] + geom[2]['Dp_rab_out_i'][n]) / 2)) + 0.005
            delta_y = ((geom[2]['h_rab_inl_i'][n] + geom[2]['h_rab_out_i'][n]) / 2) * 0.005
            F_gap = pi * D_y * delta_y
            G_gap_rab_gas = 0.6 * F_gap * sqrt((P1**2 - P2**2) / (4 * P1 * V1))        
        
        G2 = G1 - G_gap_rab_gas

        betta2 = degrees(asin((G2 * sqrt(T2w_)) / (m_(T2w_, fuel, sch) * P1w_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

        # Точка 2w_
        ################################################################  
        I2w_ = I2 + (W2**2) / 2e3
        T2w_ = T(I2w_, fuel)
        P2w_= P1w_ * sigmma_rk
        V2w_ = sch[1]['R_g'] * 1e3 * T2w_ / P2w_ 
        S2w_ = S(T2w_, fuel)  
        point2w_ = {'I2w_':I2w_, 'T2w_':T2w_, 'S2w_':S2w_, 'P2w_':P2w_, 'V2w_':V2w_}    
        ################################################################

        W2u = W2 * cos(radians(betta2))
        W2a = W2 * sin(radians(betta2))
        C2u = W2u - U2
        C2 = sqrt(C2u**2 + W2a**2)

        # Точка 2s
        ################################################################       
        I2s = I(T0_, fuel) - distr[0]['Hs_st_i'][n]
        T2s = T(I2s, fuel) 
        P2s = P2w_ * (T2 / T2w_)**(k(T2w_, fuel)/(k(T2w_, fuel) - 1))
        V2s = sch[1]['R_g'] * 1e3 * T2s / P2s 
        S2s = S(T2s, fuel)  
        point2s = {'I2s':I2s, 'T2s':T2s, 'S2s':S2s, 'P2s':P2s, 'V2s':V2s}
        ################################################################
    
        # Точка 2*
        ################################################################    
        I2_ = I2 + (C2**2) / 2e3 
        T2_ = T(I2_, fuel) 
        lambda_C2 = lamda(C2, T2_, fuel, sch)
        PI_lambda_C2 = PI_lamda(lambda_C2, T2_, fuel)
        tau_lambda_C2 = tau_lamda(lambda_C2, T2_, fuel) 
        q_lambda_C2 = q_lamda(lambda_C2, T2_, fuel) 
        P2_ = P2 / PI_lambda_C2
        V2_ = sch[1]['R_g'] * 1e3 * T2_ / P2_ 
        S2_ = S(T2_, fuel)  
        point2_ = {'I2_':I2_, 'T2_':T2_, 'S2_':S2_, 'P2_':P2_, 'V2_':V2_}
        ################################################################

        # Точка 2s_
        ################################################################       
        I2s_ = I(T1, fuel) - Hs_rk
        T2s_ = T(I2s_, fuel) 
        P2s_ = P2w_ * (T2 / T2w_)**(k(T2w_, fuel)/(k(T2w_, fuel) - 1))
        V2s_ = sch[1]['R_g'] * 1e3 * T2s_ / P2s_ 
        S2s_ = S(T2s_, fuel)  
        point2s_ = {'I2s_':I2s_, 'T2s_':T2s_, 'S2s_':S2s_, 'P2s_':P2s_, 'V2s_':V2s_}
        ################################################################
    
        if  W2u >= U2:
            alpha2 = degrees(asin(W2a / C2))
        else:
            alpha2 = 180 - degrees(asin(W2a / C2))
        
        C2u = C2 * cos(radians(alpha2)) 
        C2a = C2 * sin(radians(alpha2)) 

        a2c = sqrt(k(T2, fuel) * sch[1]['R_g'] * 1e3 * T2)
        M2c = C2 / a2c
        a2w = sqrt(k(T2, fuel) * sch[1]['R_g'] * 1e3 * T2) 
        M2w = W2 / a2w

        profile_rab = profiling(D_k = geom[2]['Dk_rab_out_i'][n] + geom[2]['h_rab_out_i'][n], width = geom[2]['S_rab_i'][n], height = geom[2]['h_rab_out_i'][n], 
                                alpha_0 = betta1, alpha_1 = betta2, M_c1 = M2w, value_1 = value1_rab, 
                                value_2 = value2_rab, value_3 = value3_rab, method_2 = 'rab')

        psi_rab = losses(fuel, temperature = T2, pressure = P2, velocity_inlet = W1, 
        velocity_outlet = W2, R = sch[1]['R_g'] * 1e3, M_inlet = M1w, M_outlet = M2w,
        lamda_velocity = lambda_W2, height = geom[2]['h_rab_out_i'][n], D_k = geom[2]['Dk_rab_out_i'][n], 
        pitch = profile_rab['t_rab'] * 1e-3, width = profile_rab['B_rab'] * 1e-3, chord = profile_rab['b_rab'] * 1e-3, 
        te_radius = profile_rab['r_out_rab'] * 1e-3, C_max = profile_rab['Cmax_rab'] * 1e-3, 
        a_inlet = profile_rab['a_inl_rab'] * 1e-3, a_outlet = profile_rab['a_out_rab'] * 1e-3, 
        blade_inlet_angle = profile_rab['betta0sc_rab'], inlet_angle = betta1, outlet_angle = betta2, 
        design_inc = profile_rab['betta0sc_rab'] - betta1, 
        coef = coef_rab, B = B_rab, SorU = SorU_rab, nu = Mu(T2), ks = ks_rab, method_1 = method_losses_rab, method_2 = 'rab')
        
        psi = psi_rab[0]  
    H_rk = I1 - I2
    Ls_st_ = I(T0_, fuel) - I(T2_, fuel)

    delta_L_lake_sopl = (G_gap_sopl_gas / consumption) * Ls_st_
    if method_bandage == False:
        delta_L_lake_rab = ksi_leaks(((geom[2]['Dsr_rab_inl_i'][n] + geom[2]['Dsr_rab_out_i'][n]) / 2) / geom[2]['h_rab_out_i'][n], distr[0]['ro_st_i'][n]) * Ls_st_ * (geom[2]['delta_r_i'][n] / geom[2]['h_rab_out_i'][n])
    else:
        delta_L_lake_rab = (G_gap_rab_gas / G1) * Ls_st_
    delta_L_tr = (1.701 / G1) * (((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2)**2) * (((pi * ((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2) * distr[0]['periodicity']) / 100)**3) * (((1 / V1) + (1 / V2)) / 2)
    
    I2_sh_ = I2_ + delta_L_lake_sopl + delta_L_lake_rab + delta_L_tr
    T2_sh_ = T(I2_sh_, fuel)
    I2_sh = I2 + delta_L_lake_sopl + delta_L_lake_rab + delta_L_tr
    T2_sh = T(I2_sh, fuel)  
    P2_sh_ = P2 * (T2_sh_ / T2_sh)**(k(T2_sh, fuel)/(k(T2_sh, fuel) - 1))
    T2s_sh_ =  T0_ * (P2_sh_ / P0_)**((k(T0_, fuel) - 1) / k(T0_, fuel))
    I2s_sh_ = I(T2s_sh_, fuel)

    L_st_ = I0_ - I2_sh_
    etta_st = (I0_ - I2_sh_) / (I0_ - I2s_sh_)
    Lu = (C1u * U1 + C2u * U2) / 1e3   
    N_i = Lu * G2 * 1e3
    C_fict = sqrt(distr[0]['Hs_st_i'][n]*2e3)
    # ro = ((W2s**2 - W1**2) - (U2**2 - U1**2)) / C_fict**2

    U_Cfict = (U1 + U2) / (2 * C_fict)
    delta_L_sopl = I(T1, fuel) - I(T1s, fuel)
    delta_L_rab = I(T2, fuel) - I(T1 * (P2 / P1)**((k(T1, fuel) - 1) / k(T1, fuel)), fuel)
    delta_Lvs = C2**2 / 2e3
    G_cold_sopl = 0
    G_cold_rab = 0
    parametrs = {'Hs_st_':distr[0]['Hs_st_i'][n], 'ro_st':distr[0]['ro_st_i'][n], 'Hs_sa':Hs_sa, 'H_sa':H_sa, 'Hs_rk':Hs_rk, 'H_rk':H_rk, 'T0_':T0_, 'P0_':P0_, 
                'T1s':T1s, 'T1':T1, 'P1':P1, 'T1_':T1_, 'P1_':P1_, 'T1w_':T1w_, 'P1w_':P1w_, 
                'T2s':T2s, 'T2':T2, 'P2':P2, 'T2_':T2_, 'P2_':P2_, 'T2w_':T2w_, 'P2w_':P2w_, 
                'T2s_':T2s_, 'C1s':C1s, 'C1':C1, 'C1u':C1u, 'C1a':C1a, 'W1':W1, 'W1u':W1u, 'W1a':W1a, 'U1':U1, 'alpha1':alpha1, 'betta1':betta1,
                'C2':C2, 'C2u':C2u, 'C2a':C2a, 'W2s':W2s, 'W2':W2, 'W2u':W2u, 'W2a':W2a,  'U2':U2, 'alpha2':alpha2,'betta2':betta2,
                'fi':fi_sopl[0], 'psi':psi_rab[0], 'Y_s_sopl':fi_sopl[1], 'Y_s_rab':psi_rab[1], 'Y_p_sopl':fi_sopl[2],'Y_ p_rab':psi_rab[2],
                'Y_sec_sopl':fi_sopl[3], 'Y_sec_rab':psi_rab[3], 'Y_tl_sopl':fi_sopl[4], 'Y_tl_rab':psi_rab[4],
                'Y_te_sopl':fi_sopl[5], 'Y_te_rab':psi_rab[5], 'Y_cl_sopl':fi_sopl[6], 'Y_cl_rab':psi_rab[6], 'M1c':M1c, 'M2c':M2c, 'M1w':M1w, 'M2w':M2w, 'G0':consumption,
                'G_g_sopl_gas':G_gap_sopl_gas, 'G_cold_sopl':G_cold_sopl, 'G1':G1, 'G_g_rab_gas':G_gap_rab_gas, 'G_cold_rab':G_cold_rab, 'G2':G2, 'Ls_st_':Ls_st_, 'dL_lake_sopl':delta_L_lake_sopl, 'dL_lake_rab':delta_L_lake_rab, 
                'dL_tr':delta_L_tr, 'dL_sopl':delta_L_sopl, 'dL_rab':delta_L_rab, 'dL_vs':delta_Lvs, 'L_st_':L_st_, 'Lu':Lu, 'C_fict':C_fict, 'U_Cfict':U_Cfict, 
                'etta_st':etta_st, 'N_i':N_i}


    points = [point0_, point1s, point1, point1_, point1w_, point2s, point2, point2_, point2s_, point2w_]
    velocity_triangle = [C1, W1, U1, alpha1, betta1, C2, W2, U2, alpha2, betta2]

    return parametrs, points, velocity_triangle, n, G2, profile_sopl, profile_rab

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
#     coef_sopl = 0.005, coef_rab = 0.005, B_sopl = 0.25, B_rab = 0.25, SorU_sopl = 0, SorU_rab = 1, ks_sopl = 1e-6, ks_rab = 1e-6, method_losses_sopl = 'ANM', method_losses_rab = 'ANM', method_bandage = False, n = 0)

# print(stage_table(stages[0], method = "parameters"))
# print(stage_table(stages[5], method = "profile sopl"))
# print(stage_table(stages[6], method = "profile rab"))
# hs_plot(point0_ = stages[1][0], point1s = stages[1][1], point1 = stages[1][2], point1_ = stages[1][3], point1w_ = stages[1][4], 
#         point2s = stages[1][5], point2 = stages[1][6], point2_ = stages[1][7], point2s_ = stages[1][8], point2w_ = stages[1][9], i = stages[3], method = 'notcold')
# velocity_triangle_plot(C_1 = stages[2][0], W_1 = stages[2][1], U_1 = stages[2][2], alpha_1 = stages[2][3], betta_1 = stages[2][4], C_2 = stages[2][5], W_2 = stages[2][6], U_2 = stages[2][7], alpha_2 = stages[2][8], betta_2 = stages[2][9], i = stages[3])

    # Lu = U2 * (C1u + C2u) / 1e3
    # Lu = (C1**2 - C2**2 + W2**2 - W1**2) / 2000
    # Lu_2 = U2 * (W1u + W2u)/ 1000
    # Lu_3 = U2 * (C1u + C2u) / 1000
#     etta_u = Lu / (distr[0]['Hs_st_i'][i] - 0.8 * ((C2)**2 / 2e3))
#     etta_u_ = 2 * ( U1 * C1 * cos(radians(alpha1)) + U2 * C2 * cos(radians(alpha2))) / ( Cfict**2 - 0.8 * C2**2)

