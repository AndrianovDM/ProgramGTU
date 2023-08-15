from distribution import *
from stagePlot import *
from stageTable import *
from stageFunction import *
from profiling import *
from Losses import *
from math import *


@st.cache
def stagecold(fuel, sch, geom, distr, consumption = None, value1_sopl = None, value1_rab = None, 
        value2_sopl = None, value2_rab = None, value3_sopl = None, value3_rab = None,
        coef_sopl = None, coef_rab = None, B_sopl = None, B_rab = None,
        SorU_sopl = None, SorU_rab = None, ks_sopl = None, ks_rab = None,  
        g_standard_sopl = None, g_standard_rab = None,
        T_metal_sopl = None, T_metal_rab = None, 
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

    # Точка 2s
    ################################################################    
    I2s = I(T0_, fuel) - distr[0]['Hs_st_i'][n]
    T2s = T(I2s, fuel)
    P2s = P0_ * (T2s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
    V2s = sch[1]['R_g'] * 1e3 * T2s / P2s 
    S2s = S(T2s, fuel)  
    point2s = {'I2s':I2s, 'T2s':T2s, 'S2s':S2s, 'P2s':P2s, 'V2s':V2s}
    ################################################################

    I1 = I(T0_, fuel) - (C1**2 /2e3)
    T1 = T(I1, fuel)
    P1 = P0_ * (T1s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
    V1 = sch[1]['R_g'] * 1e3 * T1 / P1   
    S1 = S(T1, fuel)   
    density1 = P1 / (sch[1]['R_g'] * T1 * 1e3)

    if n == 0:
        G_gap_sopl_gas = 0
    else:
        D_y_sopl = 0.5 * (((geom[2]['Dk_sopl_inl_i'][n] + geom[2]['Dk_sopl_out_i'][n]) / 2) + ((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2)) + 0.005
        delta_y_sopl = ((geom[2]['h_sopl_inl_i'][n] + geom[2]['h_sopl_out_i'][n]) / 4) * 0.005
        F_gap_sopl = pi * D_y_sopl * delta_y_sopl
        G_gap_sopl_gas = 0.6 * F_gap_sopl * sqrt((P0_**2 - P1**2) / (4 * P0_ * V0_))

    G_sopl_gas = consumption - G_gap_sopl_gas

    sigmma_sa = PI_lambda_C1s / PI_lambda_C1
    alpha1 = degrees(asin((G_sopl_gas * sqrt(T0_)) / (m_(T0_, fuel, sch) * P0_ * q_lambda_C1 * geom[2]['F_1_out_i'][n] * sigmma_sa))) 

    T_out_compressor_air_ = sch[0]['t_b_'] + 273.15
    coef_S_sopl_cold_ = 0.5
    coef_t_sopl_cold_ = 0.5

    Re_sopl = (C1 * geom[2]['S_sopl_i'][n] * density1) / (coef_S_sopl_cold_ * Mu(T1))
    Sgs = (sin(radians(distr[0]['alpha_0_i'][n])) / sin(radians(alpha1))) * sqrt(((2 * coef_S_sopl_cold_) / (coef_t_sopl_cold_ * sin(radians(distr[0]['alpha_0_i'][n] + alpha1)) * (cos(radians((distr[0]['alpha_0_i'][n] - alpha1)/2))**2))) - 1)
    T_inl_sopl_air_ = T_out_compressor_air_ + 25
    g_sopl = (g_standard_sopl * 25 * 2.4 * (T0_ / T_inl_sopl_air_)**0.25) / ((Re_sopl**0.34) * (Sgs**0.58) * coef_t_sopl_cold_ * sin(radians(alpha1)))

    T_metal_sopl = T_metal_sopl + 273.15
    T_max_sopl = 0.985 * (distr[1]['T0_i_'][n] + 50)
    if n == 0:
        G_cold_sopl_ = ((T_max_sopl - T_metal_sopl) / (T_metal_sopl - T_inl_sopl_air_)) * g_sopl * 1
    else:
        G_cold_sopl_ = ((T_max_sopl - T_metal_sopl) / (T_metal_sopl - T_inl_sopl_air_)) * g_sopl * 1.1
      
    G_cold_sopl = sch[2]['G_k'] * G_cold_sopl_
    G1_sm = consumption + G_cold_sopl   

    coef_alpha_sopl = (0.28 * 2.4 ) / ((Re_sopl**0.34) * (Sgs**0.58) * coef_t_sopl_cold_ * sin(radians(alpha1)))
    q_cold_sopl_ = coef_alpha_sopl * ((T_max_sopl - T_metal_sopl) / (distr[0]['Hs_st_i'][n] / Cp(T0_, fuel)))
    coef_tau_sopl = q_cold_sopl_ / (1 - distr[0]['ro_st_i'][n])

    n_1 = 1 / (1 - (((k(T0_, fuel) - 1) / k(T1s, fuel)) * (fi**2 + coef_tau_sopl)))

    # Точка 1s cold
    ################################################################
    T1s_cold = T0_ *((P1 / P0_)**((n_1 - 1) / n_1))
    I1s_cold = I(T1s_cold, fuel)
    P1s_cold = P1
    V1s_cold = sch[1]['R_g'] * 1e3 * T1s_cold / P1s_cold   
    S1s_cold = S(T1, fuel)   
    point1s_cold = {'I1s_cold':I1s_cold, 'T1s_cold':T1s_cold, 'S1s_cold':S1s_cold, 'P1s_cold':P1s_cold, 'V1s_cold':V1s_cold}
    ################################################################

    H_sa_cold = (fi**2) * ((I0_ - I1s_cold) - (q_cold_sopl_ * distr[0]['Hs_st_i'][n]))
    C1_gas = sqrt(2e3 * H_sa_cold)
    alpha1_air = alpha1

    ksi_a_sopl = (1 + (G_cold_sopl_ * 0.6 * (sin(radians(alpha1))/sin(radians(alpha1))))) / (1 + G_cold_sopl_)
    ksi_u_sopl = (1 + (G_cold_sopl_ * 0.6 * (sin(radians(alpha1))/sin(radians(alpha1))))) / (1 + G_cold_sopl_)  
    _G_cold_sopl_ = G_cold_sopl_ * (sch[2]['G_k'] / G_sopl_gas)

    T_out_sopl_air_ = T_inl_sopl_air_ + ((q_cold_sopl_ * distr[0]['Hs_st_i'][n]) / (fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_sopl_air_) * _G_cold_sopl_))
    I_out_sopl_air_ = fuel_inter(air.get('h'), air.get('temperature_K'), T_out_sopl_air_)

    ksi_sopl = (1 + (((T_out_sopl_air_ / T1s_cold)**1.25) * G_cold_sopl_)) / (1 + G_cold_sopl_)

    C1_sm = C1_gas * sqrt((ksi_a_sopl * sin(radians(alpha1)))**2 + (ksi_u_sopl * cos(radians(alpha1)))**2)
    alpha1_sm =  degrees(atan((ksi_a_sopl / ksi_u_sopl) * tan(radians(alpha1))))

    W1_sm = sqrt(C1_sm**2 + U1**2 - 2 * U1 * C1_sm * cos(radians(alpha1_sm)))
    betta1_sm = degrees(asin((C1_sm * sin(radians(alpha1_sm))) / W1_sm))
    W1a_sm = W1_sm * sin(radians(betta1_sm))
    W1u_sm = W1_sm * cos(radians(betta1_sm))

    # Точка 1 sm
    ################################################################
    I1_sm = (I(T0_ * ksi_sopl, fuel) - (C1_sm**2 / 2e3))
    T1_sm = T(I1_sm, fuel)
    P1_sm = P1
    V1_sm = sch[1]['R_g'] * 1e3 * T1_sm / P1_sm   
    S1_sm = S(T1_sm, fuel)   
    point1_sm = {'I1_sm':I1_sm, 'T1_sm':T1_sm, 'S1_sm':S1_sm, 'P1_sm':P1_sm, 'V1_sm':V1_sm}
    ################################################################

    # Точка 1* sm
    ################################################################
    T1_sm_ = T0_ * ksi_sopl
    I1_sm_ = I(T0_ * ksi_sopl, fuel)
    P1_sm_ = P1 * (T1_sm_ / T1_sm)**(k(T1_sm, fuel) / (k(T1_sm, fuel) - 1))
    V1_sm_ = sch[1]['R_g'] * 1e3 * T1_sm_ / P1_sm_   
    S1_sm_ = S(T1_sm_, fuel)   
    point1_sm_ = {'I1_sm_':I1_sm_, 'T1_sm_':T1_sm_, 'S1_sm_':S1_sm_, 'P1_sm_':P1_sm_, 'V1_sm_':V1_sm_}
    ################################################################

    # Точка 1w* sm
    ################################################################
    I1w_sm_ = I1_sm + ((W1_sm**2) / 2e3)
    T1w_sm_ = T(I1w_sm_, fuel)
    P1w_sm_ = P1 * (T1w_sm_ / T1_sm)**(k(T1_sm, fuel) / (k(T1_sm, fuel) - 1))
    V1w_sm_ = sch[1]['R_g'] * 1e3 * T1w_sm_ / P1w_sm_   
    S1w_sm_ = S(T1w_sm_, fuel)   
    point1w_sm_ = {'I1w_sm_':I1w_sm_, 'T1w_sm_':T1w_sm_, 'S1w_sm_':S1w_sm_, 'P1w_sm_':P1w_sm_, 'V1w_sm_':V1w_sm_}
    ################################################################

    a1c = sqrt(k(T1_sm, fuel) * sch[1]['R_g'] * 1e3 * T1_sm)
    M1c = C1_sm / a1c
    a1w = sqrt(k(T1_sm, fuel) * sch[1]['R_g'] * 1e3 * T1_sm)
    M1w = W1_sm / a1w
    a0c = sqrt(k(T0_, fuel) * sch[1]['R_g'] * 1e3 * T0_)
    M0c = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))) / a1c

    profile_sopl = profiling(D_k = geom[2]['Dk_sopl_out_i'][n] + geom[2]['h_sopl_out_i'][n], width = geom[2]['S_sopl_i'][n], height = geom[2]['h_sopl_out_i'][n], 
    alpha_0 = distr[0]['alpha_0_i'][n], alpha_1 = alpha1_sm, M_c1 = M1c, 
    value_1 = value1_sopl, value_2 = value2_sopl, value_3 = value3_sopl, method_2 = 'sopl')

    fi_sopl = losses(fuel, temperature = T1_sm, pressure = P1_sm, velocity_inlet = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))), 
    velocity_outlet = C1_sm, R = sch[1]['R_g'] * 1e3, M_inlet = M0c, M_outlet = M1c,
    lamda_velocity = lamda(C1_sm, T0_, fuel, sch), height = geom[2]['h_sopl_out_i'][n], D_k = geom[2]['Dk_sopl_out_i'][n], 
    pitch = profile_sopl['t_sopl'] * 1e-3, width = profile_sopl['B_sopl'] * 1e-3, chord = profile_sopl['b_sopl'] * 1e-3, 
    te_radius = profile_sopl['r_out_sopl'] * 1e-3, C_max = profile_sopl['Cmax_sopl'] * 1e-3, 
    a_inlet = profile_sopl['a_inl_sopl'] * 1e-3, a_outlet = profile_sopl['a_inl_sopl'] * 1e-3, 
    blade_inlet_angle = profile_sopl['alpha0sc_sopl'], inlet_angle = distr[0]['alpha_0_i'][n], outlet_angle = alpha1_sm, 
    design_inc = profile_sopl['alpha0sc_sopl'] - distr[0]['alpha_0_i'][n], 
    coef = coef_sopl, B = B_sopl, SorU = SorU_sopl, nu = Mu(T1_sm), ks = ks_sopl, method_1 = method_losses_sopl, method_2 = 'sopl')

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

        # Точка 2s
        ################################################################    
        I2s = I(T0_, fuel) - distr[0]['Hs_st_i'][n]
        T2s = T(I2s, fuel)
        P2s = P0_ * (T2s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
        V2s = sch[1]['R_g'] * 1e3 * T2s / P2s 
        S2s = S(T2s, fuel)  
        point2s = {'I2s':I2s, 'T2s':T2s, 'S2s':S2s, 'P2s':P2s, 'V2s':V2s}
        ################################################################

        I1 = I(T0_, fuel) - (C1**2 /2e3)
        T1 = T(I1, fuel)
        P1 = P0_ * (T1s / T0_)**(k(T0_, fuel) / (k(T0_, fuel) - 1))
        V1 = sch[1]['R_g'] * 1e3 * T1 / P1   
        S1 = S(T1, fuel)   
        density1 = P1 / (sch[1]['R_g'] * T1 * 1e3)

        if n == 0:
            G_gap_sopl_gas = 0
        else:
            D_y_sopl = 0.5 * (((geom[2]['Dk_sopl_inl_i'][n] + geom[2]['Dk_sopl_out_i'][n]) / 2) + ((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2)) + 0.005
            delta_y_sopl = ((geom[2]['h_sopl_inl_i'][n] + geom[2]['h_sopl_out_i'][n]) / 4) * 0.005
            F_gap_sopl = pi * D_y_sopl * delta_y_sopl
            G_gap_sopl_gas = 0.6 * F_gap_sopl * sqrt((P0_**2 - P1**2) / (4 * P0_ * V0_))

        G_sopl_gas = consumption - G_gap_sopl_gas

        sigmma_sa = PI_lambda_C1s / PI_lambda_C1
        alpha1 = degrees(asin((G_sopl_gas * sqrt(T0_)) / (m_(T0_, fuel, sch) * P0_ * q_lambda_C1 * geom[2]['F_1_out_i'][n] * sigmma_sa))) 

        T_out_compressor_air_ = sch[0]['t_b_'] + 273.15
        coef_S_sopl_cold_ = (profile_sopl['B_sopl'] * 1e-3) / (profile_sopl['b_sopl'] * 1e-3)
        coef_t_sopl_cold_ = (profile_sopl['t_sopl'] * 1e-3) / (profile_sopl['b_sopl'] * 1e-3)

        Re_sopl = (C1 * geom[2]['S_sopl_i'][n] * density1) / (coef_S_sopl_cold_ * Mu(T1))
        Sgs = (sin(radians(distr[0]['alpha_0_i'][n])) / sin(radians(alpha1))) * sqrt(((2 * coef_S_sopl_cold_) / (coef_t_sopl_cold_ * sin(radians(distr[0]['alpha_0_i'][n] + alpha1)) * (cos(radians((distr[0]['alpha_0_i'][n] - alpha1)/2))**2))) - 1)
        T_inl_sopl_air_ = T_out_compressor_air_ + 25
        g_sopl = (g_standard_sopl * 25 * 2.4 * (T0_ / T_inl_sopl_air_)**0.25) / ((Re_sopl**0.34) * (Sgs**0.58) * coef_t_sopl_cold_ * sin(radians(alpha1)))

        T_metal_sopl = T_metal_sopl + 273.15
        T_max_sopl = 0.985 * (distr[1]['T0_i_'][n] + 50)
        if n == 0:
            G_cold_sopl_ = ((T_max_sopl - T_metal_sopl) / (T_metal_sopl - T_inl_sopl_air_)) * g_sopl * 1
        else:
            G_cold_sopl_ = ((T_max_sopl - T_metal_sopl) / (T_metal_sopl - T_inl_sopl_air_)) * g_sopl * 1.1
        
        G_cold_sopl = sch[2]['G_k'] * G_cold_sopl_
        G1_sm = consumption + G_cold_sopl   

        coef_alpha_sopl = (0.28 * 2.4 ) / ((Re_sopl**0.34) * (Sgs**0.58) * coef_t_sopl_cold_ * sin(radians(alpha1)))
        q_cold_sopl_ = coef_alpha_sopl * ((T_max_sopl - T_metal_sopl) / (distr[0]['Hs_st_i'][n] / Cp(T0_, fuel)))
        coef_tau_sopl = q_cold_sopl_ / (1 - distr[0]['ro_st_i'][n])

        n_1 = 1 / (1 - (((k(T0_, fuel) - 1) / k(T1s, fuel)) * (fi**2 + coef_tau_sopl)))

        # Точка 1s cold
        ################################################################
        T1s_cold = T0_ *((P1 / P0_)**((n_1 - 1) / n_1))
        I1s_cold = I(T1s_cold, fuel)
        P1s_cold = P1
        V1s_cold = sch[1]['R_g'] * 1e3 * T1s_cold / P1s_cold   
        S1s_cold = S(T1, fuel)   
        point1s_cold = {'I1s_cold':I1s_cold, 'T1s_cold':T1s_cold, 'S1s_cold':S1s_cold, 'P1s_cold':P1s_cold, 'V1s_cold':V1s_cold}
        ################################################################

        H_sa_cold = (fi**2) * ((I0_ - I1s_cold) - (q_cold_sopl_ * distr[0]['Hs_st_i'][n]))
        C1_gas = sqrt(2e3 * H_sa_cold)
        alpha1_air = alpha1

        ksi_a_sopl = (1 + (G_cold_sopl_ * 0.6 * (sin(radians(alpha1))/sin(radians(alpha1))))) / (1 + G_cold_sopl_)
        ksi_u_sopl = (1 + (G_cold_sopl_ * 0.6 * (sin(radians(alpha1))/sin(radians(alpha1))))) / (1 + G_cold_sopl_)  
        _G_cold_sopl_ = G_cold_sopl_ * (sch[2]['G_k'] / G_sopl_gas)

        T_out_sopl_air_ = T_inl_sopl_air_ + ((q_cold_sopl_ * distr[0]['Hs_st_i'][n]) / (fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_sopl_air_) * _G_cold_sopl_))
        I_out_sopl_air_ = fuel_inter(air.get('h'), air.get('temperature_K'), T_out_sopl_air_)

        ksi_sopl = (1 + (((T_out_sopl_air_ / T1s_cold)**1.25) * G_cold_sopl_)) / (1 + G_cold_sopl_)

        C1_sm = C1_gas * sqrt((ksi_a_sopl * sin(radians(alpha1)))**2 + (ksi_u_sopl * cos(radians(alpha1)))**2)
        alpha1_sm =  degrees(atan((ksi_a_sopl / ksi_u_sopl) * tan(radians(alpha1))))

        W1_sm = sqrt(C1_sm**2 + U1**2 - 2 * U1 * C1_sm * cos(radians(alpha1_sm)))
        betta1_sm = degrees(asin((C1_sm * sin(radians(alpha1_sm))) / W1_sm))
        W1a_sm = W1_sm * sin(radians(betta1_sm))
        W1u_sm = W1_sm * cos(radians(betta1_sm))

        # Точка 1 sm
        ################################################################
        I1_sm = (I(T0_ * ksi_sopl, fuel) - (C1_sm**2 / 2e3))
        T1_sm = T(I1_sm, fuel)
        P1_sm = P1
        V1_sm = sch[1]['R_g'] * 1e3 * T1_sm / P1_sm   
        S1_sm = S(T1_sm, fuel)   
        point1_sm = {'I1_sm':I1_sm, 'T1_sm':T1_sm, 'S1_sm':S1_sm, 'P1_sm':P1_sm, 'V1_sm':V1_sm}
        ################################################################

        # Точка 1* sm
        ################################################################
        T1_sm_ = T0_ * ksi_sopl
        I1_sm_ = I(T0_ * ksi_sopl, fuel)
        P1_sm_ = P1 * (T1_sm_ / T1_sm)**(k(T1_sm, fuel) / (k(T1_sm, fuel) - 1))
        V1_sm_ = sch[1]['R_g'] * 1e3 * T1_sm_ / P1_sm_   
        S1_sm_ = S(T1_sm_, fuel)   
        point1_sm_ = {'I1_sm_':I1_sm_, 'T1_sm_':T1_sm_, 'S1_sm_':S1_sm_, 'P1_sm_':P1_sm_, 'V1_sm_':V1_sm_}
        ################################################################

        # Точка 1w* sm
        ################################################################
        I1w_sm_ = I1_sm + ((W1_sm**2) / 2e3)
        T1w_sm_ = T(I1w_sm_, fuel)
        P1w_sm_ = P1 * (T1w_sm_ / T1_sm)**(k(T1_sm, fuel) / (k(T1_sm, fuel) - 1))
        V1w_sm_ = sch[1]['R_g'] * 1e3 * T1w_sm_ / P1w_sm_   
        S1w_sm_ = S(T1w_sm_, fuel)   
        point1w_sm_ = {'I1w_sm_':I1w_sm_, 'T1w_sm_':T1w_sm_, 'S1w_sm_':S1w_sm_, 'P1w_sm_':P1w_sm_, 'V1w_sm_':V1w_sm_}
        ################################################################

        a1c = sqrt(k(T1_sm, fuel) * sch[1]['R_g'] * 1e3 * T1_sm)
        M1c = C1_sm / a1c
        a1w = sqrt(k(T1_sm, fuel) * sch[1]['R_g'] * 1e3 * T1_sm)
        M1w = W1_sm / a1w
        a0c = sqrt(k(T0_, fuel) * sch[1]['R_g'] * 1e3 * T0_)
        M0c = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))) / a1c

        profile_sopl = profiling(D_k = geom[2]['Dk_sopl_out_i'][n] + geom[2]['h_sopl_out_i'][n], width = geom[2]['S_sopl_i'][n], height = geom[2]['h_sopl_out_i'][n], 
        alpha_0 = distr[0]['alpha_0_i'][n], alpha_1 = alpha1_sm, M_c1 = M1c, 
        value_1 = value1_sopl, value_2 = value2_sopl, value_3 = value3_sopl, method_2 = 'sopl')

        fi_sopl = losses(fuel, temperature = T1_sm, pressure = P1_sm, velocity_inlet = (distr[0]['C0a_i'][n] / sin(radians(distr[0]['alpha_0_i'][n]))), 
        velocity_outlet = C1_sm, R = sch[1]['R_g'] * 1e3, M_inlet = M0c, M_outlet = M1c,
        lamda_velocity = lamda(C1_sm, T0_, fuel, sch), height = geom[2]['h_sopl_out_i'][n], D_k = geom[2]['Dk_sopl_out_i'][n], 
        pitch = profile_sopl['t_sopl'] * 1e-3, width = profile_sopl['B_sopl'] * 1e-3, chord = profile_sopl['b_sopl'] * 1e-3, 
        te_radius = profile_sopl['r_out_sopl'] * 1e-3, C_max = profile_sopl['Cmax_sopl'] * 1e-3, 
        a_inlet = profile_sopl['a_inl_sopl'] * 1e-3, a_outlet = profile_sopl['a_inl_sopl'] * 1e-3, 
        blade_inlet_angle = profile_sopl['alpha0sc_sopl'], inlet_angle = distr[0]['alpha_0_i'][n], outlet_angle = alpha1_sm, 
        design_inc = profile_sopl['alpha0sc_sopl'] - distr[0]['alpha_0_i'][n], 
        coef = coef_sopl, B = B_sopl, SorU = SorU_sopl, nu = Mu(T1_sm), ks = ks_sopl, method_1 = method_losses_sopl, method_2 = 'sopl')

        fi = fi_sopl[0]

    psi = 1
    W2s = sqrt(2e3 * Hs_rk + W1_sm**2)
    W2 = W2s * psi

    I2s_ = (I1 - Hs_rk)
    T2s_ = T(I2s_, fuel)

    I2_ = I2s_ + (((1 / psi**2) - 1) * (W2**2 / 2e3))
    T2_ = T(I2_, fuel)

    density2_ = P2s / (sch[1]['R_g'] * 1e3 * T2_)

    lambda_W2s = lamda(W2s, T2_, fuel, sch)
    lambda_W2 =  lamda(W2, T2_, fuel, sch)

    PI_lambda_W2s = PI_lamda(lambda_W2s, T2_, fuel)
    PI_lambda_W2 = PI_lamda(lambda_W2, T2_, fuel)
    q_lambda_W2s = q_lamda(lambda_W2s, T2_, fuel)
    q_lambda_W2 = q_lamda(lambda_W2, T2_, fuel)
    sigmma_rk = PI_lambda_W2s / PI_lambda_W2

    betta2_ = degrees(asin((G1_sm * sqrt(T2_)) / (m_(T2_, fuel, sch) * P1w_sm_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

    coef_t_rab_cold_ = 0.5
    coef_S_rab_cold_ = 0.5

    if method_bandage == False:
        G_gap_rab_gas = (0.5 * 0.006 * G1_sm) / (coef_t_rab_cold_ * sin(radians(betta2_)))
    else:
        D_y = 0.5 * (((geom[2]['Dp_sopl_inl_i'][n] + geom[2]['Dp_sopl_out_i'][n]) / 2) + ((geom[2]['Dp_rab_inl_i'][n] + geom[2]['Dp_rab_out_i'][n]) / 2)) + 0.005
        delta_y = ((geom[2]['h_rab_inl_i'][n] + geom[2]['h_rab_out_i'][n]) / 2) * 0.005
        F_gap = pi * D_y * delta_y
        G_gap_rab_gas = 0.6 * F_gap * sqrt((P1_sm**2 - P2s**2) / (4 * P1_sm * V1_sm))        
    
    G_rab_gas = G1_sm - G_gap_rab_gas
    betta2 = degrees(asin((G_rab_gas * sqrt(T2_)) / (m_(T2_, fuel, sch) * P1w_sm_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

    Re_rab = (W2 * geom[2]['S_rab_i'][n] * density2_) / (coef_S_rab_cold_ * Mu(T2_))

    Sgr = (sin(radians(betta1_sm)) / sin(radians(betta2))) * sqrt(((2 * coef_S_rab_cold_) / (coef_t_rab_cold_ * sin(radians(betta1_sm + betta2)) * (cos(radians((betta1_sm - betta2)/2))**2))) - 1)
    tetta = ((geom[2]['Dsr_rab_inl_i'][n] + geom[2]['Dsr_rab_out_i'][n]) / 2) / geom[2]['h_rab_out_i'][n]
    Su_rab = U2 / (W2 * tetta)

    T_inl_rab_air_ = T_out_compressor_air_ + 40
    g_rab = (g_standard_rab * 25 * (1 + (0.8 * Su_rab**0.42)) * 2.55 * ((T0_ / T_inl_rab_air_)**0.25 )) / ((Re_rab**0.34) * (Sgr**0.58) * coef_t_rab_cold_ * sin(radians(betta2))) 
    I2_compressor_ = fuel_inter(air.get('h'), air.get('temperature_K'), T_inl_rab_air_)
    I_inl_rab_air_w_ = I2_compressor_ + 40 * fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_) - (((U2**2) * (((tetta - 1) / tetta)**2)) / 2e3)
    T_inl_rab_air_w_ = T(I_inl_rab_air_w_, fuel)

    T_metal_rab = T_metal_rab + 273.15
    T_rab_air = 950

    G_cold_rab_ = ((T_metal_rab - T_rab_air) / (T_rab_air - T_inl_rab_air_w_)) * g_rab * 1.1
    G_cold_rab = sch[2]['G_k'] * round(G_cold_rab_, 3)
    G2_sm = G1_sm + G_cold_rab

    coef_alpha_rab = (1.12 * 10**(-2)) * (g_rab / g_standard_rab) * ((T_inl_rab_air_w_ / T0_)**0.25)
    q_cold_rab_ = coef_alpha_rab * ((T_metal_rab - T_rab_air) * Cp(T0_, fuel) / distr[0]['Hs_st_i'][n])
    coef_tau_rab = q_cold_rab_ / (distr[0]['ro_st_i'][n] + (W1_sm**2 / (2e3 * distr[0]['Hs_st_i'][n])))
    n_2 = 1 / (1 - (((k(T2_, fuel) - 1) / k(T2_, fuel)) * (psi**2 + coef_tau_rab)))

   # Точка 2 cold
    ################################################################
    T2_cold = T1w_sm_ * ((P2s / P1w_sm_)**((n_2 - 1) / n_2))
    P2s_cold = P2s
    I2_cold = I(T2_cold, fuel)
    ################################################################

    Hs_rk_cold = (I1_sm - I2_cold) - (q_cold_rab_ * distr[0]['Hs_st_i'][n])
    W2_gas = psi * sqrt(W1_sm**2 + 2e3 * Hs_rk_cold)
    C2_gas = sqrt(W2_gas**2 + U2**2 - 2 * W2_gas * U2 * cos(radians(betta2)))

    # Точка 2* gas
    ################################################################
    I2_gas_ = I2_cold + ((1 / psi**2) - 1) * (W2_gas**2 / 2e3) + (C2_gas**2 / 2e3)
    T2_gas_ = T(I2_gas_, fuel)
    ################################################################

    _G_cold_rab_ = G_cold_rab_ * (sch[2]['G_k'] / sch[2]['G_t_notcool'])

    delta_T_air_rab = (((q_cold_rab_ * distr[0]['Hs_st_i'][n]) / (_G_cold_rab_ *  fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_w_))) + ((U2**2) * (2 - (((tetta - 1) / tetta)**2)) /(2e3 *  fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_w_))) - ((0.7 * W2_gas * U2 * cos(radians(betta2))) / (1e3 *  fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_w_))))
    T_out_rab_air_w_ = T_inl_rab_air_w_ + delta_T_air_rab
    I_out_rab_air_w_ = I(T_out_rab_air_w_, fuel)

    betta2_air = betta2
    ksi_a_rab = (1 + (G_cold_rab_ * 0.7 * (sin(radians(betta2_air))/sin(radians(betta2))))) / (1 + G_cold_rab_)
    ksi_u_rab = (1 + (G_cold_rab_ * 0.7 * (sin(radians(betta2_air))/sin(radians(betta2))))) / (1 + G_cold_rab_)  
    ksi_rab = (1 + (((T_out_rab_air_w_ / T2_gas_)**1.25) * G_cold_rab_)) / (1 + G_cold_rab_)

    W2_sm = W2_gas * sqrt((ksi_a_rab  * sin(radians(betta2)))**2 + (ksi_u_rab * cos(radians(betta2)))**2)
    betta2_sm =  degrees(atan((ksi_a_rab / ksi_u_rab) * tan(radians(betta2))))
    W2u_sm = W2_sm * cos(radians(betta2_sm))
    W2a_sm = W2_sm * sin(radians(betta2_sm))
    C2_sm = sqrt(W2_sm**2 + U2**2 - (2 * U2 * W2_sm * cos(radians(betta2_sm))))

    # Точка 2 sm
    ################################################################
    I2_sm = I(T2_gas_ * ksi_rab, fuel) - (C2_sm**2 / 2e3)
    T2_sm = T(I2_sm, fuel)
    P2_sm = P2s
    V2_sm = sch[1]['R_g'] * 1e3 * T2_sm / P2_sm   
    S2_sm = S(T2_sm, fuel)   
    point2_sm = {'I2_sm':I2_sm, 'T2_sm':T2_sm, 'S2_sm':S2_sm, 'P2_sm':P2_sm, 'V2_sm':V2_sm}
    ################################################################

    # Точка 2 sm*
    ################################################################
    T2_sm_ = T2_gas_ * ksi_rab
    I2_sm_ = I(T2_sm_, fuel)
    P2_sm_ = P2s * (T2_sm_ / T2_sm)**(k(T2_sm, fuel) / (k(T2_sm, fuel) - 1))
    V2_sm_ = sch[1]['R_g'] * 1e3 * T2_sm / P2_sm   
    S2_sm_ = S(T2_sm, fuel)   
    point2_sm_ = {'I2_sm_':I2_sm_, 'T2_sm_':T2_sm_, 'S2_sm_':S2_sm_, 'P2_sm_':P2_sm_, 'V2_sm_':V2_sm_}
    ################################################################

    if  W2_sm * cos(radians(betta2)) >= U2:
        alpha2_sm = degrees(asin((W2_sm * sin(radians(betta2_sm))) / C2_sm))
    else:
        alpha2_sm = 180 - degrees(asin((W2_sm * sin(radians(betta2_sm))) / C2_sm))

    # Точка 2w* sm
    ################################################################  
    I2w_sm_ = I2_sm + (W2_sm**2) / 2e3
    T2w_sm_ = T(I2w_sm_, fuel)
    P2w_sm_ = P2_sm * (T2w_sm_ / T2_sm)**(k(T2_sm, fuel) / (k(T2_sm, fuel) - 1))
    V2w_sm_ = sch[1]['R_g'] * 1e3 * T2w_sm_ / P2w_sm_ 
    S2w_sm_ = S(I2w_sm_, fuel)  
    point2w_sm_ = {'I2w_sm_':I2w_sm_, 'T2w_sm_':T2w_sm_, 'S2w_sm_':S2w_sm_, 'P2w_sm_':P2w_sm_, 'V2w_sm_':V2w_sm_}    
    ################################################################

    # Точка 2s sm_
    ################################################################       
    I2s_sm_ = I(T1_sm, fuel) - Hs_rk
    T2s_sm_ = T(I2s_sm_, fuel) 
    P2s_sm_ = P2s
    V2s_sm_ = sch[1]['R_g'] * 1e3 * T2s_sm_ / P2s_sm_ 
    S2s_sm_ = S(T2s_sm_, fuel)  
    point2s_sm_ = {'I2s_sm_':I2s_sm_, 'T2s_sm_':T2s_sm_, 'S2s_sm_':S2s_sm_, 'P2s_sm_':P2s_sm_, 'V2s_sm_':V2s_sm_}
    ################################################################

    C1u_sm = C1_sm * cos(radians(alpha1_sm))
    C1a_sm = C1_sm * sin(radians(alpha1_sm))
    C2u_sm = C2_sm * cos(radians(alpha2_sm))
    C2a_sm = C2_sm * sin(radians(alpha2_sm))

    a2c = sqrt(k(T2_sm, fuel) * sch[1]['R_g'] * 1e3 * T2_sm)
    M2c = C2_sm / a2c
    a2w = sqrt(k(T2_sm, fuel) * sch[1]['R_g'] * 1e3 * T2_sm)
    M2w = W2_sm / a2w

    profile_rab = profiling(D_k = geom[2]['Dk_rab_out_i'][n] + geom[2]['h_rab_out_i'][n], width = geom[2]['S_rab_i'][n], height = geom[2]['h_rab_out_i'][n], 
    alpha_0 = betta1_sm, alpha_1 = betta2_sm, M_c1 = M2w, value_1 = value1_rab, 
    value_2 = value2_rab, value_3 = value3_rab, method_2 = 'rab')

    psi_rab = losses(fuel, temperature = T2_sm, pressure = P2_sm, velocity_inlet = W1_sm, 
    velocity_outlet = W2_sm, R = sch[1]['R_g'] * 1e3, M_inlet = M1w, M_outlet = M2w,
    lamda_velocity = lamda(W2_sm, T2_sm_, fuel, sch), height = geom[2]['h_rab_out_i'][n], D_k = geom[2]['Dk_rab_out_i'][n], 
    pitch = profile_rab['t_rab'] * 1e-3, width = profile_rab['B_rab'] * 1e-3, chord = profile_rab['b_rab'] * 1e-3, 
    te_radius = profile_rab['r_out_rab'] * 1e-3, C_max = profile_rab['Cmax_rab'] * 1e-3, 
    a_inlet = profile_rab['a_inl_rab'] * 1e-3, a_outlet = profile_rab['a_out_rab'] * 1e-3, 
    blade_inlet_angle = profile_rab['betta0sc_rab'], inlet_angle = betta1_sm, outlet_angle = betta2_sm, 
    design_inc = profile_rab['betta0sc_rab'] - betta1_sm, 
    coef = coef_rab, B = B_rab, SorU = SorU_rab, nu = Mu(T2_sm), ks = ks_rab, method_1 = method_losses_rab, method_2 = 'rab')

    while abs(psi - psi_rab[0]) >= 10e-6:
        
        psi = psi_rab[0]     
        W2 = W2s * psi

        I2s_ = (I1 - Hs_rk)
        T2s_ = T(I2s_, fuel)

        I2_ = I2s_ + (((1 / psi**2) - 1) * (W2**2 / 2e3))
        T2_ = T(I2_, fuel)

        density2_ = P2s / (sch[1]['R_g'] * 1e3 * T2_)

        lambda_W2s = lamda(W2s, T2_, fuel, sch)
        lambda_W2 =  lamda(W2, T2_, fuel, sch)

        PI_lambda_W2s = PI_lamda(lambda_W2s, T2_, fuel)
        PI_lambda_W2 = PI_lamda(lambda_W2, T2_, fuel)
        q_lambda_W2s = q_lamda(lambda_W2s, T2_, fuel)
        q_lambda_W2 = q_lamda(lambda_W2, T2_, fuel)
        sigmma_rk = PI_lambda_W2s / PI_lambda_W2

        betta2_ = degrees(asin((G1_sm * sqrt(T2_)) / (m_(T2_, fuel, sch) * P1w_sm_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

        coef_t_rab_cold_ = (profile_rab['B_rab'] * 1e-3) / (profile_rab['b_rab'] * 1e-3)
        coef_S_rab_cold_ = (profile_rab['t_rab'] * 1e-3) / (profile_rab['b_rab'] * 1e-3)

        if method_bandage == False:
            G_gap_rab_gas = (0.5 * 0.006 * G1_sm) / (coef_t_rab_cold_ * sin(radians(betta2_)))
        else:
            D_y = 0.5 * (((geom[2]['Dp_sopl_inl_i'][n] + geom[2]['Dp_sopl_out_i'][n]) / 2) + ((geom[2]['Dp_rab_inl_i'][n] + geom[2]['Dp_rab_out_i'][n]) / 2)) + 0.005
            delta_y = ((geom[2]['h_rab_inl_i'][n] + geom[2]['h_rab_out_i'][n]) / 2) * 0.005
            F_gap = pi * D_y * delta_y
            G_gap_rab_gas = 0.6 * F_gap * sqrt((P1_sm**2 - P2s**2) / (4 * P1_sm * V1_sm))        
        
        G_rab_gas = G1_sm - G_gap_rab_gas
        betta2 = degrees(asin((G_rab_gas * sqrt(T2_)) / (m_(T2_, fuel, sch) * P1w_sm_ * geom[2]['F_2_out_i'][n] * q_lambda_W2 * sigmma_rk))) 

        Re_rab = (W2 * geom[2]['S_rab_i'][n] * density2_) / (coef_S_rab_cold_ * Mu(T2_))

        Sgr = (sin(radians(betta1_sm)) / sin(radians(betta2))) * sqrt(((2 * coef_S_rab_cold_) / (coef_t_rab_cold_ * sin(radians(betta1_sm + betta2)) * (cos(radians((betta1_sm - betta2)/2))**2))) - 1)
        tetta = ((geom[2]['Dsr_rab_inl_i'][n] + geom[2]['Dsr_rab_out_i'][n]) / 2) / geom[2]['h_rab_out_i'][n]
        Su_rab = U2 / (W2 * tetta)

        T_inl_rab_air_ = T_out_compressor_air_ + 40
        g_rab = (g_standard_rab * 25 * (1 + (0.8 * Su_rab**0.42)) * 2.55 * ((T0_ / T_inl_rab_air_)**0.25 )) / ((Re_rab**0.34) * (Sgr**0.58) * coef_t_rab_cold_ * sin(radians(betta2))) 
        I2_compressor_ = fuel_inter(air.get('h'), air.get('temperature_K'), T_inl_rab_air_)
        I_inl_rab_air_w_ = I2_compressor_ + 40 * fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_) - (((U2**2) * (((tetta - 1) / tetta)**2)) / 2e3)
        T_inl_rab_air_w_ = T(I_inl_rab_air_w_, fuel)

        T_metal_rab = T_metal_rab + 273.15
        T_rab_air = 950

        G_cold_rab_ = ((T_metal_rab - T_rab_air) / (T_rab_air - T_inl_rab_air_w_)) * g_rab * 1.1
        G_cold_rab = sch[2]['G_k'] * round(G_cold_rab_, 3)
        G2_sm = G1_sm + G_cold_rab

        coef_alpha_rab = (1.12 * 10**(-2)) * (g_rab / g_standard_rab) * ((T_inl_rab_air_w_ / T0_)**0.25)
        q_cold_rab_ = coef_alpha_rab * ((T_metal_rab - T_rab_air) * Cp(T0_, fuel) / distr[0]['Hs_st_i'][n])
        coef_tau_rab = q_cold_rab_ / (distr[0]['ro_st_i'][n] + (W1_sm**2 / (2e3 * distr[0]['Hs_st_i'][n])))
        n_2 = 1 / (1 - (((k(T2_, fuel) - 1) / k(T2_, fuel)) * (psi**2 + coef_tau_rab)))

    # Точка 2 cold
        ################################################################
        T2_cold = T1w_sm_ * ((P2s / P1w_sm_)**((n_2 - 1) / n_2))
        P2s_cold = P2s
        I2_cold = I(T2_cold, fuel)
        ################################################################

        Hs_rk_cold = (I1_sm - I2_cold) - (q_cold_rab_ * distr[0]['Hs_st_i'][n])
        W2_gas = psi * sqrt(W1_sm**2 + 2e3 * Hs_rk_cold)
        C2_gas = sqrt(W2_gas**2 + U2**2 - 2 * W2_gas * U2 * cos(radians(betta2)))

        # Точка 2* gas
        ################################################################
        I2_gas_ = I2_cold + ((1 / psi**2) - 1) * (W2_gas**2 / 2e3) + (C2_gas**2 / 2e3)
        T2_gas_ = T(I2_gas_, fuel)
        ################################################################

        _G_cold_rab_ = G_cold_rab_ * (sch[2]['G_k'] / sch[2]['G_t_notcool'])

        delta_T_air_rab = (((q_cold_rab_ * distr[0]['Hs_st_i'][n]) / (_G_cold_rab_ *  fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_w_))) + ((U2**2) * (2 - (((tetta - 1) / tetta)**2)) /(2e3 *  fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_w_))) - ((0.7 * W2_gas * U2 * cos(radians(betta2))) / (1e3 *  fuel_inter(air.get('Cp'), air.get('temperature_K'), T_inl_rab_air_w_))))
        T_out_rab_air_w_ = T_inl_rab_air_w_ + delta_T_air_rab
        I_out_rab_air_w_ = I(T_out_rab_air_w_, fuel)

        betta2_air = betta2
        ksi_a_rab = (1 + (G_cold_rab_ * 0.7 * (sin(radians(betta2_air))/sin(radians(betta2))))) / (1 + G_cold_rab_)
        ksi_u_rab = (1 + (G_cold_rab_ * 0.7 * (sin(radians(betta2_air))/sin(radians(betta2))))) / (1 + G_cold_rab_)  
        ksi_rab = (1 + (((T_out_rab_air_w_ / T2_gas_)**1.25) * G_cold_rab_)) / (1 + G_cold_rab_)

        W2_sm = W2_gas * sqrt((ksi_a_rab  * sin(radians(betta2)))**2 + (ksi_u_rab * cos(radians(betta2)))**2)
        betta2_sm =  degrees(atan((ksi_a_rab / ksi_u_rab) * tan(radians(betta2))))
        W2u_sm = W2_sm * cos(radians(betta2_sm))
        W2a_sm = W2_sm * sin(radians(betta2_sm))
        C2_sm = sqrt(W2_sm**2 + U2**2 - (2 * U2 * W2_sm * cos(radians(betta2_sm))))

        # Точка 2 sm
        ################################################################
        I2_sm = I(T2_gas_ * ksi_rab, fuel) - (C2_sm**2 / 2e3)
        T2_sm = T(I2_sm, fuel)
        P2_sm = P2s
        V2_sm = sch[1]['R_g'] * 1e3 * T2_sm / P2_sm   
        S2_sm = S(T2_sm, fuel)   
        point2_sm = {'I2_sm':I2_sm, 'T2_sm':T2_sm, 'S2_sm':S2_sm, 'P2_sm':P2_sm, 'V2_sm':V2_sm}
        ################################################################

        # Точка 2 sm*
        ################################################################
        T2_sm_ = T2_gas_ * ksi_rab
        I2_sm_ = I(T2_sm_, fuel)
        P2_sm_ = P2s * (T2_sm_ / T2_sm)**(k(T2_sm, fuel) / (k(T2_sm, fuel) - 1))
        V2_sm_ = sch[1]['R_g'] * 1e3 * T2_sm / P2_sm   
        S2_sm_ = S(T2_sm, fuel)   
        point2_sm_ = {'I2_sm_':I2_sm_, 'T2_sm_':T2_sm_, 'S2_sm_':S2_sm_, 'P2_sm_':P2_sm_, 'V2_sm_':V2_sm_}
        ################################################################

        if  W2_sm * cos(radians(betta2)) >= U2:
            alpha2_sm = degrees(asin((W2_sm * sin(radians(betta2_sm))) / C2_sm))
        else:
            alpha2_sm = 180 - degrees(asin((W2_sm * sin(radians(betta2_sm))) / C2_sm))

        # Точка 2w* sm
        ################################################################  
        I2w_sm_ = I2_sm + (W2_sm**2) / 2e3
        T2w_sm_ = T(I2w_sm_, fuel)
        P2w_sm_ = P2_sm * (T2w_sm_ / T2_sm)**(k(T2_sm, fuel) / (k(T2_sm, fuel) - 1))
        V2w_sm_ = sch[1]['R_g'] * 1e3 * T2w_sm_ / P2w_sm_ 
        S2w_sm_ = S(I2w_sm_, fuel)  
        point2w_sm_ = {'I2w_sm_':I2w_sm_, 'T2w_sm_':T2w_sm_, 'S2w_sm_':S2w_sm_, 'P2w_sm_':P2w_sm_, 'V2w_sm_':V2w_sm_}    
        ################################################################

        # Точка 2s sm_
        ################################################################       
        I2s_sm_ = I(T1_sm, fuel) - Hs_rk
        T2s_sm_ = T(I2s_sm_, fuel) 
        P2s_sm_ = P2s
        V2s_sm_ = sch[1]['R_g'] * 1e3 * T2s_sm_ / P2s_sm_ 
        S2s_sm_ = S(T2s_sm_, fuel)  
        point2s_sm_ = {'I2s_sm_':I2s_sm_, 'T2s_sm_':T2s_sm_, 'S2s_sm_':S2s_sm_, 'P2s_sm_':P2s_sm_, 'V2s_sm_':V2s_sm_}
        ################################################################

        C1u_sm = C1_sm * cos(radians(alpha1_sm))
        C1a_sm = C1_sm * sin(radians(alpha1_sm))
        C2u_sm = C2_sm * cos(radians(alpha2_sm))
        C2a_sm = C2_sm * sin(radians(alpha2_sm))

        a2c = sqrt(k(T2_sm, fuel) * sch[1]['R_g'] * 1e3 * T2_sm)
        M2c = C2_sm / a2c
        a2w = sqrt(k(T2_sm, fuel) * sch[1]['R_g'] * 1e3 * T2_sm)
        M2w = W2_sm / a2w

        profile_rab = profiling(D_k = geom[2]['Dk_rab_out_i'][n] + geom[2]['h_rab_out_i'][n], width = geom[2]['S_rab_i'][n], height = geom[2]['h_rab_out_i'][n], 
        alpha_0 = betta1_sm, alpha_1 = betta2_sm, M_c1 = M2w, value_1 = value1_rab, 
        value_2 = value2_rab, value_3 = value3_rab, method_2 = 'rab')

        psi_rab = losses(fuel, temperature = T2_sm, pressure = P2_sm, velocity_inlet = W1_sm, 
        velocity_outlet = W2_sm, R = sch[1]['R_g'] * 1e3, M_inlet = M1w, M_outlet = M2w,
        lamda_velocity = lamda(W2_sm, T2_sm_, fuel, sch), height = geom[2]['h_rab_out_i'][n], D_k = geom[2]['Dk_rab_out_i'][n], 
        pitch = profile_rab['t_rab'] * 1e-3, width = profile_rab['B_rab'] * 1e-3, chord = profile_rab['b_rab'] * 1e-3, 
        te_radius = profile_rab['r_out_rab'] * 1e-3, C_max = profile_rab['Cmax_rab'] * 1e-3, 
        a_inlet = profile_rab['a_inl_rab'] * 1e-3, a_outlet = profile_rab['a_out_rab'] * 1e-3, 
        blade_inlet_angle = profile_rab['betta0sc_rab'], inlet_angle = betta1_sm, outlet_angle = betta2_sm, 
        design_inc = profile_rab['betta0sc_rab'] - betta1_sm, 
        coef = coef_rab, B = B_rab, SorU = SorU_rab, nu = Mu(T2_sm), ks = ks_rab, method_1 = method_losses_rab, method_2 = 'rab')
        
        psi = psi_rab[0] 

    Ls_st_ = I(T0_, fuel) - I(T2_sm_, fuel)
    Lu_sm = (G1_sm * U1 * C1u_sm + G2_sm * U2 * C2u_sm) / (sch[2]['G_t_notcool'] * 1e3)
    
    delta_L_lake_sopl = (G_gap_sopl_gas * Lu_sm) / sch[2]['G_t_notcool']
    if method_bandage == False:
        delta_L_lake_rab = ksi_leaks(((geom[2]['Dsr_rab_inl_i'][n] + geom[2]['Dsr_rab_out_i'][n]) / 2) / geom[2]['h_rab_out_i'][n], distr[0]['ro_st_i'][n]) * Lu_sm * (geom[2]['delta_r_i'][n] / geom[2]['h_rab_out_i'][n])
    else:
        delta_L_lake_rab = (G_gap_rab_gas* Lu_sm) / sch[2]['G_t_notcool']
    delta_L_tr = (1.701 / G1_sm) * (((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2)**2) * (((pi * ((geom[2]['Dk_rab_inl_i'][n] + geom[2]['Dk_rab_out_i'][n]) / 2) * distr[0]['periodicity']) / 100)**3) * (((1 / V1_sm) + (1 / V2_sm)) / 2)
    delta_L_n = 0.6 * (((U1 + U2) / 2)**2) /1e3
    delta_Lvs = ((C2_sm**2 )/ 2e3) * (G2_sm / sch[2]['G_t_notcool'])
    delta_L_sopl = I(T1, fuel) - I(T1s, fuel)
    delta_L_rab = I(T2_sm, fuel) - I(T1_sm * (P2_sm / P1_sm)**((k(T1_sm, fuel) - 1) / k(T1_sm, fuel)), fuel)

    L_st_ = Lu_sm - delta_L_lake_sopl - delta_L_lake_rab - delta_L_n * round(_G_cold_rab_, 3) - delta_L_tr


    lambda_out_air_sopl = lamda(C1_gas * 0.6, T_out_sopl_air_, fuel, sch)
    PI_lambda_out_air_sopl = PI_lamda(lambda_out_air_sopl, T2_, fuel)
    P_out_air_sopl_ = P1 / PI_lambda_out_air_sopl

    T2_air_sopl_ = T_out_sopl_air_ * (P2_sm_ / P_out_air_sopl_)**((k(T_out_sopl_air_, fuel) - 1) / (k(T_out_sopl_air_, fuel)))
    I2_air_sopl_ = I(T2_air_sopl_, fuel)
    Hs_air_sopl_ = I_out_sopl_air_ - I2_air_sopl_

    lambda_out_air_rab = lamda(W2_gas * 0.7, T_out_rab_air_w_, fuel, sch)
    PI_lambda_out_air_rab = PI_lamda(lambda_out_air_rab, T_out_rab_air_w_, fuel)
    P_out_air_rab_ = P2s / PI_lambda_out_air_rab 

    Hs_air_rab = I_out_sopl_air_ - I2_air_sopl_

    T2_air_rab_ = T_out_rab_air_w_ * (P2_sm_ / P_out_air_rab_)**((k(T_out_rab_air_w_, fuel) - 1) / (k(T_out_rab_air_w_, fuel)))
    I2_air_rab_ = I(T2_air_rab_, fuel)
    Hs_air_rab_ = I_out_rab_air_w_ - I2_air_rab_
    
    sumH_st = _G_cold_sopl_ * Hs_air_sopl_ + _G_cold_rab_ * Hs_air_rab_

    C_fict = sqrt(distr[0]['Hs_st_i'][n]*2e3)
    U_Cfict = (U1 + U2) / (2 * C_fict)

    etta_st = L_st_ / (distr[0]['Hs_st_i'][n] - delta_Lvs + sumH_st)
    N_i = Lu_sm * G2_sm * 1e3
    
    parametrs = {'Hs_st_':distr[0]['Hs_st_i'][n], 'ro_st':distr[0]['ro_st_i'][n],'Hs_sa':Hs_sa, 'H_sa':H_sa_cold, 'Hs_rk':Hs_rk, 'H_rk':Hs_rk_cold, 
    'T0_':T0_, 'P0_':P0_, 'T1s':T1s, 'T1':T1_sm, 'P1':P1_sm, 'T1_':T1_sm_, 'P1_':P1_sm_, 'T1w_':T1w_sm_, 'P1w_':P1w_sm_, 'T2s':T2s, 'T2':T2_sm, 'P2':P2_sm,
    'T2_':T2_sm_, 'P2_':P2_sm_, 'T2w_':T2w_sm_, 'P2w_':P2w_sm_, 'T2s_':T2s_sm_, 'C1s':C1s, 'C1':C1_sm, 'C1u':C1u_sm, 'C1a':C1a_sm, 'W1':W1_sm, 'W1u':W1u_sm, 'W1a':W1a_sm,'U1':U1,
    'alpha1':alpha1_sm, 'betta1':betta1_sm, 'C2':C2_sm, 'C2u':C2u_sm, 'C2a':C2a_sm,'W2s':W2s, 'W2':W2_sm, 'W2u':W2u_sm, 'W2a':W2a_sm, 'U2':U2, 'alpha2':alpha2_sm, 'betta2':betta2_sm,
    'fi':fi_sopl[0], 'psi':psi_rab[0], 'Y_s_sopl':fi_sopl[1], 'Y_s_rab':psi_rab[1], 'Y_p_sopl':fi_sopl[2],'Y_ p_rab':psi_rab[2], 'Y_sec_sopl':fi_sopl[3], 'Y_sec_rab':psi_rab[3], 'Y_tl_sopl':fi_sopl[4], 'Y_tl_rab':psi_rab[4],
    'Y_te_sopl':fi_sopl[5], 'Y_te_rab':psi_rab[5], 'Y_cl_sopl':fi_sopl[6], 'Y_cl_rab':psi_rab[6], 'M1c':M1c, 'M2c':M2c, 'M1w':M1w, 'M2w':M2w,
    'G0':consumption, 'G_g_sopl_gas':G_gap_sopl_gas, 'G_cold_sopl':G_cold_sopl, 'G1':G1_sm, 'G_g_rab_gas':G_gap_rab_gas, 'G_cold_rab':G_cold_rab, 'G2':G2_sm,
    'Ls_st_':Ls_st_, 'dL_lake_sopl':delta_L_lake_sopl, 'dL_lake_rab':delta_L_lake_rab,
    'dL_tr':delta_L_tr, 'dL_sopl':delta_L_sopl, 'dL_rab':delta_L_rab, 'dL_vs':delta_Lvs, 'L_st_':L_st_, 'Lu':Lu_sm, 'C_fict':C_fict, 'U_Cfict':U_Cfict, 'etta_st':etta_st, 'N_i':N_i}

    points = [point0_, point1s, point1_sm, point1_sm_, point1w_sm_, point2s, point2_sm, point2_sm_, point2s_sm_, point2w_sm_]
    velocity_triangle = [C1_sm, W1_sm, U1, alpha1_sm, betta1_sm, C2_sm, W2_sm, U2, alpha2_sm, betta2_sm] 
    
    # termod_param_sopl = {'T1s_cold':T1s_cold, 'P1s_cold':P1s_cold, 'T_inl_sopl_air_':T_inl_sopl_air_, 'T_out_sopl_air_':T_out_sopl_air_,'P_out_air_sopl_':P_out_air_sopl_,  'q_cold_sopl_':q_cold_sopl_, 'g_sopl':g_sopl}
    # termod_param_rab =  {'T2_cold':T2_cold, 'P2s_cold': P2s_cold, 'T_inl_rab_air_':T_inl_rab_air_, 'T_out_rab_air_w_':T_out_rab_air_w_, 'P_out_air_rab_':P_out_air_rab_,  'q_cold_rab_':q_cold_rab_, 'g_rab':g_rab}
    return parametrs, points, velocity_triangle, n, G2_sm, profile_sopl, profile_rab
   
    # termod_param_sopl = {'Hs_sa':Hs_sa, 'H_sa_cold':H_sa_cold, 'T0_':T0_, 'P0_':P0_, 'T1s':T1s, 'T1_sm':T1_sm, 'P1_sm':P1_sm, 'T1_sm_':T1_sm_, 'P1_sm_':P1_sm_, 'T1w_sm_':T1w_sm_, 'P1w_sm_':P1w_sm_, 'T1s_cold':T1s_cold, 'P1s_cold':P1s_cold, 'T_inl_sopl_air_':T_inl_sopl_air_, 'T_out_sopl_air_':T_out_sopl_air_,'P_out_air_sopl_':P_out_air_sopl_, 'sigmma_sa': sigmma_sa, 'G_gap_sopl_gas':G_gap_sopl_gas, 'G_cold_sopl':G_cold_sopl, 'G1_sm':G1_sm, 'q_cold_sopl_':q_cold_sopl_, 'g_sopl':g_sopl}
    # kinematic_param_sopl = {'C1s':C1s, 'fi':fi, 'C1':C1, 'C1_sm':C1_sm, 'C1u_sm':C1u_sm, 'C1a_sm':C1a_sm, 'W1_sm':W1_sm, 'W1u_sm':W1u_sm, 'W1a_sm':W1a_sm, 'U1':U1, 'M1c':M1c, 'M1w':M1w, 'alpha1_sm':alpha1_sm, 'betta1_sm':betta1_sm}
    # termod_param_rab =  {'Hs_rk':Hs_rk, 'Hs_rk_cold':Hs_rk_cold, 'T2s':T2s, 'T2_sm':T2_sm, 'P2_sm':P2_sm, 'T2_sm_':T2_sm_, 'P2_sm_':P2_sm_, 'T2w_sm_':T2w_sm_, 'P2w_sm_':P2w_sm_, 'T2s_sm_':T2s_sm_, 'T2_cold':T2_cold, 'P2s_cold': P2s_cold, 'T_inl_rab_air_':T_inl_rab_air_, 'T_out_rab_air_w_':T_out_rab_air_w_, 'P_out_air_rab_':P_out_air_rab_, 'sigmma_rk': sigmma_rk, 'G_gap_rab_gas':G_gap_rab_gas, 'G_cold_rab':G_cold_rab, 'G2_sm':G2_sm, 'q_cold_rab_':q_cold_rab_, 'g_rab':g_rab}
    # kinematic_param_rab = {'W2s':W2s, 'psi':psi, 'W2':W2, 'W2_sm':W2_sm, 'W2u_sm':W2u_sm, 'W2a_sm':W2a_sm, 'C2_sm':C2_sm, 'C2u_sm':C2u_sm, 'C2a_sm':C2a_sm, 'U2':U2, 'M2c':M2c, 'M2w':M2w, 'alpha2_sm':alpha2_sm, 'betta2_sm':betta2_sm}
    # return termod_param_sopl, kinematic_param_sopl, termod_param_rab, kinematic_param_rab, points, n


# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# schemegtu = scheme(fuel = fuel, ele_power = 153*10e5, t_c = 1060,
#                                 t_a = 15, P_a = 100.5*10e2, epsilon = 10.9, coefficient_pressure_loss = 0.97,
#                                 etta_c_c = 0.995, etta_mch = 0.995, etta_e_g = 0.982,
#                                 etta_is_t = 0.89, etta_is_k = 0.87, leak_rate = 0.005, t_w = 850, z = 4, method = 'cold')

# geometrygtu = geometry(fuel = fuel, sch = schemegtu, number_of_steps = 4, axial_speed_input = 180, axial_speed_outlet = 300,
#              D_sr__h_2_z = 3.5, K_s = 0.065, K_r = 0.06, radial_clearance = 0.0008, method = 'root')

# parametergtu = parameters(fuel = fuel, sch = schemegtu, geom = geometrygtu, alpha_02_i = [65, 80, 85, 90], delta_H = [0.5, 0.55, 0.45], periodicity = 50)


# stages = stagecold(fuel = fuel, sch = schemegtu, geom = geometrygtu, distr = parametergtu, consumption = schemegtu[2]['G_t_notcool'],
#         value1_sopl = 2, value1_rab = 2, 
#         value2_sopl = 0.79, value2_rab = 1, value3_sopl = 0.001455, value3_rab = 0.00133,
#         coef_sopl = 0.005, coef_rab = 0.005, B_sopl = 0.25, B_rab = 0.25,
#         SorU_sopl = 0, SorU_rab = 1, ks_sopl = 1e-6, ks_rab = 1e-6,  
#         g_standard_sopl = 0.032, g_standard_rab = 0.017,
#         T_metal_sopl = 900, T_metal_rab = 600, 
#         method_losses_sopl = 'ANM', method_losses_rab = 'ANM', method_bandage = False, n = 0)

# print(stage_cold_table(stages[0], method = "parameters"))

# hs_plot(point0_ = stages[1][0], point1s = stages[1][1], point1 = stages[1][2], point1_ = stages[1][3], point1w_ = stages[1][4], 
#         point2s = stages[1][5], point2 = stages[1][6], point2_ = stages[1][7], point2s_ = stages[1][8], point2w_ = stages[1][9], i = stages[3], method = 'cold')

# velocity_triangle_plot(C_1 = stages[2][0], W_1 = stages[2][1], U_1 = stages[2][2], alpha_1 = stages[2][3], betta_1 = stages[2][4],
#                         C_2 = stages[2][5], W_2 = stages[2][6], U_2 = stages[2][7], alpha_2 = stages[2][8], betta_2 = stages[2][9], i = stages[3])


    # coef_S_sopl_cold_ = 0.6 # относительная ширина сопловой решетки 
    # coef_S_rab_cold_ = 0.85 # относительная ширина сопловой решетки 

    # coef_t_sopl_cold_ = 0.88 # относительный шаг сопловой решетки 
    # coef_t_rab_cold_ = 0.79  # относительный шаг рабочей решетки 

    # g_standard_sopl = 0.032
    # g_standard_rab = 0.017
    
    # K_gap_sopl = 1
    # K_gap_rab = 1.1