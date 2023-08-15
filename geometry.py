from scheme import *
from geometryPlot import *
from geometryTable import *

# @st.cache
def geometry(fuel, sch, number_of_steps, axial_speed_input, axial_speed_outlet, 
             D_sr__h_2_z, K_s, K_r, radial_clearance, method = None):
        
        '''

        method - метод расчета геометрии канала

        root - постоянство корневого диаметра
        
        medium - постоянство среднего диаметра

        top- постоянство периферийного диаметра
    

        '''  
        
        lambda_0 = (axial_speed_input / sqrt((2 * fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) * sch[1]['R_g'] * (sch[1]['t_c'] + 273.15) * 1000) / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) + 1)))
        lambda_z = (axial_speed_outlet / sqrt((2 * fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) * sch[1]['R_g'] * (sch[1]['t_d'] + 273.15) * 1000) / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) + 1)))
        q_0 = (lambda_0 * ((1 - ((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) - 1) / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) + 1)) * lambda_0**2)** (1 / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) - 1))) * ((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) + 1) / 2)**(1 / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) - 1)))
        q_z = (lambda_z * ((1 - ((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) - 1) / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) + 1)) * lambda_z**2)** (1 / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) - 1))) * ((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) + 1) / 2)**(1 / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) - 1)))
        F_0 = (sch[2]['G_t'] * sqrt (sch[1]['t_c'] + 273.15) / (sch[1]['P_c'] * q_0 * sqrt((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) / sch[1]['R_g'] * 10**(-3)) * (2 / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) + 1))**((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) + 1) / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_c']) - 1)))))
        F2_z = (sch[2]['G_t'] * sqrt(sch[1]['t_d'] + 273.15) / (sch[1]['P_d'] * q_z * sqrt((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) / sch[1]['R_g'] * 10**(-3)) * (2 / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) + 1))**((fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) + 1) / (fuel_inter(fuel.get('k'), fuel.get('temperature_C'), sch[1]['t_d']) - 1))))) 

        height_2_z = (sqrt(F2_z / (pi * D_sr__h_2_z)))
        D_sr_z = (height_2_z * D_sr__h_2_z)
        D_k_z = (D_sr_z - height_2_z)
        D_p_z = (D_sr_z + height_2_z) 
        height_0 = ((sqrt(D_k_z**2 + 4 * F_0 / pi) - D_k_z) / 2)

        if method == 'root':
            D_k_1 = D_k_z
            D_sr_1 = (D_k_1 + height_0)  
            D_p_1 = (D_k_1 + 2 * height_0)

        if method == 'medium':
            D_k_1 = (D_sr_z - height_0)
            D_sr_1 = D_sr_z
            D_p_1 = (D_sr_z + height_0)

        if method == 'top':
            D_p_1 = (D_p_z)
            D_sr_1 = (D_p_1 - height_0)
            D_k_1 = (D_p_z - 2 * height_0)

        n_steps = [(i) / (number_of_steps) for i in range(number_of_steps + 1)]
        n_steps_dell = np.delete(n_steps, number_of_steps - number_of_steps)
        delta_F_t_i = [np.round((2.94 * (n_steps[i]**1.753) * np.exp(-1.0693 * n_steps[i])),4) for i in range(number_of_steps + 1)]
        delta_F_t_i = np.delete(delta_F_t_i, number_of_steps)
        delta_F_t_i = np.append(delta_F_t_i, 1)
        delta_F_tr_i = np.delete(delta_F_t_i, number_of_steps - number_of_steps)
        F_2_i = [(delta_F_tr_i[i] * (F2_z - F_0) + F_0) for i in range(number_of_steps)]

        if method == 'root':
            height_i = [((sqrt(D_k_z**2 + 4 * F_2_i[i] / pi) - D_k_z) / 2) for i in range(number_of_steps)]
            D_k_i = [D_k_z for i in range(number_of_steps)]
            D_sr_i = [(D_k_i[i] + height_i[i]) for i in range(number_of_steps)]
            D_p_i = [(D_k_i[i] + 2 * height_i[i]) for i in range(number_of_steps)]
            d_k_t = [D_k_i[i] / 2 for i in range(number_of_steps)]
            d_k_t.append(D_k_z / 2)
            d_vt_t = np.insert([(D_k_i[i] + 2 * height_i[i]) / 2 for i in range(number_of_steps)], number_of_steps - number_of_steps, (D_k_z + 2 * height_0) / 2)
        
        if method == 'medium':
            height_i = [((sqrt(D_k_z**2 + 4 * F_2_i[i] / pi) - D_k_z) / 2) for i in range(number_of_steps)]
            D_sr_i = [D_sr_z for i in range(number_of_steps)]
            D_k_i = [D_sr_i[i] - height_i[i] for i in range(number_of_steps)]
            D_p_i = [(D_k_i[i] + 2 * height_i[i]) for i in range(number_of_steps)]  
            d_vt_t = [(D_sr_i[i] + height_i[i]) / 2 for i in range(number_of_steps)]
            d_vt_t = np.insert(d_vt_t, number_of_steps - number_of_steps, D_p_1 / 2)
            d_k_t = [(D_sr_i[i] - height_i[i]) / 2 for i in range(number_of_steps)]
            d_k_t = np.insert(d_k_t, number_of_steps - number_of_steps, D_k_1 / 2)

        if method == 'top':  
            height_i = [((sqrt(D_k_z**2 + 4 * F_2_i[i] / pi) - D_k_z) / 2) for i in range(number_of_steps)]
            D_p_i = [D_p_z for i in range(number_of_steps)]
            D_sr_i = [(D_p_i[i] - height_i[i]) for i in range(number_of_steps)]
            D_k_i = [D_sr_i[i] - height_i[i] for i in range(number_of_steps)]
            d_vt_t = [D_p_i[i] / 2 for i in range(number_of_steps)]
            d_vt_t.append(D_p_z / 2)
            d_k_t = np.insert([(D_p_i[i] - 2 * height_i[i]) / 2 for i in range(number_of_steps)], number_of_steps - number_of_steps, (D_p_1 - 2 * height_0) / 2)

        S_sopl_i = [(K_s * D_sr_i[i]) for i in range(number_of_steps)]
        S_rab_i = [(K_r * D_sr_i[i]) for i in range(number_of_steps)]
        delta_s_i = [0.3 * S_rab_i[i] for i in range(number_of_steps)]
        delta_r_i = [radial_clearance for i in range(number_of_steps)]

        L = [0]
        j = 1
        for i in range(number_of_steps):
            if i < number_of_steps-1:
                L.append(L[j-1] + S_sopl_i[i] + 2 * delta_s_i[i] + S_rab_i[i])
                j = j + 1
            else:
                L.append(L[j-1] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i])
                j = j + 1         
        fp_vt, residuals, rank, sv, rcond = np.polyfit(L, d_vt_t, 5, full = True)
        f_vt = sp.poly1d(fp_vt)
        fp_k, residuals, rank, sv, rcond = np.polyfit(L, d_k_t, 5, full = True)
        f_k = sp.poly1d(fp_k)
        l = np.arange(L[1], L[number_of_steps], 0.001)
        D_vt = f_vt(l)
        D_k = f_k(l)

        F_1_inl_i = [pi * (f_vt(L[i]) * 2)**2 / 4 - pi * (f_k(L[i]) * 2)**2 / 4 for i in range(number_of_steps)]            
        F_1_out_i = [pi * (f_vt(L[i] + S_sopl_i[i]) * 2)**2 / 4 - pi * (f_k(L[i] + S_sopl_i[i]) * 2)**2 / 4 for i in range(number_of_steps)]          
        F_2_inl_i = [pi * (f_vt(L[i] + S_sopl_i[i] + delta_s_i[i]) * 2)**2 / 4 - pi * ((f_k(L[i] + S_sopl_i[i] + delta_s_i[i]))*2)**2 / 4 for i in range(number_of_steps)]
        F_2_out_i = [pi * (f_vt(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]) * 2)**2 / 4 - pi * ((f_k(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]))*2)**2 / 4 for i in range(number_of_steps)]        
        Dk_sopl_inl_i = [f_k(L[i]) * 2 for i in range(number_of_steps)]
        Dk_sopl_out_i = [f_k(L[i] + S_sopl_i[i]) * 2 for i in range(number_of_steps)]   
        Dk_rab_inl_i = [f_k(L[i] + S_sopl_i[i] + delta_s_i[i]) * 2 for i in range(number_of_steps)]
        Dk_rab_out_i = [f_k(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]) * 2 for i in range(number_of_steps)]
        Dp_sopl_inl_i = [f_vt(L[i]) * 2 for i in range(number_of_steps)]
        Dp_sopl_out_i = [f_vt(L[i] + S_sopl_i[i]) * 2 for i in range(number_of_steps)]   
        Dp_rab_inl_i =  [f_vt(L[i] + S_sopl_i[i] + delta_s_i[i]) * 2 for i in range(number_of_steps)]
        Dp_rab_out_i = [f_vt(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]) * 2 for i in range(number_of_steps)]
        Dsr_sopl_inl_i = [(Dk_sopl_inl_i[i] + Dp_sopl_inl_i[i]) / 2 for i in range(number_of_steps)]
        Dsr_sopl_out_i = [(Dk_sopl_out_i[i] + Dp_sopl_out_i[i]) / 2 for i in range(number_of_steps)]
        Dsr_rab_inl_i = [(Dk_rab_inl_i[i] + Dp_rab_inl_i[i]) / 2 for i in range(number_of_steps)]
        Dsr_rab_out_i = [(Dk_rab_out_i[i] + Dp_rab_out_i[i]) / 2 for i in range(number_of_steps)]
        h_sopl_inl_i = [(Dp_sopl_inl_i[i] - Dk_sopl_inl_i[i]) / 2 for i in range(number_of_steps)]
        h_sopl_out_i = [(Dp_sopl_out_i[i] - Dk_sopl_out_i[i]) / 2 for i in range(number_of_steps)]
        h_rab_inl_i = [(Dp_rab_inl_i[i] - Dk_rab_inl_i[i]) / 2 for i in range(number_of_steps)]
        h_rab_out_i = [(Dp_rab_out_i[i] - Dk_rab_out_i[i]) / 2 for i in range(number_of_steps)]
        gamma = [degrees(atan(((height_2_z + delta_r_i[i] - height_0) / ((number_of_steps * ((S_sopl_i[i] + S_rab_i[i] + (2 * number_of_steps - 1) * delta_s_i[i]))))))) for i in range(number_of_steps)]      
        d_i = [(Dsr_rab_out_i[i] / h_rab_out_i[i]) for i in range(number_of_steps)] 

        return ({'number_of_steps':number_of_steps, 'axial_speed_input':axial_speed_input, 'axial_speed_outlet':axial_speed_outlet,
                    'D_sr__h_2_z':D_sr__h_2_z, 'K_s': K_s, 'K_r': K_r, 'radial_clearance': radial_clearance},     
                    {'lambda_0': lambda_0, 'lambda_z':lambda_z, 'q_0': q_0, 'q_z':q_z, 
                    'F_0':F_0, 'F2_z':F2_z, 'D_k_1': D_k_1, 'D_k_z':D_k_z,
                    'D_sr_1': D_sr_1, 'D_sr_z':D_sr_z, 'D_p_1':D_p_1, 'D_p_z':D_p_z,
                    'height_0':height_0, 'height_2_z':height_2_z},                    
                    {'n_steps_dell':n_steps_dell,'F_1_inl_i':F_1_inl_i,'F_1_out_i':F_1_out_i,'F_2_inl_i':F_2_inl_i, 'F_2_out_i':F_2_out_i,
                    'Dk_sopl_inl_i':Dk_sopl_inl_i, 'Dk_sopl_out_i':Dk_sopl_out_i, 'Dk_rab_inl_i': Dk_rab_inl_i, 'Dk_rab_out_i': Dk_rab_out_i,
                    'Dsr_sopl_inl_i':Dsr_sopl_inl_i,'Dsr_sopl_out_i':Dsr_sopl_out_i, 'Dsr_rab_inl_i': Dsr_rab_inl_i, 'Dsr_rab_out_i': Dsr_rab_out_i,
                    'Dp_sopl_inl_i':Dp_sopl_inl_i, 'Dp_sopl_out_i':Dp_sopl_out_i, 'Dp_rab_inl_i': Dp_rab_inl_i, 'Dp_rab_out_i': Dp_rab_out_i,
                    'h_sopl_inl_i':h_sopl_inl_i, 'h_sopl_out_i':h_sopl_out_i, 'h_rab_inl_i': h_rab_inl_i, 'h_rab_out_i': h_rab_out_i,
                    'd_i':d_i, 'S_sopl_i':S_sopl_i, 'S_rab_i':S_rab_i, 'delta_s_i':delta_s_i, 'delta_r_i':delta_r_i,
                    'gamma':gamma, 'd_vt_t':d_vt_t, 'd_k_t':d_k_t, 'method':method},
                    [n_steps_dell, F_1_inl_i, F_1_out_i, F_2_inl_i, F_2_out_i,
                    Dk_sopl_inl_i, Dk_sopl_out_i, Dk_rab_inl_i, Dk_rab_out_i,
                    Dsr_sopl_inl_i, Dsr_sopl_out_i, Dsr_rab_inl_i, Dsr_rab_out_i,
                    Dp_sopl_inl_i, Dp_sopl_out_i, Dp_rab_inl_i, Dp_rab_out_i,
                    h_sopl_inl_i, h_sopl_out_i, h_rab_inl_i, h_rab_out_i,
                    d_i, S_sopl_i, S_rab_i, delta_s_i, delta_r_i, gamma])

# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# schemegtu = scheme(fuel = fuel, ele_power = 153*10e5, t_c = 1060,
#                                 t_a = 15, P_a = 100.5*10e2, epsilon = 10.9, coefficient_pressure_loss = 0.97,
#                                 etta_c_c = 0.995, etta_mch = 0.995, etta_e_g = 0.982,
#                                 etta_is_t = 0.89, etta_is_k = 0.87, leak_rate = 0.005, t_w = 850, z = 4, method = 'notcold')

# geometrygtu = geometry(fuel = fuel, sch = schemegtu, number_of_steps = 4, axial_speed_input = 180, axial_speed_outlet = 252,
#              D_sr__h_2_z = 3.5, K_s = 0.065, K_r = 0.06, radial_clearance = 0.0008, method = 'root')

# print(geometry_table(geometrygtu[0], method='inlet'))
# print(geometry_table(geometrygtu[1], method='geometry'))
# print(geometry_table(geometrygtu[3], method = 'geometrystep'))
# flowpath_plot(geometrygtu)




