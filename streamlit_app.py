from streamlitFunction import *
from stagecold import *
from stageSection import *
from readme import *

st.set_option('deprecation.showPyplotGlobalUse', False)
st.set_page_config(
        page_title = "Program GTU",
        page_icon = Image.open("icon.ico"),
        layout = "wide",
        initial_sidebar_state = "expanded")

if check_password():
    
    panel_global = st.sidebar.radio('Этапы расчета:', ["Руководство пользователя", "I. - Этап расчета ГТУ", "II. - Этап расчета ступеней ГТУ", "II. - Этап расчета по сечениям", "IV. - Этап профилирование"])
    if panel_global == "Руководство пользователя":
        st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins><em>Руководство пользователя</em></ins></h1>", unsafe_allow_html=True)
        readme()

    if panel_global == "I. - Этап расчета ГТУ":

        panel_1 = st.sidebar.radio('Этапы расчета проточной части:', ["1. - Расчет состава топлива", 
                                                        "2. - Расчет схемы ГТУ", 
                                                        "3. - Расчет геометрии ГТУ", 
                                                        "4. - Расчет параметров по ступеням"])
        
        if panel_1 == "1. - Расчет состава топлива":
            st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет состава топлива</ins> </h1>", unsafe_allow_html=True)
            st.header('Ввод исходных данных')
            st.session_state.fuel = None

            with st.form(key ='my_form_1'):
                _CO_, _H_2_, _CH_2_, = st.columns(3)
                _H_2S_, _CO_2_, _O_2_= st.columns(3)
                _CH_4_, _C_2H_4_, _C_2H_6_, = st.columns(3)
                _C_3H_8_, _C_4H_10_ = st.columns(2)

                H_2S_ = count(counter = 'H_2S_', column = _H_2S_, name = "H_2S, %", format = "%f", ke = 'count1')
                CO_2_ = count(counter = 'CO_2_', column = _CO_2_, name = "CO_2, %", format = "%f", ke = 'count2')
                O_2_ = count(counter = 'O_2_', column = _O_2_, name = "O_2, %", format = "%f", ke = 'count3')
                CO_ = count(counter = 'CO_', column = _CO_, name = "CO, %", format = "%f", ke = 'count4')
                H_2_ = count(counter = 'H_2_', column = _H_2_, name = "H_2, %", format = "%f", ke = 'count5')
                CH_2_ = count(counter = 'CH_2_', column = _CH_2_, name = "CH_2, %", format = "%f", ke = 'count6')
                CH_4_ = count(counter = 'CH_4_', column = _CH_4_, name = "CH_4, %", format = "%f", ke = 'count7')
                C_2H_4_ = count(counter = 'C_2H_4_', column = _C_2H_4_, name = "C_2H_4, %", format = "%f", ke = 'count8')
                C_2H_6_ = count(counter = 'C_2H_6_', column = _C_2H_6_, name = "C_2H_6, %", format = "%f", ke = 'count9')
                C_3H_8_ = count(counter = 'C_3H_8_', column = _C_3H_8_, name = "C_3H_8, %", format = "%f", ke = 'count10')
                C_4H_10_ = count(counter = 'C_4H_10_', column = _C_4H_10_, name = "C_4H_10, %", format = "%f", ke = 'count11')
                st.write('H_2S = ', H_2S_, '; ', 'CO_2 = ',CO_2_, '; ', 'O_2 = ', O_2_, '; ', 'CO = ',CO_, '; ', 'H_2 = ', H_2_, '; ', 'CH_2 = ', CH_2_)
                st.write('C_2H_4 = ',C_2H_4_, '; ', 'CH_4 = ',CH_4_,'; ','C_2H_6 = ', C_2H_6_, '; ','C_3H_8 = ', C_3H_8_, '; ','C_4H_10 = ', C_4H_10_)

                if st.form_submit_button('Расчет'):
                    fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                st.session_state.fuel = fuel

                if st.form_submit_button('Показать таблицу параметров'): 
                    st.header('Результаты расчета')
                    st.table(fuel_table(temperature_C = st.session_state.fuel.get('temperature_C'), temperature_K = st.session_state.fuel.get('temperature_K'), u = st.session_state.fuel.get('u'), h = st.session_state.fuel.get('h'), s_0 = st.session_state.fuel.get('s_0'), Cp = st.session_state.fuel.get('Cp'), Cv = st.session_state.fuel.get('Cv'), k = st.session_state.fuel.get('k')))

                if st.form_submit_button('Показать графики параметров'):
                    st.pyplot(fuel_plot(x = st.session_state.fuel.get('temperature_C'), y = st.session_state.fuel.get('h'), method = 'enthalpy'))
                    st.pyplot(fuel_plot(x = st.session_state.fuel.get('temperature_C'), y = st.session_state.fuel.get('s_0'), method = 'entropy'))
                    st.pyplot(fuel_plot(x = st.session_state.fuel.get('temperature_C'), y = st.session_state.fuel.get('Cp'), method = 'capacity'))
                    st.pyplot(fuel_plot(x = st.session_state.fuel.get('temperature_C'), y = st.session_state.fuel.get('k'), method = 'isentrope'))

                if st.form_submit_button('Сохранить расчеты в Excel '):
                    Save_to_file_stage(fuel_table(temperature_C = st.session_state.fuel.get('temperature_C'), temperature_K = st.session_state.fuel.get('temperature_K'), u = st.session_state.fuel.get('u'), h = st.session_state.fuel.get('h'), s_0 = st.session_state.fuel.get('s_0'), Cp = st.session_state.fuel.get('Cp'), Cv = st.session_state.fuel.get('Cv'), k = st.session_state.fuel.get('k')), name = 'fuel_parametrs', extension ='.xlsx')

        if panel_1 == "2. - Расчет схемы ГТУ":
            st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет схемы ГТУ</ins></h1>", unsafe_allow_html=True)
            st.image(Image.open('SchemeGTU.png'))
            select_1 = st.selectbox('Выберите тип схемы ГТУ:', ('Расчет схемы без охлаждения', 'Расчет схемы с охлаждением'))
            st.header('Ввод исходных данных схемы ГТУ')
            st.session_state.count34 = None
            st.session_state.schemegtu = None
            
            if select_1 ==  'Расчет схемы без охлаждения':
                st.session_state.count34 = 'notcold'
                st.session_state.t_w_ = 0.0
                st.session_state.z_ = 0.0
                with st.form(key = 'my_form_2'):
                    _t_a_, _P_a_,  _etta_is_k_= st.columns(3)
                    _t_c_, _epsilon_, _etta_is_t_ = st.columns(3)
                    _ele_power_, _pressure_loss_, _leak_rate_ = st.columns(3)
                    _etta_c_c_, _etta_mch_, _etta_e_g_ =st.columns(3)

                    t_a_ = count(counter = 't_a_', column = _t_a_, name = "Температура на входе в комп., ℃", format = "%f", ke = 'count12')
                    P_a_ = count(counter = 'P_a_', column = _P_a_, name = "Давление на входе, Па", format = "%f", ke = 'count13')
                    etta_is_k_ = count(counter = 'etta_is_k_', column = _etta_is_k_, name = "КПД компрессора, -", format = "%f", ke = 'count14')
                    t_c_ = count(counter = 't_c_', column = _t_c_, name = "Температура на входе в турб., ℃", format = "%f", ke = 'count15')
                    epsilon_ = count(counter = 'epsilon_', column = _epsilon_, name = "Степень сжатия в компрессоре, -", format = "%f", ke = 'count16')
                    etta_is_t_ = count(counter = 'etta_is_t_', column = _etta_is_t_, name = "КПД турбины, -", format = "%f", ke = 'count17')  
                    ele_power_ = count(counter = 'ele_power_', column = _ele_power_, name = "Электрическая мощность ГТУ, Вт", format = "%f", ke = 'count18')
                    pressure_loss_ = count(counter = 'pressure_loss_', column = _pressure_loss_, name = "Коэффициент потерь давления, -", format = "%f", ke = 'count19')
                    leak_rate_ = count(counter = 'leak_rate_', column = _leak_rate_, name = "Коэффициент утечек, -", format = "%f", ke = 'count20')
                    etta_c_c_ = count(counter = 'etta_c_c_', column = _etta_c_c_, name = "Коэффициент теплоты топлива, -", format = "%f", ke = 'count21')
                    etta_mch_ = count(counter = 'etta_mch_', column = _etta_mch_, name = "Механический КПД турбины, -", format = "%f", ke = 'count22')
                    etta_e_g_ = count(counter = 'etta_e_g_', column = _etta_e_g_, name = "КПД генератора, -", format = "%f", ke = 'count23')


                    st.write('t_a = ',t_a_, '; ','P_a = ',P_a_, '; ','etta_is_k = ',etta_is_k_, '; ', 't_c = ',t_c_, '; ', 'epsilon = ',epsilon_, '; ', 'etta_is_t = ',etta_is_t_)  
                    st.write('ele_power = ',ele_power_,'; ','pressure_loss = ',pressure_loss_,'; ', 'leak_rate = ',leak_rate_)
                    st.write( 'etta_c_c = ',etta_c_c_,'; ', 'etta_mch = ',etta_mch_,'; ', 'etta_e_g = ',etta_e_g_)

                    if st.form_submit_button('Расчет '):
                        fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                        st.session_state.fuel = fuel
                        schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                            t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                            etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                            etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                        st.session_state.schemegtu = schemegtu

                        st.pyplot(braiton_plot(point_a = st.session_state.schemegtu[3], point_b_ = st.session_state.schemegtu[4], point_b = st.session_state.schemegtu[5],point_c = st.session_state.schemegtu[6], point_d_ = st.session_state.schemegtu[7], point_d = st.session_state.schemegtu[8]))
                        st.header('Результаты расчета компрессора')
                        st.table(scheme_table(st.session_state.schemegtu[0], method = 'compressor'))
                        st.header('Результаты расчета турбины')
                        st.table(scheme_table(st.session_state.schemegtu[1], method = 'turbine')) 
                        st.header('Экономические показатели ГТУ')
                        st.table(scheme_table(st.session_state.schemegtu[2], method = 'scheme')) 

                    if st.form_submit_button('Сохранить расчеты в Excel'):
                        fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                        st.session_state.fuel = fuel
                        schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                            t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                            etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                            etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                        st.session_state.schemegtu = schemegtu

                        Save_to_file_stage(scheme_table(schemegtu[0], method = 'compressor'), name = 'SchemeGTU(compressor)', extension ='.xlsx')
                        Save_to_file_stage(scheme_table(schemegtu[1], method = 'turbine'), name = 'SchemeGTU(turbine)', extension ='.xlsx')
                        Save_to_file_stage(scheme_table(schemegtu[2], method = 'scheme'), name = 'SchemeGTU(not cold)', extension ='.xlsx')
            
            if select_1 ==  'Расчет схемы с охлаждением':
                st.session_state.count34 = 'cold'
                with st.form(key = 'my_form_3'):
                    _t_a_, _P_a_,  _etta_is_k_= st.columns(3)
                    _t_c_, _epsilon_, _etta_is_t_ = st.columns(3)
                    _ele_power_, _pressure_loss_, _leak_rate_ = st.columns(3)
                    _etta_c_c_, _etta_mch_, _etta_e_g_ =st.columns(3)
                    _t_w_, _z_ = st.columns(2)
                    t_a_ = count(counter = 't_a_', column = _t_a_, name = "Температура на входе в комп., ℃", format = "%f", ke = 'count12')
                    P_a_ = count(counter = 'P_a_', column = _P_a_, name = "Давление на входе, Па", format = "%f", ke = 'count13')
                    etta_is_k_ = count(counter = 'etta_is_k_', column = _etta_is_k_, name = "КПД компрессора, -", format = "%f", ke = 'count14')
                    t_c_ = count(counter = 't_c_', column = _t_c_, name = "Температура на входе в турб., ℃", format = "%f", ke = 'count15')
                    epsilon_ = count(counter = 'epsilon_', column = _epsilon_, name = "Степень сжатия в компрессоре, -", format = "%f", ke = 'count16')
                    etta_is_t_ = count(counter = 'etta_is_t_', column = _etta_is_t_, name = "КПД турбины, -", format = "%f", ke = 'count17')  
                    ele_power_ = count(counter = 'ele_power_', column = _ele_power_, name = "Электрическая мощность ГТУ, Вт", format = "%f", ke = 'count18')
                    pressure_loss_ = count(counter = 'pressure_loss_', column = _pressure_loss_, name = "Коэффициент потерь давления, -", format = "%f", ke = 'count19')
                    leak_rate_ = count(counter = 'leak_rate_', column = _leak_rate_, name = "Коэффициент утечек, -", format = "%f", ke = 'count20')
                    etta_c_c_ = count(counter = 'etta_c_c_', column = _etta_c_c_, name = "Коэффициент теплоты топлива, -", format = "%f", ke = 'count21')
                    etta_mch_ = count(counter = 'etta_mch_', column = _etta_mch_, name = "Механический КПД турбины, -", format = "%f", ke = 'count22')
                    etta_e_g_ = count(counter = 'etta_e_g_', column = _etta_e_g_, name = "КПД генератора, -", format = "%f", ke = 'count23')
                    t_w_ = count(counter = 't_w_', column = _t_w_, name = "Допустимая температура металла лопаток, ℃", format = "%f", ke = 'count24')
                    z_ = count(counter = 'z_', column = _z_, name = "Число ступеней газовой турбины, шт", format = "%f", ke = 'count25')

                    st.write('t_a = ',t_a_, '; ','P_a = ',P_a_, '; ','etta_is_k = ',etta_is_k_, '; ', 't_c = ',t_c_, '; ', 'epsilon = ',epsilon_, '; ', 'etta_is_t = ',etta_is_t_)  
                    st.write('ele_power = ',ele_power_,'; ','pressure_loss = ',pressure_loss_,'; ', 'leak_rate = ',leak_rate_)
                    st.write( 'etta_c_c = ',etta_c_c_,'; ', 'etta_mch = ',etta_mch_,'; ', 'etta_e_g = ',etta_e_g_)
                    st.write( 't_w = ',t_w_,'; ', 'z = ',z_)
                    
                    if st.form_submit_button('Расчет '):
                        fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                        st.session_state.fuel = fuel
                        schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                            t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                            etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                            etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                        st.session_state.schemegtu = schemegtu

                        st.pyplot(braiton_plot(point_a = st.session_state.schemegtu[3], point_b_ = st.session_state.schemegtu[4], point_b = st.session_state.schemegtu[5],point_c = st.session_state.schemegtu[6], point_d_ = st.session_state.schemegtu[7], point_d = st.session_state.schemegtu[8]))
                        st.header('Результаты расчета компрессора')
                        st.table(scheme_table(st.session_state.schemegtu[0], method = 'compressor'))
                        st.header('Результаты расчета турбины')
                        st.table(scheme_table(st.session_state.schemegtu[1], method = 'turbine')) 
                        st.header('Экономические показатели ГТУ')
                        st.table(scheme_table(st.session_state.schemegtu[2], method = 'schemecool')) 

                    if st.form_submit_button('Сохранить расчеты в Excel'):
                        fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                        st.session_state.fuel = fuel
                        schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                            t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                            etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                            etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                        st.session_state.schemegtu = schemegtu
                        Save_to_file_stage(scheme_table(schemegtu[0], method = 'compressor'), name = 'SchemeGTU(compressor)', extension ='.xlsx')
                        Save_to_file_stage(scheme_table(schemegtu[1], method = 'turbine'), name = 'SchemeGTU(turbine)', extension ='.xlsx')
                        Save_to_file_stage(scheme_table(schemegtu[2], method = 'schemecool'), name = 'SchemeGTU(cold)', extension ='.xlsx')

        if panel_1 == "3. - Расчет геометрии ГТУ":
            st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет геометрии ГТУ</ins></h1>", unsafe_allow_html=True)
            st.header('Ввод исходных данных геометрии ГТУ')
            st.session_state.count33 = None
            st.session_state.geometrygtu = None
            with st.form(key = 'my_form_4'): 
                
                _K_s_, _K_r_, _radial_clearance_ = st.columns(3)
                _D_sr__h_2_z_, _number_of_steps_ = st.columns(2)
                _axial_speed_input_, _axial_speed_outlet_ = st.columns(2)

                K_s_ = count(counter = 'K_s_', column = _K_s_, name = "Коэффициент зазора сопловой, -", format = "%f", ke = 'count26')
                K_r_ = count(counter = 'K_r_', column = _K_r_, name = "Коэффициент зазора рабочей, -", format = "%f", ke = 'count27')
                radial_clearance_ = count(counter = 'radial_clearance_', column = _radial_clearance_, name = "Радиальный зазор ступени, м", format = "%f", ke = 'count28')
                D_sr__h_2_z_ = count(counter = 'D_sr__h_2_z_', column = _D_sr__h_2_z_, name = "Относительная высота лопатки последней ступени, -", format = "%f", ke = 'count29')
                number_of_steps_ = count(counter = 'number_of_steps_', column = _number_of_steps_, name = "Колличество ступеней в турбине, -", format = "%g", ke = 'count30')
                axial_speed_input_ = count(counter = 'axial_speed_input_', column = _axial_speed_input_, name = "Окружная скорость на входе в турбину, м/c", format = "%f", ke = 'count31')
                axial_speed_outlet_ = count(counter = 'axial_speed_outlet_', column = _axial_speed_outlet_, name = "Окружная скорость на выходе из турбины, м/c", format = "%f", ke = 'count32')
                st.write('K_s = ',K_s_, '; ','K_r = ',K_r_, '; ','radial_clearance = ',radial_clearance_,'; ','D_sr__h_2_z = ',D_sr__h_2_z_, '; ', 'n= ',number_of_steps_)  
                st.write('speed_input = ',axial_speed_input_,'; ', 'speed_outlet_ = ',axial_speed_outlet_,)
                
                type_1 = st.radio('Выберите тип геометрии ГТУ:', ('Dк = const','Dср = const','Dп = const'))

                if type_1 == 'Dк = const':
                    st.session_state.count33 = 'root'
                
                if type_1 == 'Dср = const':
                    st.session_state.count33 = 'medium'

                if type_1 == 'Dп = const':
                    st.session_state.count33 = 'top'

                if st.form_submit_button('Расчет  '):
                    fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                    st.session_state.fuel = fuel
                    schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                    st.session_state.schemegtu = schemegtu            
                    geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                        axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                        K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                    st.session_state.geometrygtu = geometrygtu

                    st.pyplot(flowpath_plot(st.session_state.geometrygtu))
                    st.header('Введенные параметеры')
                    st.table(geometry_table(st.session_state.geometrygtu[0], method = 'inlet'))
                    st.header('Геометрия проточной части')
                    st.table(geometry_table(st.session_state.geometrygtu[1], method = 'geometry')) 
                    st.header('Разбивка проточной части')
                    st.table(geometry_table(st.session_state.geometrygtu[3], method = 'geometrystep')) 

                if st.form_submit_button('Сохранить расчеты в Excel'):
                    
                    fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                    st.session_state.fuel = fuel
                    schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                    st.session_state.schemegtu = schemegtu
                    geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                    axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                    K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                    st.session_state.geometrygtu = geometrygtu    
                    Save_to_file_stage(geometry_table(st.session_state.geometrygtu[0], method = 'inlet'), name = 'GeometryGTU(Inlet)', extension ='.xlsx')
                    Save_to_file_stage(geometry_table(st.session_state.geometrygtu[1], method = 'geometry'), name = 'GeometryGTU(Geometry)', extension ='.xlsx')
                    Save_to_file_stage(geometry_table(st.session_state.geometrygtu[3], method = 'geometrystep'), name = 'GeometryGTU(GeometryStep)', extension ='.xlsx')

        if panel_1 == "4. - Расчет параметров по ступеням":
            st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Разбивка параметров по ступеням ГТУ</ins></h1>", unsafe_allow_html=True)
            st.header('Ввод исходных данных i-ой ступени ГТУ')

            st.session_state.parametergtu = None
            st.session_state.alpha_ = 0.0
            st.session_state.deltaH_ = 0.0 

            with st.form(key = 'my_form_5'): 
                param_1 = countRAH(int(st.session_state.number_of_steps_))
                st.session_state.alpha_ = param_1[0]
                st.session_state.deltaH_ = param_1[1]
                _periodicity_, _, = st.columns(2)
                periodicity_ = count(counter = 'periodicity_', column = _periodicity_, name = "Частота вращения ротора, 1/c", format = "%f", ke = 'count35')
                if st.form_submit_button('Расчет'):
                    fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                    st.session_state.fuel = fuel
                    schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                    st.session_state.schemegtu = schemegtu            
                    geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                        axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                        K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                    st.session_state.geometrygtu = geometrygtu
                    parametergtu = parameters(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, alpha_02_i = st.session_state.alpha_, delta_H = st.session_state.deltaH_, periodicity = st.session_state.periodicity_)
                    st.session_state.parametergtu = parametergtu            
                    st.header('Основные кимематические параметры i-ой ступен')
                    st.table(distribution_table(st.session_state.parametergtu[2], method = 'kinematics'))
                    st.header('Основные термодинамические параметры i-ой ступен')
                    st.table(distribution_table(st.session_state.parametergtu[3], method = 'termod'))

                if st.form_submit_button('Показать графики параметров'):
                    fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                    st.session_state.fuel = fuel
                    schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                    st.session_state.schemegtu = schemegtu            
                    geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                        axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                        K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                    st.session_state.geometrygtu = geometrygtu
                    
                    parametergtu = parameters(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, alpha_02_i = st.session_state.alpha_, delta_H = st.session_state.deltaH_, periodicity = st.session_state.periodicity_)
                    st.session_state.parametergtu = parametergtu
                    st.pyplot(etta_plot_3D(array = data_etta)) 
                    st.pyplot(etta_plot_2D([c_85, c_86, c_87, c_88, c_89, c_90, c_91, c_92, c_93, c_94], 
                    [mu_85, mu_86, mu_87, mu_88, mu_89, mu_90, mu_91, mu_92, mu_93, mu_94], parametergtu[0]['C2a_i_'], parametergtu[0]['mu_st_i'], parametergtu[0]['etta_st_i']))
                    
                    st.pyplot(ro_plot_3D(array = data_ro_05, method = '05'))
                    st.pyplot(ro_plot_2D([alpha_ro_01_05, alpha_ro_02_05, alpha_ro_03_05, alpha_ro_04_05, 
                                alpha_ro_05_05, alpha_ro_06_05, alpha_ro_07_05],
                                [Y_ro_01_05, Y_ro_02_05, Y_ro_03_05, Y_ro_04_05, 
                                Y_ro_05_05, Y_ro_06_05, Y_ro_07_05], method = '05'))

                    st.pyplot(ro_plot_3D(array = data_ro_08, method = '08'))
                    st.pyplot(ro_plot_2D([alpha_ro_01_08, alpha_ro_02_08, alpha_ro_03_08, alpha_ro_04_08, 
                                alpha_ro_05_08, alpha_ro_06_08, alpha_ro_07_08],
                                [Y_ro_01_08, Y_ro_02_08, Y_ro_03_08, Y_ro_04_08, 
                                Y_ro_05_08, Y_ro_06_08, Y_ro_07_08], method = '08'))
                    
                    st.pyplot(ro_plot_3D(array = data_ro_11, method = '11'))
                    st.pyplot(ro_plot_2D([alpha_ro_01_11, alpha_ro_02_11, alpha_ro_03_11, alpha_ro_04_11, 
                                alpha_ro_05_11, alpha_ro_06_11, alpha_ro_07_11],
                                [Y_ro_01_11, Y_ro_02_11, Y_ro_03_11, Y_ro_04_11, 
                                Y_ro_05_11, Y_ro_06_11, Y_ro_07_11], method = '11'))

                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'ro'))
                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'alpha'))
                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'velocity'))
                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'temperature'))
                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'pressure'))
                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'etta'))
                    st.pyplot(distribution_plot(st.session_state.parametergtu, 'heatdrop'))

                if st.form_submit_button('Сохранить расчеты в Excel'):
                    fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                    st.session_state.fuel = fuel
                    schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                    st.session_state.schemegtu = schemegtu
                    geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                    axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                    K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                    st.session_state.geometrygtu = geometrygtu    
                    parametergtu = parameters(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, alpha_02_i = st.session_state.alpha_, delta_H = st.session_state.deltaH_, periodicity = st.session_state.periodicity_)
                    st.session_state.parametergtu = parametergtu            

                    Save_to_file_stage(distribution_table(st.session_state.parametergtu[2], method = 'kinematics'), name = 'ParameterGTU(kinematics)', extension ='.xlsx')
                    Save_to_file_stage(distribution_table(st.session_state.parametergtu[3], method = 'termod'), name = 'ParameterGTU(termod)', extension ='.xlsx')
    
    if panel_global == "II. - Этап расчета ступеней ГТУ":
        if st.session_state.number_of_steps_ == 0.0:
            st.header('Необходимо выполнить: "1. - Этап расчета ГТУ"')
        else:
            stage_str = [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))]
            stage_str.append('Параметры расчета по ступеням')
            panel_2 = st.sidebar.radio('Этапы расчета ступени:', stage_str)
            
            num = [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))]
            for i in range(int(st.session_state.number_of_steps_)):

                if panel_2 == num[i]:
                    st.session_state.stages = None
                    st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет ступени ГТУ на среднем диаметре</ins></h1>", unsafe_allow_html=True)

                    st.header(f'Ввод исходных данных ступени №{i+1}')
                    with st.form(key = 'my_form_6'): 
                        select_2 = st.selectbox('Выбор типа ступени:', ('Неохлаждаемая ступень', 'Охлаждаемая ступень'))   
                        if st.form_submit_button('Выбрать'):
                            if i == 0:
                                st.session_state.consumption = []
                                st.session_state.paramstage = []

                            st.session_state.value1_sopl_, st.session_state.value1_rab_ = 0.0, 0.0
                            st.session_state.value2_sopl_, st.session_state.value2_rab_ = 0.0, 0.0
                            st.session_state.value3_sopl_, st.session_state.value3_rab_ = 0.0, 0.0
                            st.session_state.method_bandage = 0.0

                            st.session_state.g_standard_sopl, st.session_state.g_standard_rab = 0.0, 0.0
                            st.session_state.T_metal_sopl, st.session_state.T_metal_rab = 0.0, 0.0
                            
                        if select_2 == 'Неохлаждаемая ступень':
                            st.header('Параметры сопловой:') 
                            param_2 = countSRVVV(i, 'sopl')   
                            st.session_state.value1_sopl_ = param_2[0]
                            st.session_state.value2_sopl_ = param_2[1]
                            st.session_state.value3_sopl_ = param_2[2]
                            losses_sopl = st.selectbox('Модель потерь для сопловой решетки', options = ('Ainley and Mathieson losses', 'Craig and Cox losses', 'CIAM losses', 'Denton losses','Soderberg losses'))
                            
                            if st.form_submit_button('Выбрать '):
                                st.session_state.method_losses_sopl = True
                                st.session_state.coef_sopl_ = 0.0
                                st.session_state.B_sopl_ = 0.0
                                st.session_state.ks_sopl_ = 0.0
                                st.session_state.sorU_sopl_ = 0.0

                            if losses_sopl == 'Ainley and Mathieson losses':
                                param_3 = countLosses(i,'ANM', 'sopl')   
                                st.session_state.coef_sopl_ = param_3[0]
                                st.session_state.B_sopl_ = param_3[1]
                                st.session_state.method_losses_sopl = 'ANM'

                            if losses_sopl == 'Craig and Cox losses':
                                param_3 = countLosses(i,'CAC', 'sopl')   
                                st.session_state.ks_sopl_ = param_3[0]
                                st.session_state.method_losses_sopl = 'CAC'

                            if losses_sopl == 'CIAM losses':
                                st.session_state.method_losses_sopl = 'CIAM'

                            if losses_sopl == 'Denton losses':
                                param_3 = countLosses(i,'DN', 'sopl')   
                                st.session_state.sorU_sopl_ = param_3[0]
                                st.session_state.method_losses_sopl = 'DN'

                            if losses_sopl == 'Soderberg losses':
                                st.session_state.method_losses_sopl = 'SDB'                  

                            st.header('Параметры рабочей:') 
                            param_4 = countSRVVV(i, 'rab') 
                            st.session_state.value1_rab_ = param_4[0]
                            st.session_state.value2_rab_ = param_4[1]
                            st.session_state.value3_rab_ = param_4[2]
                            losses_rab = st.selectbox('Модель потерь для рабочей решетки', options = ('Ainley and Mathieson losses', 'Craig and Cox losses', 'CIAM losses', 'Denton losses','Soderberg losses'))
                            select_3 = st.radio('Наличие бандажа:', ('Отсутствует', 'Присутствует'))                   
                            
                            if st.form_submit_button('Выбрать  '):
                                if select_3 == 'Отсутствует':
                                    st.session_state.method_bandage = False
                                if  select_3 == 'Присутствует':
                                    st.session_state.method_bandage = True

                                st.session_state.method_losses_rab = True
                                st.session_state.coef_rab_ = 0.0
                                st.session_state.B_rab_ = 0.0
                                st.session_state.ks_rab_ = 0.0
                                st.session_state.sorU_rab_ = 0.0  

                            if losses_rab == 'Ainley and Mathieson losses':
                                param_4 = countLosses(i,'ANM', 'rab')   
                                st.session_state.coef_rab_ = param_4[0]
                                st.session_state.B_rab_ = param_4[1]
                                st.session_state.method_losses_rab = 'ANM'

                            if losses_rab == 'Craig and Cox losses':
                                param_4 = countLosses(i,'CAC', 'rab')   
                                st.session_state.ks_rab_ = param_4[0]
                                st.session_state.method_losses_rab = 'CAC'
                            
                            if losses_rab == 'CIAM losses':
                                st.session_state.method_losses_rab = 'CIAM'

                            if losses_rab == 'Denton losses':
                                param_4 = countLosses(i,'DN', 'rab')   
                                st.session_state.sorU_rab_ = param_4[0]
                                st.session_state.method_losses_rab = 'DN'

                            if losses_rab == 'Soderberg losses':
                                st.session_state.method_losses_rab = 'SDB' 

                            if st.form_submit_button('Расчет'):

                                fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                                st.session_state.fuel = fuel
                                schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                            t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                            etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                            etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                                st.session_state.schemegtu = schemegtu            
                                geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                                    axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                                    K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                                st.session_state.geometrygtu = geometrygtu
                                parametergtu = parameters(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, alpha_02_i = st.session_state.alpha_, delta_H = st.session_state.deltaH_, periodicity = st.session_state.periodicity_)
                                st.session_state.parametergtu = parametergtu
                                
                                if i == 0:
                                    stages = stage(fuel = st.session_state.fuel , sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, distr = st.session_state.parametergtu, consumption = st.session_state.schemegtu[2]['G_t'], 
                                        value1_sopl = st.session_state.value1_sopl_, value1_rab = st.session_state.value1_rab_ , value2_sopl = st.session_state.value2_sopl_, value2_rab = st.session_state.value2_rab_ , value3_sopl = st.session_state.value3_sopl_, value3_rab = st.session_state.value3_rab_ , 
                                        coef_sopl = st.session_state.coef_sopl_, coef_rab = st.session_state.coef_rab_, B_sopl = st.session_state.B_sopl_, B_rab = st.session_state.B_rab_, SorU_sopl = st.session_state.sorU_sopl_, SorU_rab = st.session_state.sorU_rab_, 
                                        ks_sopl = st.session_state.ks_sopl_, ks_rab = st.session_state.ks_rab_, method_losses_sopl = st.session_state.method_losses_sopl, method_losses_rab = st.session_state.method_losses_rab, method_bandage = st.session_state.method_bandage, n = i)
                                    st.session_state.stages = stages 
                                    
                                    st.session_state.consumption = [] 
                                    st.session_state.consumption.insert(0, st.session_state.stages[4])  
                                    st.write(st.session_state.consumption)  

                                    st.session_state.paramstage = []
                                    st.session_state.paramstage.insert(0, st.session_state.stages[0])
                                    st.write(st.session_state.paramstage)                            

                                else:
                                    stages = stage(fuel = st.session_state.fuel , sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, distr = st.session_state.parametergtu, consumption = st.session_state.consumption[i-1], 
                                            value1_sopl = st.session_state.value1_sopl_, value1_rab = st.session_state.value1_rab_ , value2_sopl = st.session_state.value2_sopl_, value2_rab = st.session_state.value2_rab_ , value3_sopl = st.session_state.value3_sopl_, value3_rab = st.session_state.value3_rab_ , 
                                            coef_sopl = st.session_state.coef_sopl_, coef_rab = st.session_state.coef_rab_, B_sopl = st.session_state.B_sopl_, B_rab = st.session_state.B_rab_, SorU_sopl = st.session_state.sorU_sopl_, SorU_rab = st.session_state.sorU_rab_, 
                                            ks_sopl = st.session_state.ks_sopl_, ks_rab = st.session_state.ks_rab_, method_losses_sopl = st.session_state.method_losses_sopl, method_losses_rab = st.session_state.method_losses_rab, method_bandage =  st.session_state.method_bandage, n = i)
                                    st.session_state.stages = stages
                                    
                                    if len(st.session_state.consumption) == i: 
                                        st.session_state.consumption.insert(i, st.session_state.stages[4])
                                        st.session_state.paramstage.insert(i, st.session_state.stages[0])
                                    else: 
                                        st.session_state.consumption[i] = st.session_state.stages[4]
                                        st.session_state.paramstage[i] = st.session_state.stages[0]
                                    
                                    st.write(st.session_state.consumption)    
                                    st.write(st.session_state.paramstage) 

                                st.header('Термодинамические и кинематические параметры')
                                st.table(stage_table(st.session_state.stages[0], method = 'parameters'))
                                st.header('Геометрические параметры профиля сопловой')
                                st.table(stage_table(st.session_state.stages[5], method = "profile sopl"))
                                st.header('Геометрические параметры профиля рабочей')
                                st.table(stage_table(st.session_state.stages[6], method = "profile rab"))
                                st.pyplot(hs_plot(point0_ = st.session_state.stages[1][0], point1s = st.session_state.stages[1][1], point1 = st.session_state.stages[1][2], point1_ = st.session_state.stages[1][3], point1w_ = st.session_state.stages[1][4], 
                                point2s = st.session_state.stages[1][5], point2 = st.session_state.stages[1][6], point2_ = st.session_state.stages[1][7], point2s_ = st.session_state.stages[1][8], point2w_ = st.session_state.stages[1][9], i = st.session_state.stages[3], method = 'notcold')) 
                                
                                st.pyplot(velocity_triangle_plot(C_1 = st.session_state.stages[2][0], W_1 = st.session_state.stages[2][1], U_1 = st.session_state.stages[2][2], alpha_1 = st.session_state.stages[2][3], betta_1 = st.session_state.stages[2][4], C_2 = st.session_state.stages[2][5], W_2 = st.session_state.stages[2][6], U_2 = st.session_state.stages[2][7], alpha_2 = st.session_state.stages[2][8], betta_2 = st.session_state.stages[2][9], i = st.session_state.stages[3]))
                
                        if select_2 == 'Охлаждаемая ступень':
                            st.header('Параметры сопловой:')
                            param_5 = countSTGKT(i, 'sopl')   
                            st.session_state.g_standard_sopl = param_5[0]
                            st.session_state.T_metal_sopl = param_5[1]
                            param_2 = countSRVVV(i, 'sopl')   
                            st.session_state.value1_sopl_ = param_2[0]
                            st.session_state.value2_sopl_ = param_2[1]
                            st.session_state.value3_sopl_ = param_2[2]
                            losses_sopl = st.selectbox('Модель потерь для сопловой решетки', options = ('Ainley and Mathieson losses', 'Craig and Cox losses', 'CIAM losses', 'Denton losses','Soderberg losses'))
                            
                            if st.form_submit_button('Выбрать '):
                                st.session_state.method_losses_sopl = True
                                st.session_state.coef_sopl_ = 0.0
                                st.session_state.B_sopl_ = 0.0
                                st.session_state.ks_sopl_ = 0.0
                                st.session_state.sorU_sopl_ = 0.0

                            if losses_sopl == 'Ainley and Mathieson losses':
                                param_3 = countLosses(i,'ANM', 'sopl')   
                                st.session_state.coef_sopl_ = param_3[0]
                                st.session_state.B_sopl_ = param_3[1]
                                st.session_state.method_losses_sopl = 'ANM'

                            if losses_sopl == 'Craig and Cox losses':
                                param_3 = countLosses(i,'CAC', 'sopl')   
                                st.session_state.ks_sopl_ = param_3[0]
                                st.session_state.method_losses_sopl = 'CAC'

                            if losses_sopl == 'CIAM losses':
                                st.session_state.method_losses_sopl = 'CIAM'

                            if losses_sopl == 'Denton losses':
                                param_3 = countLosses(i,'DN', 'sopl')   
                                st.session_state.sorU_sopl_ = param_3[0]
                                st.session_state.method_losses_sopl = 'DN'

                            if losses_sopl == 'Soderberg losses':
                                st.session_state.method_losses_sopl = 'SDB'

                            st.header('Параметры рабочей:') 
                            param_6 = countSTGKT(i, 'rab')   
                            st.session_state.g_standard_rab = param_6[0]
                            st.session_state.T_metal_rab = param_6[1]   
                            param_4 = countSRVVV(i, 'rab') 
                            st.session_state.value1_rab_ = param_4[0]
                            st.session_state.value2_rab_ = param_4[1]
                            st.session_state.value3_rab_ = param_4[2]
                            losses_rab = st.selectbox('Модель потерь для рабочей решетки', options = ('Ainley and Mathieson losses', 'Craig and Cox losses', 'CIAM losses', 'Denton losses','Soderberg losses'))
                            select_3 = st.radio('Наличие бандажа:', ('Отсутствует', 'Присутствует'))                   
                            
                            if st.form_submit_button('Выбрать  '):
                                if select_3 == 'Отсутствует':
                                    st.session_state.method_bandage = False
                                if  select_3 == 'Присутствует':
                                    st.session_state.method_bandage = True

                                st.session_state.method_losses_rab = True
                                st.session_state.coef_rab_ = 0.0
                                st.session_state.B_rab_ = 0.0
                                st.session_state.ks_rab_ = 0.0
                                st.session_state.sorU_rab_ = 0.0  

                            if losses_rab == 'Ainley and Mathieson losses':
                                param_4 = countLosses(i,'ANM', 'rab')   
                                st.session_state.coef_rab_ = param_4[0]
                                st.session_state.B_rab_ = param_4[1]
                                st.session_state.method_losses_rab = 'ANM'

                            if losses_rab == 'Craig and Cox losses':
                                param_4 = countLosses(i,'CAC', 'rab')   
                                st.session_state.ks_rab_ = param_4[0]
                                st.session_state.method_losses_rab = 'CAC'
                            
                            if losses_rab == 'CIAM losses':
                                st.session_state.method_losses_rab = 'CIAM'

                            if losses_rab == 'Denton losses':
                                param_4 = countLosses(i,'DN', 'rab')   
                                st.session_state.sorU_rab_ = param_4[0]
                                st.session_state.method_losses_rab = 'DN'

                            if losses_rab == 'Soderberg losses':
                                st.session_state.method_losses_rab = 'SDB' 
                            
                            if st.form_submit_button('Расчет'):

                                fuel = gas(_H_2S = st.session_state.H_2S_, _CO_2 = st.session_state.CO_2_, _O_2 = st.session_state.O_2_, _CO = st.session_state.CO_, _H_2 = st.session_state.H_2_, _CH_2 = st.session_state.CH_2_, _CH_4 = st.session_state.CH_4_, _C_2H_4 = st.session_state.C_2H_4_, _C_2H_6 = st.session_state.C_2H_6_, _C_3H_8 = st.session_state.C_3H_8_, _C_4H_10 = st.session_state.C_4H_10_, temperature_C = temperature_Cels)
                                st.session_state.fuel = fuel
                                schemegtu = scheme(fuel = st.session_state.fuel, ele_power = st.session_state.ele_power_, t_c = st.session_state.t_c_,
                                            t_a = st.session_state.t_a_, P_a = st.session_state.P_a_, epsilon = st.session_state.epsilon_, coefficient_pressure_loss = st.session_state.pressure_loss_,
                                            etta_c_c = st.session_state.etta_c_c_, etta_mch = st.session_state.etta_mch_, etta_e_g = st.session_state.etta_e_g_,
                                            etta_is_t = st.session_state.etta_is_t_, etta_is_k = st.session_state.etta_is_k_, leak_rate = st.session_state.leak_rate_, t_w = st.session_state.t_w_, z = st.session_state.z_, method = st.session_state.count34)
                                st.session_state.schemegtu = schemegtu            
                                geometrygtu = geometry(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, number_of_steps = int(st.session_state.number_of_steps_), axial_speed_input = st.session_state.axial_speed_input_,
                                    axial_speed_outlet = st.session_state.axial_speed_outlet_,D_sr__h_2_z = st.session_state.D_sr__h_2_z_, 
                                    K_s = st.session_state.K_s_, K_r = st.session_state.K_r_, radial_clearance = st.session_state.radial_clearance_, method = st.session_state.count33)
                                st.session_state.geometrygtu = geometrygtu
                                parametergtu = parameters(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, alpha_02_i = st.session_state.alpha_, delta_H = st.session_state.deltaH_, periodicity = st.session_state.periodicity_)
                                st.session_state.parametergtu = parametergtu
                                
                                if i == 0:
                                    stages = stagecold(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, distr = st.session_state.parametergtu, consumption = st.session_state.schemegtu[2]['G_t_notcool'],
                                            value1_sopl = st.session_state.value1_sopl_, value1_rab = st.session_state.value1_rab_, value2_sopl = st.session_state.value2_sopl_, value2_rab = st.session_state.value2_rab_, value3_sopl = st.session_state.value3_sopl_, value3_rab = st.session_state.value3_rab_,
                                            coef_sopl = st.session_state.coef_sopl_, coef_rab = st.session_state.coef_rab_, B_sopl = st.session_state.B_sopl_, B_rab = st.session_state.B_rab_, SorU_sopl = st.session_state.sorU_sopl_, SorU_rab = st.session_state.sorU_rab_, 
                                            ks_sopl = st.session_state.ks_sopl_, ks_rab = st.session_state.ks_rab_, g_standard_sopl = st.session_state.g_standard_sopl, g_standard_rab = st.session_state.g_standard_rab,
                                            T_metal_sopl = st.session_state.T_metal_sopl, T_metal_rab = st.session_state.T_metal_rab, method_losses_sopl = st.session_state.method_losses_sopl, method_losses_rab = st.session_state.method_losses_rab, method_bandage = st.session_state.method_bandage, n = i)
                                    st.session_state.stages = stages 
                                                            
                                    st.session_state.consumption = [] 
                                    st.session_state.consumption.insert(0, st.session_state.stages[4])  
                                    st.write(st.session_state.consumption) 

                                    st.session_state.paramstage = []
                                    st.session_state.paramstage.insert(0, st.session_state.stages[0])
                                    st.write(st.session_state.paramstage)  
                                    
                                else:
                                    stages = stagecold(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, distr = st.session_state.parametergtu, consumption = st.session_state.consumption[i-1],
                                            value1_sopl = st.session_state.value1_sopl_, value1_rab = st.session_state.value1_rab_, value2_sopl = st.session_state.value2_sopl_, value2_rab = st.session_state.value2_rab_, value3_sopl = st.session_state.value3_sopl_, value3_rab = st.session_state.value3_rab_,
                                            coef_sopl = st.session_state.coef_sopl_, coef_rab = st.session_state.coef_rab_, B_sopl = st.session_state.B_sopl_, B_rab = st.session_state.B_rab_, SorU_sopl = st.session_state.sorU_sopl_, SorU_rab = st.session_state.sorU_rab_, 
                                            ks_sopl = st.session_state.ks_sopl_, ks_rab = st.session_state.ks_rab_, g_standard_sopl = st.session_state.g_standard_sopl, g_standard_rab = st.session_state.g_standard_rab,
                                            T_metal_sopl = st.session_state.T_metal_sopl, T_metal_rab = st.session_state.T_metal_rab, method_losses_sopl = st.session_state.method_losses_sopl, method_losses_rab = st.session_state.method_losses_rab, method_bandage = st.session_state.method_bandage, n = i)
                                    st.session_state.stages = stages 

                                    if len(st.session_state.consumption) == i: 
                                        st.session_state.consumption.insert(i, st.session_state.stages[4])
                                        st.session_state.paramstage.insert(i, st.session_state.stages[0])
                                    else: 
                                        st.session_state.consumption[i] = st.session_state.stages[4]
                                        st.session_state.paramstage[i] = st.session_state.stages[0]
                                    st.write(st.session_state.consumption)                                 
                                    st.write(st.session_state.paramstage) 

                                st.header('Термодинамические и кинематические параметры')
                                st.table(stage_cold_table(st.session_state.stages[0], method = 'parameters'))
                                st.header('Геометрические параметры профиля сопловой')
                                st.table(stage_cold_table(st.session_state.stages[5], method = "profile sopl"))
                                st.header('Геометрические параметры профиля рабочей')
                                st.table(stage_cold_table(st.session_state.stages[6], method = "profile rab"))
                                st.pyplot(hs_plot(point0_ = st.session_state.stages[1][0], point1s = st.session_state.stages[1][1], point1 = st.session_state.stages[1][2], point1_ = st.session_state.stages[1][3], point1w_ = st.session_state.stages[1][4], 
                                point2s = st.session_state.stages[1][5], point2 = st.session_state.stages[1][6], point2_ = st.session_state.stages[1][7], point2s_ = st.session_state.stages[1][8], point2w_ = st.session_state.stages[1][9], i = st.session_state.stages[3], method = 'cold')) 
                                st.pyplot(velocity_triangle_plot(C_1 = st.session_state.stages[2][0], W_1 = st.session_state.stages[2][1], U_1 = st.session_state.stages[2][2], alpha_1 = st.session_state.stages[2][3], betta_1 = st.session_state.stages[2][4], C_2 = st.session_state.stages[2][5], W_2 = st.session_state.stages[2][6], U_2 = st.session_state.stages[2][7], alpha_2 = st.session_state.stages[2][8], betta_2 = st.session_state.stages[2][9], i = st.session_state.stages[3]))
                                
            if panel_2 == 'Параметры расчета по ступеням':
                st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Результаты расчета на среднем диаметре</ins></h1>", unsafe_allow_html=True)
                st.session_state.stage_dict = {}
                for n in st.session_state.paramstage[0].keys():
                    st.session_state.stage_dict[n] = list(d[n] for d in st.session_state.paramstage)
                st.session_state.stage_list = list(st.session_state.stage_dict.values()) 
                st.table(stageTable(st.session_state.stage_list))

    if panel_global == "III. - Этап расчета по сечениям":

        panel_3 = st.sidebar.radio('Этапы расчета ступени:', [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))])
        num_1 = [f'Ступень №{i+1}' for i in range(int(st.session_state.number_of_steps_))]
            
        for i in range(int(st.session_state.number_of_steps_)):
                
            if panel_3 == num_1[i]:
                st.session_state.stagessection = None
                st.session_state.method_section = True
                st.session_state.value_num_  = 0.0
                st.markdown("<h1 style='text-align: center; color: #1C2833;'><ins>Расчет ступени ГТУ по сечениям</ins></h1>", unsafe_allow_html=True)

                st.header(f'Ввод исходных данных ступени №{i+1}')
                with st.form(key = 'my_form_7'): 
                    select_4 = st.radio('Выбор закона закрутки:', ['Обратный закон закрутки: 𝑟 ∙ 𝑡𝑔(𝛼1) = 𝑐𝑜𝑛𝑠𝑡', 'Закон постоянства циркуляции: 𝐶1𝑢 ∙ 𝑟𝜑2 = 𝑐𝑜𝑛𝑠𝑡', 'Закон постоянства угла выхода: 𝛼1(𝑟) = 𝑐𝑜𝑛𝑠𝑡'])  
                    if i == 0:
                        st.session_state.data = []

                    if select_4 == 'Обратный закон закрутки: 𝑟 ∙ 𝑡𝑔(𝛼1) = 𝑐𝑜𝑛𝑠𝑡':
                        param_7 = countNum(i)
                        st.session_state.value_num_ = param_7
                        st.session_state.method_section = 'rtgconst' 

                    if select_4 == 'Закон постоянства циркуляции: 𝐶1𝑢 ∙ 𝑟𝜑2 = 𝑐𝑜𝑛𝑠𝑡':
                        param_7 = countNum(i)
                        st.session_state.value_num_ = param_7
                        st.session_state.method_section = 'C1uconst'                   

                    if select_4 == 'Закон постоянства угла выхода: 𝛼1(𝑟) = 𝑐𝑜𝑛𝑠𝑡':
                        param_7 = countNum(i)
                        st.session_state.value_num_ = param_7
                        st.session_state.method_section = 'alpha1const'                            
                        
                    if st.form_submit_button('Расчет'):                        
                        stagessection = spin_laws_stage(fuel = st.session_state.fuel, sch = st.session_state.schemegtu, geom = st.session_state.geometrygtu, 
                        distr = st.session_state.parametergtu, stg = st.session_state.stage_dict, m = int(st.session_state.value_num_), 
                        n = i, method = st.session_state.method_section)
                        st.session_state.stagessection = stagessection

                        st.table(stage_section_table(st.session_state.stagessection[1]))
                        st.pyplot(velocity_triangle_i(C_1_i = st.session_state.stagessection[2][0], W_1_i = st.session_state.stagessection[2][1], U_1_i = st.session_state.stagessection[2][2], alpha_1_i = st.session_state.stagessection[2][3], betta_1_i = st.session_state.stagessection[2][4],
                        C_2_i = st.session_state.stagessection[2][5], W_2_i = st.session_state.stagessection[2][6], U_2_i = st.session_state.stagessection[2][7], alpha_2_i = st.session_state.stagessection[2][8], betta_2_i = st.session_state.stagessection[2][9]))

    if panel_global == "IV. - Этап профилирование":
        st.caption('Находится в стадии разработки')
