import pandas as pd
import numpy as np
import streamlit as st
from PIL import Image

def check_password():
    """Returns `True` if the user had a correct password."""
    def password_entered():
        st.session_state.number_of_steps_ = 0.0
        """Checks whether a password entered by the user is correct."""
        if (st.session_state["username"] in st.secrets["passwords"] and st.session_state["password"] == st.secrets["passwords"][st.session_state["username"]]):
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # don't store username + password
            del st.session_state["username"]
        else:
            st.session_state["password_correct"] = False

    if "password_correct" not in st.session_state:
        # First run, show inputs for username + password.
        st.markdown("<h2 style='text-align: center; color: #1C2833;'><ins><em>Добро пожаловать в программу по проектированию ГТУ (version - 1.0)!</em></ins></h2>", unsafe_allow_html=True)
        col1, col2 = st.columns(2)
        col1.image(Image.open('label.png'))
        col2.markdown("<h7 style='text-align: center; color: #1C2833;'><ins><em>Авторизация пользователя</em></ins></h7>", unsafe_allow_html=True)
        
        col2.text_input("Логин", on_change=password_entered, key="username")
        col2.text_input("Пароль", type="password", on_change=password_entered, key="password")
        st.caption('Разработчик: $Андрианов$ $Дмитрий$ $Михайлович$')
        return False
    
    elif not st.session_state["password_correct"]:
        st.markdown("<h2 style='text-align: center; color: #1C2833;'><ins><em>Добро пожаловать в программу по проектированию ГТУ (version - 1.0)!</em></ins></h2>", unsafe_allow_html=True)
        col1, col2 = st.columns(2)
        col1.image(Image.open('label.png'))
        col2.markdown("<h7 style='text-align: center; color: #1C2833;'><ins><em>Авторизация пользователя</em></ins></h7>", unsafe_allow_html=True)
        # Password not correct, show input + error.
        col2.text_input("Логин", on_change=password_entered, key="username")
        col2.text_input("Пароль", type="password", on_change=password_entered, key="password")
        st.error("Пользователь неизвестен или неверный пароль😕")
        st.caption('Разработчик: $Андрианов$ $Дмитрий$ $Михайлович$')
        return False
    else:
        # Password correct.
        return True

def Save_to_file_stage(stage, name , extension ):
    with pd.ExcelWriter((name + extension)) as writer:
        stage.to_excel(writer, sheet_name='sheet1', index=False,header=False,)

def count(counter, column, name, format, ke):
    if counter not in st.session_state:
        st.session_state[counter] = 0.00
    val = column.number_input((name), format = format, value = st.session_state[counter], key = ke)
    if val:
        st.session_state[counter] = st.session_state[ke]
    return st.session_state[counter]

def countRAH(num):
    alpha_02_index, delta_H_index = [], []
    alpha_02, delta_H = [], []
    for i in range(num):
        alpha_02_index.append(f'alpha_02_{i + 1}')
        delta_H_index.append(f'delta_H_{i + 1}')

    for i in range(num):
         alpha_02_index[i], delta_H_index[i] = st.columns(2)
    for i in range(num):   
        alpha_02.append(count(counter = f'alpha_02{i}', column = alpha_02_index[i], name = f"Угол выхода ст. № {i+1}, град", format = "%f", ke = f'alpha_02_{i + 1}'))
    for i in range(num):  
        delta_H.append(count(counter = f'delta_H{i}', column = delta_H_index[i], name = f"Доля теплоперепада ст. № {i+1}, -", format = "%f", ke = f'delta_H_{i + 1}'))
    return alpha_02, delta_H

def countSTGKT(i, method_2):
    if method_2 == 'sopl':
        g_standard_sopl_index = f'g_sopl_{i + 1}'
        T_metal_sopl_index = f'T_m_sopl{i + 1}'
        g_standard_sopl_index, T_metal_sopl_index  = st.columns(2)
        g_standard_sopl = count(counter = f'g_standard_sopl{i}', column = g_standard_sopl_index, name = f"Коэфф. удельного расхода охл. СР. № {i+1}, -", format = "%f", ke = f'g_standard_sopl_{i + 1}')
        T_metal_sopl = count(counter = f'T_metal_sopl{i}', column = T_metal_sopl_index, name = f"Допустимая температура металла СР. № {i+1}, -", format = "%f", ke = f'T_metal_sopl_{i + 1}')
        return g_standard_sopl, T_metal_sopl
    
    if method_2 == 'rab':
        g_standard_rab_index = f'g_rab_{i + 1}'
        T_metal_rab_index = f'T_m_rab{i + 1}'
        g_standard_rab_index, T_metal_rab_index = st.columns(2)
        g_standard_rab = count(counter = f'g_standard_rab{i}', column = g_standard_rab_index, name = f"Коэфф. удельного расхода охл. РК. № {i+1}, -", format = "%f", ke = f'g_standard_rab_{i + 1}') 
        T_metal_rab = count(counter = f'T_metal_rab{i}', column = T_metal_rab_index, name = f"Допустимая температура металла РК № {i+1}, -", format = "%f", ke = f'T_metal_rab_{i + 1}')
        return g_standard_rab, T_metal_rab

def countSRVVV(i, method_2):

    if method_2 == 'sopl':

        value1_sopl_index_ = f'value1_sopl_{i + 1}'
        value2_sopl_index_ = f'value2_sopl_{i + 1}'
        value3_sopl_index_ = f'value3_sopl_{i + 1}'
        value1_sopl_index_, value2_sopl_index_, value3_sopl_index_ = st.columns(3)
        value1_sopl_ = count(counter = f'value1_sopl{i}', column = value1_sopl_index_, name = f"Эффективный угол СР. № {i+1}, град", format = "%f", ke = f'value1_sopl_{i + 1}')  
        value2_sopl_ = count(counter = f'value2_sopl{i}', column = value2_sopl_index_, name = f"Коэфф. толщины профиля СР. № {i+1}, -", format = "%f", ke = f'value2_sopl_{i + 1}')
        value3_sopl_ = count(counter = f'value3_sopl{i}', column = value3_sopl_index_, name = f"Коэфф. кромки СР. № {i+1}, -", format = "%f", ke = f'value3_sopl_{i + 1}')
        return value1_sopl_, value2_sopl_, value3_sopl_
    
    if method_2 == 'rab':

        value1_rab_index_ = f'value1_rab_{i + 1}'
        value2_rab_index_ = f'value2_rab_{i + 1}'
        value3_rab_index_ = f'value3_rab_{i + 1}'
        value1_rab_index_, value2_rab_index_, value3_rab_index_ = st.columns(3)
        value1_rab_ = count(counter = f'value1_rab{i}', column = value1_rab_index_, name = f"Эффективный угол РК. № {i+1}, град", format = "%f", ke = f'value1_rab_{i + 1}')  
        value2_rab_ = count(counter = f'value2_rab{i}', column = value2_rab_index_, name = f"Коэфф. толщины профиля РК. № {i+1}, -", format = "%f", ke = f'value2_rab_{i + 1}')
        value3_rab_ = count(counter = f'value3_rab{i}', column = value3_rab_index_, name = f"Коэфф. кромки РК. № {i+1}, -", format = "%f", ke = f'value3_rab_{i + 1}')
        return value1_rab_, value2_rab_, value3_rab_

def countLosses(i, method_1, method_2):

    if method_2 == 'sopl':

        if method_1 ==  'ANM':
            coef_sopl_index_ = f'coef_sopl_{i + 1}'
            B_sopl_index_ = f'B_sopl_{i + 1}'
            coef_sopl_index_, B_sopl_index_ = st.columns(2)
            coef_sopl_ = count(counter = f'coef_sopl{i}', column = coef_sopl_index_, name = f"Зазор между ротором. № {i+1}, м", format = "%f", ke = f'coef_sopl_{i + 1}')  
            B_sopl_ = count(counter = f'B_sopl{i}', column = B_sopl_index_, name = f"Зазор между бандажом. № {i+1}, м", format = "%f", ke = f'B_sopl_{i + 1}')
            return coef_sopl_, B_sopl_
        
        if method_1 ==  'CAC':
            ks_sopl_index_ = f'ks_sopl_{i + 1}'
            ks_sopl_index_, _, = st.columns(2)
            ks_sopl_ = count(counter = f'ks_sopl{i}', column = ks_sopl_index_, name = f"Коэфф. шероховатости СР. № {i+1}, -", format = "%f", ke = f'ks_sopl_{i + 1}')  
            return ks_sopl_, _ 

        if method_1 ==  'DN':
            sorU_sopl_index_ = f'sorU_sopl_{i + 1}'
            sorU_sopl_index_, _, = st.columns(2)
            sorU_sopl_ = count(counter = f'sorU_sopl{i}', column = sorU_sopl_index_, name = f"Наличие бандажа № {i+1}, -", format = "%f", ke = f'sorU_sopl_{i + 1}')  
            return sorU_sopl_, _,

    if method_2 ==  'rab':

        if method_1 ==  'ANM':
            coef_rab_index_ = f'coef_rab_{i + 1}'
            B_rab_index_ = f'B_rab_{i + 1}'
            coef_rab_index_, B_rab_index_ = st.columns(2)
            coef_rab_ = count(counter = f'coef_rab{i}', column = coef_rab_index_, name = f"Зазор между ротором. № {i+1}, м", format = "%f", ke = f'coef_rab_{i + 1}')  
            B_rab_ = count(counter = f'B_rab{i}', column = B_rab_index_, name = f"Зазор между бандажом. № {i+1}, м", format = "%f", ke = f'B_rab_{i + 1}')
            return coef_rab_, B_rab_
        
        if method_1 ==  'CAC':
            ks_rab_index_ = f'ks_rab_{i + 1}'
            ks_rab_index_, _, = st.columns(2)
            ks_rab_ = count(counter = f'ks_rab{i}', column = ks_rab_index_, name = f"Коэфф. шероховатости РК. № {i+1}, -", format = "%f", ke = f'ks_rab_{i + 1}')  
            return ks_rab_, _
        
        if method_1 ==  'DN':
            sorU_rab_index_ = f'sorU_rab_{i + 1}'
            sorU_rab_index_, _, = st.columns(2)
            sorU_rab_ = count(counter = f'sorU_rab{i}', column = sorU_rab_index_, name = f"Наличие бандажа № {i+1}, -", format = "%f", ke = f'sorU_rab_{i + 1}')  
            return sorU_rab_, _,

def countNum(i):
    value_num_index_ = f'value1_sopl_{i + 1}'
    value_num_index_, _ = st.columns(2)
    value_num_ = count(counter = f'value_num{i}', column = value_num_index_, name = f"Кол-во сечений ступени. № {i+1}, -", format = "%g", ke = f'value_num_{i + 1}')  
    return value_num_
