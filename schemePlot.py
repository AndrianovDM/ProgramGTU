import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import LinearLocator
from scipy.interpolate import make_interp_spline

def braiton_plot(point_a, point_b_ =0, point_b =0, point_c =0, point_d_=0, point_d=0):
    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize=(15,10)) # параметры окна
    ax = plt.axes()

    ax.yaxis.set_major_locator(LinearLocator(15)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(15))

    ha, hb, hb_, hc, hd, hd_ = point_a['h_a'], point_b['h_b'], point_b_['h_b_'],point_c['h_c'], point_d['h_d'], point_d_['h_d_']
    ta, tb, tb_, tc, td, td_  = point_a['t_a'], point_b['t_b'], point_b_['t_b_'], point_c['t_c'], point_d['t_d'], point_d_['t_d_']
    Pa, Pb = point_a['P_a'], point_b['P_b']

    plt.tick_params(axis='both', which='major', labelsize=15, 
                    direction = 'inout', length=10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')

    fig.suptitle('Цикл Брайтона для простой схемы ГТУ', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('S, кДж/(кгК)', fontsize=20, loc = 'center')
    plt.ylabel('h, кДж/кг',fontsize=20,loc = 'center')

    x_1 = np.array([point_b_['S_b_'], point_b['S_b'], point_c['S_c']])
    y_1 = np.array([point_b_['h_b_'], point_b['h_b'], point_c['h_c']])
    xnew_1 = np.linspace (x_1. min (), x_1. max (), 100 )
    spl_1 = make_interp_spline (x_1, y_1, k= 2 )
    y_smooth_1 = spl_1(xnew_1)
 
    x_2 = np.array([point_a['S_a'], point_d_['S_d_'], point_d['S_d']])
    y_2 = np.array([point_a['h_a'], point_d_['h_d_'], point_d['h_d']])
    xnew_2 = np.linspace(x_2. min (), x_2. max (), 300 )
    spl_2 = make_interp_spline (x_2, y_2, k= 2 )
    y_smooth_2 = spl_2(xnew_2)

    Hk, = plt.plot([point_a['S_a'], point_b_['S_b_']], [point_a['h_a'], point_b_['h_b_']], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')  
    Ht, = plt.plot([point_c['S_c'], point_d_['S_d_']], [point_c['h_c'], point_d_['h_d_']], color='r',marker='o',ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')  
    Hk_, = plt.plot([point_a['S_a'], point_b['S_b']], [point_a['h_a'], point_b['h_b']], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=2, linestyle='--')
    Ht_, = plt.plot([point_c['S_c'], point_d['S_d']], [point_c['h_c'], point_d['h_d']], color='r',marker='o',ms = 8, markerfacecolor='w',linewidth=2, linestyle='--')
    
    plt.plot (xnew_1, y_smooth_1, color='r',linewidth=3, linestyle='-')
    plt.plot (xnew_2, y_smooth_2, color='b',linewidth=3, linestyle='-')

    plt.annotate(f'$T_a = {round((ta+273.15),2)} К$', xy=(point_a['S_a'], point_a['h_a']), xytext=(point_a['S_a'] - 0.005, point_a['h_a']- 40), fontsize = 15)
    plt.annotate(f'$T_b^\prime = {round((tb+273.15),2)} К$', xy=(point_b['S_b'], point_b['h_b']), xytext=(point_b['S_b'] + 0.0001, point_b['h_b']- 40), fontsize = 15)
    plt.annotate(f'$T_b = {round((tb_+273.15),2)} К$', xy=(point_b_['S_b_'], point_b_['h_b_']), xytext=(point_b_['S_b_'] - 0.005, point_b_['h_b_']+ 25), fontsize = 15)
    
    plt.annotate(f'$T_c = {round((tc+273.15),2)} К$', xy=(point_c['S_c'], point_c['h_c']), xytext=(point_c['S_c'] - 0.025, point_c['h_c']+ 25), fontsize = 15)
    plt.annotate(f'$T_d^\prime = {round((td+273.15),2)} К$', xy=(point_d['S_d'], point_d['h_d']), xytext=(point_d['S_d'] - 0.025, point_d['h_d'] - 60), fontsize = 15)
    plt.annotate(f'$T_d = {round((td_+273.15),2)} К$', xy=(point_d_['S_d_'], point_d_['h_d_']), xytext=(point_d_['S_d_'] - 0.025, point_d_['h_d_'] - 70), fontsize = 15)
    
    ax.legend((Hk, Hk_, Ht, Ht_, Hk, Ht), [r'$H_к(теорет.) = ' f'{round((hb_ - ha), 2)} ' ' кДж/кг$', 
                                   r'$H_к(действит.) = ' f'{round((hb - ha), 2)} ' ' кДж/кг$',
                                   r'$H_т(теорет.) = ' f'{round((hc - hd_), 2)} ' ' кДж/кг$',
                                   r'$H_т(действит.) = ' f'{round((hc - hd), 2)} ' ' кДж/кг$',
                                   r'$P_a = P_d= ' f'{round(Pa/10e2, 2)} ' ' кПа$',
                                   r'$P_b = P_c= ' f'{round(Pb/10e2, 3)} ' ' кПа$'],
                                    fontsize = 15, frameon=True, framealpha=True)
    
    plt.show()