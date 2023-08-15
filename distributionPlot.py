import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib import cm
from matplotlib.ticker import MaxNLocator

def ro_plot_2D(alpha, Y, method):
    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize=(15,10)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis='both', which='major', labelsize=15, 
                    direction = 'inout', length=10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    
    # if method == '05':
    #     fig.suptitle('$Зависимость$ $параметра$ $нагружености$ $Y^*_{ст}$ $от$' r' $\alpha_2,$' ' $при:$' '$0 <$' '$C_{2a}$' '$\leq 0.65$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    # if method == '08':
    #     fig.suptitle('$Зависимость$ $параметра$ $нагружености$ $Y^*_{ст}$ $от$' r' $\alpha_2,$' ' $при:$' '$0.65 <$' '$C_{2a}$' '$\leq 0.95$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    # if method == '11':
    #     fig.suptitle('$Зависимость$ $параметра$ $нагружености$ $Y^*_{ст}$ $от$' r' $\alpha_2,$' ' $при:$' '$0.95 <$' '$C_{2a}$' '$\leq 1.1$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    
    plt.xlabel(r'$\alpha_2, {\circ}$', fontsize=20, loc = 'center')
    plt.ylabel(r'$Y^*_{ст}, -$', fontsize=20,loc = 'center')
    ro_01, = plt.plot(alpha[0], Y[0], color='navy',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
    ro_02, = plt.plot(alpha[1], Y[1], color='blue',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
    ro_03, = plt.plot(alpha[2], Y[2], color='aqua',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
    ro_04, = plt.plot(alpha[3], Y[3], color='lawngreen',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
    ro_05, = plt.plot(alpha[4], Y[4], color='yellow',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
    ro_06, = plt.plot(alpha[5], Y[5], color='orange',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
    ro_07, = plt.plot(alpha[6], Y[6], color='maroon',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 

    ax.legend((ro_01, ro_02, ro_03, ro_04, ro_05, ro_06, ro_07), 
                  [r'$\rho$ = 0.1',r'$\rho$ = 0.2',r'$\rho$ = 0.3',r'$\rho$ = 0.4',
                   r'$\rho$ = 0.5',r'$\rho$ = 0.6',r'$\rho$ = 0.7'], fontsize = 15, frameon = True, framealpha = True)
    plt.show()

def ro_plot_3D(array, method):
    fig = plt.figure(figsize=(10, 7))
    ax_3d = fig.add_subplot(projection ='3d')

    if method == '05':
        fig.suptitle('Зависимость параметра нагружености $Y^*_{ст}$ от' r' $\alpha_2,$' ' при: ' '$0 <$' '$C_{2a}$' '$\leq 0.65$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    if method == '08':
        fig.suptitle('Зависимость параметра нагружености $Y^*_{ст}$ от' r' $\alpha_2,$' ' при: ' '$0.65 <$' '$C_{2a}$' '$\leq 0.95$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    if method == '11':
        fig.suptitle('Зависимость параметра нагружености $Y^*_{ст}$ от' r' $\alpha_2,$' ' при: ' '$0.95 <$' '$C_{2a}$' '$\leq 1.1$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    
    surf = ax_3d.scatter(array[:,0], array[:,1], array[:,2], c = array[:,2], cmap = cm.jet, linewidth=0.5)
    fig.colorbar(surf)
    ax_3d.xaxis.set_major_locator(MaxNLocator(5))
    ax_3d.yaxis.set_major_locator(MaxNLocator(6))
    ax_3d.zaxis.set_major_locator(MaxNLocator(6))
    fig.tight_layout()
    ax_3d.set_xlabel(r'$Угол$ $выхода$ $\alpha_2, {\circ}$')
    ax_3d.set_ylabel(r'$Нагруженость ступени$ $Y^*_{ст}, -$')
    ax_3d.set_zlabel(r'$Степень$ $реактивности$ $\rho_{ст}$')
    plt.show() 

def etta_plot_2D(C_2a, mu, val_1, val_2, val_3):

    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize=(15,10)) # параметры окна
    ax = plt.axes()
    plt.xlim((0.35, 1.3))
    plt.ylim((0.5, 3.5))
    plt.tick_params(axis='both', which='major', labelsize=15, 
                    direction = 'inout', length=10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    # fig.suptitle(r'$КПД$ $по$ $заторможенным$ $параметрам$ $ступеней$ $газовых$ $турбин$', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel(r'$Коэффициент$ $расхода$ $\overline {C_a}, -$', fontsize=20, loc = 'center')
    plt.ylabel(r'$Коэффициент$ $нагрузки$ $\mu_{ст}, -$', fontsize=20,loc = 'center')

    etta_85 = plt.scatter(C_2a[0], mu[0], s = 40, c = 'indigo') 
    etta_86 = plt.scatter(C_2a[1], mu[1], s = 40, c = 'blue') 
    etta_87 = plt.scatter(C_2a[2], mu[2], s = 40, c = 'cornflowerblue') 
    etta_88 = plt.scatter(C_2a[3], mu[3], s = 40, c = 'aqua') 
    etta_89 = plt.scatter(C_2a[4], mu[4], s = 40, c = 'lime')
    etta_90 = plt.scatter(C_2a[5], mu[5], s = 40, c = 'lawngreen')
    etta_91 = plt.scatter(C_2a[6], mu[6], s = 40, c = 'yellow')
    etta_92 = plt.scatter(C_2a[7], mu[7], s = 40, c = 'orange')
    etta_93 = plt.scatter(C_2a[8], mu[8], s = 40, c = 'r')  
    etta_94 = plt.scatter(C_2a[9], mu[9], s = 40, c = 'maroon') 
    for i in range(len(val_1)):
        etta = plt.scatter(val_1[i], val_2[i], marker ="s", s = 60, c = 'r', edgecolor ="blue")
        plt.plot([val_1[i], val_1[i]], [0, val_2[i]], color = 'black', linewidth = 2, linestyle='--')
        plt.plot([val_1[i], 0], [val_2[i], val_2[i]], color = 'black', linewidth = 2, linestyle='--')
    
    for X, Y, Z in zip(val_1, val_2, np.round(val_3,4)):
        ax.annotate('{}'.format(Z), xy=(X, Y), xytext=(-5, 5), ha='left',
                textcoords = 'offset points', size = 15)

    ax.legend((etta_85, etta_86, etta_87, etta_88, etta_89, etta_90, etta_91, etta_92, etta_93, etta_94), 
                  [r'$\eta$ = 0.85', r'$\eta$ = 0.86', r'$\eta$ = 0.87', r'$\eta$ = 0.88',
                   r'$\eta$ = 0.89', r'$\eta$ = 0.90', r'$\eta$ = 0.91', r'$\eta$ = 0.92', r'$\eta$ = 0.93', r'$\eta$ = 0.94'], fontsize = 15, frameon = True, framealpha = True)
    plt.show()

def etta_plot_3D(array):
    fig = plt.figure(figsize=(10, 7))
    fig.suptitle(r'КПД по заторможенным параметрам ступеней газовых турбин', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    ax_3d = fig.add_subplot(projection='3d')
    surf = ax_3d.scatter(array[:,0], array[:,1], array[:,2], c = array[:,2], cmap = cm.jet, linewidth = 0.5)
    # surf = ax_3d.plot_trisurf(x, y, z, cmap=cm.jet, linewidth=0)
    fig.colorbar(surf)
    ax_3d.xaxis.set_major_locator(MaxNLocator(5))
    ax_3d.yaxis.set_major_locator(MaxNLocator(6))
    ax_3d.zaxis.set_major_locator(MaxNLocator(6))
    fig.tight_layout()
    ax_3d.set_xlabel(r'$Коэффициент$ $расхода$ $\overline {C_a}, -$')
    ax_3d.set_ylabel(r'$Коэффициент$ $нагрузки$ $\mu_{ст}, -$')
    ax_3d.set_zlabel(r'$\eta_{ст}$')
    plt.show()

def distribution_plot(arr, method):
    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize=(15,5)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis='both', which='major', labelsize=15, 
                    direction = 'inout', length=10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    
    if method == 'temperature':
        fig.suptitle('Распределение температур i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('T, К',fontsize=20,loc = 'center')
        T0_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['T0_i_'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        T2_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['T2_i_'], color='r',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        T2_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['T2_i'], color='black',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        T2s_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['T2s_i'], color='g',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        T2s_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['T2s_i_'], color='#800080',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        ax.legend((T0_i_, T2_i_, T2_i, T2s_i, T2s_i_), [r'$\overline {T_{0i}}$ - $Температура$ $на$ $входе$ $СА$',
                                                        r'$\overline {T_{2i}}$ - $Температура$ $на$ $выходе$ $из$ $РК$', 
                                                        r'$T_{2i}$ - $Статическая$ $температура$ $на$ $выходе$ $из$ $РК$', 
                                                        r'$T_{2si}$ - $Статическая$ $температура$ $на$ $выходе$ $из$ $РК$', 
                                                        r'$\overline {T_{2si}}$ - $Изонтропийная$ $температура$ $на$ $выходе$ $из$ $СА$'],
                                        fontsize = 10, frameon=True, framealpha=True)      
    
    if method == 'pressure':
        fig.suptitle('Распределение давлений i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('P, Па',fontsize=20,loc = 'center')
        P0_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['P0_i_'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        P2_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['P2_i_'], color='r',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        P2_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['P2_i'], color='black',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 

        ax.legend((P0_i_, P2_i_, P2_i,), [r'$\overline {P_{0i}}$ - $Давление$ $на$ $входе$ $в$ $СА$',
                                          r'$\overline {P_{2i}}$ - $Давление$ $на$ $выходе$ $из$ $РК$', 
                                          r'$P_{2i}$ - $Давление$ $статическое$ $на$ $выходе$ $из$ $РК$'],
                                        fontsize = 10, frameon=True, framealpha=True)   
    
    if method == 'alpha':
        fig.suptitle(r'Распределение углов $\alpha$ i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('apha, град',fontsize=20,loc = 'center')
        alpha_0_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['alpha_0_i'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        alpha_02_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['alpha_02_i'], color='r',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        ax.legend((alpha_0_i, alpha_02_i,), [r'$\alpha_{0i}$ - $Угол$ $входа$ $в$ $ступень$',
                                             r'$\alpha_{2i}$ - $Угол$ $выхода$ $из$ $ступени$'],
                                        fontsize = 10, frameon=True, framealpha=True)
    
    if method == 'velocity':    
        fig.suptitle(r'Распределение скорости $C_{2a}$ i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('C2a, м/c',fontsize=20,loc = 'center')
        C2a_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['C2a_i'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        ax.legend([C2a_i], [r'$C_{2ai}$ - $Осевая$ $составляющая$ $скорости$ $на$ $выходе$ $из$ $РК$'], fontsize = 10, frameon=True, framealpha=True)   
    
    if method == 'heatdrop':
        fig.suptitle(r'Распределение теплоперепадов i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('Lst, Hst, кДж/кг',fontsize=20,loc = 'center')
        Ls_st_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['Ls_st_i_'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        L_st_i_, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['L_st_i_'], color='g',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        Hs_st_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['Hs_st_i'], color='#800080',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 

        ax.legend([Ls_st_i_, L_st_i_, Hs_st_i], [r'$\overline {L}_{sti}^\prime$ - $Изоэнтропическая$ $полезная$ $работа$ $расширения$ $газа$',
                                          r'$\overline {L_{sti}}$ - $Полезная$ $работа$ $расширения$ $газа$',
                                          r'$H_{s_{sti}}$ - $Изоэнтропийный$ $теплоперепад$ $ступени$'],
                                        fontsize = 10, frameon=False, framealpha=False)
    
    if method == 'ro':    
        fig.suptitle(r'Распределение степени реактивности $\rho$ i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('rho, -', fontsize=20, loc = 'center')
        ro_st_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['ro_st_i'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        ax.legend([ro_st_i], [r'$\rho$ - $Степень$ $реактивности$'], fontsize = 10, frameon=True, framealpha=True)        
    
    if method == 'etta':    
        fig.suptitle(r'Распределение КПД  $\eta$ i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('etta, -', fontsize=20, loc = 'center')
        etta_st_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[0]['etta_st_i'], color='r',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        ax.legend([etta_st_i], [r'$\eta$ - $КПД$ $ступени$'], fontsize = 10, frameon=True, framealpha=True)       
    
    if method == 'gdinamyc':    
        fig.suptitle(r'Распределение газодинамических функций i-ой ступени', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('n - номер ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('lamda,q, pi,tau, -', fontsize=20, loc = 'center')
        lamda_c2a_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['lamda_c2a_i'], color='r',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')        
        q_lamda_c2a_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['q_lamda_c2a_i'], color='b',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        PI_lamda_c2a_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['PI_lamda_c2a_i'], color='g',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        tau_lamda_c2a_i, = plt.plot(np.linspace(1, arr[0]['number_of_steps'], arr[0]['number_of_steps']), arr[1]['tau_lamda_c2a_i'], color='black',marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-')
        ax.legend([lamda_c2a_i ,q_lamda_c2a_i, PI_lamda_c2a_i, tau_lamda_c2a_i], [r'$\lambda$',r'$q(\lambda)$', r'$\pi(\lambda)$',r'$\tau(\lambda)$'], fontsize = 10, frameon=True, framealpha=True)       
    plt.show()

# ro_plot_2D([alpha_ro_01_05, alpha_ro_02_05, alpha_ro_03_05, alpha_ro_04_05, 
#               alpha_ro_05_05, alpha_ro_06_05, alpha_ro_07_05],
#              [Y_ro_01_05, Y_ro_02_05, Y_ro_03_05, Y_ro_04_05, 
#               Y_ro_05_05, Y_ro_06_05, Y_ro_07_05], method = '05')

# ro_plot_3D(array = data_ro_05, method = '05')

# ro_plot_2D([alpha_ro_01_08, alpha_ro_02_08, alpha_ro_03_08, alpha_ro_04_08, 
#               alpha_ro_05_08, alpha_ro_06_08, alpha_ro_07_08],
#              [Y_ro_01_08, Y_ro_02_08, Y_ro_03_08, Y_ro_04_08, 
#               Y_ro_05_08, Y_ro_06_08, Y_ro_07_08], method = '08')
# ro_plot_3D(array = data_ro_05, method = '05')

# ro_plot_2D([alpha_ro_01_11, alpha_ro_02_11, alpha_ro_03_11, alpha_ro_04_11, 
#               alpha_ro_05_11, alpha_ro_06_11, alpha_ro_07_11],
#              [Y_ro_01_11, Y_ro_02_11, Y_ro_03_11, Y_ro_04_11, 
#               Y_ro_05_11, Y_ro_06_11, Y_ro_07_11], method = '11')
# ro_plot_3D(array = data_ro_05, method = '11')

# etta_plot_2D([c_85, c_86, c_87, c_88, c_89, c_90, c_91, etta_92, etta_93, etta_94], 
#              [mu_85, mu_86, mu_87, mu_88, mu_89, mu_90, mu_91, mu_92, mu_93, mu_94])
# etta_plot_3D(array = data_etta)