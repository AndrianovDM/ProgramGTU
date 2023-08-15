import matplotlib.pyplot as plt
import matplotlib as mpl

def fuel_plot(x, y, method):
    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize=(15,5)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis='both', which='major', labelsize=15, 
                    direction = 'inout', length=10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    
    if method == 'enthalpy':
        fig.suptitle('Зависимость энтальпии от температуры', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('t, град', fontsize=20, loc = 'center')
        plt.ylabel('h, кДж/кг',fontsize=20,loc = 'center')
        h = plt.plot(x, y, color='b', marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        ax.legend(h, [r'$h$ - $Энтальпия$ $топлива$'],
                                        fontsize = 15, frameon=True, framealpha=True)  

    if method == 'entropy':
        fig.suptitle('Зависимость энтропии от температуры', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('t, град', fontsize=20, loc = 'center')
        plt.ylabel('S, кДж/(кг*К)',fontsize=20,loc = 'center')
        s = plt.plot(x, y, color='r', marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        ax.legend(s, [r'$S$ - $Энтропия$ $топлива$'], fontsize = 15, frameon=True, framealpha=True)      
    
    if method == 'capacity':
        fig.suptitle('Зависимость изобарной теплоемкости от температуры', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('t, град', fontsize=20, loc = 'center')
        plt.ylabel('Cp, -',fontsize=20,loc = 'center')
        Cp = plt.plot(x, y, color='g', marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        ax.legend(Cp, [r'$Cp$ - $Изобарная$ $теплоемкость$ $топлива$'], fontsize = 15, frameon=True, framealpha=True)    

    if method == 'isentrope':
        fig.suptitle('Зависимость показателя изоэнтропы от температуры', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('t, град', fontsize=20, loc = 'center')
        plt.ylabel('k, -',fontsize=20,loc = 'center')
        Cp = plt.plot(x, y, color='#800080', marker='o', ms = 8, markerfacecolor='w',linewidth=3, linestyle='-') 
        ax.legend(Cp, [r'$Cp$ - $Показателя$ $изоэнтропы$ $топлива$'], fontsize = 15, frameon=True, framealpha=True)    
    
    plt.show()