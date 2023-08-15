import pandas as pd

def stage_table(array , method = True):

    if method == 'parameters':
        data_1 = {'Изоэнтропийный теплоперепад в ступени':'кДж/кг',
                'Степень реактивности ступени':'-',
                'Изоэнтропийный теплоперепад сопловой решетки':'кДж/кг',
                'Действительный теплоперепад сопловой решетки':'кДж/кг',
                'Изоэнтропийный теплоперепад рабочей решетки':'кДж/кг',
                'Действительный теплоперепад рабочей решетки':'кДж/кг',
                'Температура заторможенного потока на входе в СА':'К',
                'Давление заторможенного потока на входе в СА':'Па',
                'Изоэнтропийная температура на выходе потока из СА':'К',
                'Статическая температура на выходе потока из СА':'К',
                'Статическое давление рабочего тела на выходе из СА':'Па',
                'Температура заторможенного потока на выходе из СА':'К',
                'Давление заторможенного потока на выходе из СА':'Па',
                'Температура заторможенного потока в СА в относительной СК':'К',
                'Давление заторможенного потока в СА в относительной СК':'Па',
                'Изоэнтропийная температура на выходе потока из РК':'К',
                'Статическая температура на выходе потока из РК':'К',
                'Статическое давление рабочего тела на выходе из РК':'Па',
                'Температура заторможенного потока на выходе из РК':'К',
                'Давление заторможенного потока на выходе из РК':'Па',
                'Температура заторможенного потока в РК в относительной СК':'К',
                'Давление заторможенного потока в РК в относительной СК':'Па',
                'Изоэнтропийная температура потока в РК в относительной СК':'К',
                'Абсолютная (изоэнтропическая) скорость истечения из СА':'м/с',
                'Абсолютная скорость истечения из СА':'м/с',
                'Окружная составляющая абсолютной скорости истечения из СА':'м/с',
                'Осевая составляющая абсолютной скорости истечения из СА':'м/с',
                'Относительная скорость потока на выходе из СА':'м/с',
                'Окружная составляющая относительной скорости истечения из СА':'м/с',
                'Осевая составляющая относительной скорости истечения из СА':'м/с',
                'Окружная скорость потока на выходе из СА':'м/с',
                'Абсолютный угол выхода потока из СА':'град', 
                'Относительный угол выхода потока из СА':'град',  
                'Абсолютная скорость потока на выходе из РК':'м/с',
                'Окружная составляющая абсолютной скорости истечения из РК':'м/с',
                'Осевая составляющая абсолютной скорости истечения из РК':'м/с',               
                'Относительная (изоэнтропическая) скорость истечения из РК':'м/с',
                'Относительная скорость истечения из РК':'м/с',
                'Окружная составляющая относительной скорости истечения из РК':'м/с',
                'Осевая составляющая относительной скорости истечения из РК':'м/с',
                'Окружная скорость потока на выходе из РК':'м/с',
                'Абсолютный угол выхода потока из РК':'град', 
                'Относительный угол выхода потока из РК':'град',
                'Коэффициент потери скорости в СА':'-',
                'Коэффициент потери скорости в РК':'-',
                'Коэффициент суммарных потерь энергии в СА':'-',
                'Коэффициент суммарных потерь энергии в РК':'-',
                'Коэффициент профильныех потерь энергии в СА':'-',
                'Коэффициент профильныех потерь энергии в РК':'-',
                'Коэффициент вторичных потерь энергии в СА':'-',
                'Коэффициент вторичных потерь энергии в РК':'-',
                'Коэффициент потерь энергии связанные с утечеками в СА':'-',
                'Коэффициент потерь энергии связанные с утечеками в РК':'-',
                'Коэффициент кромочных потерь энергии в СА':'-',
                'Коэффициент кромочных потерь энергии в РК':'-',
                'Коэффициент потерь энергии связанные с утечеками через переферию в СА':'-',
                'Коэффициент потерь энергии связанные с утечеками через переферию в РК':'-',
                'Число Маха потока на выходе из СА абсолютной скорости':'-', 
                'Число Маха потока на выходе из РК абсолютной скорости':'-', 
                'Число Маха потока на выходе из СА относительной скорости':'-', 
                'Число Маха потока на выходе из РК относительной скорости':'-', 
                'Расход газа на входе в СА':'кг/c',
                'Утечки газа в СА':'кг/c',
                'Расход охлаждающего воздуха в СА':'кг/c',
                'Расход газа на входе в РК':'кг/c',
                'Утечки газа в РК':'кг/c',
                'Расход охлаждающего воздуха в РК':'кг/c',
                'Расход газа на выходе из РК':'кг/c', 
                'Изоэнтропийная работа газа от параметров полного торможения':'кДж/кг',
                'Потери энергии от утечек в СА':'кДж/кг',
                'Потери энергии от утечек в РК':'кДж/кг', 
                'Потери энергии связанные с трением':'кДж/кг', 
                'Основные потери энергии в СА':'кДж/кг',
                'Основные потери энергии в РК':'кДж/кг',
                'Потери энергии с выходной скоростью потока':'кДж/кг',                
                'Работа газа от параметров полного торможения':'кДж/кг', 
                'Окружная работа ступени':'кДж/кг', 
                'Фиктивная скорость потока':'м/c',
                'Отношение скоростей U/Сф':'-',
                'КПД ступени':'-',
                'Мощность ступени':'Вт'
                }
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    if method == 'profile sopl':
        data_1 = {'Ширина профиля СА':'мм', 
                  'Хорда профиля СА': 'мм', 
                  'Шаг решетки СА':'мм',
                  'Радиус входной кромки СА':'мм', 
                  'Радиус выходной кромки СА': 'мм',
                  'Горло решетки на входе СА':'мм',
                  'Горло решетки на выходе СА':'мм', 
                  'Координата центра СА':'мм', 
                  'Толщина профиля СА':'мм',
                  'Корневой диаметр СА':'мм', 
                  'Высота лопатки СА':'мм', 
                  'Угол установки СА':'град',
                  'Лопаточный угол на входе СА':'град',
                  'Лопаточный угол на выходе СА':'град',
                  'Угол заострения на вх. СА':'град', 
                  'Угол заострения на вых. СА':'град',
                  'Угол отгиба СА':'град',
                  'Количество лопаток СА':'шт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    if method == 'profile rab':
        data_1 = {'Ширина профиля РК':'мм', 
                  'Хорда профиля РК': 'мм', 
                  'Шаг решетки РК':'мм',
                  'Радиус входной кромки РК':'мм', 
                  'Радиус выходной кромки РК': 'мм',
                  'Горло решетки на входе РК':'мм',
                  'Горло решетки на выходе РК':'мм', 
                  'Координата центра РК':'мм',
                  'Толщина профиля РК':'мм',
                  'Корневой диаметр РК':'мм', 
                  'Высота лопатки РК':'мм', 
                  'Угол установки РК':'град',
                  'Лопаточный угол на входе РК':'град',
                  'Лопаточный угол на выходе РК':'град',
                  'Угол заострения на вх. РК':'град', 
                  'Угол заострения на вых. РК':'град',
                  'Угол отгиба РК':'град',
                  'Количество лопаток РК':'шт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  
    
    return df_merged

def stage_cold_table(array , method = True):
    
    if method == 'parameters':
        data_1 = {'Изоэнтропийный теплоперепад в ступени':'кДж/кг',
        'Степень реактивности ступени':'-',
        'Изоэнтропийный теплоперепад сопловой решетки':'кДж/кг',
        'Дейстивтельный теплоперепад в охлождаемой сопловой решетки':'кДж/кг',
        'Изоэнтропийный теплоперепад рабочей решетки':'кДж/кг',
        'Дейстивтельный теплоперепад в охлождаемой рабочей решетки':'кДж/кг',
        'Температура заторможенного потока на входе в СА':'К',
        'Давление заторможенного потока на входе в СА':'Па',
        'Изоэнтропийная температура на выходе потока из СА':'К',
        'Статическая температура смешения на выходе потока из СА':'К',
        'Статическое давление смешения рабочего тела на выходе из СА':'Па',
        'Температура смешения заторможенного потока на выходе из СА':'К',
        'Давление смешения заторможенного потока на выходе из СА':'Па',
        'Температура смешения заторможенного потока в СА в относительной СК':'К',
        'Давление смешения заторможенного потока в СА в относительной СК':'Па',
        'Изоэнтропийная температура на выходе потока из РК':'К',
        'Статическая температура смешения на выходе потока из РК':'К',
        'Статическое давление смешения рабочего тела на выходе из РК':'Па',
        'Температура смешения заторможенного потока на выходе из РК':'К',
        'Давление смешения заторможенного потока на выходе из РК':'Па',
        'Температура смешения заторможенного потока в РК в относительной СК':'К',
        'Давление смешения заторможенного потока в РК в относительной СК':'Па',
        'Изоэнтропийная температура смешения потока в РК в относительной СК':'К',
        'Абсолютная (изоэнтропическая) скорость истечения из СА':'м/с',
        'Абсолютная скорость (смешения) истечения из СА':'м/с',
        'Окружная составляющая абсолютной скорости (смешения) истечения из СА':'м/с',
        'Осевая составляющая абсолютной скорости (смешения) истечения из СА':'м/с',
        'Относительная скорость (смешения) потока на выходе из СА':'м/с',
        'Окружная составляющая относительной скорости (смешения) истечения из СА':'м/с',
        'Осевая составляющая относительной скорости (смешения) истечения из СА':'м/с',
        'Окружная скорость потока на выходе из СА':'м/с', 
        'Абсолютный угол выхода потока (смешения) из СА':'град', 
        'Относительный угол выхода потока (смешения) из СА':'град',
        'Абсолютная скорость (смешения) потока на выходе из РК':'м/с',
        'Окружная составляющая абсолютной скорости (смешения) истечения из РК':'м/с',
        'Осевая составляющая абсолютной скорости (смешения) истечения из РК':'м/с',
        'Относительная (изоэнтропическая) скорость выхода из РК':'м/с',
        'Относительная скорость (смешения) выхода из РК':'м/с',
        'Окружная составляющая относительной скорости (смешения) выхода из РК':'м/с',
        'Осевая составляющая относительной скорости (смешения) выхода из РК':'м/с',
        'Окружная скорость потока на выходе из РК':'м/с', 
        'Абсолютный угол выхода потока (смешения) из РК':'град', 
        'Относительный угол выхода потока (смешения) из РК':'град',
        'Коэффициент потери скорости в СА':'-',
        'Коэффициент потери скорости в РК':'-',
        'Коэффициент суммарных потерь энергии в СА':'-',
        'Коэффициент суммарных потерь энергии в РК':'-',
        'Коэффициент профильныех потерь энергии в СА':'-',
        'Коэффициент профильныех потерь энергии в РК':'-',
        'Коэффициент вторичных потерь энергии в СА':'-',
        'Коэффициент вторичных потерь энергии в РК':'-',
        'Коэффициент потерь энергии связанные с утечеками в СА':'-',
        'Коэффициент потерь энергии связанные с утечеками в РК':'-',
        'Коэффициент кромочных потерь энергии в СА':'-',
        'Коэффициент кромочных потерь энергии в РК':'-',
        'Коэффициент потерь энергии связанные с утечеками через переферию в СА':'-',
        'Коэффициент потерь энергии связанные с утечеками через переферию в РК':'-',
        'Число Маха потока на выходе из СА абсолютной скорости':'-', 
        'Число Маха потока на выходе из РК абсолютной скорости':'-', 
        'Число Маха потока на выходе из СА относительной скорости':'-', 
        'Число Маха потока на выходе из РК относительной скорости':'-',
        'Расход газа на входе в СА':'кг/c',
        'Утечки газа в СА':'кг/c',
        'Расход охлаждающего воздуха в СА':'кг/c',
        'Расход газа на входе в РК':'кг/c',
        'Утечки газа в РК':'кг/c',
        'Расход охлаждающего воздуха в РК':'кг/c',
        'Расход газа на выходе из РК':'кг/c',
        'Изоэнтропийная работа газа от параметров полного торможения':'кДж/кг',
        'Потери энергии от утечек в СА':'кДж/кг',
        'Потери энергии от утечек в РК':'кДж/кг', 
        'Потери энергии связанные с трением':'кДж/кг', 
        'Основные потери энергии в СА':'кДж/кг',
        'Основные потери энергии в РК':'кДж/кг',
        'Потери энергии с выходной скоростью потока':'кДж/кг', 
        'Работа газа от параметров полного торможения':'кДж/кг', 
        'Окружная работа ступени':'кДж/кг', 
        'Фиктивная скорость потока':'м/c',
        'Отношение скоростей U/Сф':'-',
        'КПД ступени':'-',
        'Мощность ступени':'Вт'}
        
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    if method == 'profile sopl':
        data_1 = {'Ширина профиля СА':'мм', 
                  'Хорда профиля СА': 'мм', 
                  'Шаг решетки СА':'мм',
                  'Радиус входной кромки СА':'мм', 
                  'Радиус выходной кромки СА': 'мм',
                  'Горло решетки на входе СА':'мм',
                  'Горло решетки на выходе СА':'мм', 
                  'Координата центра СА':'мм', 
                  'Толщина профиля СА':'мм',
                  'Корневой диаметр СА':'мм', 
                  'Высота лопатки СА':'мм', 
                  'Угол установки СА':'град',
                  'Лопаточный угол на входе СА':'град',
                  'Лопаточный угол на выходе СА':'град',
                  'Угол заострения на вх. СА':'град', 
                  'Угол заострения на вых. СА':'град',
                  'Угол отгиба СА':'град',
                  'Количество лопаток СА':'шт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    if method == 'profile rab':
        data_1 = {'Ширина профиля РК':'мм', 
                  'Хорда профиля РК': 'мм', 
                  'Шаг решетки РК':'мм',
                  'Радиус входной кромки РК':'мм', 
                  'Радиус выходной кромки РК': 'мм',
                  'Горло решетки на входе РК':'мм',
                  'Горло решетки на выходе РК':'мм', 
                  'Координата центра РК':'мм',
                  'Толщина профиля РК':'мм',
                  'Корневой диаметр РК':'мм', 
                  'Высота лопатки РК':'мм', 
                  'Угол установки РК':'град',
                  'Лопаточный угол на входе РК':'град',
                  'Лопаточный угол на выходе РК':'град',
                  'Угол заострения на вх. РК':'град', 
                  'Угол заострения на вых. РК':'град',
                  'Угол отгиба РК':'град',
                  'Количество лопаток РК':'шт'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  
    
    if method == 'termod_sopl':
        data_1 = {'Изоэнтропийный теплоперепад сопловой решетки':'кДж/кг',
                'Дейстивтельный теплоперепад в охлождаемой сопловой решетки':'кДж/кг',
                'Температура заторможенного потока на входе в СА':'К',
                'Полное давление потока на входе в СА':'Па',
                'Изоэнтропийная температура на выходе потока из СА':'К',
                'Статическая температура смешения на выходе потока из СА':'К',
                'Статическое давление смешения рабочего тела на выходе из СА':'Па',
                'Температура смешения заторможенного потока на выходе из СА':'К',
                'Давление смешения заторможенного потока на выходе из СА':'Па',
                'Температура смешения заторможенного потока в СА в относительной СК':'К',
                'Давление смешения заторможенного потока в СА в относительной СК':'Па',
                'Теоретическая температура газа за охлаждаемой СА':'К',
                'Теоретическое давление газа за охлаждаемой СА':'Па', 
                'Температура заторможенного воздуха на входе в СА':'К',
                'Температура заторможенного воздуха на выходе из СА':'К',
                'Давление заторможенного воздуха на выходе из СА':'Па',
                'Коэффициент востановления полного давления в СА':'-',
                'Утчечка газа через лабиринтное уплотнение':'кг/с',
                'Потребный расход охлаждающего воздуха':'кг/с',
                'Расход рабочего тела за сопловой решеткой':'кг/с',
                'Количество тепла, отводимого от газа в систему охлажденися СА':'-',
                'Коэффициент удельного расхода охладителя в СА':'-'}
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  
    
    if method == 'termod_rab':
        data_1 = {'Изоэнтропийный теплоперепад рабочей решетки':'кДж/кг',
                'Дейстивтельный теплоперепад в охлождаемой рабочей решетки':'кДж/кг',
                'Изоэнтропийная температура на выходе потока из РК':'К',
                'Статическая температура смешения на выходе потока из РК':'К',
                'Статическое давление смешения рабочего тела на выходе из РК':'Па',
                'Температура смешения заторможенного потока на выходе из РК':'К',
                'Давление смешения заторможенного потока на выходе из РК':'Па',
                'Температура смешения заторможенного потока в РК в относительной СК':'К',
                'Давление смешения заторможенного потока в РК в относительной СК':'Па',
                'Изоэнтропийная температура смешения потока в РК в относительной СК':'К',
                
                'Теоретическая температура газа за охлаждаемой РК':'К',
                'Теоретическое давление газа за охлаждаемой РК':'Па',
                'Температура заторможенного воздуха на входе в РК':'К',
                'Температура заторможенного воздуха на выходе из РК':'К',
                'Давление заторможенного воздуха на выходе из РК':'Па',
                
                'Коэффициент востановления полного давления в РК':'-',
                'Утчечка газа через лабиринтное уплотнение':'кг/с',
                'Потребный расход охлаждающего воздуха':'кг/с',
                'Расход рабочего тела за рабочей решеткой':'кг/с',
                'Количество тепла, отводимого от газа в систему охлажденися РК':'-',
                'Коэффициент удельного расхода охладителя в РК':'-'}  
        df_merged = pd.DataFrame(list(data_1.items()),columns=['Параметры','Размерность']).join(pd.DataFrame((array.items()),columns=['Обозначение', 'Значение']), rsuffix='_right')  

    return df_merged

def stageTable(array):
    data_1 = {'Изоэнтропийный теплоперепад в ступени':'Hs_st_',
                'Степень реактивности ступени':'ro_st',
                'Изоэнтропийный теплоперепад сопловой решетки':'Hs_sa',
                'Действительный теплоперепад сопловой решетки':'H_sa',
                'Изоэнтропийный теплоперепад рабочей решетки':'Hs_rk',
                'Действительный теплоперепад рабочей решетки':'H_rk',
                'Температура заторможенного потока на входе в СА':'T0_',
                'Давление заторможенного потока на входе в СА':'P0_',
                'Изоэнтропийная температура на выходе потока из СА':'T1s',
                'Статическая температура на выходе потока из СА':'T1',
                'Статическое давление рабочего тела на выходе из СА':'P1',
                'Температура заторможенного потока на выходе из СА':'T1_',
                'Давление заторможенного потока на выходе из СА':'P1_',
                'Температура заторможенного потока в СА в относительной СК':'T1w_',
                'Давление заторможенного потока в СА в относительной СК':'P1w_',
                'Изоэнтропийная температура на выходе потока из РК':'T2s',
                'Статическая температура на выходе потока из РК':'T2',
                'Статическое давление рабочего тела на выходе из РК':'P2',
                'Температура заторможенного потока на выходе из РК':'T2_',
                'Давление заторможенного потока на выходе из РК':'P2_',
                'Температура заторможенного потока в РК в относительной СК':'T2w_',
                'Давление заторможенного потока в РК в относительной СК':'P2w_',
                'Изоэнтропийная температура потока в РК в относительной СК':'T2s_',
                'Абсолютная (изоэнтропическая) скорость истечения из СА':'C1s',
                'Абсолютная скорость истечения из СА':'C1',
                'Окружная составляющая абсолютной скорости истечения из СА':'C1u',
                'Осевая составляющая абсолютной скорости истечения из СА':'C1a',
                'Относительная скорость потока на выходе из СА':'W1',
                'Окружная составляющая относительной скорости истечения из СА':'W1u',
                'Осевая составляющая относительной скорости истечения из СА':'W1a',
                'Окружная скорость потока на выходе из СА':'U1',
                'Абсолютный угол выхода потока из СА':'alpha1', 
                'Относительный угол выхода потока из СА':'betta1', 
                'Абсолютная скорость потока на выходе из РК':'C2',
                'Окружная составляющая абсолютной скорости истечения из РК':'C2u',
                'Осевая составляющая абсолютной скорости истечения из РК':'C2a',  
                'Относительная (изоэнтропическая) скорость истечения из РК':'W2s',
                'Относительная скорость истечения из РК':'W2',
                'Окружная составляющая относительной скорости истечения из РК':'W2u',
                'Осевая составляющая относительной скорости истечения из РК':'W2a',
                'Окружная скорость потока на выходе из РК':'U2',
                'Абсолютный угол выхода потока из РК':'alpha2', 
                'Относительный угол выхода потока из РК':'betta2',
                'Коэффициент потери скорости в СА':'fi',
                'Коэффициент потери скорости в РК':'psi',
                'Коэффициент суммарных потерь энергии в СА':'Y_s_sopl',
                'Коэффициент суммарных потерь энергии в РК':'Y_s_rab',
                'Коэффициент профильныех потерь энергии в СА':'Y_p_sopl',
                'Коэффициент профильныех потерь энергии в РК':'Y_p_rab',
                'Коэффициент вторичных потерь энергии в СА':'Y_sec_sopl',
                'Коэффициент вторичных потерь энергии в РК':'Y_sec_rab',
                'Коэффициент потерь энергии связанные с утечеками в СА':'Y_tl_sopl',
                'Коэффициент потерь энергии связанные с утечеками в РК':'Y_tl_rab',
                'Коэффициент кромочных потерь энергии в СА':'Y_te_sopl',
                'Коэффициент кромочных потерь энергии в РК':'Y_te_rab',
                'Коэффициент потерь энергии связанные с утечеками через переферию в СА':'Y_cl_sopl',
                'Коэффициент потерь энергии связанные с утечеками через переферию в РК':'Y_cl_rab',
                'Число Маха потока на выходе из СА абсолютной скорости':'M1c', 
                'Число Маха потока на выходе из РК абсолютной скорости':'M2c', 
                'Число Маха потока на выходе из СА относительной скорости':'M1w', 
                'Число Маха потока на выходе из РК относительной скорости':'M2w', 
                'Расход газа на входе в СА':'G0',
                'Утечки газа в СА':'G_g_sopl_gas',
                'Расход охлаждающего воздуха в СА':'G_cold_sopl',
                'Расход газа на входе в РК':'G1',
                'Утечки газа в РК':'G_g_rab_gas',
                'Расход охлаждающего воздуха в РК':'G_cold_rab',
                'Расход газа на выходе из РК':'G2', 
                'Изоэнтропийная работа газа от параметров полного торможения':'Ls_st_',
                'Потери энергии от утечек в СА':'dL_lake_sopl',
                'Потери энергии от утечек в РК':'dL_lake_rab',
                'Потери энергии связанные с трением':'dL_tr', 
                'Основные потери энергии в СА':'dL_sopl',
                'Основные потери энергии в РК':'dL_rab',
                'Потери энергии с выходной скоростью потока':'dL_vs',                
                'Работа газа от параметров полного торможения':'L_st_', 
                'Окружная работа ступени':'Lu', 
                'Фиктивная скорость потока':'C_fict',
                'Отношение скоростей U/Сф':'U_Cfict',
                'КПД ступени':'etta_st',
                'Мощность ступени':'N_i'}


    df_1 = pd.DataFrame(list(data_1.items()),columns=['Параметры','Обозначение']) 
    df_2_ = pd.DataFrame(['кДж/кг', '-', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'К', 'Па', 'К', 'К',
                      'Па', 'К', 'Па', 'К', 'Па', 'К', 'К', 'Па', 'К', 'Па', 'К', 'Па', 'К',
                      'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'град', 'град', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'град', 'град',
                      '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-', '-',
                      'кг/с', 'к/с', 'к/с', 'к/с,', 'кг/с', 'к/с', 'к/с', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 'кДж/кг', 
                      'м/с', '-', '-', 'Вт'], columns = ['Разм'])
    df_2 = pd.DataFrame((array))
    df12_merged = df_1.join(df_2_, rsuffix='_right')
    df_merged = df12_merged.join(df_2, rsuffix='_right')

    return df_merged

def stage_section_table(array):
    data_1 = {'Радиус сечения СА': 'radius_sopl_i', 
            'Радиус сечения РК': 'radius_rab_i',
            'Безразмерный радиус сечения СА': 'r_sopl_i_', 
            'Безразмерный радиус сечения РК': 'r_rab_i_',
            'Относительный радиус сечения СА': 'rela_sopl_i', 
            'Относительный радиус сечения РК': 'rela_rab_i_',
            'Температура заторможенного потока на входе в СА':'T0_i_',
            'Давление заторможенного потока на входе в СА':'P0_i_',
            'Статическая температура на выходе потока из СА':'T1_i',
            'Статическое давление рабочего тела на выходе из СА':'P1_i',
            'Статическая температура на выходе потока из РК':'T2_i',
            'Статическое давление рабочего тела на выходе из РК':'P2_i',
            'Изоэнтропийный теплоперепад в сечении СА':'Hs_sa_i',
            'Изоэнтропийный теплоперепад в сечении РК':'Hs_rab_i',
            'Термическая степень реактивности в сечении': 'ro_term_i',
            'Кинематическая степень реактивности в сечении': 'ro_k_i',

            'Абсолютная скорость истечения из СА':'C_1_i',
            'Окружная составляющая абсолютной скорости истечения из СА':'C_1u_i',
            'Осевая составляющая абсолютной скорости истечения из СА':'C_1a_i',
            'Относительная скорость потока на выходе из СА':'W_1_i',
            'Окружная составляющая относительной скорости истечения из СА':'W_1u_i',
            'Осевая составляющая относительной скорости истечения из СА':'W_1a_i',
            'Окружная скорость потока на выходе из СА':'U_1_i',
            'Абсолютный угол выхода потока из СА':'alpha_1_i', 
            'Относительный угол выхода потока из СА':'betta_1_i',

            'Абсолютная скорость потока на выходе из РК':'C_2_i',
            'Окружная составляющая абсолютной скорости истечения из РК':'C_2u_i',
            'Осевая составляющая абсолютной скорости истечения из РК':'C_2a_i', 
            'Относительная скорость истечения из РК':'W_2_i',
            'Окружная составляющая относительной скорости истечения из РК':'W_2u_i',
            'Осевая составляющая относительной скорости истечения из РК':'W_2a_i',
            'Окружная скорость потока на выходе из РК':'U_2_i',
            'Абсолютный угол выхода потока из РК':'alpha_2_i', 
            'Относительный угол выхода потока из РК':'betta_2_i',

            'Число Маха потока на выходе из СА абсолютной скорости':'M1c_i', 
            'Число Маха потока на выходе из РК абсолютной скорости':'M2c_i', 
            'Число Маха потока на выходе из СА относительной скорости':'M1w_i', 
            'Число Маха потока на выходе из РК относительной скорости':'M2w_i'}

    df_1 = pd.DataFrame(list(data_1.items()),columns=['Параметры','Обозначение']) 
    df_2_ = pd.DataFrame(['м', 'м', '-', '-', '-', '-', 'К', 'Па', 'К', 'Па', 'К', 'Па', 'кДж/кг', 'кДж/кг', '-', '-', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'град', 'град', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'м/с', 'град', 'град', '-', '-', '-', '-'], columns = ['Разм'])
    df_2 = pd.DataFrame((array))
    df12_merged = df_1.join(df_2_, rsuffix ='_right')
    df_merged = df12_merged.join(df_2, rsuffix ='_right')
    return df_merged
