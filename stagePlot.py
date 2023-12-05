import numpy as np
from math import *
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from scipy.interpolate import make_interp_spline
from matplotlib.patches import Arc
from matplotlib.transforms import IdentityTransform, TransformedBbox, Bbox

class AngleAnnotation(Arc):
    """
    Draws an arc between two vectors which appears circular in display space.
    """
    def __init__(self, xy, p1, p2, size=75, unit="points", ax=None,
                 text="", textposition="inside", text_kw=None, **kwargs):
        """
        Parameters
        ----------
        xy, p1, p2 : tuple or array of two floats
            Center position and two points. Angle annotation is drawn between
            the two vectors connecting *p1* and *p2* with *xy*, respectively.
            Units are data coordinates.

        size : float
            Diameter of the angle annotation in units specified by *unit*.

        unit : str
            One of the following strings to specify the unit of *size*:

            * "pixels": pixels
            * "points": points, use points instead of pixels to not have a
              dependence on the DPI
            * "axes width", "axes height": relative units of Axes width, height
            * "axes min", "axes max": minimum or maximum of relative Axes
              width, height

        ax : `matplotlib.axes.Axes`
            The Axes to add the angle annotation to.

        text : str
            The text to mark the angle with.

        textposition : {"inside", "outside", "edge"}
            Whether to show the text in- or outside the arc. "edge" can be used
            for custom positions anchored at the arc's edge.

        text_kw : dict
            Dictionary of arguments passed to the Annotation.

        **kwargs
            Further parameters are passed to `matplotlib.patches.Arc`. Use this
            to specify, color, linewidth etc. of the arc.

        """
        self.ax = ax or plt.gca()
        self._xydata = xy  # in data coordinates
        self.vec1 = p1
        self.vec2 = p2
        self.size = size
        self.unit = unit
        self.textposition = textposition

        super().__init__(self._xydata, size, size, angle=0.0,
                         theta1=self.theta1, theta2=self.theta2, **kwargs)

        self.set_transform(IdentityTransform())
        self.ax.add_patch(self)

        self.kw = dict(ha="center", va="center",
                       xycoords=IdentityTransform(),
                       xytext=(0, 0), textcoords="offset points",
                       annotation_clip=True)
        self.kw.update(text_kw or {})
        self.text = ax.annotate(text, xy=self._center, **self.kw)

    def get_size(self):
        factor = 1.
        if self.unit == "points":
            factor = self.ax.figure.dpi / 72.
        elif self.unit[:4] == "axes":
            b = TransformedBbox(Bbox.unit(), self.ax.transAxes)
            dic = {"max": max(b.width, b.height),
                   "min": min(b.width, b.height),
                   "width": b.width, "height": b.height}
            factor = dic[self.unit[5:]]
        return self.size * factor

    def set_size(self, size):
        self.size = size

    def get_center_in_pixels(self):
        """return center in pixels"""
        return self.ax.transData.transform(self._xydata)

    def set_center(self, xy):
        """set center in data coordinates"""
        self._xydata = xy

    def get_theta(self, vec):
        vec_in_pixels = self.ax.transData.transform(vec) - self._center
        return np.rad2deg(np.arctan2(vec_in_pixels[1], vec_in_pixels[0]))

    def get_theta1(self):
        return self.get_theta(self.vec1)

    def get_theta2(self):
        return self.get_theta(self.vec2)

    def set_theta(self, angle):
        pass

    # Redefine attributes of the Arc to always give values in pixel space
    _center = property(get_center_in_pixels, set_center)
    theta1 = property(get_theta1, set_theta)
    theta2 = property(get_theta2, set_theta)
    width = property(get_size, set_size)
    height = property(get_size, set_size)

    # The following two methods are needed to update the text position.
    def draw(self, renderer):
        self.update_text()
        super().draw(renderer)

    def update_text(self):
        c = self._center
        s = self.get_size()
        angle_span = (self.theta2 - self.theta1) % 360
        angle = np.deg2rad(self.theta1 + angle_span / 2)
        r = s / 2
        if self.textposition == "inside":
            r = s / np.interp(angle_span, [60, 90, 135, 180],
                                          [3.3, 3.5, 3.8, 4])
        self.text.xy = c + r * np.array([np.cos(angle), np.sin(angle)])
        if self.textposition == "outside":
            def R90(a, r, w, h):
                if a < np.arctan(h/2/(r+w/2)):
                    return np.sqrt((r+w/2)**2 + (np.tan(a)*(r+w/2))**2)
                else:
                    c = np.sqrt((w/2)**2+(h/2)**2)
                    T = np.arcsin(c * np.cos(np.pi/2 - a + np.arcsin(h/2/c))/r)
                    xy = r * np.array([np.cos(a + T), np.sin(a + T)])
                    xy += np.array([w/2, h/2])
                    return np.sqrt(np.sum(xy**2))

            def R(a, r, w, h):
                aa = (a % (np.pi/4))*((a % (np.pi/2)) <= np.pi/4) + \
                     (np.pi/4 - (a % (np.pi/4)))*((a % (np.pi/2)) >= np.pi/4)
                return R90(aa, r, *[w, h][::int(np.sign(np.cos(2*a)))])

            bbox = self.text.get_window_extent()
            X = R(angle, r, bbox.width, bbox.height)
            trans = self.ax.figure.dpi_scale_trans.inverted()
            offs = trans.transform(((X-s/2), 0))[0] * 72
            self.text.set_position([offs*np.cos(angle), offs*np.sin(angle)])

def interp_plot(x, y, color):
    xnew_1 = np.linspace (x. min (), x. max (), 100 )
    spl_1 = make_interp_spline (x, y, k = 2)
    y_smooth_1 = spl_1(xnew_1)
    return plt.plot (xnew_1, y_smooth_1, color, linewidth = 2, linestyle='--')
    
def hs_plot(point0_, point1s, point1, point1_, point1w_, point2s, point2, point2_, point2s_, point2w_, i, method):
    
    if method == 'notcold':
        h0_, h1s, h1, h1_, h1w_, h2s, h2, h2_, h2s_, h2w_ =  point0_['I0_'], point1s['I1s'], point1['I1'], point1_['I1_'], point1w_['I1w_'], point2s['I2s'], point2['I2'], point2_['I2_'], point2s_['I2s_'], point2w_['I2w_']
        s0_, s1s, s1, s1_, s1w_, s2s, s2, s2_, s2s_, s2w_  = point0_['S0_'], point1s['S1s'], point1['S1'], point1_['S1_'], point1w_['S1w_'], point2s['S2s'], point2['S2'], point2_['S2_'], point2s_['S2s_'], point2w_['S2w_']
        t0_, t1s, t1, t1_, t1w_, t2s, t2, t2_, t2s_, t2w_ =  point0_['T0_'], point1s['T1s'], point1['T1'], point1_['T1_'], point1w_['T1w_'], point2s['T2s'], point2['T2'], point2_['T2_'], point2s_['T2s_'], point2w_['T2w_']
        p0_, p1s, p1, p1_, p1w_, p2s, p2, p2_, p2s_, p2w_  = point0_['P0_'], point1s['P1s'], point1['P1'], point1_['P1_'], point1w_['P1w_'], point2s['P2s'], point2['P2'], point2_['P2_'], point2s_['P2s_'], point2w_['P2w_']
        v0_, v1s, v1, v1_, v1w_, v2s, v2, v2_, v2s_, v2w_ =  point0_['V0_'], point1s['V1s'], point1['V1'], point1_['V1_'], point1w_['V1w_'], point2s['V2s'], point2['V2'], point2_['V2_'], point2s_['V2s_'], point2w_['V2w_']
    
    if method == 'cold':     
        h0_, h1s, h1, h1_, h1w_, h2s, h2, h2_, h2s_, h2w_ =  point0_['I0_'], point1s['I1s'], point1['I1_sm'], point1_['I1_sm_'], point1w_['I1w_sm_'], point2s['I2s'], point2['I2_sm'], point2_['I2_sm_'], point2s_['I2s_sm_'], point2w_['I2w_sm_']
        s0_, s1s, s1, s1_, s1w_, s2s, s2, s2_, s2s_, s2w_  = point0_['S0_'], point1s['S1s'], point1['S1_sm'], point1_['S1_sm_'], point1w_['S1w_sm_'], point2s['S2s'], point2['S2_sm'], point2_['S2_sm_'], point2s_['S2s_sm_'], point2w_['S2w_sm_']
        t0_, t1s, t1, t1_, t1w_, t2s, t2, t2_, t2s_, t2w_ =  point0_['T0_'], point1s['T1s'], point1['T1_sm'], point1_['T1_sm_'], point1w_['T1w_sm_'], point2s['T2s'], point2['T2_sm'], point2_['T2_sm_'], point2s_['T2s_sm_'], point2w_['T2w_sm_']
        p0_, p1s, p1, p1_, p1w_, p2s, p2, p2_, p2s_, p2w_  = point0_['P0_'], point1s['P1s'], point1['P1_sm'], point1_['P1_sm_'], point1w_['P1w_sm_'], point2s['P2s'], point2['P2_sm'], point2_['P2_sm_'], point2s_['P2s_sm_'], point2w_['P2w_sm_']
        v0_, v1s, v1, v1_, v1w_, v2s, v2, v2_, v2s_, v2w_ =  point0_['V0_'], point1s['V1s'], point1['V1_sm'], point1_['V1_sm_'], point1w_['V1w_sm_'], point2s['V2s'], point2['V2_sm'], point2_['V2_sm_'], point2s_['V2s_sm_'], point2w_['V2w_sm_']
            
    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 10)) # параметры окна
    ax = plt.axes()
    plt.xlim((s2s - 0.02, s1 + 0.06))
    plt.ylim((h2s - 20, h0_ + 20))
    ax.yaxis.set_major_locator(LinearLocator(15)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(15))

    fig.suptitle('HS диаграмма ступени ГТУ', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('S, кДж/(кгК)', fontsize=20, loc = 'center')
    plt.ylabel('h, кДж/кг',fontsize=20,loc = 'center')

    pinch_0_ = plt.scatter(s2s, h0_, s = 40, c = 'blue') 
    pinch_1s = plt.scatter(s2s, h1s, s = 40, c = 'blue') 
    pinch_2s = plt.scatter(s2s, h2s, s = 40, c = 'red') 
    pinch_1_ = plt.scatter(s2, h1_, s = 40, c = 'blue')
    pinch_1w_ = plt.scatter(s2, h1w_, s = 40, c = 'blue')
    pinch_1 = plt.scatter(s2, h1, s = 40, c = 'blue')
    pinch_2s_ = plt.scatter(s2, h2s_, s = 40, c = 'red')
    pinch_2w_ = plt.scatter(s1, h2w_, s = 40, c = 'red')
    pinch_2_ = plt.scatter(s1, h2_, s = 40, c = 'red')
    pinch_2 = plt.scatter(s1, h2, s = 40, c = 'red')

    scat_0_1s = plt.plot([s2s, s2s], [h0_, h1s], color = 'b', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
    scat_0_1 = plt.plot([s2s, s2], [h0_, h1], color = 'b', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
    scat_1s_2s = plt.plot([s2s, s2s], [h1s, h2s], color = 'r', marker ='o', ms = 8, markerfacecolor ='w', linewidth = 3, linestyle = '-')  
    scat_1_2 = plt.plot([s2, s1], [h1, h2], color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
    scat_1_1w_ = plt.plot([s2, s2], [h1, h1w_], color = 'b', marker = 'o', ms = 8, markerfacecolor ='w', linewidth = 2, linestyle = '--')  
    scat_1_1_ = plt.plot([s2, s2], [h1, h1_], color = 'b', marker = 'o', ms = 8, markerfacecolor ='w', linewidth = 2, linestyle = '--')  
    scat_2_2_ = plt.plot([s1, s1], [h2, h2_], color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 2, linestyle = '--')  
    scat_1_2s_ = plt.plot([s2, s2], [h1, h2s_], color = 'r', marker = 'o', ms = 8, markerfacecolor ='w', linewidth = 2, linestyle = '--')  
    scat_2_2w_ = plt.plot([s1, s1], [h2, h2w_], color = 'r', marker = 'o', ms = 8, markerfacecolor ='w', linewidth = 2, linestyle = '--')  
    
    size_0_, = plt.plot([s2s, s2s - 0.015], [h0_, h0_], color = 'black', linewidth = 1, linestyle = '-')  
    size_2s, = plt.plot([s2s, s2s - 0.015], [h2s, h2s], color = 'black', linewidth = 1, linestyle = '-')  
    ax.arrow(*[s2s - 0.015, h0_, 0, -(h0_ - h2s)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    ax.arrow(*[s2s - 0.015, h2s, 0, (h0_ - h2s)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')

    size_1s, = plt.plot([s2s, s2s - 0.008], [h1s, h1s], color = 'black', linewidth = 1, linestyle = '-')  
    ax.arrow(*[s2s - 0.008, h0_, 0, -(h0_ - h1s)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    ax.arrow(*[s2s - 0.008, h1s, 0, (h0_ - h1s)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    
    ax.arrow(*[s2s - 0.008, h1s, 0, -(h1s - h2s)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    ax.arrow(*[s2s - 0.008, h2s, 0, (h1s - h2s)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    
    size_1_, = plt.plot([s2, s2 + ((s1 - s2) / 2) + 0.015], [h1_, h1_], color = 'black', linewidth = 1, linestyle = '-')  
    size_2_, = plt.plot([s1, s1 - (s1 - s2) / 2 + 0.015], [h2_, h2_], color = 'black', linewidth = 1, linestyle = '-')  
    ax.arrow(*[s2 + ((s1 - s2) / 2) + 0.015, h1_, 0, -(h1_ - h2_)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    ax.arrow(*[s2 + ((s1 - s2) / 2) + 0.015, h2_, 0, (h1_ - h2_)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')

    size_1, = plt.plot([s2, s2 + ((s1 - s2) / 2) + 0.008], [h1_, h1_], color = 'black', linewidth = 1, linestyle = '-')  
    size_2, = plt.plot([s1, s1 - ((s1 - s2) / 2) + 0.008], [h2, h2], color = 'black', linewidth = 1, linestyle = '-')  
    ax.arrow(*[s2 + ((s1 - s2) / 2) + 0.008, h1_, 0, -(h1_ - h2)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')
    ax.arrow(*[s2 + ((s1 - s2) / 2) + 0.008, h2, 0, (h1_ - h2)], width = 0.00001, length_includes_head = True, head_length = 10, head_width = 0.001, fc = 'black', ec = 'black')


    P0_ = interp_plot(np.array([s2s - 0.005, s2s, s2s + 0.005]), np.array([h0_ - 1, h0_, h0_ + 4 - i]), 'blue')
    P1 = interp_plot(np.array([s2s - 0.005, s2s, s2, s2 + 0.005]), np.array([h1s - 1, h1s, h1, h1 + 4 - i]), 'blue')
    P1_ = interp_plot(np.array([s2 - 0.005, s2, s2 + 0.005]), np.array([h1_ - 1, h1_, h1_ + 4 - i]), 'blue')
    P1w_ = interp_plot(np.array([s2 - 0.005, s2, s2 + 0.005]), np.array([h1w_ - 1, h1w_, h1w_ + 4 - i]), 'blue')
    P2_ = interp_plot(np.array([s1 - 0.005, s1, s1 + 0.005]), np.array([h2_ - 1, h2_, h2_ + 4 - i]), 'red')
    P2w_ = interp_plot(np.array([s1 - 0.005, s1, s1 + 0.005]), np.array([h2w_ - 1, h2w_, h2w_ + 4 - i ]), 'red')
    P2 = interp_plot(np.array([s2s - 0.005, s2s, s2]), np.array([h2s - 1, h2s, h2s_]), 'red')
    P_2 = interp_plot(np.array([s2, s1, s1 + 0.005]), np.array([h2s_, h2, h2 + 4 - i]), 'red')
   
    ax.annotate('{}'.format(r'$\overline{0}$'), xy = (s2s, h0_), xytext = (-10, 7), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$1_{s}$'), xy = (s2s, h1s), xytext = (-25, -25), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$1$'), xy = (s2, h1), xytext = (-15, -20), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$\overline{1}_w$'), xy = (s2, h1w_), xytext = (5, -25), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$\overline{1}$'), xy = (s2, h1_), xytext = (5, -25), ha = 'left', textcoords = 'offset points', size = 20)

    ax.annotate('{}'.format(r'$2_{s}$'), xy = (s2s, h2s), xytext = (-10, -25), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$2^\prime_{s}$'), xy = (s2, h2s_), xytext = (-10, -25), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$2$'), xy = (s1, h2), xytext = (-5, -25), ha = 'left', textcoords = 'offset points', size = 20)
    ax.annotate('{}'.format(r'$\overline{2}_w$'), xy = (s1, h2w_), xytext = (-30, -25), ha = 'left', textcoords = 'offset points', size = 20)    
    ax.annotate('{}'.format(r'$\overline{2}$'), xy = (s1, h2_), xytext = (-25, -25), ha = 'left', textcoords = 'offset points', size = 20)

    ax.annotate(r'$ \overline{P_0} = $' f'{np.round(p0_ / 1000, 2)} кПа', xy = (s2s + 0.005, h0_+ 4 - i), xytext = (-5, 5), ha = 'left', textcoords = 'offset points', size = 12)
    ax.annotate(r'$ \overline{P_1} = $' f'{np.round(p1_ / 1000, 2)} кПа', xy = (s2 + 0.005, h1_ + 4 - i), xytext = (-5, -20), ha = 'left', textcoords = 'offset points', size = 12)
    ax.annotate(r'$ \overline{P_{1w}} = $' f'{np.round(p1w_ / 1000, 2)} кПа', xy = (s2 + 0.005, h1w_ + 4 - i), xytext = (0, -15), ha = 'left', textcoords = 'offset points', size = 12)
    ax.annotate(r'$ {P_1} = $' f'{np.round(p1 / 1000, 2)} кПа', xy = (s2 + 0.005, h1 + 4 - i), xytext = (0, -10), ha = 'left', textcoords = 'offset points', size = 12)
    ax.annotate(r'$ \overline{P_2} = $' f'{np.round(p2_ / 1000, 2)} кПа', xy = (s1 + 0.005, h2_ + 4 - i), xytext = (0, -15), ha = 'left', textcoords = 'offset points', size = 12)
    ax.annotate(r'$ {P_2} = $' f'{np.round(p2 / 1000, 2)} кПа', xy = (s1 + 0.005, h2 + 4 - i), xytext = (0, -15), ha = 'left', textcoords = 'offset points', size = 12)
    ax.annotate(r'$ \overline{P_{2w}} = $' f'{np.round(p2w_ / 1000, 2)} кПа', xy = (s1 + 0.005, h2w_ + 4 - i), xytext = (0, -15), ha = 'left', textcoords = 'offset points', size = 12)
    
    ax.annotate(r'$ H_{sСТ} = $' f'{np.round((h0_ - h2s), 2)} кДж/кг', xy = (s2s - 0.015, h0_ - ((h0_ - h2s) / 2) ), xytext = (-3, -60), ha = 'right', textcoords = 'offset points', size = 12, rotation = 90)
    ax.annotate(r'$ H_{sСА} = $' f'{np.round((h0_ - h1s), 2)} кДж/кг', xy = (s2s - 0.008, h0_ - ((h0_ - h1s) / 2) ), xytext = (-3, -60), ha = 'right', textcoords = 'offset points', size = 12, rotation = 90)
    ax.annotate(r'$ H_{sРК} = $' f'{np.round((h1s - h2s), 2)} кДж/кг', xy = (s2s - 0.008, h1s - ((h1s - h2s) / 2) ), xytext = (-3, -40), ha = 'right', textcoords = 'offset points', size = 12, rotation = 90)
    ax.annotate(r'$ L_{СТ} = $' f'{np.round((h1_ - h2), 2)} кДж/кг', xy = (s2 + ((s1 - s2) / 2) + 0.008, h1_ - ((h1_ - h2) / 2) ), xytext = (-3, -40), ha = 'right', textcoords = 'offset points', size = 12, rotation = 90)
    ax.annotate(r'$ L_{СТ}^* = $' f'{np.round((h1_ - h2_), 2)} кДж/кг', xy = (s2 + ((s1 - s2) / 2) + 0.015, h1_ - ((h1_ - h2_) / 2) ), xytext = (-3, -40), ha = 'right', textcoords = 'offset points', size = 12, rotation = 90)
    
    ax.legend((pinch_0_, pinch_1s, pinch_2s, pinch_1_, pinch_1w_, pinch_1, pinch_2s_, pinch_2w_, pinch_2_, pinch_2), 
    [r'Точка $\overline {0}: $ $\overline {T}_0$ =' f'{np.round(t0_, 2)} К', 
    r'Точка $1_s: $  $T_{1s}$ =' f'{np.round(t1s, 2)} К',
    r'Точка $2_s: $  $T_{2s}$ =' f'{np.round(t2s, 2)} К',
    r'Точка $\overline {1}: $ $\overline {T}_1$ =' f'{np.round(t1_, 2)} К',
    r'Точка $\overline {1}_w: $ $\overline {T}_{1w}$ =' f'{np.round(t1w_, 2)} К',
    r'Точка $1: $  $T_1$ =' f'{np.round(t1, 2)} К',
    r'Точка $2^\prime_{s}: $  $T^\prime_{s}$ =' f'{np.round(t2s_, 2)} К',
    r'Точка $\overline {2}_w: $ $\overline {T}_{2w}$ =' f'{np.round(t2w_, 2)} К',
    r'Точка $\overline {2}: $ $\overline {T}_2$ =' f'{np.round(t2_, 2)} К',
    r'Точка $2: $  $T_2$ =' f'{np.round(t2, 2)} К'], ncol = 2, fontsize = 12, frameon = True, framealpha = True)

    plt.tick_params(axis = 'both', which = 'major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    plt.show()

def velocity_triangle_plot(C_1, W_1, U_1, alpha_1, betta_1, C_2, W_2, U_2, alpha_2, betta_2, i):

    # plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 10)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis ='both', which='major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')

    ax.yaxis.set_major_locator(LinearLocator(10)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(10))

    fig.suptitle(f'Треугольник скоростей на среднем диаметре ступени №{i + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('X, м/c', fontsize = 20, loc = 'center')
    plt.ylabel('Y, м/c', fontsize = 20, loc = 'center')

    ax.arrow(*[0, 0, -C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1))], width = 1.5, length_includes_head = True, head_length = 30, head_width = 8, fc = 'blue', ec = 'blue')
    ax.arrow(*[0, 0, -(C_1 * cos(radians(alpha_1)) - U_1), -C_1 * sin(radians(alpha_1))], width = 1.5, length_includes_head = True, head_length = 30, head_width = 8, fc = 'blue', ec = 'blue')
    ax.arrow(*[-(C_1 * cos(radians(alpha_1)) - U_1), -C_1 * sin(radians(alpha_1)), -U_1, 0], width = 1.5, length_includes_head = True, head_length = 30, head_width = 8, fc = 'blue', ec = 'blue')
    
    ax.arrow(*[0, 0, W_2 * cos(radians(betta_2)) - U_2, -W_2 * sin(radians(betta_2))], width = 1.5, length_includes_head = True, head_length = 30, head_width = 8, fc = 'red', ec = 'red')
    ax.arrow(*[0, 0, W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2))], width = 1.5, length_includes_head = True, head_length = 30, head_width = 8, fc = 'red', ec = 'red')
    ax.arrow(*[ W_2 * cos(radians(betta_2)), -C_2 * sin(radians(alpha_2)),-U_2,0], width = 1.5, length_includes_head = True, head_length = 30, head_width = 8, fc = 'red', ec = 'red')

    c1 = plt.scatter(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1)), s = 1, c = 'blue') 
    w1 = plt.scatter(-W_1 * cos(radians(betta_1)), -W_1 * sin(radians(betta_1)), s = 1, c = 'blue') 
    u1 = plt.scatter(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1)), s = 1, c = 'blue') 
    al1 = plt.scatter(-W_1 * cos(radians(betta_1)), -W_1 * sin(radians(betta_1)), s = 1, c = 'blue') 
    bet1 = plt.scatter(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1)), s = 1, c = 'blue') 

    c2 = plt.scatter(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2)), s = 1, c = 'red') 
    w2 = plt.scatter(W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2)), s = 1, c = 'red') 
    u2 = plt.scatter(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2)), s = 1, c = 'red') 
    al2 = plt.scatter(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2)), s = 1, c = 'red') 
    bet2 = plt.scatter(W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2)), s = 1, c = 'red')


    p0 = [(-C_1 * cos(radians(alpha_1)), 0), (0, 0)]
    p1 = [(0, W_2 * cos(radians(betta_2))), (0, 0)]
    point, = ax.plot(*(0, 0), color = 'black', marker = 'o', ms = 8, markerfacecolor = 'w')
    line0, = plt.plot(*p0, color = 'black', linewidth = 2, linestyle = '-')
    line1, = plt.plot(*p1, color = 'black', linewidth = 2, linestyle = '-')

    p2 = [(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1))), (0, 0)]
    p3 = [(-W_1 * cos(radians(betta_1)), -W_1 * sin(radians(betta_1))), (0, 0)]    

    p4 = [(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2))), (0, 0)]
    p5 = [(W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2))), (0, 0)]

    alpha1 = AngleAnnotation((0, 0), p0[0], p2[0], ax = ax, size = 150, text = r"$\alpha_1$", color = "black", textposition = "pixels", text_kw = dict(fontsize = 20, color = "black"))
    betta1 = AngleAnnotation((0, 0), p0[0], p3[0], ax = ax, size = 500, text = r"$\beta_1$", color = "black", textposition = "inside", text_kw = dict(fontsize = 20, color = "black"))
    alpha2 = AngleAnnotation((0, 0), p4[0], p1[1], ax = ax, size = 400, text = r"$\alpha_2$", color = "black", textposition = "inside", text_kw = dict(fontsize = 20, color = "black"))
    betta2 = AngleAnnotation((0, 0), p5[0], p1[1], ax = ax, size = 150, text = r"$\beta_2$", color = "black", textposition = "pixels", text_kw = dict(fontsize = 20, color = "black"))

    ax.annotate(r'$C_1$', xy = (-C_1 * cos(radians(alpha_1)) / 2,  -C_1 * sin(radians(alpha_1)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$W_1$', xy = (-W_1 * cos(radians(betta_1)) / 2,  -W_1 * sin(radians(betta_1)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$U_1$', xy = (((-C_1 * cos(radians(alpha_1))) + (-W_1 * cos(radians(betta_1)))) / 2, ((-C_1 * sin(radians(alpha_1))) + (-W_1 * sin(radians(betta_1)))) / 2 ), xytext = (0, 10), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)

    ax.annotate(r'$C_2$', xy = (C_2 * cos(radians(alpha_2)) / 2,  -C_2 * sin(radians(alpha_2)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$W_2$', xy = (W_2 * cos(radians(betta_2)) / 2,  -W_2 * sin(radians(betta_2)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$U_2$', xy = (((C_2 * cos(radians(alpha_2))) + (W_2 * cos(radians(betta_2)))) / 2, ((-C_2 * sin(radians(alpha_2))) + (-W_2 * sin(radians(betta_2)))) / 2 ), xytext = (0, 10), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)

    leg1 = ax.legend((c1, w1, u1, al1, bet1), 
    [r'$C_1 = $ ' f'{np.round(C_1, 2)} м/c',
    r'$W_1 = $ ' f'{np.round(W_1, 2)} м/c',
    r'$U_1 = $ ' f'{np.round(U_1, 2)} м/c',
    r'$\alpha_1 = $ ' f'{np.round(alpha_1, 2)} град',
    r'$\beta_1 = $ ' f'{np.round(betta_1, 2)} град'], loc = 2,  fontsize = 12, frameon = True, framealpha = True)

    plt.gca().add_artist(leg1)

    leg2 = ax.legend((c2, w2, u2, al2, bet2), 
    [r'$C_2 = $ ' f'{np.round(C_2, 2)} м/c',
    r'$W_2 = $ ' f'{np.round(W_2, 2)} м/c',
    r'$U_2 = $ ' f'{np.round(U_2, 2)} м/c',
    r'$\alpha_2 = $ ' f'{np.round(alpha_2, 2)} град',
    r'$\beta_2 = $ ' f'{np.round(betta_2, 2)} град'], loc = 1,  fontsize = 12, frameon = True, framealpha = True)

    fig.set_figwidth(15)
    fig.set_figheight(8)
    plt.show() 
