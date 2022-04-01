import os

import numpy as np
from PyQt5.QtGui import QPainter
from numpy import sin, cos, sqrt, deg2rad, zeros, linspace
import pyqtgraph as pg
import pandas as pd
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
from functools import partial
from scipy.integrate import solve_ivp

import sys


class fluid:
    def __init__(self):
        self.g = 9.81
        self.ro = 1


class air(fluid):
    def __init__(self):
        super().__init__()
        self.ro = 1.225
        """- h -| - t -         - g -                       - p -           - ρ -       - μ -
                Temperature| Acceleration of Gravity| Absolute Pressure|  Density| Dynamic Viscosity
          (m)          (deg_C)	|       (m/s2)      |    (10^4 N/m2)	    |(kg/m3)     (10^-5 N s/m2)
        
            
        """
        self.decity_mat = [[-1000, 21.50, 9.810, 11.39, 1.347, 1.821],
                           [0, 15.00, 9.807, 10.13, 1.225, 1.789],
                           [1000, 8.50, 9.804, 8.988, 1.112, 1.758],
                           [2000, 2.00, 9.801, 7.950, 1.007, 1.726],
                           [3000, -4.49, 9.797, 7.012, 0.9093, 1.694],
                           [4000, -10.98, 9.794, 6.166, 0.8194, 1.661],
                           [5000, -17.47, 9.791, 5.405, 0.7364, 1.628],
                           [6000, -23.96, 9.788, 4.722, 0.6601, 1.595],
                           [7000, -30.45, 9.785, 4.111, 0.5900, 1.561],
                           [8000, -36.94, 9.782, 3.565, 0.5258, 1.527],
                           [9000, -43.42, 9.779, 3.080, 0.4671, 1.493],
                           [10000, -49.90, 9.776, 2.650, 0.4135, 1.458],
                           [15000, -56.50, 9.761, 1.211, 0.1948, 1.422],
                           [20000, -56.50, 9.745, 0.5529, 0.08891, 1.422],
                           [25000, -51.60, 9.730, 0.2549, 0.04008, 1.448],
                           [30000, -46.64, 9.715, 0.1197, 0.01841, 1.475],
                           [40000, -22.80, 9.684, 0.0287, 0.00399, 61.601],
                           [50000, -2.5, 9.654, 0.007978, 0.001027, 1.704],
                           [60000, -26.13, 9.624, 0.002196, 0.00030, 97, 1.584],
                           [70000, -53.57, 9.594, 0.00052, 0.00008283, 1.438],
                           [80000, -74.51, 9.564, 0.00011, 0.00001846, 1.321]]


class bullet:
    def __init__(self, d=0.27 * 0.0254, m=0.01, v0=400, cb=None, cd=None):

        self.conv = 703.0696
        self.mass = m
        self.v0 = v0
        self.d = d

        self.cb = cb
        self.cd = cd

        self.form = 1
        self.ro = 10
        self.cd_ef = 1
        if cb:
            self.calc_drag()

    def calc_bal(self):
        self.cb = self.mass / (self.d ** 2 * self.form)

    def calc_bal_cd(self):
        self.cb = self.mass * 4 / (self.cd * self.d ** 2 * np.pi) / self.conv

    def calc_drag(self):
        self.cd_ef = 0.5 / (self.cb * self.conv)  # =0.5*(cd*A/m) since this turns fotce into a

    def print_stats(self):
        pass

    def calc_form(self, l):
        n = l ** 2 + 0.25
        self.form = 2 / n * np.sqrt(4 - 1 / n)

    def calc_drag_v(self, v):
        self.cd = 8 / (self.ro * np.pi * v ** 2 * self.d ** 2)

    def calc_self(self):
        if self.cb:
            self.calc_drag()  # todo else


class tikka(bullet):
    def __init__(self):
        super().__init__(0.27 * 0.0254, 0.01,400,0.4)
        self.calc_drag()


class scope(QWidget):
    def __init__(self):
        super().__init__()
        self.t = 'mil'
        self.dyz = np.zeros(2)
        self.ang = np.zeros(2)
        self.clicks = np.zeros(2)
        self._lable_set()

    def _lable_set(self):
        self.layout = QGridLayout()
        self.setLayout(self.layout)
        self.click_lab = []
        for n, i in enumerate(['H', 'V']):
            lab = QLabel(i)
            jj = QLineEdit()
            jj.setReadOnly(True)
            jj.setText('0')
            self.click_lab.append(jj)
            self.layout.addWidget(jj, n, 1)
            self.layout.addWidget(lab, n, 0)

    def set_v(self, disp):
        self.dyz = disp[1:]
        ang = np.arctan(self.dyz / disp[0])
        self.ang = np.rad2deg(ang)
        if self.t == 'mil':  # ie millrad
            self.clicks = ang * 1000
        else:
            self.clicks = self.ang * 60
        self.clicks = self.clicks.astype(int)
        for i in range(2):
            self.click_lab[i].setText(str(self.clicks[i]))
        #


class Arrow(QWidget):
    def __init__(self):
        super().__init__()
        # self.setMinimumSize(15,15)
        self.s1 = np.array([30, 30])
        self.s2 = self.s1 // 2
        self.setFixedSize(*self.s1)
        self.line = [0, 0]

    def paintEvent(self, event):
        p = QPainter(self)
        p.setPen(Qt.blue)
        points = self.s2 + self.line
        p.drawLine(*self.s2, *points)

    # def resizeEvent(self, event):
    #     self.scale


class wind(QWidget):
    def __init__(self, par, v=10, v_ang=0, ang=0):
        super().__init__()
        self.par = par
        self.w_vect = np.zeros(3)
        self.ang = ang
        self.v = v
        self.v_ang = v_ang
        self._set_layout()

    def _set_layout(self):
        self.lay = QVBoxLayout()
        l1 = QHBoxLayout()
        self.scale_factor = 3
        self.arrow = Arrow()
        l2 = QHBoxLayout()
        l2.addWidget(self.arrow)

        self.ang_rot = QDial()
        self.ang_rot.setRange(0, 360)
        self.ang_rot.setValue(270 + self.ang)
        self.ang_rot.valueChanged.connect(self.set_ang)

        self.v_slide = QSlider()
        self.v_slide.setRange(0, 50)
        self.v_slide.setValue(self.v)
        self.v_slide.valueChanged.connect(self.set_v)

        self.angle_slide = QSlider()
        self.angle_slide.setRange(-10, 10)
        self.angle_slide.setValue(self.v_ang)
        self.angle_slide.valueChanged.connect(self.set_ang)
        l1.addWidget(self.ang_rot)
        l1.addWidget(self.angle_slide)
        l1.addWidget(self.v_slide)
        self.lay.addLayout(l2)
        self.lay.addLayout(l1)
        self.setLayout(self.lay)
        self.paint_v()

    def set_ang(self):
        ang_v = 270 - self.ang_rot.value()
        self.ang = deg2rad(ang_v % 360)
        self.rs()

    def set_v(self):
        self.v = self.v_slide.value()
        self.rs()

    def set_v_ang(self):
        self.v_ang = deg2rad(self.angle_slide.value())
        self.rs()

    def rs(self):
        self.set_eq()
        self.paint_v()
        self.par.upd(self.w_vect)

    def paint_v(self):
        if self.v != 0:
            self.arrow.line = self.v * self.scale_factor * np.array((cos(self.ang), -sin(self.ang)))
            self.arrow.update()

    def rot(self):
        return np.array([cos(self.v_ang) * cos(self.ang), cos(self.v_ang) * sin(self.ang), sin(self.v_ang)])

    def set_eq(self):
        self.w_vect = self.rot() * self.v

    def r_v(self):
        self.v = np.random.rand(0, 30)

    def r_a(self):
        self.ang = np.random.rand(0, 30)


class dragOb(pg.GraphicsLayoutWidget):
    def __init__(self, parent=None, th=5, v=300):
        super().__init__(parent)

        self.v0 = v  # km / h

        self.wind = np.zeros(3)  # par hold wind object
        self.theta = th

        self.max_n = [10000, 1000]
        self.lims = [0, 100]

        d_key = ['vmag', 'vx', 'vy', 'vz',
                 't', 'x', 'y', 'z', 'inv',
                 'ang', 'E', 'E_norm', 'in ang']
        self.data = {i: 0 for i in d_key}
        self.conversions = {'v': {'m/s': 3.6, 'mph': 1.6 * 3.6, 'kph': 1 / 3.6},
                            'ang': {'rad': np.rad2deg, 'deg': np.deg2rad}}

        self._set_plot()
        self._set_const()
        self.lim_su()
        self.update()

    def _set_plot(self):
        # when not image
        self.xyz = {'xz': self.addPlot(1, 0, colspan=2), 'xy': self.addPlot(0, 1), 'yz': self.addPlot(0, 0)}
        self.plots = {}

        for i, j in self.xyz.items():
            self.plots[i] = j.plot(pen='c', width=3)
        self.t = np.array([])
        self.pos = np.array([])
        self.vel = np.array([])
        # self.z = []
        self.eq_vals = np.array([])

    def lim_su(self):
        nn = [['xMin', 'xMax'], ['yMin', 'yMax']]
        jj = [[0, 20], [0, 10], [10, 20]]
        for n, i in enumerate(self.xyz.values()):
            if n == 2:
                i.invertX(True)
            for ii in range(2):
                if jj[n][ii] == 0:
                    k = {nn[ii][0]: 0}
                else:
                    k = {nn[ii][0]: -jj[n][ii], nn[ii][1]: jj[n][ii]}
                i.setLimits(**k)

    def _set_const(self):
        print('seting const')
        self.fluid = air()
        self.bullet = tikka()

        self.rad = np.deg2rad(self.theta)
        self.z0 = np.zeros(4)

    def set_init(self):
        self.rad = np.deg2rad(self.theta)
        self.z0 = np.array([0, 0, 0, np.cos(self.rad), 0, np.sin(self.rad)]) * self.v0

    def _set_no_air(self):
        self.no_air = np.array((np.cos(self.rad), sin(self.rad))).reshape(2, 1) * self.v0 * self.t

    def solve_x(self):

        def f_x(t, z):  # x,y,z; vx,vy, vz:::vx,vy, vz; ax, ay, az,,
            rel_z = z[3:] - self.wind  # for wind only affect the acelleration w
            l2 = -self.fluid.ro * rel_z * np.linalg.norm(rel_z) * self.bullet.cd_ef  # cdA/m, a=f/m
            l2[2] -= self.fluid.g

            return np.concatenate([z[3:], l2])

        eq_vals = solve_ivp(f_x, self.lims, self.z0)
        self.t = eq_vals.t
        xv = eq_vals.y  # .reshape((2,3, -1))
        self.pos = xv[:3]  # , 0]
        self.vel = xv[3:]

    def set_sc(self):
        self.setXRange(0, self.lims[0], padding=0.005)
        self.setYRange(0, self.lims[1], padding=0)

    def update(self, w=None):
        if w is not None:
            self.wind = w
        self.set_init()
        self.solve_x()

        self.plots['xz'].setData(self.pos[0], self.pos[2])
        self.plots['xy'].setData(self.pos[0], self.pos[1])
        self.plots['yz'].setData(self.pos[2], -1 * self.pos[1])
        self.calc_ground()

    def e_calc(self, v):
        return v ** 2 * 0.5 * self.bullet.mass

    def set_outputs(self, ic, ax=2, d=0):
        def interp(x_p, x_v, y_v):
            return (y_v[1] - y_v[0]) / (x_v[1] - x_v[0]) * (x_p - x_v[0]) + y_v[0]

        interp_x = self.pos[ax, ic - 1:ic + 1]

        for n, x in enumerate('xyz'):  # range(3):
            interp_y = self.vel[n, ic - 1:ic + 1]
            interp_y_p = self.pos[n, ic - 1:ic + 1]
            cor_dat = interp(d, interp_x, interp_y)
            self.data['v' + x] = cor_dat
            self.data[x] = interp(d, interp_x, interp_y_p)
        # print(f'interpx {interp_x}, interp{interp_y_p}')

        self.data['vmag'] = np.linalg.norm(self.vel)
        self.data['t'] = self.t[ic]

        self.data['inv'] = self.v0
        self.data['in ang'] = self.theta

        self.data['ang'] = np.arctan(self.vel[2, ic] / np.linalg.norm(self.vel[:2, ic]))
        self.data['E'] = self.e_calc(self.data['vmag'])
        self.data['Emag'] = self.e_calc(self.data['vz'])
        self.data['E_norm'] = self.data['E'] * sin(self.data['ang'])
        return np.array((self.data['x'], self.data['y'], self.data['z']))

    def calc_ground(self):
        ic = np.where(self.pos[2] < 0)[0][0]
        # print('Ground at: ', ic)
        self.set_outputs(ic)
        # self.data[x] = y[incept]
        # excel
        # save
        pass

    def calc_target(self, d):
        incept = np.where(self.pos[0] > d)[0][0]
        # center based on line sigt dist to find height, hdist
        return self.set_outputs(incept, 0, d)

    def set_target(self, d, ang):
        # r = hyp
        dist = np.array(cos(ang), sin(ang)) * d
        # todo

        """ resolve for initial angle for now just use dist x, todo also set for zero range,
        mode: edit val, mode predefined val
        todo temp wind, altitude, dirrection, min output
        todo negative ange
        todo for air0
        for when target at angle: ie solve theta = x,y, solve zeroped d, adjust
        # todo bullet v0# todo barrel len# todo draw center todo rep
        # todo gtype# todo mass grain, d in"""
        self.calc_target(dist[0])

    def save_vals(self):
        dirs = os.getcwd()
        def_name = dirs + '/' + 'Bullet Angle Drag.csv'
        # noinspection PyArgumentList
        f = QFileDialog.getSaveFileName(self, 'Hello', def_name, 'CSV (*.csv)')
        f_0 = f[0]
        if f_0 != '':
            dirn = os.path.dirname(f_0)

            if not os.path.exists(dirn):
                os.makedirs(dirn)

            pd_d = []
            t = np.linspace(1, 90, 89)
            # th = np.deg2rad(t)
            w = np.zeros(3)

            for i in t:
                self.theta = i
                self.update(w)
                # print(self.data)
                pd_d.append({i: j for i, j in self.data.items()})

            df = pd.DataFrame.from_records(pd_d)

            df.to_csv(f[0])


class Window(QMainWindow):
    # noinspection PyArgumentList
    def __init__(self):
        super().__init__()

        # self.setAcceptDrops(True)
        self.defult = [5, 300, 50]
        self.setWindowTitle('QMainWindow')
        self.drag = dragOb(v=self.defult[1], th=self.defult[0])
        self.setCentralWidget(self.drag)
        self.new_d = None
        self.running = False
        self.ls = []

        self._set_tool()
        self._set_data()
        self._set_tar()
        self._set_tool_bar()

    def _set_tool_bar(self):
        self.tb = QToolBar()
        self.addToolBar(self.tb)
        mo = ['Bullet perameters defined', 'Solve Mode', 'full control', 'Save', 'Load', 'New Bull']
        self.at = {}
        for i in mo:
            j = QAction(i)
            self.at[i] = j
            self.tb.addAction(j)
            j.triggered.connect(partial(self.mode_set, i))

    def mode_set(self, mo):
        ls = []
        if mo == 'Save':  # todo else grey
            self.drag.save_vals()

        elif mo == 'Bullet perameters defined':
            ls = ['v']
        elif mo == 'Solve Mode':
            ls = ['v', 'ang']
        elif mo == 'Load':
            self.load_j()
            ls = ['v']
        elif mo == 'New Bull':
            ls = ['v']
            self.set_con()
        self.set_grey(ls)

    def set_grey(self, ls):  # todo perm list to hold affected so wont work, but reset next
        print('mode, ls: ', ls)

        for i in self.ls:
            self.vals[i].setReadOnly(False)
            self.vals[i].setStyleSheet("background-color: lightgrey; }")
        self.ls = []
        for i in ls:
            self.vals[i].setReadOnly(True)
            self.vals[i].setStyleSheet('background-color:  grey; }')
            self.ls.append(i)

    # noinspection PyArgumentList
    def _set_tool(self):
        tool_wig = QWidget()
        tool_dock = QDockWidget('Tools')
        tool_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
        tool_dock.setFeatures(QDockWidget.DockWidgetMovable)
        tool_dock.setWidget(tool_wig)
        self.addDockWidget(Qt.BottomDockWidgetArea, tool_dock)

        self.layout = QGridLayout()
        tool_wig.setLayout(self.layout)

        self.input_v = ['ang', 'v', 'tar']
        # self.rd = ['deg', 'rad']

        self.r_v = {}
        self.vals = {}
        """check if rad"""
        # for i in range(2):
        #     j = QRadioButton(self.rd[i])
        #     self.layout.addWidget(j, 0, i)
        #     j.clicked.connect(partial(self.rad, self.rd[i]))
        #     self.r_v[self.rd[i]] = j

        for n, i in enumerate(self.input_v):
            val = QLineEdit()
            lab = QLabel(i)
            val.setText(str(self.defult[n]))
            self.vals[i] = val
            self.layout.addWidget(lab, 1, n)
            self.layout.addWidget(val, 2, n)
            val.editingFinished.connect(partial(self.change_v, i))

        self.d = QDial()
        self.d.setMinimum(0)
        self.d.setMaximum(90)
        self.s = QSlider()
        self.s.setMaximum(1000)
        self.s.setMinimum(0)

        self.d.setValue(self.defult[0])
        self.s.setValue(self.defult[1])
        self.s.valueChanged.connect(self.set_s)
        self.d.valueChanged.connect(self.set_d)

        self.layout.addWidget(self.d, 3, 0)
        self.layout.addWidget(self.s, 3, 1)

    def _set_data(self):
        # noinspection PyArgumentList
        self.data_w = QWidget()
        self.data_dock = QDockWidget('data')
        self.data_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
        self.data_dock.setFeatures(QDockWidget.DockWidgetMovable)
        self.data_dock.setWidget(self.data_w)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.data_dock)

        self.wind = wind(self, v=self.defult[2])
        self.w_dock = QDockWidget('Wind')
        self.w_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
        self.w_dock.setFeatures(QDockWidget.DockWidgetMovable)
        self.w_dock.setWidget(self.wind)
        self.addDockWidget(Qt.RightDockWidgetArea, self.w_dock)

        d_lay = QGridLayout()
        self.data_w.setLayout(d_lay)

        self.data_v = {}
        n = 0
        ni = 0
        for i, j in self.drag.data.items():
            m = QLabel(i)
            k = QLineEdit(str(j))
            self.data_v[i] = k
            k.setReadOnly(True)
            # noinspection PyArgumentList
            d_lay.addWidget(m, n, ni)
            # noinspection PyArgumentList
            d_lay.addWidget(k, n + 1, ni)
            ni += 1
            if ni >= 4:
                ni = 0
                n += 2

    def _set_tar(self):
        self.tar = scope()
        tar_dock = QDockWidget('Target')
        tar_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)
        tar_dock.setFeatures(QDockWidget.DockWidgetMovable)
        tar_dock.setWidget(self.tar)
        self.addDockWidget(Qt.BottomDockWidgetArea, tar_dock)

    def reset_vals(self):
        for i, j in self.drag.data.items():
            jj = round(j, 2)
            self.data_v[i].setText(str(jj))

    def change_v(self, i):
        val = int(self.vals[i].text())
        if i == 'ang':
            self.d.setValue(val)
        elif i == 'v':
            self.s.setValue(val)
        else:
            dis = self.drag.calc_target(val)
            self.tar.set_v(dis)

    def set_d(self):
        # print('connected, d')
        v = self.d.value()
        self.drag.theta = v
        self.vals['ang'].setText(str(v))
        self.upd()
        pass

    def set_s(self):
        # print('connected, s')
        v = self.s.value()
        self.drag.v0 = v
        self.vals['v'].setText(str(v))
        self.upd()
        pass

    def set_wind(self):
        self.upd(self.wind.w_vect)

    def upd(self, w=None):
        self.drag.update(w)
        self.reset_vals()
        dis = self.drag.calc_target(int(self.vals['tar'].text()))
        self.tar.set_v(dis)

    def save_load(self):
        ls = []
        ls_2 = []

        ra = np.linspace(1, 90)
        v = 300
        self.drag.v0 = v
        for i in ra:
            self.drag.theta = i
            self.drag.update()
            va = self.drag.data
            ls.append(va)
            ls_2.append(va['vmag'])
        # setData(ra, ls_2)
        # self.on_save(ls)

    def load_j(self):
        pa = os.getcwd()+'/Bull_Config.json'
        if os.path.exists(pa):
            df = pd.DataFrame.from_json(pa)
            di = df.to_dict()
            self.new_d = bullet()
            self.new_d.vals = di

            self.init_bull()
            # self.bu(d)
        else:
            self.set_con()

    def set_con(self):
        w = wizardGun(self)
        if w.exec_():
            self.init_bull()

    def init_bull(self):
        self.drag.bullet = self.new_d
        if self.new_d.v0 > self.s.maximum():
            self.s.setMaximum(self.new_d.v0)
        self.s.setValue(int(self.new_d.v0))
        self.new_d.calc_self()
        self.drag.lim_su()
        self.drag.update()

class wizardGun(QWizard):
    def __init__(self, parent):
        # noinspection PyArgumentList
        super().__init__(parent)
        self.par = parent
        self.la = 1
        self.pa = {}
        self.finished.connect(self.fn)
        print('Wizard')
        plist =['Intro', 'Bullet Basic', 'Bullet Advanced', 'Conclusion']
        for n, i in enumerate(plist):
            j = QWizardPage(self)
            j.setTitle(i)
            self.setPage(n, j)
            self.pa[i] = j

        self._set_intro()
        self._set_bul_bas('Bullet Basic')
        # self._set_bul_bas('Bullet Advanced')
        # self._detail()
        self._conc()

    def _set_intro(self):
        lay = QGridLayout()
        self.but = []
        self.pa['Intro'].setLayout(lay)
        lay.addWidget(QLabel('Basic or advanced'),0,0,1,2)
        for n, i in enumerate(['Basic', 'Advanced']):
            ii = QRadioButton(i)
            lay.addWidget(ii,1,n)
            if n == 0:
                ii.click()
            ii.clicked.connect(partial(self.layer, n + 1))
            self.but.append(ii)

    def _conc(self):
        pass

    def layer(self,i):
        if self.but[i].isChecked(): # assuming just got checked
            self.la = i
        else:
            self.but[(i+1) % 2].click()

    def nextId(self):
        print(self.currentId())
        if self.currentId() == 0:
            return 1
        elif self.currentId() == 3:
            return -1
        else:
            return 3

    # noinspection PyArgumentList
    def _set_bul_bas(self,pag):

        lay = QGridLayout()
        self.pa[pag].setLayout(lay)
        self.basic_per = {}
        self.bull_proper = {}
        n = 0
        for i in ['Bullet Mass', 'Muzzel V', 'Balisic Coef']:
            lay.addWidget(QLabel(i), n, 0)
            j = QLineEdit()
            j.editingFinished.connect(partial(self.bul_set, i))
            self.pa[pag].registerField(i+'*', j)
            lay.addWidget(j, n, 1)
            self.basic_per[i] = j
            n += 1

    def bul_set(self, i):
        fac = None
        conversion = {'mm': 0.001, 'cm': 0.01, 'm/s': 1, 'km/h': 1 / 3.6, 'kg/m**2': 1,
                      'm': 1, 'in': 0.0254, 'gr': 0.1,  'psi': 706, 'lb': 0.8, 'g': 0.001, 'ft/s': 0.3}
        v = self.basic_per[i].text()

        for jj, j in conversion.items():
            if jj in v:
                fac = j
                v = v.replace(jj, '')

        f = float(v)
        if i == 'Bullet D':

            if f < 1:
                de = 'in'
            else:
                de = 'mm'
        elif i == 'Bullet Mass':
            if f > 10:
                de = 'gr'
            else:
                de = 'g'
        elif i == 'Muzzel V':
            de = 'm/s'
        else:
            if f < 50:
                de = 'psi'
            else:
                de = 'kg/m**2'

        fac_2 = fac if fac else conversion[de]
        val = f * fac_2
        self.bull_proper[i] = val

    def _detail(self):
        self.det = QWizardPage(self)
        self.addPage(self.det)
        lay = QGridLayout()
        self.det.setLayout(lay)

    def fn(self):
        bul = bullet(m=self.bull_proper['Bullet Mass'],v0=self.bull_proper['Muzzel V'])
        if self.hasVisitedPage(2):

            bul.d = self.bull_proper['d']
            bul.form = self.bull_proper['form']
        else:
            bul.cb = self.bull_proper['Balisic Coef']
            bul.calc_drag()

        self.par.new_d = bul
        self.save_config()

    def save_config(self):
        df = pd.DataFrame.from_dict({i:[j] for i,j in self.bull_proper.items()})
        df.to_json('Bull_Config.json')
        pass

if __name__ == "__main__":
    app = QApplication(sys.argv)
    audio_app = Window()
    audio_app.show()
    sys.exit(app.exec_())

    # todo add val condiition, angle zero
    # todo gen bullet gen gun # todo fix units
    # todo change mode
    # todo _detail bul geometry
    # todo run bull func in bul
    # todo later when advanced: self.la,[nextId]
    # len()
    # d, m known
    # disp image?

