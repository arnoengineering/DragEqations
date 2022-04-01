    # def upd(self):
    #     # fft
    #     self.line_image()
    #     data_np = self.mic.read()
    #     data_np = data_np * self.input_scale
    #     freq_fft = np.abs(rfft(data_np)/self.chunk)
    #     # print(freq_fft.max())
    #     # self.save_data()
    #     psd = np.log10(freq_fft) * 20
    #     print(psd.max())
    #     print(data_np)
    #     # roll down one and replace leading edge with new data
    #     for i in self.view_ls:
    #         if i != 'img':
    #             print('window: ', i)
    #             self.windows[i].update(psd)
    #     self.windows['tools'].inp.update_vals([self.mic.m, psd.max()])
    #     self.windows['start'].waveform_plot.setData(self.windows['start'].x, data_np)
"""
    # def explicit_euler(self, func_mat, dt, z0): # zout, l_f, i=10000
    # 
    #     zout = np.array(z0) # z0 array
    #     i = 0
    #     while i < self.max_n[0]:
    #         zout[i + 1] = zout[i]+f(zout[i]) * dt
    #         if zout[i + 1, 2] <= 0:  #  break when y < 0
    #             break
    #         i += 1
    # 
    #     del_f = (zout[-1, 1] - zout[-2, 1]) / abs(zout[-1, 2] - zout[-2, 2]) * zout[-2, 2]
    #     l_f = del_f + zout[-2, 1]
    #     # tout = np.linspace[0, dt * (length(zout) - 1), length(zout) - 1)
    #    return l_f, zout, i
    # def drag_eq(self, z):
    # 
    # 
    # def main_s(self):
    #     L0 = 0
    #     Lth0 = [0, 0]
    #     i = 0
    # 
    #     while i < self.max_n[1]:
    #         zv, L, i = self.explicit_euler(drag, t, z0) # calculates
    #         er = 1 - L0 / L
    #         L0 = L # calculates
    #         if er <= er_tol: # bbrakes loop if target reached
    #             break
    #         i = i + 1
    #         t = t / 2
plot(zv(:, 1), zv(:, 2), 'b')
title('troectory of stone')
xlabel('Position x (m)')
ylabel('Position y (m)') 

# plot

plot(zv(:, 1), zv(:, 2), 'b')
fplot(norm_eqx, norm_eqy, [0, t_m], 'r')

xlabel('Position x (m)')
ylabel('Position y (m)')
legend('damped', 'undaped')


c = c2
theta_lim = [0, 0, 90] # flaceholder

# loops untiler is < tol,, for each use domain to calc max l by looping
# through, then use max values for next itter
while i2 < maxn:
    sub_j = 10 # 10 subintervals

    Lth_c = zeros((sub_j, 2)) # reset
    th_dom = linspace(min(theta_lim), max(theta_lim), sub_j) # range

    for j in range(sub_j):
        th = th_dom(j) # get
        rad = deg2rad(th) # rad

        z02 = np.array([0, 0, cos(rad), sin(rad)]) * U2 # initial
        zv, L, _ = explicit_euler(drag, t2, z02) # solve
        Lth_c[j] = [th, L] # appends


    Lmax, ind = np.max(Lth_c, 3, 1) # matrix  # todo nlargest, rocet dyn

    theta_lim = Lth_c[ind[:,2], 1] # uses

    Lth = [Lth_c[ind[1,2], 1], Lmax[1, 2]] # same

    er = np.abs(Lth - Lth0) # vector

    Lth0 = Lth # set
    if all(er < er_tol2): # since er is a vector can use all
        break
"""