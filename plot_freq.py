# -*- coding: utf-8 -*-
"""
Universidade Estadual de Campinas
IM 342 - Análise de Máquinas Rotativas
Pedro Lucas - Ra 263117

Resposta em Frequência do Sistema
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt  
import vibration_toolbox as vbt

de = 12E-3
le = 850E-3
dd = 100E-3
ld = 10E-3
a = 2/3*le
b = le-a
rho = 7850
vol_d = (np.pi*((dd/2)**2))*ld
m = rho*vol_d
E = 200E9
I = np.pi*(de**4)/64
Id = m*((dd/2)**2)/4
Ip = m*((dd/2)**2)/2

k1_n = (a**3)+(b**3)
k1_d = (a**3)*(b**3)
k2_n = le*(a-b)
k2_d = (a**2)*(b**2)
k3_n = le
k3_d = a*b

k1 = 3*E*I*(k1_n/k1_d)
k2 = 3*E*I*(k2_n/k2_d)
k3 = 3*E*I*(k3_n/k3_d)

f1 = []
f2 = []
f3 = []
f4 = []

W = 0
alpha = 2E-4
freq_hz = 10000
freq_ref = round(freq_hz/0.16) #[rad/s]

#for W in range(0,freq_ref,6):    
M = np.array([[m,0,0,0],
              [0,m,0,0],
              [0,0,Id,0],
              [0,0,0,Id]])

G = np.array([[0,0,0,0],
              [0,0,0,0],
              [0,0,0,(-Ip*W)],
              [0,0,(Ip*W),0]])

K = np.array([[k1,0,0,-k2],
              [0,k1,k2,0],
              [0,k2,k3,0],
              [-k2,0,0,k3]])

damp = alpha*K
C = damp + G

sys = vbt.VibeSystem(M,C,K)

def plot_freq_resp(self, modes=None, ax0=None, ax1=None, **kwargs):
    
    if ax0 is None or ax1 is None:
        fig, ax = plt.subplots(2)
        if ax0 is not None:
            _, ax1 = ax
        if ax1 is not None:
            ax0, _ = ax
        else:
            ax0, ax1 = ax

    omega, magdb, phase = self.freq_response(modes=modes)

    ax0.plot(omega, magdb[0,0], **kwargs, color='r', label='H11')
    ax1.plot(omega, phase[0,0], **kwargs, color='r')
    
    # ax0.plot(omega, magdb[0,2], **kwargs, color='g', label='H13')
    # ax1.plot(omega, phase[0,2], **kwargs, color='g')    
    
    ax0.plot(omega, magdb[0,3], **kwargs, color='b', label='H14')
    ax1.plot(omega, phase[0,3], **kwargs, color='b')
    
    ax0.plot(omega, magdb[1,2], **kwargs, color='orange', label='H23')
    ax1.plot(omega, phase[1,2], **kwargs, color='orange')

    # ax0.plot(omega, magdb[1,3], **kwargs, color='m', label='H24')
    # ax1.plot(omega, phase[1,3], **kwargs, color='m')
     
    ax0.plot(omega, magdb[2,2], **kwargs, color='c', label='H33')
    ax1.plot(omega, phase[2,2], **kwargs, color='c')

    # ax0.plot(omega, magdb[1,0], **kwargs, color='r', label='H20')
    # ax1.plot(omega, phase[1,0], **kwargs, color='r')    
    
    # ax0.plot(omega, magdb[2,0], **kwargs, color='y', label='H30')
    # ax1.plot(omega, phase[2,0], **kwargs, color='y')

   
    for ax in [ax0, ax1]:
        ax.set_xlim(0,max(omega))
        ax.yaxis.set_major_locator(
            mpl.ticker.MaxNLocator(prune='lower'))
        ax.yaxis.set_major_locator(
            mpl.ticker.MaxNLocator(prune='upper'))
    
    ax0.set_title('Resposta em Frequência do Sistema Amortecido')
    ax0.set_ylabel('Magnitude $(dB)$')
    legend = ax0.legend(loc='upper right', shadow = True)
    legend.get_frame().set_facecolor('white')
    ax1.set_ylabel('Ângulo de Fase $(°)$')
    ax1.set_xlabel('Frequência (rad/s)')
        
    return ax0, ax1

plot_freq_resp(sys)
# plt.savefig('freq_resp_amortecido1.pdf',dpi=600)



