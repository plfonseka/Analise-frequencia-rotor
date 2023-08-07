# -*- coding: utf-8 -*-
"""
Universidade Estadual de Campinas
IM 342 - Análise de Máquinas Rotativas
Pedro Lucas - Ra 263117

Diagrama dos Fatores de Amortecimento
"""
import numpy as np
import vibration_toolbox as vbt
import matplotlib.pyplot as plt

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

for W in range(0,freq_ref,6):    
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
    freq = sys.damping_ratio
    
    f1.append(freq[0])
    f2.append(freq[1])
    f3.append(freq[2])
    f4.append(freq[3])
    
omega = np.linspace(0,10000,np.ceil(freq_ref/6))
   
fig1, ax1 = plt.subplots()
ax1.set_xlabel('Rotação [Hz]')
ax1.set_ylabel('Fator de Amortecimento')
ax1.plot(omega,f1,'m',label='GDL1')
ax1.plot(omega,f2,'orange',label='GDL2')
ax1.plot(omega,f3,'g',label='GDL3')
ax1.plot(omega,f4,'c',label='GDL4')
# ax1.plot(omega,omega,'b--')
# ax1.plot(omega,2*omega,'r--')
ax1.margins(x=0.0)
ax1.legend(loc='upper right')
# plt.axis('image')
# plt.xticks(range(0,5001,100))
# plt.yticks(range(0,10001,100))
# plt.savefig('fator_amortecimento.png',dpi=600)