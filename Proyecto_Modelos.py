# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 23:05:40 2020

@author: Valentina Hoyos
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

####### Parámetros ##############
# Paso de host intermedio a humanos
b1=1
b2=0.1

# Vida de Taenias adultas
mu1=1/5 # Una taenia adulta no vive más de 5 años en el humano
mu2=1/5 # Para esto ver el paper "Peru.pdf" página 6
mu3=1/5

# Competencia en humanos
d1=2
d2=5
d3=0.4
d4=1
d5=2
d6=2

# Paso de huevos (feces) a host intermedio
alpha1=0.5
alpha2=0.5
alpha3=1

# Competencia en cerdo
gamma1=0.5
gamma2=1


tsim=np.linspace(0,100)

def model(z,t):
    Hs=b1*z[3]-mu1*z[0]-d1*z[0]*z[1]-d2*z[0]*z[2]
    Ha=b1*z[4]-mu2*z[1]-d3*z[0]*z[1]-d4*z[1]*z[2]
    Hb=b2*z[5]-mu3*z[2]-d5*z[2]*z[0]-d6*z[1]*z[2]
    Cs=alpha1*z[6]-b1*z[3]-gamma1*z[3]*z[4]
    Ca=alpha2*z[7]-b1*z[4]-gamma2*z[3]*z[4]
    Vb=alpha3*z[8]-b2*z[5]
    Es=mu1*z[0]-alpha1*z[6]
    Ea=mu2*z[1]-alpha2*z[7]
    Eb=mu3*z[2]-alpha3*z[8]
    
    return [Hs, Ha,Hb,Cs,Ca,Vb,Es,Ea,Eb]
# HS(0), HA(1), HB(2),Cs(3),Ca(4), Vb(5), Es(6), Ea(7), Eb(8)
# Condiciones iniciales
z=[2,1,1,1,1,1,2,1,1]
sol=odeint(model,z,tsim)

Hs=sol[:,0]
Ha=sol[:,1]
Hb=sol[:,2]
Cs=sol[:,3]
Ca=sol[:,4]
Vb=sol[:,5]
Es=sol[:,6]
Ea=sol[:,7]
Eb=sol[:,8]

plt.plot(tsim,Hs,color='red',label='Hs')
plt.legend(loc='upper right')
plt.plot(tsim,Ha,color='green',label='Ha')
plt.legend(loc='upper right')
plt.plot(tsim,Hb,color='blue',label='Ha')
plt.legend(loc='upper right')
plt.xlabel('Tiempo')
plt.ylabel('Individuos')
plt.title('Modelo de taeniasis en el humano')
plt.show()

plt.plot(tsim,Cs,color='blue',label='Cs')
plt.legend(loc='upper right')
plt.plot(tsim,Ca,color='red',label='Ca')
plt.legend(loc='upper right')
plt.xlabel('Tiempo')
plt.ylabel('Individuos')
plt.title('Modelo de taeniasis en el cerdo')
plt.show()

plt.plot(tsim,Vb,color='blue',label='Vb')
plt.legend(loc='upper right')
plt.xlabel('Tiempo')
plt.ylabel('Individuos')
plt.title('Modelo de taeniasis en la vaca')
plt.show()

plt.plot(tsim,Es,color='blue',label='Es')
plt.legend(loc='upper right')
plt.plot(tsim,Ea,color='red',label='Ea')
plt.legend(loc='upper right')
plt.plot(tsim,Eb,color='green',label='Eb')
plt.legend(loc='upper right')
plt.xlabel('Tiempo')
plt.ylabel('Individuos')
plt.title('Modelo de taeniasis en el ambiente')
plt.show()
