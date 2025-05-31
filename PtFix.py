"""
Resolução de problema em SEP utilizando o Método de Ponto Fixo

- Implementado para sistemas em dois barramentos;

João Vítor Drumond
2025/1
"""

import numpy as np
from math import acos, sin, cos

"""Dados do Sistema"""
Vbase = 127         #Tensão de base para p.u. 
Va = 127 + 0*1j     #Tensão da fonte (barra A)
Zs = 0.05 + 0.1*1j  #Impedância da fonte
Zab = 0.05 + 0.1*1j #Impedância da LT (entre as barras A e B)

"""Dados da Carga"""
Sl = 1e3  #Potência aparente [VA]
fp = 0.95 #Fator de potência

"""Equivalência de fontes"""
Ia = Va/Zs      #Corrente da fonte
Ya = 1/Zs       #Admitância equivalente da fonte
Yab = 1/Zab     #Admitância da LT

th_l = acos(fp) #Ângulo do fator de potência da carga
Sl_d = Sl*cos(th_l) + Sl*1j*sin(th_l) #Potência da carga em função das componentes ativa e reativa

"""Dados da carga"""
Vb = Va                     #Assume-se inicialmente
Zl = (Va**2)/np.conj(Sl_d)  #Impedância da carga
Yl = 1/Zl                   #Admitância da carga

"""Matriz Admitâncias"""
Y_sys = np.array([[(Ya + Yab), -Yab],
                  [-Yab, (Yl+Yab)]])

"""Condição Inicial"""
I_comp = 0 #Corrente de compensação na barra B (de carga)
#Injeções de corrente
I_bus = np.array([[Ia],
                  [I_comp]])
#Tensões iniciais
V_bus = np.matmul(np.linalg.inv(Y_sys),I_bus)
Va = V_bus[0][0]                #Novo valor da tensão na barra A
Vb = V_bus[1][0]                #Novo valor da tensão na carga
I_l = Yl*Vb                     #Corrente linear de carga
Ib = np.conj(Sl_d)/np.conj(Vb)  #Corrente na carga
I_comp = I_l-Ib                 #Novo valor da corrente de compensação

"""Laço"""
err = 1e-7 #Erro projetado
i = 0      #Contagem de iterações
while 1:
    i = i+1

    #Injeções de corrente
    I_bus = np.array([[Ia],
                  [I_comp]])
    
    #Tensões nos barramentos
    V_bus = np.matmul(np.linalg.inv(Y_sys),I_bus)
    Van = V_bus[0][0]
    Vbn = V_bus[1][0]

    #Corrente linear de carga
    I_l = Yl*Vbn

    #Corrente na carga
    Ib = np.conj(Sl_d)/np.conj(Vbn)

    #Corrente de compensação
    I_comp = I_l-Ib

    #Erro
    mat_err = np.array([[abs(Van-Va)/Vbase],
                        [abs(Vbn-Vb)/Vbase]])
    
    if (mat_err[0][0] < err) and (mat_err[1][0] < err):
        break

    Va = Van
    Vb = Vbn

print('******** Saídas ********\n')
print(f'Ia    = ({round(abs(Ia),2)}| {round(np.rad2deg(np.angle(Ia)),5)}°) A')
print(f'Ib    = ({round(abs(Ib),2)}| {round(np.rad2deg(np.angle(Ib)),5)}°) A')
print(f'Icmp  = ({round(abs(I_comp),2)}| {round(np.rad2deg(np.angle(I_comp)),5)})° A')
print(f'Va    = ({round(abs(Van),2)}| {round(np.rad2deg(np.angle(Van)),5)}°) V')
print(f'Vb    = ({round(abs(Vbn),2)}| {round(np.rad2deg(np.angle(Vbn)),5)}°) V\n')
print(f'Erro    = {mat_err}\n')
print(f'Total de iterações: {i}')