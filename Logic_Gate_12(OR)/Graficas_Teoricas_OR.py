"""
Este codigo esta orientado a hacer las curvas teoricas de costo metabolico
de funcionamiento de los FFL bajo logica OR
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd
from tqdm import tqdm
#%%
Kpx = 200            #Tasa de creacion de proteina X
Kpy = 300            #Tasa de creacion de proteina Y
Kpz = 100            #Tasa de creacion de proteina Z

gammamx = 1/5        #Tasa de degradacion de ARNmX
gammamy = 1/7        #Tasa de degradacion de ARNmY
gammamz = 1/10        #Tasa de degradacion de ARNmZ

muX     =1/20            #Tasa de degradacion de proteina X
muY     =1/40            #Tasa de degradacion de proteina Y
muZ     =1/30            #Tasa de degradacion de proteina Z

My = 10
Mz = 25
valores_posibles_Kx = [1,2,3,4,5,6,7,8,9,10]

valor_Y_estacionario = (Kpy/muY)*My
valor_Z_estacionario = (Kpz/muZ)*Mz

# %%
Hill = 0
Valor = []
for Kx in valores_posibles_Kx:
    
    Mx = Kx/gammamx
    
    valor_X_estacionario = (Kpx/muX)*Mx


    Kxy  = valor_X_estacionario        #Coeficiente de interaccion proteina X con ARNmY
    Kxz  = valor_X_estacionario         #Coeficiente de interaccion proteina X con ARNmZ
    Kyz  = 2*valor_Y_estacionario    

    Ky = (My*gammamy)*(((valor_X_estacionario**Hill) + (Kxy**Hill))/(valor_X_estacionario**Hill)) 
    print(Ky)
    #print(Ky)
    resultado = ((Kpy*Ky)/(gammamy*muY))*((valor_X_estacionario**Hill)/(valor_X_estacionario**Hill + Kxy**Hill))
    Valor.append(resultado)
#Valor
# %%
print(Ky)
# %%
(Kpy*Ky)/(gammamy*muY)

# %%
Hill = 3
(valor_X_estacionario**Hill)/(valor_X_estacionario**Hill + Kxy**Hill)
# %%
