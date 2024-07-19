#%%
import numpy as np
from tqdm import tqdm
from numba import jit,njit
import pandas as pd
import json
import matplotlib.pyplot as plt
#%%

def Hill_Activation(cantidad, sensitivity, expresion_level, hill):
    return expresion_level*((cantidad**hill)/(cantidad**hill + sensitivity**hill))

def Hill_Represion(cantidad, sensitivity, expresion_level, hill):
    return expresion_level*((sensitivity**hill)/(cantidad**hill + sensitivity**hill))
#%%
Kpx = 200            #Tasa de creacion de proteina X
Kpy = 300            #Tasa de creacion de proteina Y

gammamx = 1/5        #Tasa de degradacion de ARNmX
gammamy = 1/7        #Tasa de degradacion de ARNmY

muX     =1/20            #Tasa de degradacion de proteina X
muY     =1/40            #Tasa de degradacion de proteina Y



Kx = 4
Ky = 10
Hill = 2

Mx = Kx/gammamx
valor_X_estacionario = (Kpx/muX)*Mx

Kxy  = valor_X_estacionario/2       #Coeficiente de interaccion proteina X con ARNmY
Kxz  = valor_X_estacionario         #Coeficiente de interaccion proteina X con ARNmZ

#My = (Ky/gammamy)*((Kxy**Hill)/(valor_X_estacionario**Hill + Kxy**Hill))
#valor_Y_estacionario = (Kpy/muY)*My


#Dinamica ARNmx
@njit
def funcion_creacion_ARNmX():
    return Kx

@njit
def funcion_degradacion_ARNmX(cantidad_mX):
    return gammamx*cantidad_mX

#Dinamica X
@njit
def funcion_creacion_X(cantidad_mX):
    return Kpx*cantidad_mX

@njit
def funcion_degradacion_X(cantidad_X):
    return muX * cantidad_X 


#Dinamica ARNmy
@njit
def funcion_creacion_ARNmY(cantidad_X):
    return Ky*((Kxy**Hill)/(cantidad_X**Hill + Kxy**Hill))

@njit
def funcion_degradacion_ARNmY(cantidad_mY):
    return gammamy*cantidad_mY


#Dinamica Y
@njit
def funcion_creacion_Y(cantidad_mY):
    return Kpy*cantidad_mY

@njit
def funcion_degradacion_Y(cantidad_Y):
    return muY * cantidad_Y 

@njit
def modelo_constitutivo(cantidad_mX, cantidad_mY, cantidad_X, cantidad_Y):

    propensidad_creacion_ARNmX = funcion_creacion_ARNmX()
    propensidad_creacion_ARNmY = funcion_creacion_ARNmY(cantidad_X)

    propensidad_creacion_proteinaX = funcion_creacion_X(cantidad_mX)
    propensidad_creacion_proteinaY = funcion_creacion_Y(cantidad_mY)

    propensidad_degradacion_ARNmX = funcion_degradacion_ARNmX(cantidad_mX)
    propensidad_degradacion_ARNmY = funcion_degradacion_ARNmY(cantidad_mY)

    propensidad_degradacion_proteinaX = funcion_degradacion_X(cantidad_X)
    propensidad_degradacion_proteinaY = funcion_degradacion_Y(cantidad_Y)

    return propensidad_creacion_ARNmX, propensidad_creacion_ARNmY, propensidad_creacion_proteinaX, propensidad_creacion_proteinaY, propensidad_degradacion_ARNmX, propensidad_degradacion_ARNmY, propensidad_degradacion_proteinaX, propensidad_degradacion_proteinaY

@njit('f8[:](f8[:],f8)')
def Gillespie(trp0,tmax):
    """
    Esta funcion se emplea solamente para hacer la evolución de un paso individual en la celula. Evoluciona no un paso temporal, 
    pero si temporalmente la cantidad de veces que pueda evolucionar antes del tmax en una corrida
    """
    
    t,ARNmX, ARNmY, proteinaX, proteinaY =trp0 

    while t < tmax:
        s_1, s_2, s_3, s_4, s_5, s_6, s_7, s_8 = modelo_constitutivo(ARNmX, ARNmY, proteinaX, proteinaY)
        S_T = s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 

        τ = (-1/S_T)*np.log(np.random.rand())
        x = np.random.rand()

        if x <= (s_1)/S_T:
            ARNmX += 1

        elif x<= (s_1 + s_2)/S_T:
            ARNmY += 1
        
        elif x <= (s_1 + s_2 + s_3)/S_T :
            proteinaX+=1
        
        elif x <= (s_1 + s_2 + s_3 + s_4)/S_T :
            proteinaY+= 1
        
        elif x <= (s_1 + s_2 + s_3 + s_4 + s_5)/S_T :
            ARNmX-= 1
        
        elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6)/S_T :
            ARNmY-=1

        elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7)/S_T :

            proteinaX-=1

        elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8)/S_T :
            proteinaY-=1

        t+=τ
    return np.array([t,ARNmX, ARNmY, proteinaX, proteinaY]) 

@njit('f8[:,:](f8[:],f8[:])')
def Estado_celula(X0,tiempos):

    X = np.zeros((len(tiempos),len(X0)))
    X[0] = X0
    
    for i in range(1,len(tiempos)):
        X[i] = Gillespie(X[i-1],tiempos[i])
    
    return X

x0 = np.array([0., 0., 0., 0., 0.])

num_cel = 100 #número de células 
celulas = np.array([Estado_celula(x0,np.arange(0.,700.,2.)) for i in tqdm(range(num_cel))])
# %%
celulas_promedio = np.mean(celulas, axis = 0)
# %%
import matplotlib.pyplot as plt
plt.plot(celulas_promedio[:,1], color = "green", label = "X constitutivo")
plt.plot(celulas_promedio[:,2], color = "red", label = "Y Activacion")
plt.legend()
# %%
plt.plot(celulas_promedio[:,2])
# %%
print(My)
# %%
celulas_promedio[:,2][-1]
# %%
