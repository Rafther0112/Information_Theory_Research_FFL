#%%
#Importe de librerias
import numpy as np
from tqdm import tqdm
from numba import jit,njit
import pandas as pd
import json
from Logic_Gate_Function import logic_gate_function_7, Hill_Activation, Hill_Represion
#___________________________________________________________________________________________________
#Parametros de simulacion

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
#___________________________________________________________________________________________________
diccionario_global_FFL_I3 = {}
valores_posibles_Hill = [1,2,3, 4]
valores_posibles_Kx = [1, 2,3,4,5,6,7,8, 9, 10]

for Hill in valores_posibles_Hill:
    distribucion_proteina_X = []
    distribucion_proteina_Y = []
    distribucion_proteina_Z = []

    for Kx in tqdm(valores_posibles_Kx):
        
        Mx = Kx/gammamx
        valor_X_estacionario = (Kpx/muX)*Mx

        valor_Y_estacionario = (Kpy/muY)*My
        valor_Z_estacionario = (Kpz/muZ)*Mz

        Kxy  = valor_X_estacionario        #Coeficiente de interaccion proteina X con ARNmY
        Kxz  = valor_X_estacionario         #Coeficiente de interaccion proteina X con ARNmZ
        Kyz  = 2*valor_Y_estacionario         #Coeficiente de interaccion proteina Y con ARNmZ

        Ky = (My*gammamy)*(((valor_X_estacionario**Hill) + (Kxy**Hill))/(valor_X_estacionario**Hill))

        @njit
        def funcion_creacion_ARNmX():
            return Kx

        @njit
        def funcion_creacion_ARNmY(cantidad_X):
            return Ky*((cantidad_X**Hill)/(cantidad_X**Hill + Kxy**Hill))

        @njit
        def funcion_creacion_ARNmZ(cantidad_X, cantidad_Y):

            ARNmZ_interaction_X = Hill_Represion(cantidad_X, Kxy, Ky, Hill)
            ARNmZ_interaction_Y = Hill_Activation(cantidad_Y, Kxy, Ky, Hill)
            K_parameters = [1,1,1,1,1]
            retorno = logic_gate_function_7(ARNmZ_interaction_X, ARNmZ_interaction_Y, K_parameters, Ky)
            #NECESITO REVISAR EL CUARTO PARAMETRO DE ESTA FUNCION PARA SABER COMO HACER QUE CUADRE BIEN LA ESCALA QUE QUEREMOS
            return retorno

        @njit
        def funcion_creacion_X(cantidad_mX):
            return Kpx*cantidad_mX

        @njit
        def funcion_creacion_Y(cantidad_mY):
            return Kpy*cantidad_mY

        @njit
        def funcion_creacion_Z(cantidad_mZ):
            return Kpz*cantidad_mZ

        @njit
        def funcion_degradacion_ARNmX(cantidad_mX):
            return gammamx*cantidad_mX

        @njit
        def funcion_degradacion_ARNmY(cantidad_mY):
            return gammamy*cantidad_mY

        @njit
        def funcion_degradacion_ARNmZ(cantidad_mZ):
            return gammamz*cantidad_mZ

        @njit
        def funcion_degradacion_X(cantidad_X):
            return muX * cantidad_X 

        @njit
        def funcion_degradacion_Y(cantidad_Y):
            return muY * cantidad_Y 

        @njit
        def funcion_degradacion_Z(cantidad_Z):
            return muZ * cantidad_Z 
        
        @njit
        def modelo_constitutivo(cantidad_mX, cantidad_mY, cantidad_mZ, cantidad_X, cantidad_Y,cantidad_Z):

            propensidad_creacion_ARNmX = funcion_creacion_ARNmX()
            propensidad_creacion_ARNmY = funcion_creacion_ARNmY(cantidad_X)
            propensidad_creacion_ARNmZ = funcion_creacion_ARNmZ(cantidad_X, cantidad_Y)

            propensidad_creacion_proteinaX = funcion_creacion_X(cantidad_mX)
            propensidad_creacion_proteinaY = funcion_creacion_Y(cantidad_mY)
            propensidad_creacion_proteinaZ = funcion_creacion_Z(cantidad_mZ)

            propensidad_degradacion_ARNmX = funcion_degradacion_ARNmX(cantidad_mX)
            propensidad_degradacion_ARNmY = funcion_degradacion_ARNmY(cantidad_mY)
            propensidad_degradacion_ARNmZ = funcion_degradacion_ARNmZ(cantidad_mZ)

            propensidad_degradacion_proteinaX = funcion_degradacion_X(cantidad_X)
            propensidad_degradacion_proteinaY = funcion_degradacion_Y(cantidad_Y)
            propensidad_degradacion_proteinaZ = funcion_degradacion_Z(cantidad_Z)

            return propensidad_creacion_ARNmX, propensidad_creacion_ARNmY, propensidad_creacion_ARNmZ, propensidad_creacion_proteinaX, propensidad_creacion_proteinaY, propensidad_creacion_proteinaZ, propensidad_degradacion_ARNmX, propensidad_degradacion_ARNmY, propensidad_degradacion_ARNmZ, propensidad_degradacion_proteinaX, propensidad_degradacion_proteinaY, propensidad_degradacion_proteinaZ
        
        @njit('f8[:](f8[:],f8)')
        def Gillespie(trp0,tmax):

            t,ARNmX, ARNmY, ARNmZ, proteinaX, proteinaY, proteinaZ =trp0 

            while t < tmax:
                s_1, s_2, s_3, s_4, s_5, s_6, s_7, s_8, s_9, s_10, s_11, s_12 = modelo_constitutivo(ARNmX, ARNmY, ARNmZ, proteinaX, proteinaY, proteinaZ)
                S_T = s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10 + s_11 + s_12

                τ = (-1/S_T)*np.log(np.random.rand())
                x = np.random.rand()

                if x <= (s_1)/S_T:
                    ARNmX += 1

                elif x<= (s_1 + s_2)/S_T:
                    ARNmY += 1
                
                elif x <= (s_1 + s_2 + s_3)/S_T :
                    ARNmZ+=1
                
                elif x <= (s_1 + s_2 + s_3 + s_4)/S_T :
                    proteinaX+=1
                
                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5)/S_T :
                    proteinaY+= 1

                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6)/S_T :
                    proteinaZ += 1
                
                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7)/S_T :
                    ARNmX-= 1
                
                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8)/S_T :
                    ARNmY-=1

                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9)/S_T :
                
                    ARNmZ-= 1

                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10)/S_T :

                    proteinaX-=1

                elif x <= (s_1 + s_2 + s_3 + s_4 + s_5 + s_6 + s_7 + s_8 + s_9 + s_10 + s_11)/S_T :
                    proteinaY-=1

                else: 
                    proteinaZ-=1

                t+=τ
            return np.array([t,ARNmX, ARNmY, ARNmZ, proteinaX, proteinaY, proteinaZ]) 
        
        @njit('f8[:,:](f8[:],f8[:])')
        def Estado_celula(X0,tiempos):

            X = np.zeros((len(tiempos),len(X0)))
            X[0] = X0
            
            for i in range(1,len(tiempos)):
                X[i] = Gillespie(X[i-1],tiempos[i])
            
            return X
        
        x0 = np.array([0., 0., 0., 0., 0., 0., 0.])

        num_cel = 10000 #número de células 
        celulas = np.array([Estado_celula(x0,np.arange(0.,700.,2.)) for i in tqdm(range(num_cel))])

        distribuciones_propias_X = celulas[:,0:,4]
        distribuciones_propias_Y = celulas[:,0:,5]
        distribuciones_propias_Z = celulas[:,0:,6]

        distribucion_proteina_X.append(distribuciones_propias_X)
        distribucion_proteina_Y.append(distribuciones_propias_Y)
        distribucion_proteina_Z.append(distribuciones_propias_Z)

#    diccionario_global_FFL_I3[f"Coeficiente_Hill_{Hill}"] = [distribucion_proteina_X, distribucion_proteina_Y, distribucion_proteina_Z]
#    np.save('Simulacion_FFL_I3_AND_final.npy', diccionario_global_FFL_I3)



