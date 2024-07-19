#%%
#Importe de librerias

import numpy as np
from tqdm import tqdm

from numba import jit,njit
import pandas as pd

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

diccionario_global_FFL_C2 = {}
#valores_posibles_Hill = [1,2,3, 4]
#valores_posibles_Kx = [1, 2,3,4,5,6,7,8, 9, 10]

valores_posibles_Hill = [4]
valores_posibles_Kx = [10]

for Hill in valores_posibles_Hill:
    distribucion_proteina_X = []
    distribucion_proteina_Y = []
    distribucion_proteina_Z = []

    for Kx in tqdm(valores_posibles_Kx):

        Mx = Kx/gammamx
        valor_X_estacionario = (Kpx/muX)*Mx

        Kxy  = valor_X_estacionario/2        #Coeficiente de interaccion proteina X con ARNmY
        Kxz  = valor_X_estacionario/2   

        Ky = 3
        Kz = 25

        valor_Y_estacionario = Ky*(Kpy/(muY*gammamy))*((Kxy**Hill)/(valor_X_estacionario**Hill + Kxy**Hill))
        
        Kyz  = valor_Y_estacionario/2

        @njit
        def funcion_creacion_ARNmX():
            return Kx

        @njit
        def funcion_creacion_ARNmY(cantidad_X):
            return Ky*((Kxy**Hill)/(cantidad_X**Hill + Kxy**Hill))

        @njit
        def funcion_creacion_ARNmZ(cantidad_X, cantidad_Y):
            creacion_ARNmZ = Kz*((cantidad_X**Hill)/(cantidad_X**Hill + Kxz**Hill))*((Kyz**Hill)/(cantidad_Y**Hill + Kyz**Hill))
            return creacion_ARNmZ

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

        num_cel = 1000 #número de células 
        celulas = np.array([Estado_celula(x0,np.arange(0.,700.,2.)) for i in tqdm(range(num_cel))])

        distribuciones_propias_X = celulas[:,0:,4]
        distribuciones_propias_Y = celulas[:,0:,5]
        distribuciones_propias_Z = celulas[:,0:,6]

        distribucion_proteina_X.append(distribuciones_propias_X)
        distribucion_proteina_Y.append(distribuciones_propias_Y)
        distribucion_proteina_Z.append(distribuciones_propias_Z)

#    diccionario_global_FFL_C2[f"Coeficiente_Hill_{Hill}"] = [distribucion_proteina_X, distribucion_proteina_Y, distribucion_proteina_Z]
#    np.save('Simulacion_FFL_C2_AND_final.npy', diccionario_global_FFL_C2)
#%% 
print(valor_Y_estacionario)        
#%%
celulas_promedio = np.mean(celulas, axis=0)
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # 1 fila, 3 columnas

axs[0].plot(celulas_promedio[:,4])
axs[0].set_title('Protein X')
axs[1].plot(celulas_promedio[:,5])
axs[1].set_title('Protein Y')
axs[2].plot(celulas_promedio[:,6])
axs[2].set_title('Protein Z')
fig.suptitle('Logic Gate 2 C2', fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("Logic_Gate_2_Coherent_2.jpg", dpi = 500)
# %%

Distribucion_X = celulas[:,:,-3]
DIstribucion_Z = celulas[:,:,-1]
#%%
Informacion = np.zeros((len(Distribucion_X), len(Distribucion_X)))
for tiempo_i in np.arange(0,len(Distribucion_X)):
    for tiempo_j in np.arange(0, len(Distribucion_X)):
        data_I1 = {'X': Distribucion_X[:,tiempo_i],
                'Y': DIstribucion_Z[:,tiempo_j]}
        Cov_matrix_I1 = np.array(pd.DataFrame.cov(pd.DataFrame(data_I1)))
        Informacion[tiempo_i][tiempo_j] = (1/2)*np.log2((Cov_matrix_I1[0][0]* Cov_matrix_I1[1][1])/(Cov_matrix_I1[0][0]* Cov_matrix_I1[1][1] - (Cov_matrix_I1[0][1])**2))
#%%
plt.imshow(Informacion)
# %%
np.nanmax(Informacion)
#%%