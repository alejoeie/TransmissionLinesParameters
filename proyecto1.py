# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 17:17:47 2020

@authors: Esteban, Alejandro
"""


#Proyecto Transmi
#II 2020

import numpy as np 
import pandas as pd


#*===========================================FUNCIONES===================================================*
#Declaramos ciertas variables requeridas:
epsilon = 1.7975103582e+10
def matriz_Coeficientes(m1, m2):
    if len(m1) == len(m2) and len(m1[0]) == len(m2[0]):
        m3 = []
        for i in range(len(m1)):
            m3.append([])
            for j in range(len(m1[0])):
                m3[i].append((epsilon)*np.log(m1[i][j]/m2[i][j]))
        return m3
    else:
        return None

#Para el calculo de la matriz de impedancias

def matriz_Impedancias(mtrx1, mtrx2):
    mxF = []
    for i in range(len(mtrx1)):
        mxF.append([])
        for j in range(len(mtrx1[0])):
            z = 0+1j
            if i == j:
                mxF[i].append(Rkk + Rkk_prima + z*(4*mu*np.pi*f*np.log((Dkk_prima/dkk))))
            else:
                mxF[j].append(Rkk_prima + z*(4*mu*np.pi*f*np.log(MatDkk_prima/dkk)))
        return mxF
    

#==========================================PROCEDIMIENTO=========================================
#Paso1: Calculo del GMR y req de los haces dobles 
#Radio de conductor ACSR Condor 
r = (0.02774/2) #m

#Distancia entre conductores del haz
d = 0.3 #m 

#GMR de los haces dobles
gmr = np.sqrt(r*d)



print('El GMR de los haces dobles es: ', gmr,'m')
#Conversión pies a metros
gmr_hiloguarda = 0.00001*0.3048
print('El GMR de hilos guarda es', gmr_hiloguarda,'m')

#Calculo de req de los haces dobles
#req = np.sqrt()

#Cálculo de Dkk'
rho = 100 #ohms m
f = 60 #Hz 
Dkk_prima = 658.5*np.sqrt(rho/f) #m

print('La distancia entre cada fase y su retorno a tierra es: ', Dkk_prima,'m')

#Definimos distancias entre conductores //d77 y d88 son los gmr de los neutros
d11 = d22 = d33 = d44 = d55 = d66  = gmr #
d77=d88 = gmr_hiloguarda
d21 = d12 = d23 = d32 = d45 = d54 = d56 = d65 = 5#m
d13 = d31 = d46 = d64 = 10 #m
d14 = d41 = d25 = d52 = d36 = d63 = 8.65 # 
d15 = d51 = d42 = d24 = d62 = d26 = d35 = d53 = 9.9911 # m
d16 = d61 = d34 = d43 = 13.22 #m 

#Definimos distancias entre conductores y conductores de retorno 
D11_prima = D22_prima = D33_prima = D44_prima = D55_prima = D66_prima = D77_prima = D88_prima = D78_prima=D87_prima= Dkk_prima #m
D14_prima = D41_prima = D25_prima = D52_prima = D36_prima = D63_prima = 850.164 #m
D15_prima = D42_prima = D53_prima = D26_prima = 855.1637 #m
D43_prima = D16_prima = 860.1635 #m
D24_prima = D51_prima = D35_prima = D62_prima = 845.1643 #m
D34_prima = D61_prima = 840.1645 #m
D12_prima = D45_prima = D23_prima = D56_prima = 855.12 #m
D21_prima = D54_prima = D32_prima = D65_prima = 845.12 #m
D13_prima = D46_prima = 860.12 #m
D31_prima = D64_prima = 840.12 #m
D17_prima=D48_prima=D71_prima=D84_prima=D18_prima=D81_prima=D47_prima=D74_prima=845
D27_prima=D58_prima=D72_prima=D85_prima=D28_prima=D82_prima=D57_prima=d75_prima=840
D37_prima=D68_prima=D73_prima=D86_prima=D38_prima=D67_prima=D83_prima=D76_prima=835

#Se definen las distancias entre los conductores y los neutros
d77=d88=gmr_hiloguarda
estima = round(gmr_hiloguarda,2)
d17=d47=d71=d74=5
d27=d57=d72=d75=10
d37=d67=d73=d76=15
d18=d48=d81=d84=9.6897*f 
d28=d58=d82=d85=13
d38=d68=d83=d86=17.1432*f
d78=d87=7.95*f
#Matriz de distancias entre conductores y neutros Dkk, para k = 1,..,8.  
dkk = np.matrix([[d11, d12, d13, d14, d15, d16, d17, d18],
                 [d21, d22, d23, d24, d25, d26, d27, d28],
                 [d31, d32, d33, d34, d35, d36, d37, d38],
                 [d41, d42, d43, d44, d45, d46, d47, d48],
                 [d51, d52, d53, d54, d55, d56, d57, d58],
                 [d61, d62, d63, d64, d65, d66, d67, d68],
                 [d71, d72, d73, d74, d75, d76, d77, d78],
                 [d81, d82, d83, d84, d85, d86, d87, d88]])
print('La matriz de distancias entre conductores, considerando neutros, (Dkk) en metros es: \n', dkk)

#Matriz de distancias entre conductores y conductores de retorno  Dkk'
MatDkk_prima = np.matrix([[D11_prima, D12_prima, D13_prima, D14_prima, D15_prima, D16_prima, D17_prima,D18_prima],
                 [D21_prima, D22_prima, D23_prima, D24_prima, D25_prima, D26_prima, D27_prima, D28_prima],
                 [D31_prima, D32_prima, D33_prima, D34_prima, D35_prima, D36_prima, D37_prima, D38_prima],
                 [D41_prima, D42_prima, D43_prima, D44_prima, D45_prima, D46_prima, D47_prima, D48_prima],
                 [D51_prima, D52_prima, D53_prima, D54_prima, D55_prima, D56_prima, D57_prima, D58_prima],
                 [D61_prima, D62_prima, D63_prima, D64_prima, D65_prima, D66_prima, D67_prima, D68_prima],
                 [D71_prima, D72_prima, D73_prima, D74_prima, d75_prima, D76_prima, D77_prima, D78_prima],
                 [D81_prima, D82_prima, D83_prima, D84_prima, D85_prima, D86_prima, D78_prima, D87_prima]]) 

print('La matriz de distancias entre los conductores y conductores de retorno (Dkk´) en metros es: \n', MatDkk_prima)

#Resistencia de los conductores y conductores de retorno
Rkk = 8.742e-05 #ohms/metro
Rkk_prima = 9.869e-07*f
print('La resistencia de los conductores en ohms por metro es {} y la resistencia de los \n conductores de retorno en ohms por metro es {}'.format(Rkk, Rkk_prima))

#Calculo de Matriz de impedancias, existen dos casos: Elementos de la diagonal (if) y elementos fuera de la 
#diagonal (else)

columnas = filas = 8 
z_kk = np.zeros([filas, columnas], dtype=np.complex_)

mu = 10e-07
print(mu)
#for i in range(filas):
#    for k in range(columnas):
#        z = 0+1j
        
#        if i == k:
#            z_kk[i,k] = Rkk + Rkk_prima + z*(4*mu*np.pi*f*np.log((Dkk_prima/dkk)))
#        else:
#            z_kk[i,k] = Rkk_prima + z*(4*mu*np.pi*f*np.log(MatDkk_prima/dkk)) 
z_kkred = z_kk.round(decimals=2)
z_kkeq = matriz_Impedancias(dkk, Dkk_prima)
print('La matriz de impedancias es: \n', z_kkeq)

#Nombramos Za
#Matriz de 6x6
Za = np.zeros([6, 6], dtype=np.complex_)
#Extremos los datos de la matriz z_kk
for i in range(6):
    for j in range(6):
        Za[i, j] = z_kk[i, j]
Zared = Za.round(decimals=2)
print('Za:', Zared)

#Zb de tamaño 6x2
Zb = np.zeros([6,2], dtype = np.complex_)
for i in range(6):
    for j in range(2):
        Zb[i, j] = z_kk[i, j]
Zbred = Zb.round(decimals=2)
print('Zb:', Zb)

#Zc matriz de 2x6
Zc = np.zeros([2,6], dtype = np.complex_)
for i in range(2):
    for j in range(6):
        Zc[i, j] = z_kk[i, j]
Zcred = Zc.round(decimals=2)
print('Zc:', Zcred)

#Zd matriz de 2x2

Zd = np.zeros([2,2], dtype = np.complex_)
for i in range(2):
    for j in range(2):
        Zd[i, j] = z_kk[i, j]
Zdred = Zd.round(decimals=2)
print('Zd:', Zdred)
        
#Calculamos Zp = Za - Zb*Zd^-1*Zc
Zd_inv = np.linalg.inv(Zd)
Zp = np.zeros([6,6], dtype = np.complex_)
Zp_equiv=np.matmul(Zb,Zd_inv)
Zp_nueva = Za - np.matmul(Zp_equiv, Zc)

print('Zp \n:', Zp_nueva)
        

#Se calcula la matriz Hkm para obtener la matriz de coeficientes
constante_dielec = 8.8542e-012
h11=h44=26.6 #m
h12=h45=31.6 #m
h13=h46=36.6 #m
h21=h54=31.6 #m
h22=h55=36.6 #m
h23=h56=41.6 #m
h31=h64=36.6 #m
h32=h65=41.6 #m
h33=h66=46.6 #m

h34=h61=37.6 #m
h35=h62=42.49 #m
h36=h63=47.40 #m
h38=h67=h76=h83=52.46 #m
h37=h68=h73=h86=51.8 #m
h24=h51=32.76 #m
h25=h52=37.6 #m
h26=h53=42.49 #m
h28=h57=h75=h82=47.5 #m
h27=h58=h72=h85= 46.8 #m
h16=h43=37.6 #m
h15=h42=32.76 #m
h14=h41=27.97 #m
h18=h47=h74=h81=42.6 #m
h17=h48=h71=h84=42.2 #m 
h88= h77 = 57 #m
h78 = h87 = 57.55 #m

hkm = np.matrix([[h11, h12, h13, h14, h15, h16, h17, h18],
                [h21, h22, h23, h24, h25, h26, h27, h28],
                [h31, h32, h33, h34, h35, h36, h37, h38],
                [h41, h42, h43, h44, h45, h46, h47, h48],
                [h51, h52, h53, h54, h55, h56, h57, h58],
                [h61, h62, h63, h64, h65, h66, h67, h68],
                [h71, h72, h73, h74, h75, h76, h77, h78],
                [h81, h82, h83, h84, h85, h86, h87, h88]])
hkmred=hkm.round(decimals=2)
#print("H_km:\n",hkm)

#Se definen las distancias de los conductores aereos
a11=a22=a33=a44=a55=a66=gmr
a12=a21=a23=a32=a65=a56=a54=a45=5
a13=a31=a46=a64=10
a14=a41=a52=a25=a63=a36=8.65
a35=a53=a26=a62=a24=a42=a15=a51=9.9911
a37=a73=a86=a68=5.2117
a27=a72=a58=a85=10.2
a17=a71=a48=a84=15.2
a38=a83=a67=a76=9.794
a28=a82=a57=a75=13.15
a18=a47=a81=a74=17.32
a77=a88=gmr
a16=a61=a34=a43=13.22
a78=a87=7.95 
A_kk = np.matrix([[a11, a12, a13, a14, a15, a16, a17, a18], 
                  [a21, a22, a23, a24, a25, a26, a27, a28],
                  [a31, a32, a33, a34, a35, a36, a37, a38],
                  [a41, a42, a43, a44, a45, a46, a47, a48],
                  [a51, a52, a53, a54, a55, a56, a57, a58],
                  [a61, a62, a63, a64, a65, a66, a67, a68],
                  [a71, a72, a73, a74, a75, a76, a77, a78],
                  [a81, a82, a83, a84, a85, a86, a87, a77]],dtype = float, copy=True)

A_kkred=A_kk.round(decimals = 2)
print(A_kkred)

#Se llama la funcion que define la matriz de capacitancias

P_km = matriz_Coeficientes(hkmred, A_kkred)
if P_km == None:
    print("Error de dimension\n")
else:
    for fila in P_km:
        print("[", end=" ")
        for elemento in fila:
            print(elemento, end= " ")
        print("]")
    


#Subdividiendo la matriz en Pa, Pb, Pc y Pd
#Pa es de tamaño 6x6
P_a = np.matrix([[a11, a12, a13, a14, a15, a16],
                [a21, a22, a23, a24, a25, a26],
                [a31, a32, a33, a34, a35, a36],
                [a41, a42, a43, a44, a45, a46],
                [a51, a52, a53, a54, a55, a56],
                [a61, a62, a63, a64, a65, a66]],dtype = float, copy= True)
print("PA\n ", P_a)




#Pb con tamaño 6x2
P_b = np.matrix([[a17, a18],
                 [a27, a28],
                 [a37,a38],
                 [a47,a48],
                 [a57, a58],
                 [a67, a68]], dtype = float, copy = True)

#Pc con dimensión 2x6
P_c = np.matrix([[a71, a72, a73, a74, a75, a76],[a81, a82, a83, a84, a85, a86]], dtype = float)

#Pd con dimension 2x2
P_d = np.matrix([[a77, a78],[a87, a88]], dtype = float, copy = True)

#Calculo de matriz de capacitancias Cp = (Pa - PbPd^-1Pc^-1)
Pc_equiv = P_c.I
Pd_equiv = P_d.I
#resultado = np.matmul(P_b,Pd_equiv)
#Cp =P_a - np.matmul(resultado,Pc_equiv)
#print("Cp:\n", Cp)