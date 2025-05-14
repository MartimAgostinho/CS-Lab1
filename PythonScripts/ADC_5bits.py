import numpy as np
from sympy import *
import matplotlib.pyplot as plt 

numbits=5
res=[0, 0, 0, 0, 0]
vr=1
vcn = vr/2


C_ip   = [0, 0, 1, 1]
C_in   = [0, 0, 0, 1]
C_Bi  = [1, 2, 0, 0]
C_B   = 4/3
C_B0l = 1
res_final= np.array([ 0 for _ in range(10000)])

s_LC = sum( C_Bi ) + C_B0l
s_MC = sum( C_ip ) + sum(C_in)
i=0
#for vin in range(-5000, 5000, 1):
vin = 0.3#vin/10000 
print(vin)
vin_p= vin/2 + vcn
vin_n= - vin/2 + vcn
if(vin_p > vin_n):
    res[numbits-1] = 1
else:
    res[numbits-1] = 0
vx_p = vin_p
vx_n = vin_n
print(vx_p)
print(vx_n)
for n in range(numbits -2 ,-1, -1):
    print(res)
    negado = np.logical_not(res).astype(int)
    s_LB_n = C_Bi[n]*res[n+1]
    s_LB_p = C_Bi[n]*negado[n+1]
    if(n==2):
        s_MB_n = C_ip[n]*res[n+1]
        s_MB_p = C_ip[n]*negado[n+1]
    else:
        if(res[n+1] == 1):
            s_MB_p = -C_in[n]*res[n+1]
            s_MB_n = C_in[n]*res[n+1]
        if(res[n+1] == 0):
            s_MB_p = C_in[n]*negado[n+1]
            s_MB_n = -C_in[n]*negado[n+1]
    vx_n = vx_n + vr/2*( s_MB_n*(s_LC + C_B) + C_B*s_LB_n ) /( s_MC*C_B + s_MC*s_LC + s_LC*C_B)
    vx_p = vx_p + vr/2*( s_MB_p*(s_LC + C_B) + C_B*s_LB_p ) /( s_MC*C_B + s_MC*s_LC + s_LC*C_B)
    print(vx_p)
    print(vx_n)
    if(vx_p > vx_n):
        res[n] = 1
    else:
        res[n] = 0

print("res")
print(np.flip(res))
res = np.array(np.flip(res))
decimal = res.dot(2 ** np.arange(res.size - 1, -1, -1))
res_final[i]=int(decimal)
print(decimal)
i+=1

x = (np.array([ n for n in range(-5000,5000) ])/10000)
plt.step(x, res_final,where='mid')
#plt.plot(x, res_final, 'o', color='orange')

plt.xlabel("Vin")
plt.ylabel("Code Out")
plt.title("ADC transfer function")
plt.legend()
plt.grid(True)
plt.show()