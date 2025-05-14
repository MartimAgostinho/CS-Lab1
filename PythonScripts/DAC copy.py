import numpy as np
from sympy import *
import matplotlib.pyplot as plt  

#Sc    = sum Ci
#Scb   = sum Ci*Bi
#Sc_B  = sum CBi
#Scb_B = sum CBi*Bi

def dac_out(vr,C_B0l,C_B,C_ip, C_in,C_Bi, vi, code ):
    #Ver  C_i, C_Bi e code sao arrays
    #Ver tamanho das arrays
    
    numbits = len(code)
    code = np.flip(code)
    print("\nDEBUG "+ str(int(numbits/2)) +"\n")
    #TODO otimizar
    s_LC = sum( C_Bi ) + C_B0l
    print(s_LC)
    s_MC = sum( C_ip ) + sum(C_in)#+ (C_B*s_LC)/(C_B + s_LC) #- 1
    print(s_MC)

    s_LB = 0 
    s_MB_p = 0 
    s_MB_n = 0 

    for i in range( 0, 5):
        s_LB += C_Bi[i]*code[i]
    for i in range( 0, 7):
        s_MB_p += C_ip[i]*code[i+5]
        s_MB_n += C_in[i]*(~code[i+5])
    return vi + \
        (vr*( (s_MB_p +s_MB_n)*(s_LC + C_B) + C_B*s_LB ) )/( s_MC*C_B + s_MC*s_LC + s_LC*C_B)
#exemplo para 5 bits lsb e 7 msb

C_ip   = [1, 1, 2, 4, 8, 16, 32]
C_in   = [0, 1, 2, 4, 8, 16, 32]
C_Bi  = [1, 2, 4, 8, 16]
C_B   = 32/31
C_B0l = 1
vr    = 1
vi    = 0

res_dac = [ 0 for _ in range(4096)]

for n in range (0,4096):
    
    tmp = [ int(i) for i in bin(n)[2:] ]
    if len(tmp) < 12:
        code = np.append( [ 0 for _ in range( 12 - len(tmp)  ) ],tmp  )
    
    else:
        code = np.asarray(tmp)

    #print(tmp)
    res_dac[n] = dac_out(vr,C_B0l,C_B,C_ip, C_in,C_Bi, vi, code )
    print(str(code)+" -> "+str( res_dac[n] ))


x = [ n for n in range(0,4096) ]
plt.step(x, res_dac,where='mid')
plt.plot(x, res_dac, 'o', color='orange')

plt.xlabel("DAC OUT")
plt.ylabel("Code")
plt.title("DAC transfer function")
plt.legend()
plt.grid(True)
plt.show()
