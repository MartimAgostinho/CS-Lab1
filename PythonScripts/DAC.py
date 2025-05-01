import numpy as np
from sympy import *
import matplotlib.pyplot as plt  

import datetime

TYPE = np.float64

def num2bin(n,n_bits):
    
    tmp = [ int(i) for i in bin(n)[2:] ]
    if len(tmp) < n_bits:
        code = np.append( [ int(0) for _ in range( int(n_bits) - len(tmp)  ) ],tmp  )
    else:
        code = np.asarray(tmp,dtype=int)
    
    return np.flip(code)

def dac_TF(C_ai, C_bi,C_M0, C_Li,C_B,C_B0l,vi,vr,code_matrix = 0,n_bits = -1):

    C_abi = C_ai + C_bi
    s_LC = sum( C_Li ) + C_B0l 
    s_MC = sum( C_abi) + C_M0
    s_Cb = sum( C_bi )*s_LC

    if n_bits == -1:
        n_bits = len(C_ai) + len(C_Li) + 1

    if( code_matrix == 0 ):
        code_matrix  = [ num2bin(i,n_bits) for i in range(0,2**n_bits) ] 
    
    g  = 1/( C_B*( s_LC + s_MC ) + s_LC*s_MC )
    
    tf = np.empty(shape=(2**n_bits),dtype=TYPE)
    for j in range(0,2**n_bits ):

        code    = code_matrix[j]
        Lsb_qty = len(C_Li)
        s_LB    = 0 
        s_MB    = code[Lsb_qty]*C_M0

        for i in range( 0, Lsb_qty):
            s_LB += C_Li[i]*code[i]

        for i in range( 0, len(C_abi)):
            s_MB += C_abi[i]*code[i+Lsb_qty+1]

        tf[j] = vi + vr*g*( (C_B+s_LC)*s_MB +  s_LB*C_B - s_Cb)
    return tf


def get_Vlsb_real( dac_tf ):
    
    dac_len = len(dac_tf) - 1
    return ( dac_tf[ dac_len  ] - dac_tf[0] ) / dac_len 

def get_INL(dac_tf, Vlsb_real):

    dac_len = len(dac_tf)
    return (dac_tf - np.arange(dac_len) * Vlsb_real - dac_tf[0]) / Vlsb_real

def get_DNL(dac_tf, Vlsb_real):

    return ((dac_tf[1:] - dac_tf[:-1]) / Vlsb_real) - 1


n_bits = 4
C_ai   = np.asarray([ 2**(i-1) for i in range(1,int(n_bits/2) + 1) ])
C_bi   = np.asarray([ 2**(i-1) for i in range(1,int(n_bits/2) + 1) ])
C_Li   = np.asarray([ 2**i for i in range(0,int(n_bits/2)) ])
C_M0   = 1 
C_B    = (2**(len(C_Li)))/( (2**(len(C_Li))) - 1 )
C_B0l  = 1
vr     = 1
vi     = 0


n_bits    = len(C_ai) + len(C_Li) + 1
codes     = [ num2bin(i,n_bits) for i in range(0,2**n_bits) ]# Avoids recalculating all codes for all monte carlo iterations
neg_codes = [ (code - 1)*(-1) for code in codes ]
print(neg_codes)
start = datetime.datetime.now()

Vx_p = dac_TF(
            C_ai,
            C_bi,
            C_M0, 
            C_Li,
            C_B,
            C_B0l,
            vi,
            vr,
            code_matrix = codes,
            n_bits = n_bits)
Vx_n = dac_TF(
            C_ai,
            C_bi,
            C_M0, 
            C_Li,
            C_B,
            C_B0l,
            vi,
            vr,
            code_matrix = neg_codes,
            n_bits = n_bits)

res_dac = Vx_p - Vx_n

x = [ n for n in range(0,2**n_bits) ]
print("dac_Transfer_function time = "+ str( (datetime.datetime.now() - start).total_seconds() ) )

Vlsb_real = get_Vlsb_real(res_dac)
print("VLSB DAC")

inl_res = get_INL(res_dac, Vlsb_real)
print("INL DAC")

dnl_res = get_DNL(res_dac, Vlsb_real)
print("DNL DAC")

print(Vlsb_real)
print(inl_res)
print("\n========================DNL========================\n")
print(dnl_res)
plt.step(x, res_dac,where='mid')
#plt.plot(x, res_dac, 'o', color='orange')
plt.xlabel("DAC OUT")
plt.ylabel("Code")
plt.title("DAC transfer function")
#  plt.legend()
plt.grid(True)
plt.show()

plt.step(x, inl_res,where='mid')
#plt.plot(x, res_dac, 'o', color='orange')
plt.xlabel("INL OUT")
plt.ylabel("Code")
plt.title("INL Values")
plt.grid(True)
plt.show()

plt.step(x[:-1], dnl_res,where='mid')
#plt.plot(x, res_dac, 'o', color='orange')
plt.xlabel("DNL OUT")
plt.ylabel("Code")
plt.title("DNL Values")
plt.grid(True)
plt.show()
