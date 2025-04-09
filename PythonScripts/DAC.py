import numpy as np
from sympy import *

Vy,Vx,Vi,Vr, Cb, Cb0l, Scb_b, Scb, Sc, Sc_b, n, N, m, M = symbols( "V_y V_x V_i V_r C_B C_{B0l} S_{cb_B} S_{cB} S_c S_{c_B} n N m M" )

#print( latex(vy),latex(vx),latex(vi),latex(vr), latex(Cb), latex(Cb0l), latex(Scb_B),latex(Scb),latex(Sc),latex(Sc_b))

#Cbi = IndexedBase('C_B')
#bbi = IndexedBase('b_B')
#Ci  = IndexedBase('C')
#bi  = IndexedBase('b')
#
#Scb_b = summation( Cbi[n]*bbi[n], (n,0,N) )
#Sc_b  = summation( Cbi[n], (n,0,N))
#
#Scb = summation( Ci[m]*bi[m], (m,0,M) )
#Sc  = summation( Ci[m], (m,0,M))

Vy = (  Cb*(Vx-Vi) + (Vr*Scb_b)  )/( Cb + Cb0l + Scb )

f1 = Vi + ( ( Vy*Cb + Vr*Sc_b )/(Sc) ) - Vx

Vx = solve( f1,Vx )

print(latex(cancel(Vx[0])))