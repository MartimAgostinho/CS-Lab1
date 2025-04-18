import numpy as np
from sympy import *
import matplotlib.pyplot as plt  

Vy,Vx,Vi,Vr, Cb, Cb0l, Scb_B, Scb, Sc, Sc_B, n, N, m, M = symbols( "V_y V_x V_i V_r C_B C_{B0l} S_{cb_B} S_{cb} S_c S_{c_B} n N m M" )


#Sc    = sum Ci
#Scb   = sum Ci*Bi
#Sc_B  = sum CBi
#Scb_B = sum CBi*Bi

Vy = (  Cb*(Vi-Vx) + (Vr*Scb_B)  ) / \
     ( Cb0l - Cb + Sc_B )

eq1 = Vi - Vx + ( ( Vy*Cb + Vr*Scb )/(Cb+Sc) )

Vx = cancel( expand( solve( eq1,Vx )[0]))

print(latex( collect( Vx, [Vi, Vr] ) )  )




