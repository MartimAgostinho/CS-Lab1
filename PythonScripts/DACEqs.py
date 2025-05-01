import numpy as np
from sympy import *
import matplotlib.pyplot as plt  

Vy,Vx,Vi,Vr, Cb, Cb0l, s_LC, s_MC, s_LB, s_Cb, s_MB, n, N, m, M = symbols( 'V_y V_x V_i V_r C_B C_{B0l} \sigma_{LC} \sigma_{MC} \sigma_{LB} \sigma_{Cb} \sigma_{MB} n N m M' )

#s_LC = sum( C_Bi ) + C_B0l - C_B 
#s_MC = sum( C_i )  + C_B   #- 1
#s_LB = sum( C_Bi*bi_LSB )
#s_MB = sum( C_i*bi_MSB ) 

Vy = (Cb*Vx +Vr*s_LB) /\
     (Cb + s_LC)

eq1 =  Vx*( Cb + s_MC ) - Vr*s_MB - Vy*Cb \
     - Vi*(s_MC + (Cb*s_LC)/(Cb + s_LC) ) + Vr*s_Cb

Vx = cancel( expand( solve( eq1,Vx )[0]))

print(latex( collect( Vx, [Vi, Vr] ) )  )




