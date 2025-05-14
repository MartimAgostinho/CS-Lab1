import numpy as np
from sympy import *
import matplotlib.pyplot as plt  

Vy,Vx,Vxp1,Vxp2,Vi,Vr, Cb, Cb0l,dCLi, dCMi, s_LC, s_MC, s_LB, s_Cb, s_MB, n, N, m, M = symbols( 'V_y V_x V_{x}^{\phi_n} V_{x}^{\phi_{n+1}} V_i V_r C_B C_{B0l} {\Delta}C_{Li} {\Delta}C_{Mi} \sigma_{LC} \sigma_{MC} \sigma_{LB} \sigma_{Cb} \sigma_{MB} n N m M' )


dVy = ( Cb*(Vxp2-Vxp1) + Vr*(dCLi) )/(s_LC+Cb)
eq = Vxp1 - Vxp2 + ( Vr*( dCMi/(s_MC+Cb) ) ) + dVy * (Cb/(s_MC+Cb))

Vx = cancel( expand( solve( eq,Vxp2 )[0]))

print(latex( collect( Vx, [Vxp1, Vr] ) )  )

#s_LC = sum( C_Bi ) + C_B0l - C_B 
#s_MC = sum( C_i )  + C_B   #- 1
#s_LB = sum( C_Bi*bi_LSB )
#s_MB = sum( C_i*bi_MSB ) 

'''
Vy = (Cb*Vx +Vr*s_LB) /\
     (Cb + s_LC)

eq1 =  Vx*( Cb + s_MC ) - Vr*s_MB - Vy*Cb \
     - Vi*(s_MC + (Cb*s_LC)/(Cb + s_LC) ) + Vr*s_Cb

Vx = cancel( expand( solve( eq1,Vx )[0]))
\frac{V_{r} \left(- C_{B} {\Delta}C_{Li} - \sigma_{LC} {\Delta}C_{Mi}\right) + V_{x}^{\phi_n} \left(C_{B}^{2} - \sigma_{LC} \sigma_{MC}\right)}{C_{B}^{2} - \sigma_{LC} \sigma_{MC}}


''' 