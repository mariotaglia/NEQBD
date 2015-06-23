subroutine switchpot(k)
use system
use params
use mainm
use externalforce
implicit none

integer k
real alpha, pH, f
integer JJJ
real*8 temp

!!!!! switch external force
     if(eftype.eq.1) then ! circular motion
     ex = (abs(float(mod((k+eperiod),2*eperiod))/float(eperiod)-1.0)-0.5)*2.0 ! from 0 to 1 V shape
     ey = (abs(float(mod((k+eperiod+eperiod/2),2*eperiod))/float(eperiod)-1.0)-0.5)*2.0 ! from 0 to 1 V shape
     else if (eftype.eq.2) then ! vibration only
     ex = (abs(float(mod((k+eperiod),2*eperiod))/float(eperiod)-1.0)-0.5)*2.0 ! from 0 to 1 V shape
     ey = 0
     else if (eftype.eq.3) then ! better circular
     ex = -sin(float(k)/float(2*eperiod)*2*3.14159)
     ey = cos(float(k)/float(2*eperiod)*2*3.14159)
     endif 

!!!!!!!!! switch potential
if (switchtype.eq.0) then ! switchtype = 0, sudden jump
     if(mod(k, period).eq.0) then
     if(mod((k/period),2).eq.1) then
     rint(1,1) = r1a
     eint(1,1) = e1a
     zint(1) = z1a
     rint(2,2) = r2a
     eint(2,2) = e2a
     zint(2) = z2a
     switch = 0
     else
     rint(1,1) = r1b
     eint(1,1) = e1b
     zint(1) = z1b
     rint(2,2) = r2b
     eint(2,2) = e2b
     zint(2) = z2b
     switch = 1
     endif
     rint(1,2) = (rint(1,1)+rint(2,2))/2
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
     endif
else if (switchtype.eq.1) then ! switchtype = 1, sine function
     alpha = float(k+period)/float(period)*pi
     rint(1,1) = (r1a-r1b)*(cos(alpha)+1.0)/2.0 + r1b
     eint(1,1) = (e1a-e1b)*(cos(alpha)+1.0)/2.0 + e1b
     zint(1) = (z1a-z1b)*(cos(alpha)+1.0)/2.0 + z1b
     rint(2,2) = (r2a-r2b)*(cos(alpha)+1.0)/2.0 + r2b
     eint(2,2) = (e2a-e2b)*(cos(alpha)+1.0)/2.0 + e2b
     zint(2) = (z2a-z2b)*(cos(alpha)+1.0)/2.0 + z2b

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
else if (switchtype.eq.2) then ! switchtype = 2, triangular function, linear
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     rint(1,1) = r1a - (r1a-r1b)*alpha
     eint(1,1) = e1a - (e1a-e1b)*alpha
     zint(1) = z1a - (z1a-z1b)*alpha
     rint(2,2) = r2a - (r2a-r2b)*alpha
     eint(2,2) = e2a - (e2a-e2b)*alpha
     zint(2) = z2a - (z2a-z2b)*alpha

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
else if (switchtype.eq.3) then ! switchtype = 3, triangular function, area
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     rint(1,1) = (r1a**2 - (r1a**2-r1b**2)*alpha)**(0.5)
     eint(1,1) = (e1a**2 - (e1a**2-e1b**2)*alpha)**(0.5)
     zint(1) = (z1a**2 - (z1a**2-z1b**2)*alpha)**(0.5)
     rint(2,2) = (r2a**2 - (r2a**2-r2b**2)*alpha)**(0.5)
     eint(2,2) = (e2a**2 - (e2a**2-e2b**2)*alpha)**(0.5)
     zint(2) = (z2a**2 - (z2a**2-z2b**2)*alpha)**(0.5)

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
else if (switchtype.eq.4) then ! switchtype = 4, triangular function, volume 
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     rint(1,1) = (r1a**3 - (r1a**3-r1b**3)*alpha)**(1.0/3.0)
     eint(1,1) = (e1a**3 - (e1a**3-e1b**3)*alpha)**(1.0/3.0)
     zint(1) = (z1a**3 - (z1a**3-z1b**3)*alpha)**(1.0/3.0)
     rint(2,2) = (r2a**3 - (r2a**3-r2b**3)*alpha)**(1.0/3.0)
     eint(2,2) = (e2a**3 - (e2a**3-e2b**3)*alpha)**(1.0/3.0)
     zint(2) = (z2a**3 - (z2a**3-z2b**3)*alpha)**(1.0/3.0)

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
else if (switchtype.eq.5) then ! switchtype = 5, pH-type function
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=alpha*4.0+3.0 
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
else if (switchtype.eq.6) then ! switchtype = 6, pH-type function blunt tip, see Fig 3-ii PNAS 2014
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=1.0/(1.0+exp((-alpha+0.5)*20.0))*4.0 + 3.0
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
else if (switchtype.eq.7) then ! switchtype = 7, pH-type function sharp tip, see Fig 3-iii PNAS 2014
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=(exp(10.0*(alpha-0.5))-exp(-10.0*(alpha-0.5)))/74.5 + 5.0
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)

else if (switchtype.eq.17) then ! switchtype = 7, pH-type function sharp tip, see Fig 3-iii PNAS 2014
     temp=20.0
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=(exp(temp*(alpha-0.5))-exp(-temp*(alpha-0.5)))*2.0/(exp(temp*0.5)-exp(-temp*0.5)) + 5.0
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)


else if (switchtype.eq.27) then ! switchtype = 7, pH-type function sharp tip, see Fig 3-iii PNAS 2014
     temp=40.0
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=(exp(temp*(alpha-0.5))-exp(-temp*(alpha-0.5)))*2.0/(exp(temp*0.5)-exp(-temp*0.5)) + 5.0
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)

else if (switchtype.eq.37) then ! switchtype = 7, pH-type function sharp tip, see Fig 3-iii PNAS 2014
     temp=80.0
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=(exp(temp*(alpha-0.5))-exp(-temp*(alpha-0.5)))*2.0/(exp(temp*0.5)-exp(-temp*0.5)) + 5.0
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)


else if (switchtype.eq.8) then ! switchtype = 8, LJ blunt
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     alpha=1.0/(1.0+exp((-alpha+0.5)*20.0))
     rint(1,1) = r1a - (r1a-r1b)*alpha
     eint(1,1) = e1a - (e1a-e1b)*alpha
     zint(1) = z1a - (z1a-z1b)*alpha
     rint(2,2) = r2a - (r2a-r2b)*alpha
     eint(2,2) = e2a - (e2a-e2b)*alpha
     zint(2) = z2a - (z2a-z2b)*alpha

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)

else if (switchtype.eq.9) then ! switchtype = 9, LJ, sharp
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     alpha=(exp(10.0*(alpha-0.5))-exp(-10.0*(alpha-0.5)))/74.5/4.0 + 0.5
     rint(1,1) = r1a - (r1a-r1b)*alpha
     eint(1,1) = e1a - (e1a-e1b)*alpha
     zint(1) = z1a - (z1a-z1b)*alpha
     rint(2,2) = r2a - (r2a-r2b)*alpha
     eint(2,2) = e2a - (e2a-e2b)*alpha
     zint(2) = z2a - (z2a-z2b)*alpha

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)

else if (switchtype.eq.10) then ! switchtype = 8, LJ blunt
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     alpha=1.0/(1.0+exp((-alpha+0.2)*25.0))
     rint(1,1) = r1a - (r1a-r1b)*alpha
     eint(1,1) = e1a - (e1a-e1b)*alpha
     zint(1) = z1a - (z1a-z1b)*alpha
     rint(2,2) = r2a - (r2a-r2b)*alpha
     eint(2,2) = e2a - (e2a-e2b)*alpha
     zint(2) = z2a - (z2a-z2b)*alpha

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)

else if (switchtype.eq.20) then ! switchtype = 2, triangular function with maximum at tau/N

     JJJ = 5
     temp = float(mod(k,2*period))/2.0/float(period)*float(JJJ)

     if(temp.le.1.0) then
          alpha = float(mod(k,period*2/JJJ))/float(period)/2.0*float(JJJ)
     else
          alpha = 1.0 - (float(mod(k-period*2/JJJ,2*period))/float(period)/2.0/(1.0-1.0/float(JJJ)))
     endif

     pH=alpha*4.0+3.0
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a
     eint(2,2) = e2a
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
 
 else if (switchtype.eq.30) then ! switchtype = 30, pH-type function up to 5

     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=alpha*2.0+3.0 
     f=1.0/(1.0 + 10**(5.0-pH))
     rint(1,1) = r1a 
     eint(1,1) = e1a
     zint(1) = z1a - (z1a-z1b)*f
     rint(2,2) = r2a 
     eint(2,2) = e2a 
     zint(2) = z2a - (z2a-z2b)*f
     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)

else
     alpha = float(iniMCstep+period)/float(period)*pi
     rint(1,1) = (r1a-r1b)*(cos(alpha)+1.0)/2.0 + r1b
     eint(1,1) = (e1a-e1b)*(cos(alpha)+1.0)/2.0 + e1b
     zint(1) = (z1a-z1b)*(cos(alpha)+1.0)/2.0 + z1b
     rint(2,2) = (r2a-r2b)*(cos(alpha)+1.0)/2.0 + r2b
     eint(2,2) = (e2a-e2b)*(cos(alpha)+1.0)/2.0 + e2b
     zint(2) = (z2a-z2b)*(cos(alpha)+1.0)/2.0 + z2b

     rint(1,2) = (rint(1,1)+rint(2,2))/2.0
     eint(1,2) = (eint(1,1)*eint(2,2))**(0.5)
     rint(2,1) = (rint(1,1)+rint(2,2))/2.0
     eint(2,1) = (eint(1,1)*eint(2,2))**(0.5)
endif ! other switchtype mean no switch



end subroutine

