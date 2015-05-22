subroutine switchpot(k)
use system
use params
use mainm
implicit none

integer k
real alpha, pH, f

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
else if (switchtype.eq.2) then ! switchtype = 2, triangular function
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
else if (switchtype.eq.3) then ! switchtype = 3, triangular function
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
else if (switchtype.eq.4) then ! switchtype = 3, triangular function
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
else if (switchtype.eq.6) then ! switchtype = 5, pH-type function
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
else if (switchtype.eq.7) then ! switchtype = 5, pH-type function
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
