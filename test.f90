     integer period 
     period = 100
     do k = 1, 5*period
     alpha = abs(float(mod((k+period),2*period))/float(period)-1.0) ! from 0 to 1 V shape
     pH=(exp(10.0*(alpha-0.5))-exp(-10.0*(alpha-0.5)))/74.5 + 5
     f=1.0/(1.0 + 10**(5.0-pH))
     print*, k, pH
     enddo
     end

