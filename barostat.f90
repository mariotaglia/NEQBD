subroutine berendsen(targetpres, barostat_const)
use system
use BD

real barostat_const
real eargetpres

real kai ! linear scale factor

kai = (1.0-(-pres+targetpres)/barostat_const)**(1.0/di)
!print*, kai, pres, targetpres
xlim = xlim*kai

xpos(:,:) = xpos(:,:)*kai 
end subroutine
