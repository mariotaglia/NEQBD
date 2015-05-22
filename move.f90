
subroutine BDmove
use system
use BD
implicit none
real mean, sd
real pre1, pre2
integer i, j
real rand1

pre1 = (jump**2)
pre2 = sqrt(2.0*temperature)*jump

mean = 0.0
sd = 1.0


do i = 1, Npart
do j = 1, di
call randomgauss(mean, sd, rand1)
xpos(i,j) = xpos(i, j) + pre1*forces(i,j) + pre2*rand1
xpos(i,j) = mod(xpos(i,j)+xlim, xlim)
enddo
enddo

end


