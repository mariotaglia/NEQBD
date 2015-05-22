subroutine averagepot
use system
use nlist
use mpi
use BD
implicit none
integer ll,i, j, k

interface
function LJenergy(i, j)
use system
integer i, j
real :: LJenergy
endfunction

function LJEenergy(i, j)
use system
integer i, j
real :: LJEenergy
endfunction

function Yenergy(i, j)
use system
integer i, j
real :: Yenergy
endfunction

function YAVenergy(i, j)
use system
integer i, j
real :: YAVenergy
endfunction

function tableenergy(i, j)
use system
integer i, j
real :: tableenergy
endfunction

function LJNenergy(i, j)
use system
integer i, j
real :: LJNenergy
endfunction

end interface


procedure(LJenergy), pointer :: energypointer => NULL()

select case (forcetype)
case ('LJ0') ! standard LJ... no aditional parameters need
energypointer => LJenergy
case ('LJN')
energypointer => LJNenergy
case ('YUK')
energypointer => Yenergy
case ('YAV')
energypointer => YAVenergy
case ('ATR')
energypointer => LJEenergy
case ('TAB')
energypointer => tableenergy
end select

tabpos = 0.0
tabpot = 0.0

xpos(1,1) = 0.0
xpos(1,2) = 0.0
xpos(2,1) = 0.0

do i = 1, types
do j = 1, types

ntab(i,j) = maxn

tp(1) = i
tp(2) = j


do k = 1, maxn ! divide potential in steps

tabpos(i,j,k) = float(k)/float(maxn)*cutoff
xpos(2,2) = tabpos(i,j,k)

do ll = 1, 2*period
call switchpot(ll) ! update potential
!print*, i,j, ll, zint(tp(1)), zint(tp(2)) 
tabpot(i,j,k) = tabpot(i,j,k) + energypointer(1,2)/float(2*period)
!print*, i,j,k, tabpot(i,j,k)
enddo

enddo
enddo
enddo

forcetype = 'TAB'

end

