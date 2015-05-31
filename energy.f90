subroutine penergy
use system
use nlist
use mpi
use BD
implicit none
integer l,i
integer, external :: whichcell
integer icell, k
integer jcell
integer ccc

real, external :: dist
real, external :: eenergy


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

call update_mpi_pos ! send the position from proc 0 to all others

senergy_tosend = 0.0

do ccc = 1, Nlistproc(rank+1)
jcell = listproc(rank+1,ccc)

  l = head(jcell) ! head of cell
  do while (l.ne.0) ! loop over all particles in jcell

  do k = 1, 3**di ! loop over all neighbors of jcell

  icell = cellneighbor(jcell,k) 
   
  if(icell.eq.jcell) then
  i = list(l) ! same cell, start from the next particle
  elseif(icell.lt.jcell) then
  i = head(icell) ! different cell, start from the head particle of icell
  else
  i = 0 ! only icell < jcell to not overcount interactions
  endif

  do while (i.ne.0)
  if(dist(l,i).lt.cutoff) then
  senergy_tosend = senergy_tosend + energypointer(l, i)
  endif
  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

  senergy_tosend = senergy_tosend + eenergy(l) ! external contribution to energy

l = list(l) ! next in line
enddo ! while l
enddo ! jcell

call update_mpi_energy ! sum the energies from all procs 

end

!real function energy(l)
!use system
!use nlist
!implicit none
!integer l,i
!real e, r
!real vect
!real, external :: dist
!integer ll
!integer, external :: whichcell
!integer icell, k

!energy = 0

!ll = whichcell(l)

!do k = 1, 3**di
!icell = cellneighbor(ll,k)
!i = head(icell)
!do while(i.ne.0)
!if(i.ne.l) then

!e = eint(tp(i), tp(l))
!r = rint(tp(i), tp(l))

!vect = dist(l, i)
!energy = energy + e*((r/vect)**12)
!endif
!i = list(i)
!enddo
!enddo ! k

!end

subroutine plotpot
use system
use nlist
use mpi
use BD
implicit none
integer l,i
integer, external :: whichcell
integer icell, k
integer jcell
integer ccc
real temp
real, external :: dist

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

xpos(1,1) = 0.0
xpos(1,2) = 0.0
xpos(2,1) = 0.0

open(file='POT11.dat', unit=901)
open(file='POT12.dat', unit=902)
open(file='POT22.dat', unit=903)

do i = 1, 100
xpos(2,2) = float(i)/100.0*cutoff

tp(1) = 1
tp(2) = 1
temp = energypointer(1,2)
write(901,*)xpos(2,2), temp

tp(1) = 1
tp(2) = 2
temp = energypointer(1,2)
write(902,*)xpos(2,2), temp

tp(1) = 2
tp(2) = 2
temp = energypointer(1,2)
write(903,*)xpos(2,2), temp
enddo

close(901)
close(902)
close(903)

end


