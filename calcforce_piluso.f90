
subroutine calcforce
use system
use nlist
use BD
use mpi
implicit none
integer l,i

real, external :: dist

interface

function eforce(i)
use system
integer i
real, dimension(di) :: eforce
endfunction
function LJforce(i, j)
use system
integer i, j
real, dimension(di) :: LJforce
endfunction
function HSforce(i, j)
use system
integer i, j
real, dimension(di) :: HSforce
endfunction
function Yforce(i, j)
use system
integer i, j
real, dimension(di) :: Yforce
endfunction
function YAVforce(i, j)
use system
integer i, j
real, dimension(di) :: YAVforce
endfunction
function LJEforce(i, j)
use system
integer i, j
real, dimension(di) :: LJEforce
endfunction
function LJNforce(i, j)
use system
integer i, j
real, dimension(di) :: LJNforce
endfunction
function tableforce(i, j)
use system
integer i, j
real, dimension(di) :: tableforce
endfunction
end interface

integer, external :: whichcell
integer icell, k
real fff(di)
integer jcell
integer ccc

call update_mpi_pos ! send the position from proc 0 to all others

forces_tosend = 0

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
!!! Calc force between l and i

    if(dist(l,i).lt.cutoff) then


    fff(:) = HSforce(l, i)


   ! DEBUG
   ! print*, l,i,fff

    forces_tosend(l, :) = forces_tosend(l, :) + fff(:)
    forces_tosend(i, :) = forces_tosend(i, :) - fff(:)
    endif

  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

fff(:) = eforce(l) ! external force for particle l
forces_tosend(l, :) = forces_tosend(l, :) + fff(:)

l = list(l) ! next in line
enddo ! while l
enddo ! jcell

call update_mpi_force ! send the forces from all procs to proce 0
    
!DEBUG
!stop

end


subroutine calcpres
use system
use nlist
use BD
use mpi
implicit none
integer l,i

real, external :: dist

interface
function HSforce(i, j)
use system
integer i, j
real, dimension(di) :: HSforce
endfunction
function eforce(i)
use system
integer i
real, dimension(di) :: eforce
endfunction
function LJforce(i, j)
use system
integer i, j
real, dimension(di) :: LJforce
endfunction
function Yforce(i, j)
use system
integer i, j
real, dimension(di) :: Yforce
endfunction
function YAVforce(i, j)
use system
integer i, j
real, dimension(di) :: YAVforce
endfunction
function LJEforce(i, j)
use system
integer i, j
real, dimension(di) :: LJEforce
endfunction
function LJNforce(i, j)
use system
integer i, j
real, dimension(di) :: LJNforce
endfunction

function tableforce(i, j)
use system
integer i, j
real, dimension(di) :: tableforce
endfunction


end interface

real, external :: distk
integer, external :: whichcell
integer icell, k
real fff(di)
integer jcell
integer ccc
integer kk

call update_mpi_pos ! send the position from proc 0 to all others

pres_tosend = 0

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
!!! Calc force*vector between l and i

    if(dist(l,i).lt.cutoff) then
    fff(:) = HSforce(l, i)
    do kk = 1, di 
    pres_tosend = pres_tosend + distk(l, i, kk)*fff(kk) 

! this is an expensive way to do this, note that
!
! f = f(r)*rx/r*i + f(r)*ry/r*j
! r = rx*i + ry*j
! fxr = f(r)*rx^2*i/r + f(r)*ry^2*j/r
! fxr = f(r)*r

    enddo
    endif

  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

fff = eforce(l) ! external contribution to virial
do kk = 1, di
pres_tosend = pres_tosend + xpos(l,kk)*fff(kk)	
enddo

l = list(l) ! next in line
enddo ! while l
enddo ! jcell

call update_mpi_pres ! send the forces from all procs to proce 0

pres = (Npart*temperature + (1.0/di)*pres)/(xlim**di)
end


