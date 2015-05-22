
subroutine makelist
use system
use nlist
use mpi
implicit none
integer icell, k
integer, external :: whichcell
integer flagerror_tosend, flagerror

flagerror = 0
flagerror_tosend = 0

head = 0

do k = 1, Npart
icell = whichcell(k)
if(icell.lt.0) then
print*, 'cell out of range'
print*, 'icell, k, xpos(k,1), xpos(k,2), xlim'
print*, icell, k, xpos(k,1), xpos(k,2), xlim 
flagerror_tosend = 1
CALL MPI_ALLREDUCE(flagerror_tosend, flagerror, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif
list(k) = head(icell)
head(icell) = k
enddo

CALL MPI_ALLREDUCE(flagerror_tosend, flagerror, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
if(flagerror.ne.0) then
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif

end


subroutine makecellneighbor
use system
use nlist
use mpi
implicit none
integer dx(3**di, di)
integer x(di), x0(di), i,k,j
integer xtemp, N, NN
integer icell

! cell system !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sidecell = int(xlim/cutoff)+1
Ncell = sidecell**di

allocate(list(Npart))
allocate(head(Ncell))
allocate(cellneighbor(Ncell,3**di))

do i = 1, (3**di) ! create dx
xtemp = 0

do j = 1, di
N = i-xtemp
NN = 3**j
dx(i,j) = mod(N, NN)
xtemp = xtemp + dx(i, j)
enddo

do j = 1, di
dx(i,j)=dx(i,j)/(3**(j-1))-1
enddo

enddo

do i = 0, sidecell**di-1 ! loop over all cells and get their coordinates
  xtemp = 0
  do j = 1, di
  x(j) = mod((i-xtemp), (sidecell**j))
  xtemp = xtemp+ x(j)
  enddo
  do j = 1, di
  x(j)=x(j)/(sidecell**(j-1)) 
  enddo

  do k = 1, 3**di

  do j = 1, di
  x0(j) = x(j) + dx(k,j)
  if(x0(j).ge.sidecell)x0(j)=x0(j)-sidecell
  if(x0(j).lt.0)x0(j)=x0(j)+sidecell
  x0(j) = x0(j)
  enddo
  
  icell = 0
  do j = 1, di    
  icell = icell + (x0(j))*sidecell**(j-1)
  enddo

  cellneighbor(i+1,k) = icell+1
!  print*, i, x(1), x(2), x0(1), x0(2)
!  print*, i, k, cellneighbor(i,k)
  enddo
enddo 

end

integer function whichcell(l)
use nlist
use system
implicit none
integer icell, l, i

icell = 1

do i = 1, di
icell = icell + (int((xpos(l,i))/xlim*float(sidecell)))*sidecell**(i-1)
enddo
whichcell = icell
end


