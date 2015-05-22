
subroutine reanalize
use system
use mpi
use mainm
implicit none
integer k, j
character (len = 3), dimension(Npart) :: atom, atoms
integer Nparttemp
real temp
real xp2(Npart,2), xp3(Npart,3)
LOGICAL :: file_exists
integer overdrive
integer Nperiod

overdrive = 1

atoms(1) = 'C '
atoms(2) = 'O '

INQUIRE(FILE="BDin.xyz", EXIST=file_exists) 

if(file_exists.neqv..true.) then
if(rank.eq.0)print*, 'Cannot find BDini.xyz'
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif

print*, 'Read input from file'
print*, 'and reanalize it...'

open(unit=15, file='BDin.xyz')
do while (.true.)

read(15, *, END=500) Nparttemp
if(Npart.ne.Nparttemp) then
   print*, 'Number of particles in infile is incorrect!'
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif
read(15, *, END=500) iniMCstep
if(rank.eq.0)print*, 'Read timestep', iniMCstep

do j = 1, Npart
if(di.eq.2)read(15, *, END=500)atom(j), xp2(j,1), xp2(j,2) , temp
if(di.eq.3)read(15, *, END=500)atom(j), xp3(j,1), xp3(j,2) , xp3(j,3)
enddo

do j = 1, Npart
if(di.eq.2)xpos(j,:) = xp2(j,:)
if(di.eq.3)xpos(j,:) = xp3(j,:)
enddo

Nparti=0

do j = 1, Npart
if(atom(j).eq.atoms(1)) then
tp(j) = 1
Nparti(1) = Nparti(1) + 1
elseif (atom(j).eq.atoms(2)) then
tp(j) = 2
Nparti(2) = Nparti(2) + 1
else 
   print*, 'Unrecognized particle type'
   print*, j, atom(j)
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif
enddo

overdrive=1
Nperiod = int(k/period)+1

call makelist
call analisis0(iniMCstep, Nperiod)
call switchpot(iniMCstep)
call analisis(iniMCstep, Nperiod, overdrive)

enddo ! read file

 500  continue

end

subroutine randompos
use system
!use IFPORT
implicit none
integer j, i, k
print*, 'Random positions'
do i = 1, Npart
  do j = 1, di
  xpos(i, j) = rand()*xlim
  enddo
  if(i.le.int(Npart*fa)) then
  tp(i) = 1
  Nparti(1) = Nparti(1) + 1
  else
  tp(i) = 2
  Nparti(2) = Nparti(2) + 1
  endif
enddo
end


subroutine squarepos
use system
use mpi
!use IFPORT
implicit none
integer j, i, l, k
integer side
integer tpt, randpos

if(rank.eq.0)print*, 'REGULAR INITIAL LATTICE'

!!!! 2D !!!!

if (di.eq.2) then
   side=int(sqrt(float(Npart)))
   tp(:) = 2
   Nparti = 0
  
   if((side**2).ne.Npart) then
     if(rank.eq.0)print*,'NOT A SQUARE'
     call MPI_FINALIZE(ierr) ! finaliza MPI
     stop
   endif

   do i = 1, side
   do j = 1, side
   l = i+(j-1)*side
     xpos(l, 1) = (xlim/float(side))*(float(i)-0.5)
     xpos(l, 2) = (xlim/float(side))*(float(j)-0.5)
     if(l.le.int(Npart*fa)) then
        tp(l) = 1
        Nparti(1) = Nparti(1) + 1
     endif
   enddo
   enddo

   Nparti(2) = Npart - Nparti (1)
endif

!!!! 3D !!!!

if (di.eq.3) then
   side=int((float(Npart))**(1.0/3.0))
   tp(:) = 2
   Nparti = 0
  
   if((side**3).ne.Npart) then
     if(rank.eq.0)print*,'NOT A CUBE'
     call MPI_FINALIZE(ierr) ! finaliza MPI
     stop
   endif

   do i = 1, side
   do j = 1, side
   do k = 1, side
   l = i+(j-1)*side+(k-1)*side*side
     xpos(l, 1) = (xlim/float(side))*(float(i)-0.5)
     xpos(l, 2) = (xlim/float(side))*(float(j)-0.5)
     xpos(l, 3) = (xlim/float(side))*(float(k)-0.5)
     if(l.le.int(Npart*fa)) then
        tp(l) = 1
        Nparti(1) = Nparti(1) + 1
     endif
   enddo
   enddo
   enddo

   Nparti(2) = Npart - Nparti (1)
endif


! scramble list
do l = 1, Npart
randpos = int(rand()*float(Npart))+1
tpt = tp(randpos)
tp(randpos) = tp(l)
tp(l) = tpt
enddo

end

subroutine squareideal
use system
use mpi
!use IFPORT
implicit none
integer j, i, l
integer side
integer tpt, randpos

if(di.ne.2) then
if(rank.eq.0)print*, 'SQURE IDEAL VALID FOR 2D ONLY'
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif

side=int(sqrt(float(Npart/2)))

tp(:) = 2
Nparti = 0

if((side**2).ne.Npart/2) then
if(rank.eq.0)print*,'NOT A SQUARE (x2)'
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif

do i = 1, side
do j = 1, side
l = i+(j-1)*side
  xpos(l, 1) = (xlim/float(side))*(float(i)-0.25)
  xpos(l, 2) = (xlim/float(side))*(float(j)-0.25)
  tp(l) = 1
  Nparti(1) = Nparti(1) + 1
enddo
enddo

do i = 1, side
do j = 1, side
l = i+(j-1)*side + Npart/2
  xpos(l, 1) = (xlim/float(side))*(float(i)-0.75)
  xpos(l, 2) = (xlim/float(side))*(float(j)-0.75)
  tp(l) = 2
  Nparti(2) = Nparti(2) + 1
enddo
enddo

end


subroutine readpos
use system
use mpi
use mainm
implicit none
integer k, j
character (len = 3), dimension(Npart) :: atom, atoms
integer Nparttemp
real temp
real xp2(Npart,2), xp3(Npart,3)
LOGICAL :: file_exists

atoms(1) = 'C '
atoms(2) = 'O '


INQUIRE(FILE="BDin.xyz", EXIST=file_exists) 

if(file_exists.neqv..true.) then
if(rank.eq.0)print*, 'Cannot find BDini.xyz'
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif

if(rank.eq.0)print*, 'Read input from file'

open(unit=15, file='BDin.xyz')
do while (.true.)

read(15, *, END=500) Nparttemp
if(Npart.ne.Nparttemp) then
   if(rank.eq.0)print*, 'Number of particles in infile is incorrect!'
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif
read(15, *, END=500) iniMCstep
if(rank.eq.0)print*, 'Read timestep', iniMCstep

do j = 1, Npart
if(di.eq.2)read(15, *, END=500)atom(j), xp2(j,1), xp2(j,2) , temp
if(di.eq.3)read(15, *, END=500)atom(j), xp3(j,1), xp3(j,2) , xp3(j,3)
enddo

do j = 1, Npart
if(di.eq.2)xpos(j,:) = xp2(j,:)
if(di.eq.3)xpos(j,:) = xp3(j,:)
enddo
enddo ! read file

 500  continue

do j = 1, Npart
if(atom(j).eq.atoms(1)) then
tp(j) = 1
Nparti(1) = Nparti(1) + 1
elseif (atom(j).eq.atoms(2)) then
tp(j) = 2
Nparti(2) = Nparti(2) + 1
else 
   print*, 'Unrecognized particle type'
   print*, j, atom(j)
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif
enddo

end


