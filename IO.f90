
subroutine loadconf
use system
use BD
use mainm
use mpi
use random
implicit none
character basura
integer i, j, k
character (len=30) :: filename
logical file_exists

open(file='pressure.dat', unit=999)
open(file='density.dat', unit=998)

!!!!!!!!!!!!!!!!!!!!!!! read from file !!!!!!!!!!!!!!!!!!!!!!!!!!

INQUIRE(FILE="input.txt", EXIST=file_exists) 

if(file_exists.neqv..true.) then
if(rank.eq.0)print*, 'Cannot find input.txt'
!   call printinput
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
endif


open(file='input.txt', unit=20)

read(20, *)basura
read(20, *)di
if(rank.eq.0)print*, 'Dimensions', di

read(20, *)basura
read(20, *)Npart
if(rank.eq.0)print*, 'Number of particle', Npart

read(20, *)basura
read(20, *)forcetype
if(rank.eq.0)print*, 'Type of potential ', forcetype

select case (forcetype)
case ('LJ0') ! standard LJ... no aditional parameters need
if(rank.eq.0)print*, 'Repulsive 12 - LJ' 

case ('LJN')
if(rank.eq.0)print*, 'Repulsive N (even) - LJ'
read(20, *) basura
read(20, *) NLJ, NNLJ, alfaLJ
if(rank.eq.0)print*, 'N LJ (rep):', NLJ, 'N LJ (atr)', NNLJ, 'Alpha LJ:', alfaLJ

case ('YUK')
if(rank.eq.0)print*, 'Yukawa attraction between A and B + LJ12'
read(20, *) basura
read(20, *) eexp, elen
if(rank.eq.0)print*,  'Prefactor exponetial attracion:',  eexp, ' decay: ', elen

case ('YAV')
if(rank.eq.0)print*, 'Yukawa attraction between A and B + LJ12 with averaged charged'
read(20, *) basura
read(20, *) eexp, eexpav(1,1), eexpav(2,2), eexpav(1,2) , elen
eexpav(2,1)=eexpav(1,2)
if(rank.eq.0)print*,  'Prefactor exponetial attracion:', eexp, ' * ',  eexpav, ' decay: ', elen

case ('ATR')
if(rank.eq.0)print*, 'Yukawa-type attraction between all particles'
read(20, *) basura
read(20, *) eexp, elen
if(rank.eq.0)print*,  'Prefactor exponetial attracion:',  eexp, ' decay: ', elen

case ('TAB')
if(rank.eq.0)print*, 'Tabulated potential -- STATIC'
if(rank.eq.0)print*, 'Read files TAB.txt'
tabpos = 0.0
tabpot = 0.0
tabforce = 0.0

do i = 1, types
do j = i, types
write(filename,'(A4, I3.3, A1, I3.3, A4)')'TAB.', i,'.',j,'.txt'
open(file=filename, unit=920)
k=0
do 
k=k+1
read(920, *, END=50)tabpos(i,j,k), tabpot(i,j,k)
tabpos(j,i,k)=tabpos(i,j,k)
tabpot(j,i,k)=tabpot(i,j,k)
enddo
 50 continue
 ntab(i,j) = k-1
 ntab(j,i) = k-1
 close(920)
enddo
enddo

end select

read(20, *) basura
read(20, *) jump
if(rank.eq.0)print*, '(D*dt)^0.5', jump

read(20, *) basura
read(20, *) temperature
if(rank.eq.0)print*, 'Temperature', temperature

read(20, *) basura
read(20, *) r1a, r1b, r2a, r2b
read(20, *) basura
read(20, *) e1a, e1b, e2a, e2b
read(20, *) basura
read(20, *) z1a, z1b, z2a, z2b

if(rank.eq.0)print*, r1a, ' < r1 < ', r1b
if(rank.eq.0)print*, r2a, ' < r2 < ', r2b

if(rank.eq.0)print*, e1a, ' < e1 < ', e1b
if(rank.eq.0)print*, e2a, ' < e2 < ', e2b

if(rank.eq.0)print*, z1a, ' < z1 < ', z1b
if(rank.eq.0)print*, z2a, ' < z2 < ', z2b


read(20, *) basura
read(20, *) readinput 
if(rank.eq.0)print*, 'Readinput', readinput
if(readinput.eq.1) then
if(rank.eq.0)print*, 'Reading from BDin.xyz'
endif

read(20, *) basura
read(20, *) switchtype
if(rank.eq.0)print*, 'Switchtype', switchtype

read(20, *) basura
read(20, *) interpart
if(rank.eq.0)print*, 'Average interparticle distance', interpart

density = 1.0/(interpart**di)
if(rank.eq.0)print*, 'Density', density

xlim = (Npart/density)**(1.0/di)
if(rank.eq.0)print*, 'Box side', xlim
if(rank.eq.0)print*, 'Box Area/Volume', xlim**di

read(20, *) basura
read(20, *) fa
if(rank.eq.0)print*, 'Fraction of A', fa
if(rank.eq.0)print*, 'Molecules of type A', int(fa*Npart)
if((fa.lt.0.0).or.(fa.gt.1)) then
if(rank.eq.0)print*, 'Check fa!'
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif

read(20, *) basura
read(20, *) period
if(rank.eq.0)print*, 'Switching HALF period', period

read(20, *) basura
read(20, *) period2, phase
if(rank.eq.0)print*, 'Switching HALF period for particles 2 ', period2, ',phase ', phase


read(20, *) basura
read(20, *) aveflag
if(aveflag.eq.1) then
if(rank.eq.0)print*, 'Averaged potential'
endif


read(20, *) basura
read(20, *) MCstep
if(rank.eq.0)print*, 'BD steps', MCstep

read(20, *) basura
read(20, *) saveevery
if(rank.eq.0)print*, 'Save coordinates every', saveevery

read(20, *) basura
read(20, *) offset
if(rank.eq.0)print*, 'Save coordinate offset', offset
if(offset.ge.saveevery) then
if(rank.eq.0)print*, 'Offset > Save every, set offset to 0'
offset=0
endif

read(20, *) basura
read(20, *) flaggr
if(rank.eq.0)print*, 'Save gr?', flaggr

read(20, *) basura
read(20, *) flagbo
if(rank.eq.0)print*, 'Save bo?', flagbo

if(di.gt.2) then
  if(rank.eq.0)print*, 'bo only for 2D calculations'
  flagbo = 0
endif

read(20, *) basura
read(20, *) flagclus, flagxyz, BOcutoff, ddd_cutoff
if(rank.eq.0)print*, 'Find clusters?', flagclus
if(rank.eq.0)print*, 'Use BOcutoff', BOcutoff
if(rank.eq.0)print*, 'save xyz?', flagxyz

if(di.gt.2) then 
  if(rank.eq.0)print*, 'cluster analysis only for 2D calculations'
  flagclus = 0  
  BOcutoff = 0
  flagxyz = 0
endif

read(20, *) basura
read(20, *) flagcr
if(rank.eq.0)print*, 'Save crystals?', flagcr

if(di.ne.2) then
  if(rank.eq.0)print*, 'crystal analysis valid only for 2D calculations'
  flagclus = 0  
  BOcutoff = 0
  flagxyz = 0
endif

read(20, *) basura
read(20, *) flagdis
if(rank.eq.0)print*, 'Save energy?', flagdis

read(20, *)basura
read(20, *)saveallfiles
if(rank.eq.0)print*,'Save all bo and gr files?', saveallfiles

read(20, *) basura
read(20, *) saveeveryother, saveeveryothere
if(rank.eq.0)print*, 'Save data every ', saveeveryother, ' periods'
if(rank.eq.0)print*, 'Save energy data every ', saveeveryothere, ' periods'

read(20, *) basura
read(20, *) savepoints, savepointse
if(rank.eq.0)print*, 'Number of points to average inside period', savepoints
if(rank.eq.0)print*, 'Number of points to average inside period (energy)', savepointse

read(20, *) basura
read(20, *) cutoff
if(rank.eq.0)print*, 'Cutoff', cutoff
cutoff12 = cutoff**12.0

read(20, *) basura
read(20, *) calcpresevery
if(rank.eq.0)print*, 'Calculate pressure every (and update barostat)', calcpresevery

read(20, *) basura
read(20, *) barostat, barostat_const, targetpres
if(rank.eq.0)print*, 'Use barostat? ', barostat, ' Constant ', barostat_const, ' Target pres', targetpres

read(20, *) basura
read(20, *) seed
if(rank.eq.0)print*, 'Seed', seed
close (20)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Variable initalization an allocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(Nparti(types))
allocate(xpos(Npart,di))
allocate(forces(Npart,di))
allocate(forces_tosend(Npart,di))
allocate(tp(Npart))
allocate(rtp(Npart))

Nparti = 0

end

subroutine savecurrent(k)
use system
implicit none
integer k, j
character (len=3) :: atom
write(10, *) Npart
write(10, *) k
do j = 1, Npart

   if(tp(j).eq.1) then 
   atom = ' C ' 
   else
   atom = ' O ' 
   endif

if(di.eq.2)write(10, *)atom, xpos(j,1), xpos(j,2) , 0.0, rint(tp(j),tp(j))
if(di.eq.3)write(10, *)atom, xpos(j,1), xpos(j,2) , xpos(j,3), rint(tp(j),tp(j))
enddo
   write(11, *)(rint(tp(j),tp(j)), j=1,Npart)
end


subroutine makespline
use system
use mpi
use nlist
implicit none
integer i,j, k
real xspl, yspl1, yspl2
real x(maxn), y(maxn)
character (len=30) :: filename
integer error
integer, parameter :: npoints = 100

! check that potential is long enough



do k = 1, types
do j = 1, types
if(tabpos(k,j,ntab(k,j)).lt.cutoff) then
print*, 'Tabulated potential maximum of', tabpos(k,j,ntab(k,j))
print*, 'is shorter than cutoff of', cutoff
endif
enddo
enddo

! calculate the derivative at each position using the spline function

!do i = 1, types
!do j = 1, types

!x=tabpos(i,j,:)
!y=tabpot(i,j,:)

!do k = 1, ntab(i, j) 
!xspl = tabpos(i,j,k)
!call Interpol_Akima(ntab(i,j),error,xspl,yspl1,X,Y)
!xspl = tabpos(i,j,k)+delta
!call Interpol_Akima(ntab(i,j),error,xspl,yspl2,X,Y)
!tabforce(i,j,k) = -(yspl2-yspl1)/delta
!enddo ! k
!enddo
!enddo

! make the force in the origin very large
!tabforce(:,:,0) = tabforce(:,:,1)*1e3


! save splined data
do i = 1, types
do j = 1, types
if(rank.eq.0) then
write(filename,'(A7, I3.3, A1, I3.3, A4)')'POTSPL.', i,'.',j,'.dat'
open(file=filename, unit=820)
write(filename,'(A7, I3.3, A1, I3.3, A4)')'FORSPL.', i,'.',j,'.dat'
open(file=filename, unit=830)
endif

x(:)=tabpos(i,j,:)
y(:)=tabpot(i,j,:)


do k = 1, maxn
xspl = cutoff/float(maxn-1)*float(k-1)

if((xspl.ge.(tabpos(i,j,2))).and.((xspl+delta).le.(tabpos(i,j,(ntab(i,j)-3))))) then
call Interpol_Akima(ntab(i,j),error,xspl,yspl1,X,Y)
call Interpol_Akima(ntab(i,j),error,(xspl+delta),yspl2,X,Y)
tabpotspl(i,j,k) = yspl1
tabforcespl(i,j,k) = (yspl1-yspl2)/delta
else if(xspl.lt.(tabpos(i,j,2))) then
tabpotspl(i,j,k) = 1e31
tabforcespl(i,j,k) = 1e31
else if ((xspl+delta).gt.(tabpos(i,j,(ntab(i,j)-3)))) then
tabpotspl(i,j,k) = 0.0
tabforcespl(i,j,k) = 0.0
endif 
tabposspl(i,j,k) = xspl

if (rank.eq.0) then
write(820,*)xspl,tabpotspl(i,j,k)
write(830,*)xspl,tabforcespl(i,j,k)
endif

enddo

do k = 1, maxn-1
tabpotm(i,j,k)=(tabpotspl(i,j,k+1)-tabpotspl(i,j,k))/(tabposspl(i,j,k+1)-tabposspl(i,j,k))
tabforcem(i,j,k)=(tabforcespl(i,j,k+1)-tabforcespl(i,j,k))/(tabposspl(i,j,k+1)-tabposspl(i,j,k))
enddo

close(820)
close(830)

enddo
enddo

end

