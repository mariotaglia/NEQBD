
module mpi
include 'mpif.h'
integer rank, size, ierr
integer stat(MPI_STATUS_SIZE)
integer tag, source, dest
parameter(tag = 0)
integer err
integer flagsolver
integer ier_tosend
double  precision norma_tosend
integer, allocatable :: Nlistproc(:)
integer, allocatable :: listproc(:,:)
endmodule

module random
integer seed
end module

module externalforce
integer eftype
real ef ! magnitude of external force
integer eperiod ! period of rotation of external force
real ex, ey ! current angle of external force
end module

module params
real :: pi = acos(-1.0)
endmodule

module system
integer aveflag
integer saveallfiles
real r1a, r2a, r1b, r2b
real e1a, e2a, e1b, e2b
real z1a, z2a, z1b, z2b
integer Npart
integer di
integer, parameter :: types = 2
Integer, allocatable :: Nparti(:)
real xlim
real, allocatable :: xpos(:,:)
real, allocatable :: xposCOM(:,:) ! position of COM averaged over one period
real, allocatable :: xposCOMt(:,:) ! position of COM averaged over one period
integer, allocatable :: tp(:)
integer, allocatable :: rtp(:)
integer grcounter
integer switch
real fa
real temperature
real density
real interpart
integer period
integer switchtype
real cutoff
real BOcutoff
real ddd_cutoff
real cutoff12
character*3 forcetype
integer saveCOM
integer, parameter :: maxn = 2000

real NLJ
real NNLJ
real eLJ
real alfaLJ
real eexp, elen
real eexpav(types,types)
real rint(types, types), eint(types, types), zint(types)

! tabulated potential
real tabpot(types,types,maxn) 
real tabpos(types,types,maxn) 
real tabforce(types,types,maxn) 
integer ntab(types,types)
real, parameter :: delta=1e-5 ! delta for numerical derivative

real tabpotspl(types,types,maxn) 
real tabpotm(types,types,maxn) 
real tabposspl(types,types,maxn) 
real tabforcespl(types,types,maxn) 
real tabforcem(types,types,maxn) 

endmodule

module bo
integer, allocatable :: bo2co(:), bo4co(:), bo6co(:), bo1co(:),bo3co(:)
real, allocatable :: bo4r(:,:,:)
real, allocatable :: bo6r(:,:,:)
real, allocatable :: bo2r(:,:,:)
real, allocatable :: bo1r(:,:,:)
real, allocatable :: bo3r(:,:,:)
real, allocatable ::  bo4rt(:), bo6rt(:)
real, allocatable :: bo2rt(:), bo1rt(:), bo3rt(:)
real, allocatable ::  bo4r_tosend(:,:,:)
real, allocatable ::  bo6r_tosend(:,:,:)
real, allocatable ::  bo2r_tosend(:,:,:)
real, allocatable ::  bo1r_tosend(:,:,:)
real, allocatable ::  bo3r_tosend(:,:,:)
real, allocatable ::  bo4rt_tosend(:)
real, allocatable ::  bo6rt_tosend(:)
real, allocatable ::  bo2rt_tosend(:)
real, allocatable ::  bo1rt_tosend(:)
real, allocatable ::  bo3rt_tosend(:)

complex, allocatable :: bo2aux(:)
complex, allocatable :: bo4aux(:)
complex, allocatable :: bo6aux(:)
complex, allocatable :: bo1aux(:)
complex, allocatable :: bo3aux(:)

complex, allocatable :: bo2aux_tosend(:)
complex, allocatable :: bo4aux_tosend(:)
complex, allocatable :: bo6aux_tosend(:)
complex, allocatable :: bo1aux_tosend(:)
complex, allocatable :: bo3aux_tosend(:)

integer, allocatable :: NVECINOS2(:)
integer, allocatable :: NVECINOS2_tosend(:)
integer, allocatable :: VECINOS2(:,:)
integer, allocatable :: VECINOS2_tosend(:,:)
integer, allocatable :: NVECINOS4(:)
integer, allocatable :: NVECINOS4_tosend(:)
integer, allocatable :: VECINOS4(:,:)
integer, allocatable :: VECINOS4_tosend(:,:)
integer, allocatable :: NVECINOS6(:)
integer, allocatable :: NVECINOS6_tosend(:)
integer, allocatable :: VECINOS6(:,:)
integer, allocatable :: VECINOS6_tosend(:,:)
integer, allocatable :: NVECINOS1(:)
integer, allocatable :: NVECINOS1_tosend(:)
integer, allocatable :: VECINOS1(:,:)
integer, allocatable :: VECINOS1_tosend(:,:)
integer, allocatable :: NVECINOS3(:)
integer, allocatable :: NVECINOS3_tosend(:)
integer, allocatable :: VECINOS3(:,:)
integer, allocatable :: VECINOS3_tosend(:,:)
integer, parameter :: VV = 20

integer, allocatable :: clusters(:)
integer, allocatable :: clustersize(:)

endmodule


module BD
real jump
real, allocatable :: forces(:, :)
real pres
real, allocatable :: forces_tosend(:, :)
real pres_tosend
real senergy
real senergy_tosend
end module

module grsystem
real dbin ! binning distance
integer Ngrmax
real, allocatable :: GR(:,:,:)
real, allocatable :: GR_tosend(:,:,:)
endmodule

module nlist
integer sidecell, Ncell ! number of cells
integer, allocatable :: head(:)
integer, allocatable :: list(:)
integer, allocatable :: cellneighbor(:,:)
endmodule

module mainm
integer MCstep
integer iniMCstep
integer saveevery
integer offset
integer saveeveryother
integer saveeveryothere
integer savepoints
integer savepointse
integer flaggr
integer flagdis
integer flagxyz
integer flagclus
integer flagbo
integer flagcr
integer readinput
integer barostat
real barostat_const
integer calcpresevery
real targetpres
endmodule

module analysis
real U0, U1, DU, UP, DUP
integer cperiod, ck
endmodule

program BDOSC
use system
use nlist
use BD
use random
use mpi
use mainm
use params
use analysis
!use IFPORT
implicit none
integer i, j, k
!real U0
!real expdif
!real xpos0(di)
real, external :: randomgauss
real alpha
integer Nperiod ! period number 
integer overdrive
! timming
real start, finish
integer :: timming = 0

overdrive=0

!!!!!!! INIT MPI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!! INTRO !!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(rank.eq.0)print*, 'Brownian Dynamic Simulator'
if(rank.eq.0)print*, 'Git Version: ', _VERSION 
if(rank.eq.0)print*, 'Number of processes', size
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call loadconf ! load configuration from file


call makecellneighbor ! make cell list


if(aveflag.eq.1)call averagepot ! average potential and make a tabulated one

if(forcetype.eq.'TAB') then
call makespline
if(rank.eq.0)call plotpot ! save potentials to disk
endif

call proc ! distribute list on processors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output stuff
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0) then
open(unit=10, file='BD.xyz')
open(unit=11, file='sizes.dat')
endif

if((flagclus.eq.1).and.(flagxyz.eq.1)) then
if(rank.eq.0) then
open(unit=1010, file='BDcluster.xyz')
open(unit=1012, file='cluster2.dat')
open(unit=1014, file='cluster4.dat')
open(unit=1016, file='cluster6.dat')
open(unit=1022, file='cluster_tagged_2.dat')
endif
endif

call allocatebin  ! allocates binning for g(r)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
UP = 0
DUP = 0


iniMCstep = 1
if(readinput.eq.0) then
!call randompos
call squarepos
else if(readinput.eq.1) then
call readpos
else if(readinput.eq.3) then
call readpos
MCstep=iniMCstep+period*2
else if(readinput.eq.-1) then
call squareideal
else if(readinput.eq.2) then
call reanalize
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif

call switchpot(iniMCstep-1) ! updates the potential to current time step

!!!!!!!!!!!!!!!! MAIN MC LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cperiod = saveeveryothere ! for energy disipation
xposCOM = xpos

do k = iniMCstep, MCstep

Nperiod = int((k-iniMCstep)/period)+1

if(timming.eq.1)call cpu_time(start)
call makelist

if(timming.eq.1)call cpu_time(finish)
if(timming.eq.1) then
  if(rank.eq.0)print '("Make list time = ",f6.4," seconds.")',finish-start
endif

call analisis0(k, Nperiod)
!!!!!!!!! switch potential
call switchpot(k)
!!!!!!!!! integrate COM
call intCOM(k)

call analisis(k, Nperiod, overdrive)


if(timming.eq.1)call cpu_time(start)
call calcforce ! calculate forces in the system
if(timming.eq.1)call cpu_time(finish)

!DEBUG
!do i = 1, Npart
!print*, i, forces(i,1),forces(i,2)
!enddo
!stop

if(timming.eq.1) then
 if(rank.eq.0)print '("Calcforce time = ",f6.4," seconds.")',finish-start
endif

if(timming.eq.1)call cpu_time(start)
if(rank.eq.0)call BDmove ! only processor 0 moves particles 
if(timming.eq.1)call cpu_time(finish)
if(timming.eq.1) then
  if(rank.eq.0)print '("BDmove time = ",f6.4," seconds.")',finish-start
endif

if(mod(k, calcpresevery).eq.0) then
call calcpres ! calculate pressure

if(rank.eq.0) then
!print*, k, pres
write(999,*)k, pres
write(998,*)k, Npart/(xlim**di)
endif

if(barostat.eq.1) then ! call barostat
call berendsen(targetpres, barostat_const)
endif
endif

if(timming.eq.1) then
call MPI_FINALIZE(ierr) ! finaliza MPI
stop
endif

enddo ! k

!!!! close IO stuff

close(10)
close(20)

open(unit=10, file='status')
write(10,*)'OK'
close(10)

call MPI_FINALIZE(ierr) ! finaliza MPI
end

subroutine analisis(k, Nperiod, overdrive)
use system
use BD
use mpi
use mainm
use params
use analysis
implicit none
integer k, Nperiod, overdrive

if(((flagdis.eq.1).and.(mod(Nperiod, saveeveryothere).eq.0).and.(mod(k, int(period/savepointse)).eq.0)).and.&
overdrive.eq.0) then
         call penergy
         U1 = senergy
         DU = U0-U1

         if (cperiod.eq.Nperiod) then
           DUP = DUP + DU
           ck = k
         else
            if(rank.eq.0) then
            open(unit=26, file='denergyP.dat', access='append')
            write(26,*), ck, DUP/float(savepointse)
            close(26)
            endif
            cperiod = Nperiod
            ck = k
            DUP = DU
            UP = U0
         endif 

         if(rank.eq.0) then
         print*,'Step = ', k, ' DU = ', DU

         open(unit=25, file='denergy.dat', access='append')
         write(25,*), k, DU
         close(25)
         endif
endif



!!!!!!! save coordinates and do analysis
   if(((flaggr.eq.1).and.(mod(Nperiod, saveeveryother).eq.0).and.(mod(k, int(period/savepoints)).eq.0)).or.&
overdrive.eq.1) then
         call savegr(k)
   endif
   if(((flagbo.eq.1).and.(mod(Nperiod, saveeveryother).eq.0).and.(mod(k, int(period/savepoints)).eq.0)).or.&
overdrive.eq.1) then
         call psavebo(k)
   endif

   if(((flagclus.eq.1).and.(mod(Nperiod, saveeveryother).eq.0).and.(mod(k, int(period/savepoints)).eq.0)).or.&
overdrive.eq.1) then
         call clusterbo3(k)
         if(rank.eq.0)call findcluster2(k)
   endif

   if(((flagclus.eq.2).and.(mod(Nperiod, saveeveryother).eq.0).and.(mod(k, int(period/savepoints)).eq.0)).or.&
overdrive.eq.1) then
         call clusternn(k)
         if(rank.eq.0)call findcluster2(k)
   endif



if (rank.eq.0) then
   if(((flagcr.eq.1).and.(mod(Nperiod, saveeveryother).eq.0).and.(mod(k, int(period/savepoints)).eq.0)).or.&
overdrive.eq.1) then
         call savecr(k)
   endif
   if(((mod(k, saveevery).eq.offset).or.(k.eq.MCstep)).and.overdrive.eq.0) then
         call savecurrent(k) ! save coordinates
         print*, 'Step =', k, 'd1, d2 =', rint(1,1), rint(2,2)
    endif
endif
end ! main loop

subroutine analisis0(k, Nperiod)
use system
use BD
use mpi
use mainm
use params
use analysis
implicit none
integer k, Nperiod, overdrive

!!!!!!!!!!! for energy dissipation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

if(((flagdis.eq.1).and.(mod(Nperiod, saveeveryothere).eq.0).and.(mod(k, int(period/savepointse)).eq.0))) then
         call penergy
         U0 = senergy
         
         if (cperiod.eq.Nperiod) then
           UP = UP + U0
           ck = k
         else
           if(rank.eq.0) then
            open(unit=21, file='energyP.dat', access='append')
            write(21,*), ck, UP/float(savepointse)
            close(21)
            endif
         endif 

         if(rank.eq.0) then
         print*,'Step = ', k, ' U = ', U0

         open(unit=20, file='energy.dat', access='append')
         write(20,*), k, U0
         close(20)
         endif
endif
end

subroutine intCOM(k) ! integrates COM positions over one period of external field...
use system
use externalforce
implicit none
real, external :: distkCOM
integer k, kk, i,j

do i =  1, Npart
do kk = 1, di
xposCOMt(i,kk) = xposCOMt(i,kk) + distkCOM(i,i,kk)
enddo ! kk
enddo ! i

if(mod(k, 2*eperiod).eq.0) then
 xposCOM = xposCOM + xposCOMt/float(eperiod)/2.0
do i =  1, Npart
do j = 1, di
 xposCOM = mod(xposCOM(i,j)+xlim, xlim)
enddo ! j
enddo ! i
 xposCOMt = 0.0
endif

end subroutine
