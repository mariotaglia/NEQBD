subroutine savecr(step)
use system
use grsystem
use params
use mpi
implicit none

integer step
integer j,k, i, ii

real, external :: angle
real, external :: dist
integer l

!!! crystal

integer list(4)
real listdist(4)
integer list2(4)
real ccutoff, cutoffpos
real a1, d1
integer flag(4)
integer ll
real w(4,2), v(2)
real crystal
integer NNN
real vect
integer tempi
real tempr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine number of sites that belong to a lattice
!
! only works for binary mixtures
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


if(rank.eq.0) then
open(file='crystalfrac.dat', unit=500, access='append')
endif

cutoffpos = 0.2*(max(r1a,r1b,r2a,r2b)+min(r1a,r1b,r2a,r2b))  

crystal = 0.0
do j=1, Npart

    listdist(:) = xlim
 
    do k = 1, Npart
    vect = dist(j,k)

    if((vect.lt.listdist(1)).and.(k.ne.j)) then

      listdist(1) = vect ! drops the largest one
      list(1) = k
  
      do i = 1, 3 ! sort the list
      do ii = i+1, 4
      if(listdist(i).lt.listdist(ii)) then ! swap
      tempi = list(i)
      tempr = listdist(i)
      list(i)= list(ii)
      listdist(i)= listdist(ii)
      list(ii) = tempi
      listdist(ii) = tempr
      endif
      enddo
      enddo 
    endif
    enddo

!    print*, k, listdist(1), listdist(2),listdist(3),listdist(4)
 

    if ((tp(list(1)).ne.tp(j)).and.(tp(list(2)).ne.tp(j))&
    &.and.(tp(list(3)).ne.tp(j)).and.(tp(list(4)).ne.tp(j))) then

    d1 = dist(list(1),j)
    a1 = angle(list(1),j)
   
    do k = 1,4 
    w(k,1) = d1*sin(a1) 
    w(k,2) = d1*cos(a1)
    a1 = a1 + pi/2;  
    enddo ! k
 
    flag = 0
    
    do ll = 2,4   
    do l = 2,4
    d1 = dist(list(l),j)
    a1 = angle(list(l),j)
    v(1) = d1*sin(a1) 
    v(2) = d1*cos(a1)   
    vect= (v(1)-w(ll,1))**2 + (v(2)-w(ll,2))**2
    if(vect.le.cutoffpos**2) then
        flag(ll) = 1
        list2(ll) = list(l)
        exit
    endif      
    enddo ! ll
    enddo ! l
    
    flag(1) = flag(2)*flag(3)*flag(4)
    
    if(flag(1).eq.1) then
    crystal = crystal + 1.0
    endif 

    endif ! N=4
enddo ! j

write(500, *)step, crystal/float(Npart)
print*, 'Step:', step, 'Ncrystal: ', crystal
close(500)
end

subroutine addgr
use system
use grsystem
implicit none

integer i, j 
real vect
real, external :: dist

GR = 0 ; !reset GR
do i = 1, Npart
do j = 1, Npart
vect = dist(i, j)
GR(int(vect/dbin)+1,tp(i),tp(j)) = GR(int(vect/dbin)+1,tp(i),tp(j))+1.0
enddo
enddo

end


subroutine savegr(step)
use system
use grsystem
use params
use mpi
implicit none
real rpos(Ngrmax)

integer step
integer i, j,k
real vect
character (len=30) :: filename
integer maxsv

real fshellr(types,types), fshellsum(types,types)
integer fshell(types, types), peak(types, types)

real, external :: angle
real, external :: dist
real aa
real GRtemp(Ngrmax)
real firstpeak
integer l

call paddgr

if(rank.eq.0) then
do i = 1, types
do j = 1, types
write(filename,'(A3, I3.3, A1, I3.3, A4)')'CN.', i,'.',j,'.dat'
open(file=filename, unit=300+i*types+j, access='append')
write(filename,'(A7, I3.3, A1, I3.3, A4)')'cutoff.', i,'.',j,'.dat'
open(file=filename, unit=400+i*types+j, access='append')
enddo
enddo
endif

!!!!!! save gr

maxsv = 300

do i = 1, types
GR(:,i,:) = GR(:,i,:)/Nparti(i)
enddo !i


do i = 1, Ngrmax
rpos(i) = (float(i)-0.5)*dbin
GR(i,:,:) = GR(i,:,:)/(2.0*pi*rpos(i)*dbin)
enddo !i
GR(1,:,:) = 0.0 ! remove self-interaction

if(saveallfiles.eq.1) then
do i = 1, types
do j = i, types
write(filename,'(A3, I3.3, A1, I3.3, A1, I8.8, A4)')'GR.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=30)
do k = 1,maxsv
write(30,*)rpos(k), GR(k,i,j)
enddo !k
close(30)
enddo !i
enddo !j
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! find minima of GR

GRtemp = 0
firstpeak = max(r1a,r1b,r2a,r2b)*1.2 ! maximum location of first peak
l = int(firstpeak/dbin)+1
fshell = 0

do i = 1, types
do j = 1, types
GRtemp(1:l) = GR(1:l,i,j)
peak(i,j)=maxloc(GRtemp,1)
do k = peak(i,j), Ngrmax-1
if(GR(k+1,i,j)>GR(k,i,j)) then 
fshell(i,j)=k
fshellr(i,j)=rpos(k)
exit
endif
enddo ! k
enddo ! j
enddo ! i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integrate first shells

do i = 1, types
do j = 1, types
fshellsum(i,j) = 0.0
do k = 1, fshell(i,j)
fshellsum(i,j)=fshellsum(i,j)+GR(k,i,j)*dbin*rpos(k)*2.0*pi
enddo ! k
enddo ! j
enddo ! i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! save all to disk

do i = 1, types
do j = 1, types
write(300+i*types+j, *)step, fshellsum(i,j)
close(300+i*types+j)
write(400+i*types+j, *)step, fshellr(i,j)
close(400+i*types+j)
enddo ! j
enddo ! i

end



subroutine allocatebin
use system
use grsystem
use bo
real vect
integer i

dbin = min(r1a, r2a, r1b,r2b)/5

Ngrmax = int(cutoff/dbin)

ALLOCATE(GR(Ngrmax,types,types))
ALLOCATE (bo2co(Npart))
ALLOCATE (bo4co(Npart))
ALLOCATE (bo6co(Npart))
ALLOCATE(GR_tosend(Ngrmax,types,types))

ALLOCATE(bo4r(types,types, Ngrmax))
ALLOCATE(bo6r(types, types,Ngrmax))
ALLOCATE(bo2r(types, types, Ngrmax))
ALLOCATE(bo4rt(Ngrmax))
ALLOCATE(bo6rt(Ngrmax))
ALLOCATE(bo2rt(Ngrmax))
ALLOCATE(bo4r_tosend(types,types, Ngrmax))
ALLOCATE(bo6r_tosend(types, types, Ngrmax))
ALLOCATE(bo2r_tosend(types, types, Ngrmax))
ALLOCATE(bo4rt_tosend(Ngrmax))
ALLOCATE(bo6rt_tosend(Ngrmax))
ALLOCATE(bo2rt_tosend(Ngrmax))

ALLOCATE(bo2aux(Npart))
ALLOCATE(bo4aux(Npart))
ALLOCATE(bo6aux(Npart))

ALLOCATE (NVECINOS2(Npart))
ALLOCATE (NVECINOS2_tosend(Npart))
ALLOCATE (VECINOS2(Npart, VV))
ALLOCATE (VECINOS2_tosend(Npart, VV))
ALLOCATE (NVECINOS4(Npart))
ALLOCATE (NVECINOS4_tosend(Npart))
ALLOCATE (VECINOS4(Npart, VV))
ALLOCATE (VECINOS4_tosend(Npart, VV))
ALLOCATE (NVECINOS6(Npart))
ALLOCATE (NVECINOS6_tosend(Npart))
ALLOCATE (VECINOS6(Npart, VV))
ALLOCATE (VECINOS6_tosend(Npart, VV))



GR = 0
end

subroutine savebo(step)
use system
use grsystem
use params
use mpi
implicit none
real rpos(Ngrmax)

integer step
integer i, j,k
real vect
character (len=30) :: filename
integer maxsv

real fshellr(types,types), fshellsum(types,types)
integer fshell(types, types), peak(types, types)

complex ctemp4(types), ctemp6(types), ctemp2(types)
complex ctemp4t, ctemp6t, ctemp2t
complex bo4(types, types), bo6(types, types), bo2(types,types)
complex bo4t, bo6t, bo2t
complex cnum
integer NN(types), NNt

real, external :: angle
real, external :: dist
real aa
real bo4r(types,types, Ngrmax),bo6r(types, types, Ngrmax), bo2r(types, types, Ngrmax)
real bo4rt(Ngrmax),bo6rt(Ngrmax), bo2rt(Ngrmax)
real GRtemp(Ngrmax)
real firstpeak
integer l
integer iBOmax, iBOmin
real BOmax, BOmin
integer ii

BOmax = max(r1a,r2a,r1b,r2b)*3.0
BOmin = min(r1a,r2a,r1b,r2b)

iBOmax = int(BOmax/dbin) ! number of dbins 
iBOmin = int(BOmin/dbin) ! number of dbins 

bo4r = 0.0
bo6r = 0.0
bo2r = 0.0
bo4rt = 0.0
bo6rt = 0.0
bo2rt = 0.0

do ii = iBOmin, iBOmax
rpos(ii) = (float(ii)-0.5)*dbin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate BO

cnum = cmplx(0,1)
bo4 = cmplx(0,0)
bo6 = cmplx(0,0)
bo2 = cmplx(0,0)
bo4t = cmplx(0,0)
bo6t = cmplx(0,0)
bo2t = cmplx(0,0)

do k = 1, Npart

NN = 0
NNt = 0
ctemp4 = cmplx(0,0)
ctemp6 = cmplx(0,0)
ctemp2 = cmplx(0,0)
ctemp4t = cmplx(0,0)
ctemp6t = cmplx(0,0)
ctemp2t = cmplx(0,0)

do j = 1, Npart

vect = dist(k, j)

if((vect.lt.rpos(ii)).and.(k.ne.j)) then
      NN(tp(j)) = NN(tp(j)) + 1
      NNt = NNt + 1
      aa = angle(k,j)
      ctemp4(tp(j)) = ctemp4(tp(j)) + exp(4.0*cnum*aa)
      ctemp4t = ctemp4t + exp(4.0*cnum*aa)
      ctemp6(tp(j)) = ctemp6(tp(j)) + exp(6.0*cnum*aa)
      ctemp6t = ctemp6t + exp(6.0*cnum*aa)
      ctemp2(tp(j)) = ctemp2(tp(j)) + exp(2.0*cnum*aa)
      ctemp2t = ctemp2t + exp(2.0*cnum*aa)
endif
enddo ! j

do i = 1, types
      if(NN(i).ne.0) then
      bo2(tp(k),i) = bo2(tp(k),i) + ctemp2(i)/cmplx(float(NN(i)),0.0)
      bo4(tp(k),i) = bo4(tp(k),i) + ctemp4(i)/cmplx(float(NN(i)),0.0)
      bo6(tp(k),i) = bo6(tp(k),i) + ctemp6(i)/cmplx(float(NN(i)),0.0)
      endif
enddo ! i

if(NNt.ne.0) then
bo2t = bo2t + ctemp2t/cmplx(float(NNt),0.0)
bo4t = bo4t + ctemp4t/cmplx(float(NNt),0.0)
bo6t = bo6t + ctemp6t/cmplx(float(NNt),0.0)
endif

enddo ! k

bo2r(:,:,ii) = abs(bo2(:,:)/float(Npart))
bo2rt(ii) = abs(bo2t/float(Npart))
bo4r(:,:,ii) = abs(bo4(:,:)/float(Npart))
bo4rt(ii) = abs(bo4t/float(Npart))
bo6r(:,:,ii) = abs(bo6(:,:)/float(Npart))
bo6rt(ii) = abs(bo6t/float(Npart))
enddo ! ii


if(rank.eq.0) then
do i = 1, types
do j = i, types
write(filename,'(A4, I3.3, A1, I3.3, A1, I8.8, A4)')'BO4.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=40)
write(filename,'(A4, I3.3, A1, I3.3, A1, I8.8, A4)')'BO6.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=140)
write(filename,'(A4, I3.3, A1, I3.3, A1, I8.8, A4)')'BO2.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=240)
do ii = iBOmin, iBOmax
write(40,*)rpos(ii), bo4r(i, j, ii)
write(140,*)rpos(ii), bo6r(i, j, ii)
write(240,*)rpos(ii), bo2r(i, j, ii)
!print*, rpos(ii), bo4r(i, j, ii),i,j
!print*, rpos(ii), bo6r(i, j, ii),i,j
enddo ! ii
close(40)
close(140)
close(240)
enddo !i
enddo !j

write(filename,'(A5, I8.8, A4)')'BO4t.', step,'.dat'
open(file=filename,unit=40)
write(filename,'(A5, I8.8, A4)')'BO6t.', step,'.dat'
open(file=filename,unit=140)
write(filename,'(A5, I8.8, A4)')'BO2t.', step,'.dat'
open(file=filename,unit=240)

do ii = iBOmin, iBOmax
write(40,*)rpos(ii), bo4rt(ii)
write(140,*)rpos(ii), bo6rt(ii)
write(240,*)rpos(ii), bo6rt(ii)
!print*, rpos(ii), bo4r(i, j, ii),i,j
!print*, rpos(ii), bo6r(i, j, ii),i,j
enddo ! ii
close(40)
close(140)
close(240)


do i = 1, types
do j = 1, types
write(filename,'(A4, I3.3, A1, I3.3, A4)')'BO4.', i,'.',j,'.dat'
open(file=filename, unit=40, access='append')
write(40,*)step, maxval(bo4r(i, j, :))
close(40)

write(filename,'(A4, I3.3, A1, I3.3, A4)')'BO6.', i,'.',j,'.dat'
open(file=filename, unit=140, access='append')
write(140,*)step, maxval(bo6r(i, j, :))
close(140)

write(filename,'(A4, I3.3, A1, I3.3, A4)')'BO2.', i,'.',j,'.dat'
open(file=filename, unit=240, access='append')
write(240,*)step, maxval(bo2r(i, j, :))
close(240)

enddo
enddo

write(filename,'(A5, A3)')'BO4t.', 'dat'
open(file=filename, unit=40, access='append')
write(40,*)step, maxval(bo4rt(:))
close(40)

write(filename,'(A5, A3)')'BO6t.','dat'
open(file=filename, unit=140, access='append')
write(140,*)step, maxval(bo6rt(:))
close(140)

write(filename,'(A5, A3)')'BO2t.','dat'
open(file=filename, unit=240, access='append')
write(240,*)step, maxval(bo2rt(:))
close(240)

endif ! rank
end

subroutine psavebo(step)
use system
use grsystem
use params
use mpi
use bo
use nlist
implicit none
real rpos(Ngrmax)

integer step
integer i, j,k
real vect
character (len=30) :: filename
integer maxsv

real fshellr(types,types), fshellsum(types,types)
integer fshell(types, types), peak(types, types)

complex ctemp4(types), ctemp6(types), ctemp2(types)
complex ctemp4t, ctemp6t, ctemp2t
complex bo4(types, types), bo6(types, types), bo2(types,types)
complex bo4t, bo6t, bo2t
complex cnum
integer NN(types), NNt

real, external :: angle
real, external :: dist
real aa
real GRtemp(Ngrmax)
real firstpeak
integer l
integer iBOmax, iBOmin
real BOmax, BOmin
integer ii, nnn, ccc
integer icell, jcell

BOmax = cutoff
BOmin = min(r1a,r2a,r1b,r2b)

iBOmax = int(BOmax/dbin) ! number of dbins 
iBOmin = int(BOmin/dbin) ! number of dbins 

bo4r_tosend = 0.0
bo6r_tosend = 0.0
bo2r_tosend = 0.0
bo4rt_tosend = 0.0
bo6rt_tosend = 0.0
bo2rt_tosend = 0.0

bo4r = 0.0
bo6r = 0.0
bo2r = 0.0
bo4rt = 0.0
bo6rt = 0.0
bo2rt = 0.0

do ii = iBOmin, iBOmax
rpos(ii) = (float(ii)-0.5)*dbin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate BO

cnum = cmplx(0,1)
bo4 = cmplx(0,0)
bo6 = cmplx(0,0)
bo2 = cmplx(0,0)
bo4t = cmplx(0,0)
bo6t = cmplx(0,0)
bo2t = cmplx(0,0)


do ccc = 1, Nlistproc(rank+1)
jcell = listproc(rank+1,ccc)

  l = head(jcell) ! head of cell
  do while (l.ne.0) ! loop over all particles in jcell

  NN = 0
  NNt = 0
  ctemp4 = cmplx(0,0)
  ctemp6 = cmplx(0,0)
  ctemp2 = cmplx(0,0)
  ctemp4t = cmplx(0,0)
  ctemp6t = cmplx(0,0)
  ctemp2t = cmplx(0,0)

  do k = 1, 3**di ! loop over all neighbors of jcell

  icell = cellneighbor(jcell,k)
  i = head(icell) 

  do while (i.ne.0) ! loop over all particles in icell

vect = dist(l, i)
if((i.ne.l).and.(dist(l,i).lt.cutoff).and.(vect.lt.rpos(ii))) then
      NN(tp(i)) = NN(tp(i)) + 1
      NNt = NNt + 1
      aa = angle(l,i)
      ctemp4(tp(i)) = ctemp4(tp(i)) + exp(4.0*cnum*aa)
      ctemp4t = ctemp4t + exp(4.0*cnum*aa)
      ctemp6(tp(i)) = ctemp6(tp(i)) + exp(6.0*cnum*aa)
      ctemp6t = ctemp6t + exp(6.0*cnum*aa)
      ctemp2(tp(i)) = ctemp2(tp(i)) + exp(2.0*cnum*aa)
      ctemp2t = ctemp2t + exp(2.0*cnum*aa)
endif

  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

do nnn = 1, types
      if(NN(nnn).ne.0) then
      bo2(tp(l),nnn) = bo2(tp(l),nnn) + ctemp2(nnn)/cmplx(float(NN(nnn)),0.0)
      bo4(tp(l),nnn) = bo4(tp(l),nnn) + ctemp4(nnn)/cmplx(float(NN(nnn)),0.0)
      bo6(tp(l),nnn) = bo6(tp(l),nnn) + ctemp6(nnn)/cmplx(float(NN(nnn)),0.0)
      endif
enddo ! nnn

if(NNt.ne.0) then
bo2t = bo2t + ctemp2t/cmplx(float(NNt),0.0)
bo4t = bo4t + ctemp4t/cmplx(float(NNt),0.0)
bo6t = bo6t + ctemp6t/cmplx(float(NNt),0.0)
endif


l = list(l) ! next in line
enddo ! while l
enddo ! jcell

bo2r_tosend(:,:,ii) = abs(bo2(:,:)/float(Npart))
bo2rt_tosend(ii) = abs(bo2t/float(Npart))
bo4r_tosend(:,:,ii) = abs(bo4(:,:)/float(Npart))
bo4rt_tosend(ii) = abs(bo4t/float(Npart))
bo6r_tosend(:,:,ii) = abs(bo6(:,:)/float(Npart))
bo6rt_tosend(ii) = abs(bo6t/float(Npart))
enddo ! ii

call update_mpi_bo ! sum the energies from all procs 


if(rank.eq.0) then
if(saveallfiles.eq.1) then
do i = 1, types
do j = i, types
write(filename,'(A4, I3.3, A1, I3.3, A1, I8.8, A4)')'BO4.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=40)
write(filename,'(A4, I3.3, A1, I3.3, A1, I8.8, A4)')'BO6.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=140)
write(filename,'(A4, I3.3, A1, I3.3, A1, I8.8, A4)')'BO2.', i,'.',j,'.', step,'.dat'
open(file=filename,unit=240)
do ii = iBOmin, iBOmax
write(40,*)rpos(ii), bo4r(i, j, ii)
write(140,*)rpos(ii), bo6r(i, j, ii)
write(240,*)rpos(ii), bo2r(i, j, ii)
!print*, rpos(ii), bo4r(i, j, ii),i,j
!print*, rpos(ii), bo6r(i, j, ii),i,j
enddo ! ii
close(40)
close(140)
close(240)
enddo !i
enddo !j
endif

if(saveallfiles.eq.1) then
write(filename,'(A5, I8.8, A4)')'BO4t.', step,'.dat'
open(file=filename,unit=40)
write(filename,'(A5, I8.8, A4)')'BO6t.', step,'.dat'
open(file=filename,unit=140)
write(filename,'(A5, I8.8, A4)')'BO2t.', step,'.dat'
open(file=filename,unit=240)

do ii = iBOmin, iBOmax
write(40,*)rpos(ii), bo4rt(ii)
write(140,*)rpos(ii), bo6rt(ii)
write(240,*)rpos(ii), bo6rt(ii)
!print*, rpos(ii), bo4r(i, j, ii),i,j
!print*, rpos(ii), bo6r(i, j, ii),i,j
enddo ! ii
close(40)
close(140)
close(240)
endif

do i = 1, types
do j = 1, types
write(filename,'(A4, I3.3, A1, I3.3, A4)')'BO4.', i,'.',j,'.dat'
open(file=filename, unit=40, access='append')
write(40,*)step, maxval(bo4r(i, j, :))
close(40)

write(filename,'(A4, I3.3, A1, I3.3, A4)')'BO6.', i,'.',j,'.dat'
open(file=filename, unit=140, access='append')
write(140,*)step, maxval(bo6r(i, j, :))
close(140)

write(filename,'(A4, I3.3, A1, I3.3, A4)')'BO2.', i,'.',j,'.dat'
open(file=filename, unit=240, access='append')
write(240,*)step, maxval(bo2r(i, j, :))
close(240)

enddo
enddo

write(filename,'(A5, A3)')'BO4t.', 'dat'
open(file=filename, unit=40, access='append')
write(40,*)step, maxval(bo4rt(:))
close(40)

write(filename,'(A5, A3)')'BO6t.','dat'
open(file=filename, unit=140, access='append')
write(140,*)step, maxval(bo6rt(:))
close(140)

write(filename,'(A5, A3)')'BO2t.','dat'
open(file=filename, unit=240, access='append')
write(240,*)step, maxval(bo2rt(:))
close(240)

endif
end

subroutine paddgr
use system
use grsystem
use nlist
use mpi
implicit none

integer i, j, ccc
integer icell, jcell
real vect
real, external :: dist
integer k,l
GR_tosend = 0.0
GR = 0 ; !reset GR

do ccc = 1, Nlistproc(rank+1)
jcell = listproc(rank+1,ccc)

  l = head(jcell) ! head of cell
  do while (l.ne.0) ! loop over all particles in jcell

  do k = 1, 3**di ! loop over all neighbors of jcell

  icell = cellneighbor(jcell,k)
  i = head(icell) 

  do while (i.ne.0) ! loop over all particles in icell
 
vect = dist(l, i)
if(vect.lt.cutoff) then
GR_tosend(int(vect/dbin)+1,tp(l),tp(i)) = GR_tosend(int(vect/dbin)+1,tp(l),tp(i))+1.0
endif

  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

l = list(l) ! next in line
enddo ! while l
enddo ! jcell

call update_mpi_gr ! sum the energies from all procs 

end

subroutine clusterbo(step)
use system
use grsystem
use params
use mpi
use bo
use nlist
use mainm
implicit none
real rpos(Ngrmax)

integer step
integer i, j,k
real vect
character (len=30) :: filename
integer maxsv

real fshellr(types,types), fshellsum(types,types)
integer fshell(types, types), peak(types, types)

complex ctemp4(types), ctemp6(types), ctemp2(types)
complex ctemp4t, ctemp6t, ctemp2t
complex bo4(types, types), bo6(types, types), bo2(types,types)
complex bo4t, bo6t, bo2t
complex cnum
integer NN(types), NNt

real, external :: angle
real, external :: dist
real aa
real GRtemp(Ngrmax)
real firstpeak
integer l
integer iBOmax, iBOmin
real BOmax, BOmin
integer ii, nnn, ccc
integer icell, jcell
real ddd, ddd_cutoff
real F2, F4, F6


ddd_cutoff = 0.5

bo2aux=cmplx(0,0)
bo4aux=cmplx(0,0)
bo6aux=cmplx(0,0)

cnum=cmplx(0,1)

NVECINOS2_tosend = 0.0
NVECINOS2 = 0.0
VECINOS2_tosend = 0.0
VECINOS2 = 0.0
NVECINOS4_tosend = 0.0
NVECINOS4 = 0.0
VECINOS4_tosend = 0.0
VECINOS4 = 0.0
NVECINOS6_tosend = 0.0
NVECINOS6 = 0.0
VECINOS6_tosend = 0.0
VECINOS6 = 0.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate BO


do ccc = 1, Nlistproc(rank+1)
jcell = listproc(rank+1,ccc)

  l = head(jcell) ! head of cell
  do while (l.ne.0) ! loop over all particles in jcell

  NNt = 0

  do k = 1, 3**di ! loop over all neighbors of jcell

  icell = cellneighbor(jcell,k)
  i = head(icell) 

  do while (i.ne.0) ! loop over all particles in icell

  vect = dist(l, i)

if((i.ne.l).and.(vect.lt.BOcutoff).and.(tp(i).ne.tp(l))) then
      NNt = NNt + 1
      aa = angle(l,i)
      bo2aux(l) = bo2aux(l) + exp(2.0*cnum*aa)
      bo4aux(l) = bo4aux(l) + exp(4.0*cnum*aa)
      bo6aux(l) = bo6aux(l) + exp(6.0*cnum*aa)
endif

  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

bo2aux(l) = bo2aux(l)/cmplx(NNt,0)
bo4aux(l) = bo4aux(l)/cmplx(NNt,0)
bo6aux(l) = bo6aux(l)/cmplx(NNt,0)

l = list(l) ! next in line
enddo ! while l
enddo ! jcell

!!!! DETERMINE BOND

do ccc = 1, Nlistproc(rank+1)
jcell = listproc(rank+1,ccc)

  l = head(jcell) ! head of cell
  do while (l.ne.0) ! loop over all particles in jcell
  NNt = 0

  do k = 1, 3**di ! loop over all neighbors of jcell

  icell = cellneighbor(jcell,k)
  i = head(icell) 

  do while (i.ne.0) ! loop over all particles in icell

  vect = dist(l, i)

if((i.ne.l).and.(vect.lt.BOcutoff)) then

      ddd = bo2aux(l)*conjg(bo2aux(i))/abs(bo2aux(l))/abs(bo2aux(i)) ! bond between i and l... if ddd > 0.5 and type i != l then i and l are bonded
      if((ddd.gt.ddd_cutoff).and.(tp(i).ne.tp(l))) then
      NVECINOS2_tosend(l) = NVECINOS2_tosend(l) + 1
      if(NVECINOS2_tosend(l).eq.VV+1) then
      print*, 'NVECINOS2 too small. Stop'
      stop
      endif
      VECINOS2_tosend(l, NVECINOS2_tosend(l)) = i ! number of bonds of l
      endif

      ddd = bo4aux(l)*conjg(bo4aux(i))/abs(bo4aux(l))/abs(bo4aux(i)) ! bond between i and l... if ddd > 0.5 and type i != l then i and l are bonded
      if((ddd.gt.ddd_cutoff).and.(tp(i).ne.tp(l))) then
      NVECINOS4_tosend(l) = NVECINOS4_tosend(l) + 1
      if(NVECINOS4_tosend(l).eq.VV+1) then
      print*, 'NVECINOS4 too small. Stop'
      stop
      endif
      VECINOS4_tosend(l, NVECINOS4_tosend(l)) = i ! number of bonds of l
      endif

      ddd = bo6aux(l)*conjg(bo6aux(i))/abs(bo6aux(l))/abs(bo6aux(i)) ! bond between i and l... if ddd > 0.5 and type i != l then i and l are bonded
      if((ddd.gt.ddd_cutoff).and.(tp(i).ne.tp(l))) then
      NVECINOS6_tosend(l) = NVECINOS6_tosend(l) + 1
      if(NVECINOS6_tosend(l).eq.VV+1) then
      print*, 'NVECINOS6 too small. Stop'
      stop
      endif
      VECINOS6_tosend(l, NVECINOS6_tosend(l)) = i ! number of bonds of l
      endif
endif

  i = list(i) ! next in line
  enddo ! while i
  enddo ! k

l = list(l) ! next in line
enddo ! while l
enddo ! jcell

call update_mpi_bo_vecinos 

F2 = 0.0
F4 = 0.0
F6 = 0.0
do j = 1,Npart
if(NVECINOS2(j).gt.2)F2=F2+1.0
if(NVECINOS4(j).gt.2)F4=F4+1.0
if(NVECINOS6(j).gt.2)F6=F6+1.0
enddo
F2=F2/Npart
F4=F4/Npart
F6=F6/Npart

if(rank.eq.0) then
write(filename,'(A5, A3)')'BO4f.', 'dat'
open(file=filename, unit=40, access='append')
write(40,*)step, F4
!print*, 'Step:', step, 'BO4f', F4
close(40)

write(filename,'(A5, A3)')'BO6f.','dat'
open(file=filename, unit=140, access='append')
write(140,*)step, F6
close(140)

write(filename,'(A5, A3)')'BO2f.','dat'
open(file=filename, unit=240, access='append')
write(240,*)step, F2
close(240)

if(flagxyz.eq.1)call savecluster(step)

endif

end


subroutine savecluster(k)
use system
use bo
implicit none
integer k, j
character (len=3) :: atom
write(1010, *) Npart
write(1010, *) k
do j = 1, Npart

   if(tp(j).eq.1) then
   atom = ' C '
   else
   atom = ' O '
   endif

   write(1010, *)atom, xpos(j,1), xpos(j,2) , 0.0, 0.0
enddo
   write(1011, *)(NVECINOS4(j), j=1,Npart)
end

