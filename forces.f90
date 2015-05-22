
!######################################################3
!
! After including a new force, be sure to update: 
! force, calcpres, penergy, plotpot and avepot
!



function LJforce(l, i)
use system
implicit none
real, external :: dist
real, external :: distk
real vect14, LJpre, vectkk
real, dimension(di) :: LJforce
integer kk
integer l, i

LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))*12.0
vect14 = (dist(l, i)**14)

do kk = 1, di
    vectkk = distk(l,i,kk) ! distance in one of the dimenstions only
    LJforce(kk) = LJpre*vectkk/vect14 ! force
enddo ! kk
endfunction


function LJenergy(l, i)
use system
implicit none
real, external :: dist
real LJ(types,types), vect12, LJpre
real :: LJenergy
integer l, i

LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))
vect12 = (dist(l, i)**12)
LJenergy = LJpre/vect12 - LJpre/cutoff12

endfunction

function LJEforce(l, i)
use system
implicit none
real, external :: dist
real, external :: distk
real vect14, LJpre, vectkk, vect
real, dimension(di) :: LJEforce
integer kk
integer l, i

LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))*12.0
vect = dist(l, i)
vect14 = (vect**14)

 do kk = 1, di
    vectkk = distk(l,i,kk) ! distance in one of the dimenstions only
    LJEforce(kk) = vectkk *(LJpre/vect14 - eexp*(1.0/elen/vect+1.0/vect/vect)*exp(-vect/elen)/vect)
enddo ! kk
endfunction


function LJEenergy(l, i)
use system
implicit none
real, external :: dist
real LJ(types,types), vect12, LJpre, vect
real :: LJEenergy
integer l, i

vect = dist(l, i)
LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))
vect12 = vect**12

LJEenergy = (LJpre/vect12 - eexp*(exp(-vect/elen))) - (LJpre/cutoff12-eexp*exp(-cutoff/elen))
endfunction



function Yforce(l, i)
use system
implicit none
real, external :: dist
real, external :: distk
real vect14, LJpre, vectkk, vect
real, dimension(di) :: Yforce
integer kk
integer l, i
real signo

LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))*12.0
vect = dist(l, i)
vect14 = (vect**14)

do kk = 1, di
vectkk = distk(l,i,kk) ! distance in one of the dimenstions only
Yforce(kk) = vectkk*(LJpre/vect14 + zint(tp(i))*zint(tp(l))*eexp*(1.0/elen/vect+1.0/vect/vect)*exp(-vect/elen)/vect)
enddo ! kk
endfunction


function Yenergy(l, i)
use system
implicit none
real, external :: dist
real LJ(types,types), vect12, LJpre, vect
real :: Yenergy
integer l, i
real signo

vect = dist(l, i)
LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))
vect12 = vect**12

Yenergy = LJpre/vect12 + zint(tp(i))*zint(tp(l))*eexp*exp(-vect/elen)/vect
Yenergy = Yenergy - (LJpre/cutoff12 + zint(tp(i))*zint(tp(l))*eexp*exp(-cutoff/elen)/cutoff)
endfunction

function LJNforce(l, i)
use system
implicit none
real, external :: dist
real, external :: distk
real vectN1, vectN2, LJpre1, LJpre2, vectkk
real, dimension(di) :: LJNforce
integer kk
integer l, i

    LJpre1 = rint(tp(i),tp(l))**(NNLJ)
    LJpre2 = rint(tp(i),tp(l))**(NLJ)

    LJpre1 = LJpre1*(NNLJ)*eint(tp(i),tp(l))
    LJpre2 = LJpre2*(NLJ)*eint(tp(i),tp(l))

    vectN1 = (dist(l, i)**(NNLJ+2.0))
    vectN2 = (dist(l, i)**(NLJ+2.0))

    do kk = 1, di
    vectkk = distk(l,i,kk) ! distance in one of the dimenstions only
    LJNforce(kk) = vectkk*(LJpre2/vectN2 - alfaLJ*LJpre1/vectN1) ! force
    enddo ! kk
endfunction

function LJNenergy(l, i)
use system
implicit none
real, external :: dist
real LJ(types,types), vectN1, LJpre1, vectN2, LJpre2
real :: LJNenergy
integer l, i
real cutoffN1, cutoffN2

LJpre1 = (rint(tp(i),tp(l))**NNLJ)
LJpre2 = (rint(tp(i),tp(l))**NLJ)

LJpre1 = LJpre1*eint(tp(i),tp(l))
LJpre2 = LJpre2*eint(tp(i),tp(l))

vectn1 = (dist(l, i)**nnlj)
vectn2 = (dist(l, i)**nlj)

cutoffN1 = (cutoff**nnlj)
cutoffN2 = (cutoff**nlj)

lJNenergy = LJpre2/vectN2 - alfaLJ*LJpre1/vectN1 
LJNenergy = LJNenergy - (LJpre2/cutoffN2 - alfaLJ*LJpre1/cutoffN1) 
endfunction

function tableenergy(l, i)
use system
implicit none
real, external :: dist
real :: tableenergy
integer l, i
real xspl, yspl
integer error
integer ipos
real y1,x1,m
integer tpl, tpi

tpl=tp(l)
tpi=tp(i)

xspl = dist(l,i)
ipos = int(xspl/cutoff*(maxn-1)) + 1

y1 = tabpotspl(tpl,tpi,ipos)
x1 = tabposspl(tpl,tpi,ipos)
m = tabpotm(tpl,tpi,ipos)

yspl = m*(xspl-x1) + y1
tableenergy = yspl
endfunction

function tableforce(l, i)
use system
implicit none
real, external :: dist
real, external :: distk
real vectN, LJpre, vectkk
real, dimension(di) :: tableforce
integer kk

integer l, i
real x(0:maxn), y(0:maxn)
real xspl, yspl
integer ipos
real y1,x1,m
integer error
integer tpl, tpi

tpl=tp(l)
tpi=tp(i)

xspl = dist(l,i)
ipos = int(xspl/cutoff*(maxn-1)) + 1

!DEBUG
!if(ipos.gt.maxn) then
!print*, xspl, cutoff
!stop
!endif

y1 = tabforcespl(tpl,tpi,ipos)
x1 = tabposspl(tpl,tpi,ipos)
m = tabforcem(tpl,tpi,ipos)

! DEBUG
!print*, ipos, x1, y1, m

yspl = m*(xspl-x1) + y1

do kk = 1, di
vectkk = distk(l,i,kk) ! distance in one of the dimenstions only
tableforce(kk) = yspl*vectkk/xspl ! force
enddo ! kk
endfunction


function YAVforce(l, i)
use system
implicit none
real, external :: dist
real, external :: distk
real vect14, LJpre, vectkk, vect
real, dimension(di) :: YAVforce
integer kk
integer l, i
real signo

LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))*12.0
vect = dist(l, i)
vect14 = (vect**14)

do kk = 1, di
vectkk = distk(l,i,kk) ! distance in one of the dimenstions only
YAVforce(kk) = vectkk*(LJpre/vect14 + zint(tp(i))*zint(tp(l))* &
eexp*eexpav(tp(i),tp(l))*(1.0/elen/vect+1.0/vect/vect)*exp(-vect/elen)/vect)
enddo ! kk
endfunction


function YAVenergy(l, i)
use system
implicit none
real, external :: dist
real LJ(types,types), vect12, LJpre, vect
real :: YAVenergy
integer l, i
real signo

vect = dist(l, i)
LJpre = (rint(tp(i),tp(l))**12)*eint(tp(i),tp(l))
vect12 = vect**12

YAVenergy = LJpre/vect12 + zint(tp(i))*zint(tp(l))*eexp*exp(-vect/elen)/vect
YAVenergy = YAVenergy - (LJpre/cutoff12 + zint(tp(i))*zint(tp(l))*eexp* &
eexpav(tp(i),tp(l))*exp(-cutoff/elen)/cutoff)
endfunction


