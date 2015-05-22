
real function dist(l, i)
use system
implicit none
integer i, j, l
real vect
real DX(di)

vect = 0;
do j = 1, di
    DX(j)=xpos(l, j)-xpos(i, j)
    DX(j)=DX(j)-nint(DX(j)/xlim)*xlim
    vect = vect + DX(j)**2
enddo
dist = vect**(0.5)
end


real function dist2(l, i)
use system
implicit none
integer i, j, l
real vect
real DX(di)
vect = 0;
do j = 1, di
    DX(j)=xpos(l, j)-xpos(i, j)
    DX(j)=DX(j)-nint(DX(j)/xlim)*xlim
    vect = vect + DX(j)**2
enddo
dist2 = vect
end


real function angle(l, i)
use system
implicit none
integer i,  l
real DX, DY

DX=xpos(l, 1)-xpos(i, 1)
DX=DX-nint(DX/xlim)*xlim

DY=xpos(l, 2)-xpos(i, 2)
DY=DY-nint(DX/xlim)*xlim

angle = atan2(DY, DX)
end



real function distk(l, i, j)
use system
implicit none
integer i, j, l
    distk=xpos(l, j)-xpos(i, j)
    distk=distk-nint(distk/xlim)*xlim
end

