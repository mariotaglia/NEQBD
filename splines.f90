

      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
! xai~Rs in order), and given the array y2a(1:n), which is the output from spline above,
!and given a value of x, this routine returns a cubic-spline interpolated value y.
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
! We will find the right place in the table by means of bisection.
!This is optimal if sequential calls to this routine are at random values of x. If sequential calls are in order, and closely
! spaced, one would do better to store previous values of
! klo and khi and test if they remain appropriate on the
! next call.
      khi=n
    1 if (khi-klo.gt.1) then
      k=(khi+klo)/2
      if(xa(k).gt.x)then
      khi=k
      else
      klo=k
      endif
      goto 1
      ENDIF ! klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated.
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END



      SUBROUTINE splines(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      REAL yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=500)

      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then !The lower boundary condition is set either to be natural
      y2(1)=0.
      u(1)=0.
      else !or else to have a specified first derivative.
      y2(1)=-0.5
      u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1 !This is the decomposition loop of the tridiagonal algorithm. y2 and u are used for temporary storage of the decomposed factors.  
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/p
      u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn.gt..99e30) then !The upper boundary condition is set either to be ~Snatural~T
      qn=0.
      un=0.
      else !or else to have a specified first derivative.
      qn=0.5
      un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1 ! This is the backsubstitution loop of the tridiagonal algorithm
      y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      return
      END


!********************************************************
!*          Akima spline fitting subroutine             *
!* ---------------------------------------------------- *
!* The input table is X(i), Y(i), where Y(i) is the     *
!* dependant variable. The interpolation point is xx,   *
!* which is assumed to be in the interval of the table  *
!* with at least one table value to the left, and three *
!* to the right. The interpolated returned value is yy. *
!* n is returned as an error check (n=0 implies error). *
!* It is also assumed that the X(i) are in ascending    *
!* order.                                               *
!********************************************************
Subroutine Interpol_Akima(iv,n,xx,yy,X,Y)  
  !Labels: 100,200,300
  integer i,iv,n
  real  xx,yy
  real  X (0:iv), Y (0:iv)
  real  XM(0:iv+3)
  real  Z (0:iv)

  n=1
  !Check to see if interpolation point is correct
  if (xx<X(1).or.xx>=X(iv-3)) then
    n=0 ; return
  end if
  X(0)=2.d0*X(1)-X(2)
  !Calculate Akima coefficients, a and b
  do i = 1, iv-1
    !Shift i to i+2
    XM(i+2)=(Y(i+1)-Y(i))/(X(i+1)-X(i))
  end do
  XM(iv+2)=2.d0*XM(iv+1)-XM(iv)
  XM(iv+3)=2.d0*XM(iv+2)-XM(iv+1)
  XM(2)=2.d0*XM(3)-XM(4)
  XM(1)=2.d0*XM(2)-XM(3)
  do i = 1, iv
    a=abs(XM(i+3)-XM(i+2))
    b=abs(XM(i+1)-XM(i))
    if (a+b.ne.0.d0) goto 100
    Z(i)=(XM(i+2)+XM(i+1))/2.d0
    goto 200
100 Z(i)=(a*XM(i+1)+b*XM(i+2))/(a+b)
200 end do
  !Find relevant table interval
  i=0
300 i=i+1
  if (xx>X(i)) goto 300
  i=i-1
  !Begin interpolation
  b=X(i+1)-X(i)
  a=xx-X(i)
  yy=Y(i)+Z(i)*a+(3.d0*XM(i+2)-2.d0*Z(i)-Z(i+1))*a*a/b
  yy=yy+(Z(i)+Z(i+1)-2.d0*XM(i+2))*a*a*a/(b*b)
  return
end

! End of file akima.f90
