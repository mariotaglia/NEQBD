
!****************************************************************
!**********************************************************************
        real FUNCTION RANDS (SEED)
!**********************************************************************

!-----  this is a special function for random number generation
!        on 32-bit machines that do not support long integer
!        multiplication and truncation.  the technique used is to do
!        the multiplication and addition in parts, by splitting all
!       integers in a 'high' and a 'low' part.  the algorithm is
!-----        exact, and should give machine-independent results.
!
!-----        the algorithm implemented is (following d.e. knuth):
!        seed = seed*1592653589 + 453816691
!        if (seed.lt.0) seed = seed + 1 + 2147483647
!-----        note that 1592653589 = 48603*2**15 + 30485
! 32768 = 2^15, 65536 = 2^16, 4.65661287308E-10 = 2^(-31)

INTEGER SEED, I1, I2

I2 = MOD (SEED, 32768) * 30485

I1 = MOD (SEED / 32768 * 30485, 65536) + MOD (MOD (SEED, 32768)* 48603, 65536) + I2 / 32768 + MOD (SEED / 32768, 2) * 32768 + 13849
I2 = MOD (I2, 32768) + 12659
I1 = I1 + I2 / 32768
I2 = MOD (I2, 32768)
I1 = MOD (I1, 65536)
SEED = I1 * 32768 + I2
RANDS = REAL(I1 * 32768 + I2) * 4.65661287308E-10
RETURN
END

real FUNCTION ran1(idum)
INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
REAL AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
END


real FUNCTION ran2(idum)
  INTEGER idum,idum2
  INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1)
  PARAMETER (IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791)
  PARAMETER (NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
  return
END FUNCTION ran2

subroutine randomgauss(mean, sd, rand1)
use random
use params
implicit none
REAL :: mean, sd
real, external :: RANDS
real, external :: ran1
real, external :: ran2
real rand1
real r1, r2

r1 = 0.0
r2 = 0.0

do while (r1.eq.0.0) 
r1 = ran2(seed)
enddo
do while (r2.eq.0.0)
r2 = ran2(seed)
enddo

! Now convert to normal distribution
rand1 = sd * SQRT(-2.0*LOG(r1)) * COS(2*pi*r2) + mean
end


