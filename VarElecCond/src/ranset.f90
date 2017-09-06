      SUBROUTINE RANSET(id,A,N,AMP)
!  Sets A to random values in the range +- AMP
      implicit none
      integer ISEED,N,id,I
      double precision A(N),amp, ran0

      ISEED = 1996 + 64*id
      DO I = 1,N
        A(I) = 2.0d0*(RAN0(ISEED) - 0.5d0)*AMP
      end do
      RETURN
      END
!
!  Numerical Recipes Function ran0
!
      double precision function ran0(idum)
      INTEGER idum,IA,IM,IQ,IR,MASK
      double precision  AM
      PARAMETER (IA=16807,IM=2147483647,AM=1.0d0/IM, &
                 IQ=127773,IR=2836,MASK=123459876)
      INTEGER k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ran0=AM*idum
      idum=ieor(idum,MASK)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software =$j!]Y'1,).
