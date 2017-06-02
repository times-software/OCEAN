c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine bachelet(vi,r,njrc,nr)
      implicit double precision (a-h,o-z)
      integer njrc( 4 )
      double precision vi(nr,7),r(nr)
      double precision c(6),aa(6),q(6,6),s(6,6),a(3)
      double precision, allocatable :: vcore( : )
      allocate( vcore( nr ) )
      write (6,*) 'PLEASE ENTER ZV,A1,A2,C1,C2 FOR CORE POTENTIAL.'
      read (5,*) zv,a1,a2,c1,c2
      ra1=dsqrt(a1)
      ra2=dsqrt(a2)
      do 25 i=1,nr
      vcore (i)=-zv/r(i)*(c1*errfunc(ra1*r(i))
     &                   +c2*errfunc(ra2*r(i)))
 25   continue
      write (6,*) 'YOU MAY NOW ENTER WHAT YOU WANT TO ENTER.'
      write (6,*) '1=S,2=PSO,3=PAV,4=DSO,5=DAV,6=FSO,7=FAV,0=STOP.'
 35   read (5,*) nop
      if (nop.eq.0) then 
         deallocate( vcore )
         return
      end if
      write (6,*) 'PLEASE ENTER A1,A2,A3,C1,C2,C3,C4,C5,C6.'
      read (5,*) a(1),a(2),a(3),c(1),c(2),c(3),c(4),c(5),c(6)
      pi=4.d0*datan(1.d0)
      rpi=dsqrt(pi)
      do 42 i=1,3
        do 40 k=1,3
          s(i  ,k  )=0.2500d0*rpi/(a(i)+a(k))**1.5d0
          s(i+3,k  )=0.3750d0*rpi/(a(i)+a(k))**2.5d0
          s(i  ,k+3)=0.3750d0*rpi/(a(i)+a(k))**2.5d0
          s(i+3,k+3)=0.9375d0*rpi/(a(i)+a(k))**3.5d0
 40     continue
 42   continue
      do 100  l=1  ,6
        do 45   i=l+1,6
          q(i,l)=0.d0
 45     continue
        do 75 i=1,l-1
          q(i,l)=s(i,l)
          do 55 k=1,i-1
            q(i,l)=q(i,l)-q(k,i)*q(k,l)
 55       continue
          q(i,l)=q(i,l)/q(i,i)
 75     continue
        q(l,l)=s(l,l)
        do 85 k=1,l-1
          q(l,l)=q(l,l)-q(k,i)*q(k,i)
 85     continue
        q(l,l)=dsqrt(q(l,l))
 100  continue
      call patt(q)
      do 120 i=1,6
        aa(i)=0.d0
        do 110 l=1,6
          aa(i)=aa(i)+q(i,l)*c(l)
 110    continue
 120  continue
      if (2*(nop/2).eq.nop) then
        do 235 i=1,nr
          rr=r(i)*r(i)
          vi(i,nop)=-(aa(1)+rr*aa(4))*dexp(-a(1)*rr)
     &              -(aa(2)+rr*aa(5))*dexp(-a(2)*rr)
     &              -(aa(3)+rr*aa(6))*dexp(-a(3)*rr)
 235    continue
      else
        ii=(nop+1)/2
        njrc(ii)=1
        do 240 i=1,nr
          rr=r(i)*r(i)
          vi(i,nop)=-(aa(1)+rr*aa(4))*dexp(-a(1)*rr)
     &              -(aa(2)+rr*aa(5))*dexp(-a(2)*rr)
     &              -(aa(3)+rr*aa(6))*dexp(-a(3)*rr)
     &              +vcore(i)
 240    continue
      endif
      goto 35
      end
C------------------------------------------------------------------------------
C
      SUBROUTINE PATT(Q)
      DOUBLE PRECISION Q(6,6),QI(6,6)
      DOUBLE PRECISION DELTA(6,6)
C
      ICHK=0
C
      DO 100 I=1,6
      DO 100 J=1,6
      QI(I,J)=0.D0
  100 CONTINUE
C
      QI(1,1)=1.D0/Q(1,1)
      QI(1,2)=-Q(1,2)/(Q(1,1)*Q(2,2))
      QI(1,3)=(Q(1,2)*Q(2,3)-Q(1,3)*Q(2,2))/(Q(1,1)*Q(2,2)*Q(3,3))
      QI(1,4)=(-Q(1,2)*Q(2,3)*Q(3,4)+Q(1,2)*Q(2,4)*Q(3,3)+
     & Q(1,3)*Q(2,2)*Q(3,4)-Q(1,4)*Q(2,2)*Q(3,3))/
     & (Q(1,1)*Q(2,2)*Q(3,3)*Q(4,4))
      QI(1,5)=(Q(1,2)*Q(2,3)*Q(3,4)*Q(4,5)-Q(1,2)*Q(2,3)*Q(3,5)*Q(4,4)-
     & Q(1,2)*Q(2,4)*Q(3,3)*Q(4,5)+Q(1,2)*Q(2,5)*Q(3,3)*Q(4,4)-
     & Q(1,3)*Q(2,2)*Q(3,4)*Q(4,5)+Q(1,3)*Q(2,2)*Q(3,5)*Q(4,4)+
     & Q(1,4)*Q(2,2)*Q(3,3)*Q(4,5)-Q(1,5)*Q(2,2)*Q(3,3)*Q(4,4))/
     & (Q(1,1)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5))
      QI(1,6)=(-Q(1,2)*Q(2,3)*Q(3,4)*Q(4,5)*Q(5,6)+
     & Q(1,2)*Q(2,3)*Q(3,4)*Q(4,6)*Q(5,5)+
     & Q(1,2)*Q(2,3)*Q(3,5)*Q(4,4)*Q(5,6)-
     & Q(1,2)*Q(2,3)*Q(3,6)*Q(4,4)*Q(5,5)+
     & Q(1,2)*Q(2,4)*Q(3,3)*Q(4,5)*Q(5,6)-
     & Q(1,2)*Q(2,4)*Q(3,3)*Q(4,6)*Q(5,5)-
     & Q(1,2)*Q(2,5)*Q(3,3)*Q(4,4)*Q(5,6)+
     & Q(1,2)*Q(2,6)*Q(3,3)*Q(4,4)*Q(5,5)+
     & Q(1,3)*Q(2,2)*Q(3,4)*Q(4,5)*Q(5,6)-
     & Q(1,3)*Q(2,2)*Q(3,4)*Q(4,6)*Q(5,5)-
     & Q(1,3)*Q(2,2)*Q(3,5)*Q(4,4)*Q(5,6)+
     & Q(1,3)*Q(2,2)*Q(3,6)*Q(4,4)*Q(5,5)-
     & Q(1,4)*Q(2,2)*Q(3,3)*Q(4,5)*Q(5,6)+
     & Q(1,4)*Q(2,2)*Q(3,3)*Q(4,6)*Q(5,5)+
     & Q(1,5)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,6)-
     & Q(1,6)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5))/
     & (Q(1,1)*Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5)*Q(6,6))
      QI(2,2)=1.D0/Q(2,2)
      QI(2,3)=-Q(2,3)/(Q(2,2)*Q(3,3))
      QI(2,4)=(Q(2,3)*Q(3,4)-Q(2,4)*Q(3,3))/(Q(2,2)*Q(3,3)*Q(4,4))
      QI(2,5)=(-Q(2,3)*Q(3,4)*Q(4,5)+Q(2,3)*Q(3,5)*Q(4,4)+
     & Q(2,4)*Q(3,3)*Q(4,5)-Q(2,5)*Q(3,3)*Q(4,4))/
     & (Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5))
      QI(2,6)=(Q(2,3)*Q(3,4)*Q(4,5)*Q(5,6)-Q(2,3)*Q(3,4)*Q(4,6)*Q(5,5)-
     & Q(2,3)*Q(3,5)*Q(4,4)*Q(5,6)+Q(2,3)*Q(3,6)*Q(4,4)*Q(5,5)-
     & Q(2,4)*Q(3,3)*Q(4,5)*Q(5,6)+Q(2,4)*Q(3,3)*Q(4,6)*Q(5,5)+
     & Q(2,5)*Q(3,3)*Q(4,4)*Q(5,6)-Q(2,6)*Q(3,3)*Q(4,4)*Q(5,5))/
     & (Q(2,2)*Q(3,3)*Q(4,4)*Q(5,5)*Q(6,6))
      QI(3,3)=1.D0/Q(3,3)
      QI(3,4)=-Q(3,4)/(Q(3,3)*Q(4,4))
      QI(3,5)=(Q(3,4)*Q(4,5)-Q(3,5)*Q(4,4))/(Q(3,3)*Q(4,4)*Q(5,5))
      QI(3,6)=(-Q(3,4)*Q(4,5)*Q(5,6)+Q(3,4)*Q(4,6)*Q(5,5)+
     & Q(3,5)*Q(4,4)*Q(5,6)-Q(3,6)*Q(4,4)*Q(5,5))/
     & (Q(3,3)*Q(4,4)*Q(5,5)*Q(6,6))
      QI(4,4)=1.D0/Q(4,4)
      QI(4,5)=-Q(4,5)/(Q(4,4)*Q(5,5))
      QI(4,6)=(Q(4,5)*Q(5,6)-Q(4,6)*Q(5,5))/(Q(4,4)*Q(5,5)*Q(6,6))
      QI(5,5)=1.D0/Q(5,5)
      QI(5,6)=-Q(5,6)/(Q(5,5)*Q(6,6))
      QI(6,6)=1.D0/Q(6,6)
C
C CHECK INVERSE
C
      WRITE(9,1000)
 1000 FORMAT(/4X,'QUALITY OF ANALYTIC INVERSION BY PATTNAIK ET AL'/)
      DO 150 I=1,6
      DO 200 J=1,6
      DELTA(I,J)=0.D0
      IF (I.EQ.J) DELTA(I,J)=-1.D0
C
      DO 300 ISUM=1,6
      DELTA(I,J)=DELTA(I,J)+QI(I,ISUM)*Q(ISUM,J)
  300 CONTINUE
C
  200 CONTINUE
C
      WRITE(9,2000) (DELTA(I,K),K=1,6)
 2000 FORMAT(4X,6D18.7)
C
  150 CONTINUE
C
C TRANSFER INVERSE MATRIX ONTO ORIGINAL MATRIX
C
      DO 500 L=1,6
      DO 500 M=1,6
      Q(L,M)=QI(L,M)
  500 CONTINUE
C
      RETURN
      END
