! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine dft(dl,rel,al,nr,ne,l,j,s,o,dr,r2,cq,ph,or,ra,et,x137)
  implicit none
  integer nr,ne,l(ne),s(ne),i,k,ii
  double precision j(ne),o(ne),r2(nr),dr(nr),cq(nr)
  double precision ph(nr,ne),or(nr,ne)
  double precision rel,al,ra,et,dl,den,occ,pr,xn,fx,fc,bfac, su
  double precision xn1,ux1,uc1,uxc1,xn2,ux2,uc2,uxc2,nex,ec,x137
  !
  pr=0.0001d0
  !
  if (al .gt. 0.5d0) then
     fx = 1
     fc = 1
  else
     fx = 0
     fc = 1
  end if
  !
  su = 0.0d0
  ii=0
  open( unit=99, file='vxcofr', form='formatted', status='unknown' )
  rewind 99
  do i=1,nr
     !
     xn =0.d0
     xn1=0.d0
     xn2=0.d0
     do k=1,ne
        occ=o(k)
        den=occ*ph(i,k)*ph(i,k)
!       write ( 99, * ) k, o( k ), ph( i, k ), den, xn
        if ((j(k).lt.-pr).or.(occ.gt.dble(2*l(k)+1)+pr)) then
           xn=xn+den
        else
           if (s(k).eq.1) then
              xn1=xn1+den
           else
              xn2=xn2+den
           end if
        end if
     end do
     su = su + dr( i ) * xn
     xn=xn+cq(i)
!    write ( 99, * ) cq( i ), xn
     xn1=xn1+0.5d0*xn
     xn2=xn2+0.5d0*xn
     if ((xn1+xn2)/r2(i).lt.1.d-30) then
        nex=0.d0
        ec=0.d0
        ux1=0.d0 
        ux2=0.d0 
        uc1=0.d0 
        uc2=0.d0 
     else
        call exchcorr(rel,r2(i),xn1,xn2,nex,ec,ux1,ux2,uc1,uc2,x137)
     end if
     if (ii.eq.0) bfac=14.d0/45.d0
     if ((ii.eq.1).or.(ii.eq.3)) bfac=64.d0/45.d0
     if (ii.eq.2) bfac=24.d0/45.d0
     if (ii.eq.4) bfac=28.d0/45.d0
     et=et+dl*dsqrt(r2(i))*bfac*(fc*ec*(xn1+xn2)+fx*nex)
     uxc1 = ra * ( fx*ux1 + fc*uc1 )
     uxc2 = ra * ( fx*ux2 + fc*uc2 )
     write ( 99, '(6(1x,1e15.8))' ) uxc1, uxc2, fc * ec * ( xn1 + xn2 ) + fx * nex, xn1 + xn2
     !
     do k=1,ne
        occ=o(k)
        if ((j(k).lt.-pr).or.(occ.gt.dble(2*l(k)+1)+pr)) then
           or(i,k)=or(i,k)+0.5d0*(uxc1+uxc2)
        else
           if (s(k).eq.1) then
              or(i,k)=or(i,k)+uxc1
           end if
           if (s(k).eq.2) then
              or(i,k)=or(i,k)+uxc2
           end if
        end if
     end do
     !
     ii=ii+1
     if (ii.eq.5) ii=1
  end do
! write ( 6, * ) 'total density = ', su
  !
  close( unit=99 )
  !
  return
end subroutine dft
