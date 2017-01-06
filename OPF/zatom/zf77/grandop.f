c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine grandop(nwgt,wgts,mode)
      implicit none
      integer nwgt,iz,i,j,nel,ibig,ilit
      double precision zorig,xz
      integer no(1000),nl(1000),is(1000)
      double precision xnj(1000),occ(1000),ev(1000),wgts(nwgt)
      double precision other(5),tmp(1000)
      character*1 dumcha,first,second
      character*3 mode,fn3,front
      character*4 lets,fn4
      character*10 digs
      logical four
      open(unit=77,file='tmp',form='formatted',status='unknown')
      rewind 77
      do j=1,5
        other(j)=0.d0
      end do
      do i=1,nwgt
        read (77,*) zorig,nel
        iz=zorig
        read (77,*) (tmp(j),j=1,5)
        do j=1,5
          other(j)=other(j)+wgts(i)*tmp(j)
        end do
        read (77,*) (no(j),nl(j),xnj(j),is(j),occ(j),tmp(j),j=1,nel)
        do j=1,nel
          ev(j)=ev(j)+wgts(i)*tmp(j)
        end do
      end do
      close(unit=77)
      xz=0.d0
      do j=1,nel
        xz=xz+occ(j)
      end do
      ibig=iz/10
      ilit=iz-10*ibig
      ibig=ibig+1
      ilit=ilit+1
      digs='0123456789'
      lets='spdf'
      fn4(1:1)=digs(ibig:ibig)
      fn4(2:2)=digs(ilit:ilit)
      fn3(1:1)=digs(ibig:ibig)
      fn3(2:2)=digs(ilit:ilit)
      open(unit=98,file='occinfo',form='formatted',status='unknown')
      rewind 98
      do j=1,iz-1
        read (98,'(1a1)') dumcha
      end do
      read (98,'(8x,2a1)') first,second
      close(unit=98)
      if (second.eq.' ') then
        four=.false.
      else
        four=.true.
      end if
      fn3(3:3)=first
      fn4(3:3)=first
      fn4(4:4)=second
      if (four) then
        open(unit=99,file=fn4,form='formatted',status='unknown')
        write (6,*) fn4
      else
        open(unit=99,file=fn3,form='formatted',status='unknown')
        write (6,*) fn3
      end if
      rewind 99
      write (99,'(1a7,1f16.6)') 'Etot  =',other(1)
      write (99,'(1a7,1f16.6)') 'Ekin  =',other(2)
      write (99,'(1a7,1f16.6)') 'Ecoul =',other(3)
      write (99,'(1a7,1f16.6)') 'Eenuc =',other(4)
      write (99,'(1a7,1f16.6)') 'Exc   =',other(5)
      do j=1,nel
        front(1:1)=digs(no(j)+1:no(j)+1)
        front(2:2)=lets(nl(j)+1:nl(j)+1)
        if ((mode.eq.'lda').or.(mode.eq.'scr')) then
          front(3:3)=' '
        else
          if (mode.eq.'lsd') then
            if (is(j).eq.1) then
              front(3:3)='D'
            else
              front(3:3)='u'
            end if
          else
            if (dabs(xnj(j)).gt.dble(nl(j))) then
              front(3:3)='P'
            else
              front(3:3)='M'
            end if
          end if
        end if
        write (99,'(1a3,1f17.6)') front,ev(j)
      end do
      close(unit=99)
      return
      end
