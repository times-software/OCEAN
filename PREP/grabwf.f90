      subroutine grabwf(filenum, maxband, maxnpw, kg_unshift, kg_shift, &
     &  eigen_un, eigen_sh, occ_un, occ_sh, cg_un, cg_sh, occ_max,      &
     &  unocc_max, nband, un_npw, sh_npw, noshift)

      integer :: filenum, maxband, maxnpw, nband(2),                    &
     &  iband,ii,un_npw,sh_npw,nspinor, occ_max, unocc_max
      double precision :: eigen_un(maxband), eigen_sh(maxband),         &
     &  occ_un(maxband), occ_sh(maxband), cg_un(maxband,2*maxnpw),      &
     &  cg_sh(maxband,2*maxnpw)
      integer :: kg_unshift(3,maxnpw), kg_shift(3,maxnpw)
      logical :: noshift

     
      read(filenum) un_npw,nspinor,nband(1) 
    
      read(filenum) kg_unshift(1:3,1:un_npw)
   
      read(filenum) eigen_un(1:nband(1)),occ_un(1:nband(1))
  
      do iband=1,nband(1)
       read(filenum) (cg_un(iband,ii),ii=1,2*un_npw)
      enddo
 
!      if ( nband(1) .lt. occ_max ) then
!       stop 'occ_max is less than nbands from abinit run'
!      endif
      if (.not. noshift) then
      read(filenum) sh_npw,nspinor,nband(2)

!     write(6,*) sh_npw,nspinor,nband(ikpt*2)
      read(filenum) kg_shift(1:3,1:sh_npw)

      read(filenum) eigen_sh(1:nband(2)),occ_sh(1:nband(2))

      do iband=1,nband(2)
       read(filenum) (cg_sh(iband,ii),ii=1,2*sh_npw)
      enddo

      if ( nband(2) .lt. unocc_max ) then
       stop 'unocc_max is less than nbands from abinit run'
      endif
      endif

      end 
      
