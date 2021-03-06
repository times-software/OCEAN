! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
  open( unit=99, file='sphpts', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nsphpt
  allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
  do isphpt = 1, nsphpt
     read ( 99, * ) xsph( isphpt ), ysph( isphpt ), zsph( isphpt ), wsph( isphpt )
  end do
  close( unit=99 )
  sphsu = sum( wsph( : ) )
  wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / sphsu )
  write ( 6, * ) nsphpt, ' points with weights summing to four pi '
