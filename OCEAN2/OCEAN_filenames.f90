! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_filenames

  implicit none
  private
  save

  public :: OCEAN_filenames_spectrum, OCEAN_filenames_lanc, OCEAN_filenames_energy, &
            OCEAN_filenames_ehamp, OCEAN_filenames_read_ehamp
  

  contains

  ! The file name used to store the spectrum file for a run
  subroutine OCEAN_filenames_spectrum( sys, filename, ierr )
    use OCEAN_system, only : o_system
    !
    type( o_system ), intent( in ) :: sys
    character(len=*), intent( out ) :: filename
    integer, intent( inout ) :: ierr
    !
    if( len( filename ) < 25 ) then
      ierr = 5
      return
    endif
    !

    select case ( sys%cur_run%calc_type)

      case( 'XES' )
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xesspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon

      case( 'XAS' )
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon

      case( 'RXS')
        write(filename,'(A8,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxsspct_', sys%cur_run%elname, &
            '.', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol

      case ('C2C')
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'ctcspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol

      case( 'VAL' )
        write(filename, '(A)' ) 'opcons'

      case default
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'absspct_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon

    end select
    !
  end subroutine OCEAN_filenames_spectrum

  subroutine OCEAN_filenames_lanc( sys, filename, ierr )
    use OCEAN_system, only : o_system
    !
    type( o_system ), intent( in ) :: sys
    character(len=*), intent( out ) :: filename
    integer, intent( inout ) :: ierr
    !
    if( len( filename ) < 24 ) then
      ierr = 5
      return
    endif
    !

    select case ( sys%cur_run%calc_type)

      case( 'XES' )
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'xeslanc_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon

      case( 'XAS' )
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'abslanc_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon

      case( 'RXS')
        write(filename,'(A7,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'rxlanc_', sys%cur_run%elname, &
            '.', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol

      case( 'C2C' )
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2)' ) 'ctclanc_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol

      case( 'VAL' )
        write(filename, '(A)' ) 'oplanc'

      case default
        write(filename,'(A8,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'abslanc_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon

    end select
    !
  end subroutine OCEAN_filenames_lanc

  ! The Electron--(core)-Hole energy files 
  subroutine OCEAN_filenames_energy( sys, filename, iter, hflag, ierr)

    use OCEAN_system, only : o_system
    !
    type( o_system ), intent( in ) :: sys
    character(len=*), intent( out ) :: filename
    integer, intent( in ) :: iter
    integer, intent( inout ) :: ierr
    integer, intent( in ) :: hflag(6)
    character(len=6) :: hamname
!    character(len=3000) :: hamname
!    integer :: stringLen
!    integer :: last, actual
!    write(hamname,*) hflag
!    stringLen = len(hamname)
!    last = 1
!    actual = 1
!    do while (actual < stringLen)
!    if( hamname(last:last)==' ') then
!    actual = actual + 1
!    hamname(last:last)=hamname(actual:actual)
!    hamname(actual:actual)=' '
!    else
!    last = last + 1
!    if( actual < last) &
!    actual = last
!    endif 
!    end do

    write( hamname, '(6(I1.1))' ) hflag(:)

    select case( sys%cur_run%calc_type )

      case( 'RXS' )
        if( len( filename ) < 35 ) then
          ierr = 6
          return
        endif
        write(filename, '(A6,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2,A1,I4.4,A1,A6)' ) 'ehamp_', & 
            sys%cur_run%elname, '.', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol, '.',&
            iter,'.',trim(hamname)

      case( 'VAL' )
        if( len( filename ) < 17 ) then
          ierr = 6
          return
        endif
        write(filename, '(A6,I4.4,A1,A6)' ) 'ehamp_', iter, '.', hamname
        
      case( 'XAS', 'XES' )
        if( len( filename ) < 32 ) then
          ierr = 6
          return
        endif
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4,A1,A6)' ) 'echamp_', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.',&
                 iter,'.',trim(hamname)

      case default ! Currently XAS/XES option
        if( len( filename ) < 32 ) then
          ierr = 6
          return
        endif
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4,A1,A6)' ) 'echamp_', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.',&
                 iter,'.',trim(hamname)

    end select
    !
  end subroutine OCEAN_filenames_energy

  ! The Electron--(Core)-Hole AMPlitude files 
  subroutine OCEAN_filenames_ehamp( sys, filename, iter, ierr )
    use OCEAN_system, only : o_system
    !
    type( o_system ), intent( in ) :: sys
    character(len=*), intent( out ) :: filename
    integer, intent( in ) :: iter
    integer, intent( inout ) :: ierr
    !
    select case( sys%cur_run%calc_type )

      case( 'RXS' )
        if( len( filename ) < 28 ) then
          ierr = 6
          return
        endif
        write(filename, '(A6,A2,A1,A2,A1,I2.2,A1,I5.5,A1,I2.2,A1,I4.4)' ) 'ehamp_', & 
            sys%cur_run%elname, '.', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', &
            sys%cur_run%rixs_energy, '.', sys%cur_run%rixs_pol, '.', iter

      case( 'VAL' )
        if( len( filename ) < 10 ) then
          ierr = 6
          return
        endif
        write(filename, '(A6,I4.4)' ) 'ehamp_', iter
        
      case( 'XAS', 'XES' )
        if( len( filename ) < 25 ) then
          ierr = 6
          return
        endif
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', iter

      case default ! Currently XAS/XES option
        if( len( filename ) < 25 ) then
          ierr = 6
          return
        endif
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', iter

    end select
    !
  end subroutine OCEAN_filenames_ehamp


  ! The Electron--(Core)-Hole AMPlitude files 
  subroutine OCEAN_filenames_read_ehamp( sys, filename, iter, ierr )
    use OCEAN_system, only : o_system
    !
    type( o_system ), intent( in ) :: sys
    character(len=*), intent( out ) :: filename
    integer, intent( in ) :: iter
    integer, intent( inout ) :: ierr
    !
    select case( sys%cur_run%calc_type ) 

      case( 'RXS' ) 
        if( len( filename ) < 25 ) then
          ierr = 6 
          return
        endif
        write(6,*) 'RXS!'
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', iter, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', sys%cur_run%rixs_energy
          
      case( 'C2C' )
        if( len( filename ) < 25 ) then
          ierr = 6 
          return
        endif
        write(6,*) 'C2C!'
        write( 6, * ) sys%cur_run%rixsInputCoreLevel
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, & 
              '.', iter, '_', sys%cur_run%rixsInputCoreLevel, '_', sys%cur_run%photon, '.', sys%cur_run%rixs_energy

      case default ! Currently XAS/XES option
        if( len( filename ) < 25 ) then
          ierr = 6 
          return
        endif
        write(6,*) 'DEFAULT!'
        write(filename,'(A7,A2,A1,I4.4,A1,A2,A1,I2.2,A1,I4.4)' ) 'echamp_', sys%cur_run%elname, &
              '.', iter, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon, '.', sys%cur_run%rixs_energy

    end select
    !
  end subroutine OCEAN_filenames_read_ehamp


end module OCEAN_filenames
