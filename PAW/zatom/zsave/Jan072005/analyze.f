      program analyze
      implicit none
c
c  ntot = nfit + ncheck = no. of partial waves for a given l
c  nfit = number of waves used to fit A_i coefficients for
c
c         \int_{0}^{\infinity} dr \phi_c(r) O(r) \psi_{lm,k}(r) 
c            \approx \sum_{i} A_i \psi_{lm,k}^{Pseudo} ( r_i )
c
c  nchek = no.  of partial waves used to test A_i
c          (besides those used to fit)
c
c  nrm = # of values of r_i ( and A_i )
c  nr = total # of r's at which partial waves are provided
c  lval = value of l ...
c
      integer nfit, nchek, nr, ntot, nrm, lval
c
c  mel, energy = matrix element, energy for each partial wave
c
c  fval = table of \psi_{lm}^{Pseudo} of each partial wave at each r_i
c  rval = intended value of r_i 
c  ropt = true value of r_i (i.e. a nearby one available)
c  a    = a_i
c
      double precision, allocatable :: mel( : ), ener( : )
      double precision, allocatable :: fval( :, : ), rval( : )
      double precision, allocatable :: ropt( : )
      double precision, allocatable :: a( : )
c
c
c  variables used to do least square's fitting ...
c
c  s = secular matrix
c  suse = workspace
c  sinv = inverse of s
c  t = right hand side
c
c  S( i, j ) = sum_{ k = 1, nfit } psi_k( r_i ) * psi_k( r_j )
c
c  t( i ) = sum_{ k = 1, nfit } psi_k( r_i ) * mel( k )
c
c  we solve      
c
c        S * a = t, or  a = S ** -1 * t,
c
c  because this minimizes the error in the fit.
c
c        [ mel = sum_i a_i f( r_i ) . ]
c
      double precision, allocatable :: s( :, : ), suse( :, : )
      double precision, allocatable :: sinv( :, : ), t( : )
c
c  variables used for book-keeping ... which r_i's have been found
c  based on the list of input r's ...
c
      integer ngot, ngotlast
      logical, allocatable :: got( : )
c
c  work-horse variables
c
      integer i, j, k, irm           ! loop counters
      double precision dum, r, f     ! dummy, input r, input f
      double precision mfit          ! used to output value of fit
c
      read ( 5, * ) nfit, nchek, nrm
      allocate( rval( nrm ), a( nrm ), ropt( nrm ) )
      allocate( t( nrm ), got( nrm ) )
      allocate( s( nrm, nrm ), suse( nrm, nrm ), sinv( nrm, nrm ) )
c
      read ( 5, * ) rval             ! read in the whole list of r_i's
c
      ntot = nfit + nchek
      allocate( mel( ntot ), ener( ntot ), fval( ntot, nrm ) )
c
c
c  when looping over all-electron states, want l, M.E. and energy...
c
      do i = 1, ntot
        read ( 5, '(1x,2i5,2f20.10)' ) lval, nr, mel( i ), ener( i )
        do j = 1, nr            !
          read ( 5, * ) dum     !  skip radial fcn
        end do                  !
      end do
c
c
c  when looping over pseudo states, want wave function
c
      do i = 1, ntot
        read ( 5, '(1x,1i5)' ) dum ! skip M.E., energy
c
c  reset book-keeping to identify r's
c
        do irm = 1, nrm
          got( irm ) = .false.
        end do
        ngot = 0
        ngotlast = 0
c
c
c  loop over all r's ... get r and f, record if it counts for an r_i
c  do not double-count, and do not let two r_i's be the same ...
c
        do j = 1, nr
          read ( 5, * ) dum, r, f
          f = f / r
          do irm = 1, nrm
            if ( .not. got( irm ) ) then
              if ( r .gt. rval( irm ) ) then
                got( irm ) = .true.
                ngot = ngot + 1
                fval( i, irm ) = f
                ropt( irm ) = r         ! note de facto value of r_i
              end if
            end if
          end do
c
c
c  crash if two r_i's are the same ...
c
          if ( ngot .gt. ngotlast + 1 ) stop 'bad getting'
          ngotlast = ngot  
        end do
c
      end do
c
c
c  build up s and t ...
c
      do j = 1, nrm
        do k = 1, nrm
          s( j, k ) = 0.d0
          do i = 1, nfit
            s( j, k ) = s( j, k ) + fval( i, j ) * fval( i, k )
          end do
        end do
        t( j ) = 0.d0
        do i = 1, nfit
          t( j ) = t( j ) + fval( i, j ) * mel( i )
        end do
      end do
c
c
c  find inverse of s ...
c
      call invert( nrm, nrm, s, suse, sinv )
c
c
c  get a ... 
c
      do j = 1, nrm
        a( j ) = 0.d0
        do k = 1, nrm
          a( j ) = a( j ) + sinv( j, k ) * t( k )
        end do
      end do
c
c
c  output quality of fit
c
      do i = 1, ntot
        mfit = 0.d0
        do j = 1, nrm
          mfit = mfit + a( j ) * fval( i, j )
        end do
        mfit = mfit / mel( i )
        write ( 6, '(3(2x,1e15.8))' ) mel( i ), ener( i ), mfit 
      end do
c
c
c   output nrm, 1, lval ...
c   output list of r_i's, a_i's ...
c
      open( unit=72, file='gfile', form='formatted',
     &      status='unknown' )
      rewind 72
      write ( 72, '(2x,1i5)' ) nrm, 1, lval
      do j = 1, nrm
        write ( 72, '(3(2x,1e15.8))' ) ropt( j ), a( j )
      end do
      close( unit=72 )
c
      end
