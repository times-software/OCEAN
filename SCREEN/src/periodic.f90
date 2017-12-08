module periodic 

  implicit none
  private
  save

  character(len=2), dimension(109), parameter :: elements =             &
     &    (/ 'H_', 'He', 'Li', 'Be', 'B_', 'C_', 'N_',                  &
     & 'O_', 'F_', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P_', 'S_', 'Cl', 'Ar',&
     & 'K_', 'Ca', 'Sc', 'Ti', 'V_', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',&
     & 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y_', 'Zr',&
     & 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb',&
     & 'Te', 'I_', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm',      &
     & 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',&
     & 'Ta', 'W_', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',&
     & 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U_', 'Np', 'Pu',&
     & 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db',&
     & 'Sg', 'Bh', 'Hs', 'Mt' /) 

  integer, parameter :: el_len = 109

  public :: getsymbol, getsymbol_underscore, get_atom_number

  contains

  subroutine getsymbol( zatom, satom )
    integer, intent( in ) :: zatom
    character(len=3), intent(out) :: satom

    satom = elements( zatom )
    if( satom(2:2) .eq. '_' ) satom(2:2) = ' '
!    satom(3) = ' '

  end subroutine getsymbol

  subroutine getsymbol_underscore( zatom, satom )
    integer, intent( in ) :: zatom
    character(len=3), intent(out) :: satom

    satom = elements( zatom )
!    satom(3) = ' '

  end subroutine getsymbol_underscore

  subroutine get_atom_number( zatom, satom )
    integer, intent( out ) :: zatom
    character(len=2), intent( in ) :: satom
    !
    integer :: i
    !
    zatom = -1
    do i = 1, el_len
      if( satom .eq. elements( i ) ) then
        zatom = i
        goto 111
      endif
    enddo
  
111 continue
    
  end subroutine get_atom_number

end module periodic
