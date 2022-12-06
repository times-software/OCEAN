!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdcmt.f90,v $:
! $Revision: 1.1.2.1 $
! $Author: hebhop $
! $Date: 2012/04/04 23:11:55 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE rdcmt(iUnt,Cmt)
      INTEGER iUnt, i1
      CHARACTER(300) line
      CHARACTER(4) Cmt
      CHARACTER :: TmpCmt(4), ch
      LOGICAL CmtLin


      CmtLin = .true.
      DO i1 = 1, 4
         TmpCmt(i1) = Cmt(i1:i1)
      END DO
 5    CONTINUE
      READ(iUnt,'(A1)',END=10) ch
      DO i1 = 1, 4
         IF(ch.eq.TmpCmt(i1)) then
           goto 5
         endif
      END DO
      
 10   CONTINUE
      BACKSPACE(iUnt)
      
      RETURN
      END
