      SUBROUTINE rdloss(Uloss,Dat,NPts)
      IMPLICIT NONE
c     Variables passed
      INTEGER Uloss, NPts
      DOUBLE PRECISION Dat(50000,2)

      
c     local variables
      INTEGER nterp,n,i,ierr
      DOUBLE PRECISION dum
      CHARACTER(4) comment
      CHARACTER(200) Line
      comment='#*!c'
      nterp=1
      n=0

c     Read data into Dat
 200  DO i=1,50000
c     Read past comments
         CALL rdcmt(Uloss,comment)
         READ(Uloss,*,end=250,IOSTAT=ierr) Dat(i,1), Dat(i,2)
         n=i
      END DO      
 250  BACKSPACE(Uloss)
c 250  CONTINUE
      NPts=n
      
300   CLOSE(Uloss)
      RETURN 
      END
