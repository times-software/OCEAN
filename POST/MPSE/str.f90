! FUNCTION ISTRLN (STRING)  Returns index of last non-blank
!                           character.  Returns zero if string is
!                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
!     there is a tab character here  ^

!  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

!  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
! SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*200  TMP
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
!     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

!  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

!  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

!  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

!  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
! SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END
! SUBROUTINE LOWER (STRING)  Changes A-Z to lower case.

      SUBROUTINE LOWER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 65) .OR.  (IC .GT. 90))  GOTO 10
         STRING (I:I) = CHAR (IC + 32)
   10 CONTINUE

      RETURN
      END
!***********************************************************************
!
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
!
!     Breaks string into words.  Words are seperated by one or more
!     blanks or tabs, or a comma and zero or more blanks.
!
!     ARGS        I/O      DESCRIPTION
!     ----        ---      -----------
!     S            I       CHAR*(*)  String to be broken up
!     NWORDS      I/O      Input:  Maximum number of words to get
!                          Output: Number of words found
!     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
!                          Contains words found.  WORDS(J), where J is
!                          greater then NWORDS found, are undefined on
!                          output.
!
!      Written by:  Steven Zabinsky, September 1984
!      Tab char added July 1994.
!
!**************************  Deo Soli Gloria  **************************

!  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', TAB = '	')
!     there is a tab character here               ^.

!  -- BETW    .TRUE. if between words
!     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

!  -- Maximum number of words allowed
      WORDSX = NWORDS

!  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

!  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

!  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END

!***********************************************************************
!
      SUBROUTINE BWORDS2 (S, NWORDS, WORDS)
!
!     Breaks string into words.  Words are seperated by one or more
!     blanks or tabs.
!
!     ARGS        I/O      DESCRIPTION
!     ----        ---      -----------
!     S            I       CHAR*(*)  String to be broken up
!     NWORDS      I/O      Input:  Maximum number of words to get
!                          Output: Number of words found
!     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
!                          Contains words found.  WORDS(J), where J is
!                          greater then NWORDS found, are undefined on
!                          output.
!
!      Written by:  Steven Zabinsky, September 1984
!      Tab char added July 1994.
!
!**************************  Deo Soli Gloria  **************************

!  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
!     there is a tab character here               ^.

!  -- BETW    .TRUE. if between words
!     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

!  -- Maximum number of words allowed
      WORDSX = NWORDS

!  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

!  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

!  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END


      SUBROUTINE UNTAB (STRING)
! REPLACE TABS WITH BLANKS :    TAB IS ASCII DEPENDENT
      INTEGER        ITAB , I, ILEN, ISTRLN
      PARAMETER      (ITAB = 9)
      CHARACTER*(*)  STRING, TAB*1
      EXTERNAL ISTRLN
      TAB  = CHAR(ITAB)
      ILEN = MAX(1, ISTRLN(STRING))
 10   CONTINUE 
        I = INDEX(STRING(:ILEN), TAB ) 
        IF (I .NE. 0) THEN
            STRING(I:I) = ' '
            GO TO 10
        END IF
      RETURN
! END SUBROUTINE UNTAB
      END

      logical function iscomm (line)
!     returns true if line is a comment or blank line, false otherwise
!#mn{ rewritten to allow ";*%#" as comment characters
       character*(*) line
       iscomm = ((line.eq.' ').or.(index(';*%#',line(1:1)).ne.0))
!#mn}
      return
      end
