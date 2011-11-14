C======================================================================C
      SUBROUTINE CHANGECASE(STRING, MODE)
C======================================================================C
C
C     INPUTS
C       STRING   :  A STRING WHICH IS TO BE CHANGED.
C       MODE     :  1 = UPPER CASE TO LOWER CASE
C                   2 = LOWER CASE TO UPPER CASE
C     OUTPUTS
C       STRING   :  THE STRING WHICH HAS BEEN CHANGED.
C
C----------------------------------------------------------------------C

      IMPLICIT NONE

      CHARACTER*(*) STRING
      INTEGER MODE

      INTEGER N
      INTEGER CH

      INTEGER I

C----------------------------------------------------------------------C

      N = LEN(STRING)

      DO I=1,N

         CH = ICHAR(STRING(I:I))

         IF ( ( MODE .EQ. 1 ) .AND.
     &        ( 65 .LE. CH ) .AND. ( CH .LE. 90 ) ) THEN
            CH = CH + 32
         ELSE IF ( ( MODE .EQ. 2 ) .AND.
     &           ( 97 .LE. CH ) .AND. ( CH .LE. 122 ) ) THEN
            CH = CH - 32
         END IF

         STRING(I:I) = CHAR(CH)

      END DO

      RETURN

      END

C======================================================================C

