C======================================================================C
      SUBROUTINE SNKEMATTER(MATTER, MNKE, NKE, Z2, RE, RHO, MW)
C======================================================================C
C
C     INPUTS
C       MATTER       : NAME OF MATTER.
C       MNKE         : MAXIMUM OF NKE
C
C     OUTPUTS
C       NKE          : NUMBER OF KINDS OF ELEMENTS IN MATTER.
C                      0 = I DON'T KNOW THE MATTER.
C       Z2           : Z OF EACH ELEMENT IN MATTER.
C       RE           : RATIO OF EACH ELEMENT IN MATTER.
C       RHO          : DENSITY (G/CM**3) OF MATTER.
C                      -1. = GAS
C       MW           : WEIGHT (G/MOL) OF MOLECULE IN MATTER.
C
C----------------------------------------------------------------------C

      IMPLICIT NONE

      CHARACTER*(*) MATTER
      INTEGER MNKE

      INTEGER NKE, Z2
      REAL RE, RHO, MW
      DIMENSION Z2(*),RE(*)

      INCLUDE 'DATAPATH.INC'
      CHARACTER*100 SNKE_MATTER_DAT
      PARAMETER(SNKE_MATTER_DAT=DATAPATH//'SNKE_MATTER.DAT')

      REAL M2

      CHARACTER*100 MATTER_TMP1, MATTER_TMP2

      INTEGER J, K

      REAL DUMMY, PCOEF(8)

      INCLUDE 'SNKE_MATTER.INC'


      LOGICAL FLAG_FIRSTCALL
      DATA FLAG_FIRSTCALL /.TRUE./


C----------------------------------------------------------------------C

C===== MAKE TABLE.

      IF (FLAG_FIRSTCALL) THEN
         
         OPEN(30, FILE=SNKE_MATTER_DAT, ACCESS='SEQUENTIAL',
     &        STATUS='UNKNOWN', FORM='FORMATTED')
         
         DO J=1, MAX_NUMBER_OF_TABLE
            READ(30,*) T_MATTER(J)
            
            IF(T_MATTER(J).EQ.'END') THEN
               NUMBER_OF_TABLE = J-1
               GO TO 200
            END IF
            
            READ(30,*) T_NKE(J)
            IF (T_NKE(J).GT.MAX_NUMBER_OF_NKE) THEN
               WRITE(6,*) 'ERROR : TOO LARGE ''NKE''',
     &              ' IN ''SNKE_MATTER.DAT'''
               STOP
            END IF
            
            READ(30,*) (T_Z2(K,J),K=1,T_NKE(J))
            READ(30,*) (T_RE(K,J),K=1,T_NKE(J))

            READ(30,*) T_RHO(J)
            IF ( (T_NKE(J).EQ.1) .AND. (T_RHO(J).EQ.0.) ) THEN
               CALL SCOEF(T_Z2(1,J), DUMMY, DUMMY, DUMMY,
     &              T_RHO(J), DUMMY, DUMMY, DUMMY, PCOEF)
            END IF
            
            T_MW(J) = 0.
            DO K=1, T_NKE(J)
               CALL SCOEF(T_Z2(K,J), DUMMY, DUMMY, M2, DUMMY, DUMMY,
     &              DUMMY, DUMMY, PCOEF)
               T_MW(J) = T_MW(J) + M2 * T_RE(K,J)
            END DO
            
         END DO

 200     CONTINUE

         CLOSE(30)

         FLAG_FIRSTCALL = .FALSE.

      END IF

C===== GET DATA FROM TABLE.

      MATTER_TMP1 = MATTER
      CALL CHANGECASE(MATTER_TMP1, 2)

      DO J=1, NUMBER_OF_TABLE

         MATTER_TMP2 = T_MATTER(J)
         CALL CHANGECASE(MATTER_TMP2, 2)
         
         IF(MATTER_TMP1.EQ.MATTER_TMP2) THEN

            IF (MNKE.LT.T_NKE(J)) GO TO 300

            NKE = T_NKE(J)

            DO K=1,T_NKE(J)
               Z2(K) = T_Z2(K,J)
               RE(K) = T_RE(K,J)
            END DO

            RHO = T_RHO(J)
            MW = T_MW(J)

            GO TO 500

         END IF

      END DO

 300  CONTINUE

      NKE = 0

 500  CONTINUE

      RETURN
      END

C======================================================================C
