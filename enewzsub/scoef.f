C======================================================================C
      SUBROUTINE SCOEF(Z1,MM1,M1,M2,RHO,ATRHO,VFERMI,LFCTR,PCOEF)
C
C     J. F. Ziegler, J. P. Biersack and U. Littmark
C     "The Stopping and Range of Ions in Solids"
C     Pergamon Press, New York (1985)
C
C======================================================================C
C
C       INPUTS
C         Z1     : ION ATOMIC NUMBER
C
C       OUTPUTS
C         MM1    : ATOMIC MASS OF THE MOST ABUNDANT ISOTOPE
C         M1     : ATOMIC WEIGHT OF THE MOST ABUNDANT ISOTOPE
C         M2     : ATOMIC WEIGHT OF SOLID WITH NORMAL ISOTOPES
C         RHO    : DENSITY OF SOLID IN GRAMS/CM^3
C         ATRHO  : DENSITY OF SOLID IN UNITS OF 10^22 ATOMS/CM^3
C         VFERMI : FERMI VELOCITY OF THE SOLID, IN BOHRVELOCITY
C         LFCTR  : FACTOR DETERMINING ION SCREENING LENGTH
C         PCOEF  : COEFFICIENTS FOR PROTON STOPPING CROSS-SECTION
C
C----------------------------------------------------------------------C

      IMPLICIT NONE

      INTEGER Z1
      REAL MM1, M1, M2, RHO, ATRHO, VFERMI, LFCTR, PCOEF
      DIMENSION PCOEF(8)

      INCLUDE 'DATAPATH.INC'
      CHARACTER*(*) SCOEF_DAT
      PARAMETER(SCOEF_DAT=DATAPATH//'SCOEF.DAT')

      REAL Z,A,MASS1,MASS2,DENSITY,ATOMDENSITY,FERMIVEL,
     &     FACTSL,COEFPSCS
      DIMENSION Z(92),A(92),MASS1(92),MASS2(92),
     &     DENSITY(92),ATOMDENSITY(92),FERMIVEL(92),
     &     FACTSL(92),COEFPSCS(92,8)

      COMMON /TABLESCOEF/ Z,A,MASS1,MASS2,DENSITY,ATOMDENSITY,
     &     FERMIVEL,FACTSL,COEFPSCS

      INTEGER LINENUM

      INTEGER COEFNUM

      LOGICAL FLAG_FIRSTCALL
      DATA FLAG_FIRSTCALL /.TRUE./

C----------------------------------------------------------------------C

C***** INIT *****C
      
      IF ( FLAG_FIRSTCALL .EQV. .TRUE. ) THEN
      
C***** OPEN DATA FILE *****C
         
         OPEN( UNIT=20, ERR=333, FILE=SCOEF_DAT, STATUS='OLD',
     &        ACCESS='SEQUENTIAL', FORM='FORMATTED')
         
C*****CREATE TABLE *****C
         
         DO LINENUM=1,92

            READ(20,100,ERR=998,END=997) Z(LINENUM),A(LINENUM),
     &           MASS1(LINENUM),MASS2(LINENUM),DENSITY(LINENUM),
     &           ATOMDENSITY(LINENUM), FERMIVEL(LINENUM),
     &           FACTSL(LINENUM)
 100        FORMAT(F3.0,F4.0,F9.4,F9.4,F10.6,F7.3,F8.5,F5.2)

            IF ( INT(Z(LINENUM)) .NE. LINENUM ) THEN
               WRITE(*,*) 'ERROR'
     &              //' : INVALID ORDER IN ',SCOEF_DAT,
     &              ' ; Z = ',
     &              INT(Z(LINENUM)), ' AT LINE ',LINENUM
               STOP
            END IF

         END DO
         
         DO LINENUM=1,92

            READ(20,200,ERR=998,END=997) Z(LINENUM), 
     &           (COEFPSCS(LINENUM,COEFNUM),COEFNUM=1,8)
 200        FORMAT(F3.0,F11.7,F10.7,F12.8,F8.5,F8.2,F8.5,F10.3,F10.7)

            IF ( INT(Z(LINENUM)) .NE. LINENUM ) THEN
               WRITE(*,*) 'ERROR'
     &              //' : INVALID ORDER IN ',SCOEF_DAT,
     &              ' ; Z = ',
     &              INT(Z(LINENUM)), ' AT LINE ',LINENUM+92
               STOP
            END IF

         END DO

C***** CLOSE DATA FILE *****C
         
         CLOSE(20)

C***** CHECK FLAG *****C
         
         FlAG_FIRSTCALL = .FALSE.
         
      END IF

C***** SET OUTPUTS *****C

      MM1    = A(Z1)
      M1     = MASS1(Z1)
      M2     = MASS2(Z1)
      RHO    = DENSITY(Z1)
      ATRHO  = ATOMDENSITY(Z1)
      VFERMI = FERMIVEL(Z1)
      LFCTR  = FACTSL(Z1)
      DO COEFNUM=1,8
         PCOEF(COEFNUM) = COEFPSCS(Z1,COEFNUM)
      END DO

      ATRHO = ATRHO*1E22

C***** RETURN *****C

      RETURN

C***** FILE OPEN ERROR *****C
 333  WRITE(*,*) 'CANNOT OPEN FILE ; ',SCOEF_DAT
      STOP

C***** FILE READ ERROR *****C
 998  WRITE(*,*) 'READ ERROR ; ',SCOEF_DAT
      STOP

C***** FILE IS TRUNCATED *****C
 997  WRITE(*,*) 'FILE IS TRUNCATED ; ',SCOEF_DAT
      STOP

C***** END *****C

      END

C======================================================================C
