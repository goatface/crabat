C======================================================================C
c      program ENEWZ
      SUBROUTINE ENEWZSUB(Z1,M1,E,matter1,unit_pressure,pressure,
     &     temperature,unit_thick,thick1,AFT_ENE)
C======================================================================C
C     CALCULATE ENERGY LOSS IN MATTER
C       USING ZIEGLER'S SUBROUTINE
C Modified by HAYAKAWA Seiya in 2008
C======================================================================C

      IMPLICIT NONE

      REAL AFT_ENE
      
      INTEGER MAX_NKE
      PARAMETER(MAX_NKE=10)
      INTEGER MAX_NM
      PARAMETER(MAX_NM=100)

      INTEGER Z1
      REAL M1, E

      INTEGER NM
      CHARACTER*32 MATTER(MAX_NM)
      character*32 matter1
      INTEGER NKE(MAX_NM), Z2(MAX_NKE, MAX_NM)
      REAL RE(MAX_NKE, MAX_NM), RHO(MAX_NM), MW(MAX_NM)
      REAL PRESSURE, TEMPERATURE,THICK(MAX_NM)
      real thick1
      INTEGER UNIT_PRESSURE, UNIT_TEMPERATURE, UNIT_THICK
      CHARACTER*20 C_UNIT

      REAL ENEW(0:MAX_NM)

      CHARACTER*14 COLUMN_MATTER(5)

      INTEGER Z2_TMP(MAX_NKE)
      REAL RE_TMP(MAX_NKE)
      REAL M2

      INTEGER J, K, L
      REAL DUMMY, PCOEF(8)

      REAL ENERGYNEWZ

      INCLUDE 'SNKE_MATTER.INC'
      
c      Z1=2
c      M1=4
c      E=2.5

c     layer ha 1 mai ni gentei
      NM=1
c      MATTER(1)='HE'
      
c      UNIT_PRESSURE=1
c      PRESSURE=760.0
c      TEMPERATURE=300.0
      
c      UNIT_THICK=1
c      THICK(1)=300.

C----------------------------------------------------------------------C

      CALL SNKEMATTER('MATTER', MAX_NKE, NKE(1), Z2_TMP, RE_TMP,
     &     RHO(1), MW(1))
c      write(*,*) nke(1), z2_tmp, re_tmp, rho(1), mw(1)

c      WRITE(*,*)
c      WRITE(*,'(1X, ''Input : Z, Mass(amu), Energy(MeV/amu) > '')')
c      READ(*,*) Z1, M1, E
      
c      WRITE(*,*)
c      WRITE(*,'(1X ,''Input : Number of layer > '')')
c      READ(*,*) NM
      IF (NM.GT.MAX_NM) THEN
         WRITE(*,*) 
     &        'ERROR : TOO LARGE ''NM'' ; NM = ',NM,' > ',MAX_NM
         STOP
      END IF

      DO J=1, NM
         
c         WRITE(*,*)

c         WRITE(*,'(1X, ''Input : Matter of layer#'',I2)') J

c         WRITE(*,*) '------------------------------',
c     &        ' MATTER LIST ', 
c     &        '------------------------------'
         DO K=1, NUMBER_OF_TABLE/5+1
            DO L=1,5
               IF ( ((K-1)*5+L) .LE. NUMBER_OF_TABLE ) THEN
                  COLUMN_MATTER(L) = T_MATTER(((K-1)*5+L))
               ELSE
                  COLUMN_MATTER(L) = ' '
               END IF
            END DO
c            WRITE(*,'(5(A, '' ''))') COLUMN_MATTER
         END DO
c         WRITE(*,*) '------------------------------',
c     &        '-------------', 
c     &        '------------------------------'

c         WRITE(*,'(1X, ''Input : Matter > '')')

c        READ(*,*) MATTER(J)
c     matter(1) dake tukau
c         write(*,*) matter(J), NKE(J)
        matter(1)=matter1

c        write(*,*) matter1, NKE(J)
c         write(*,*) matter(J), NKE(J)
c         write(*,*) max_nke
         CALL SNKEMATTER(MATTER(J), MAX_NKE, NKE(J), Z2_TMP, RE_TMP,
     &        RHO(J), MW(J))
c         write(*,*) nke(j), z2_tmp, re_tmp, rho(j), mw(j)
c         stop

         IF (NKE(J).EQ.0) THEN
            DO K=1, MAX_NKE
               WRITE(*,'(3X,A,A)')
     &              'Input : Z, Number of element in one molecule',
     &              ' (END = -1 -1.)'
               WRITE(*,'(5X, ''> '')') 
               READ(*,*) Z2(K,J), RE(K,J)
               IF ((Z2(K,J).EQ.-1).AND.(RE(K,J).EQ.-1.)) THEN
                  NKE(J) = K-1
                  GO TO 100
               END IF
            END DO

            WRITE(*,'(1X, ''Too many elements.'')')
            STOP

 100        CONTINUE
            
c            WRITE(*,*) 'Input : Rho(g/cm**3) ( -1. = GAS ) > '
            READ(*,*) RHO(J)

            MW(J) = 0.
            DO K=1, NKE(J)
               CALL SCOEF(Z2(K,J), DUMMY, DUMMY, M2, DUMMY, DUMMY,
     &              DUMMY, DUMMY, PCOEF)
               MW(J) = MW(J) + M2*RE(K,J)
            END DO
         ELSE 
            DO K=1,NKE(J)
               Z2(K,J) = Z2_TMP(K)
               RE(K,J) = RE_TMP(K)
            END DO
         END IF

c         WRITE(*,*)
         IF (RHO(J) .EQ. -1) THEN
c            WRITE(*,*) 'Input : Unit ( [1] torr, [2] mbar [3] atm ) > '
c            READ(*,*) UNIT_PRESSURE
            IF ( UNIT_PRESSURE .EQ. 1 ) THEN
               C_UNIT = 'torr'
            ELSE IF ( UNIT_PRESSURE .EQ. 2 ) THEN
               C_UNIT = 'mbar'
            ELSE IF ( UNIT_PRESSURE .EQ. 3 ) THEN
               C_UNIT = 'atm'
            END IF
c               WRITE(*,'(1X,A,A4,A)')
c     &              'Input : Pressure ( ', C_UNIT, ' ) of gas > '
c            READ(*,*) PRESSURE
                        
c            WRITE(*,*) 'Input : Temperature (K) > '
c            READ(*,*) TEMPERATURE
            
            CALL SRHOGAS( PRESSURE, TEMPERATURE, 
     &           UNIT_PRESSURE, 1, MW(J), RHO(J) )
            
         END IF
         
c         WRITE(*,*)
c         WRITE(*,*) 'Input : Unit ( [1] mm, [2] mg/cm**2 ) > '
c         READ(*,*) UNIT_THICK
         IF (UNIT_THICK.EQ.1) THEN
            C_UNIT = 'mm'
         ELSE
            C_UNIT = 'mg/cm**2'
         END IF
c         WRITE(*,'(1X,A,A8,A,I2,A)')
c     &        'Input : Thickness ( ', C_UNIT, ' ) of layer#',J,' > '
c         READ(*,*) THICK(J)
c     thick(1) dake tukau
         thick(j)=thick1
c         write(*,*) 'thickj',thick(j),'thick1',thick1,'rhoj',rho(j)

         IF (UNIT_THICK.EQ.1) THEN
            THICK(J) = THICK(J) / 10. * RHO(J) * 1000.
         END IF

         
      END DO
      
      ENEW(0)=E
      
      DO J=1,NM
         DO K=1,NKE(J)
            Z2_TMP(K) = Z2(K,J)
            RE_TMP(K) = RE(K,J)
         END DO
         ENEW(J) = ENERGYNEWZ(Z1,M1,NKE(J),
     &        Z2_TMP,RE_TMP,ENEW(J-1),THICK(J))
      END DO
      
c      WRITE(*,*) 
c      WRITE(*,'(1X,F10.5'' MeV/amu ('',F10.5''MeV)'')')
c     &     ENEW(0), ENEW(0)*M1


c      DO J=1,NM
c         WRITE(*,'(1X,''--- MATTER : '',A,'' '',F8.3,'' mg/cm^2 ---'')') 
c     &        MATTER(J), THICK(J)
c         WRITE(*,'(1X,F10.5,'' MeV/amu ('',F10.5,'' MeV)'')')
c     &        ENEW(J), ENEW(J)*M1
c      END DO
      
      AFT_ENE=ENEW(1)  ![MeV/amux]
c      WRITE(*,*) AFT_ENE,'[MeV] after enewz'
      
      RETURN      
      END


      
C======================================================================C

