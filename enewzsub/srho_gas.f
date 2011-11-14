C======================================================================C
      SUBROUTINE SRHOGAS( P, T, UNIT_P, UNIT_T ,MW, RHO)
C======================================================================C
C
C     INPUTS 
C       P, T     : PRESSURE AND TEMPERATURE OF GAS.
C       UNIT_P   : UNIT OF P
C                  1 = TORR
C                  2 = MBAR
C                  3 = ATM
C       UNIT_T   : UNIT OF T
C                  1 = K
C                  2 = C
C       MW       : MOLECULAR WEIGHT (G/MOL)
C
C     OUTPUTS
C       RHO      : DENSITY OF GAS (G/CM**3)
C
C----------------------------------------------------------------------C

      IMPLICIT NONE

      REAL P, T, MW
      INTEGER UNIT_P, UNIT_T
      
      REAL RHO

      REAL TORR2PA, MBAR2PA, ATM2PA
      PARAMETER(TORR2PA=133.322)
      PARAMETER(MBAR2PA=100.)
      PARAMETER(ATM2PA=101325.)

      REAL NA, KB
      PARAMETER(NA=6.0221367E23)
      PARAMETER(KB=1.380658E-23)

      REAL P_PA, T_K

C----------------------------------------------------------------------C
      
      IF ( UNIT_P .EQ. 1 ) THEN
         P_PA = P * TORR2PA
      ELSE IF ( UNIT_P .EQ. 2 ) THEN
         P_PA = P * MBAR2PA
      ELSE IF (UNIT_P .EQ. 3) THEN
         P_PA = P * ATM2PA
      ELSE
         P_PA = 0.
      END IF

      IF ( UNIT_T .EQ. 1 ) THEN
         T_K = T
      ELSE IF ( UNIT_T .EQ. 2 ) THEN
         T_K = T + 273.15
      ELSE
         T_K = 0.
      END IF

      RHO = ( MW * P_PA ) / ( KB * T_K * NA ) * 1.E-6

      RETURN

      END

C======================================================================C
