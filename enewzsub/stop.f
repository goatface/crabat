C======================================================================C
	SUBROUTINE STOP(Z1,M1,Z2,EE,SE,SN,UNITS)
C
C     J. F. Ziegler, J. P. Biersack and U. Littmark
C     "The Stopping and Range of Ions in Solids"
C     Pergamon Press, New York (1985)
C
C======================================================================C
C
C       INPUTS
C         Z1     : ION ATOMIC NUMBER
C         M1     : ION ATOMIC WEIGHT (AMU)
C         Z2     : TARGET ATOMIC NUMBER
C         EE     : ION ENERGY (KEV)
C         UNITS  : UNITS OF SE AND SN
C                    1 = EV*CM^2/1.E15
C                    2 = MEV*CM2/MG
C                    3 = EV/ANGSTROM
C                    4 = L.S.S. REDUCED UNITS
C
C       OUTPUTS
C         SE     : CALCULATED ELECTRONIC STOPPING
C         SN     : CALCULATED NUCLEAR STOPPING
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1, Z2, UNITS
	REAL M1, EE, SE, SN

	REAL MASSMA
	REAL LFCTR
	REAL E

	REAL M2,RHO,ATRHO,VFERMI,PCOEF
	DIMENSION PCOEF(8)

	REAL EPSIL

	REAL DUMMY,DUMMY8
	DIMENSION DUMMY8(8)

C----------------------------------------------------------------------C

C***** CHECK INPUT PARAMETERS *****C
    
	IF ( Z1 .GT. 92 ) GO TO 190                                           
	IF ( EE .LT. 1.E-10 ) GO TO 220                                        

C***** GET ION DATA *****C
      
	CALL SCOEF(Z1,DUMMY,MASSMA,DUMMY,DUMMY,DUMMY,DUMMY,LFCTR,DUMMY8)
	IF ( M1 .EQ. 0. ) THEN
	   M1 = MASSMA
	END IF

C***** CHECK INPUT PARAMETERS *****C

      E = EE/M1
      IF ( E .GT. 100000. ) GO TO 190

C***** GET TARGET DATA

	CALL SCOEF(Z2,DUMMY,DUMMY,M2,RHO,ATRHO,VFERMI,DUMMY,PCOEF)

C***** CALCULATE SE DEPENDING ON ION TYPE

C= PROTON =C

	IF ( Z1 .EQ. 1 ) THEN
	   CALL PSTOP(Z1,M1,Z2,M2,E,PCOEF,SE)
	   GOTO 888
	END IF

C= ALPHA =C

	IF ( Z1 .EQ. 2 ) THEN
	   CALL HESTOP(Z1,M1,Z2,M2,E,PCOEF,SE)
	   GOTO 888
	END IF

C= HEAVY ION =C

	IF ( Z1 .LE. 92 ) THEN
	   CALL HISTOP(Z1,M1,Z2,M2,E,EE,VFERMI,LFCTR,PCOEF,SE)
	   GOTO 888
	END IF

 888	CONTINUE

C***** CALCULATE SN *****C
	
	EPSIL = 32.53 * M2 * EE / 
     &             ( Z1*Z2* (M1+M2) * ( Z1**.23 + Z2**.23 ) )

	IF (EPSIL.GE.30.) THEN
	   SN = ALOG(EPSIL) / (2*EPSIL)
	ELSE 
	   DUMMY = (.01321 * EPSIL**.21226) + (.19593 * EPSIL**.5)
	   SN = .5 * ALOG( 1 + 1.1383*EPSIL ) / ( EPSIL + DUMMY )
	END IF

C***** CONVERT FROM REDUCED UNITS TO EV*CM^2/1.E15

	SN = SN * Z1*Z2 * M1 * 8.462 / ( (M1+M2) * (Z1**.23 + Z2**.23) )

C***** CONVERT TO DESIRED STOPPING UNITS

	GO TO (230,160,170,180), UNITS

C= MEV*CM^2/MG
 
 160	SE = SE * .60222/M2
	SN = SN * .60222/M2
	GO TO 230

C= EV/ANGSTROM

 170	SE = SE * ATRHO * 1.E-23
	SN = SN * ATRHO * 1.E-23
	GO TO 230

C= L.S.S. REDUCED UNITS : (DELTA EPSILON/DELTA RHO)

 180	DUMMY = ( (M1+M2) * (Z1**.6667 + Z2**.6667)**.5 ) / 
     &                ( Z1*Z2 * M1 * 8.462 )
	SE = SE * DUMMY
	SN = SN * DUMMY
	GO TO 230

C***** ERROR ; OUT OF RANGE *****C

 190	WRITE (6,210)
 210	FORMAT(1X, 'INPUT PARAMETERS EXCEED LIMIT OF THIS PROGRAM')
 220	SE = 0.
	SN = 0.

C***** RETURN *****C

 230	RETURN
	END

C======================================================================C
	SUBROUTINE PSTOP(Z1,M1,Z2,M2,E,PCOEF,SE)
C
C       PROTON ELECTRONIC STOPPING POWERS
C
C======================================================================C

	IMPLICIT NONE

	INTEGER Z1, Z2
	REAL M1, M2, E, PCOEF(8), SE

	REAL PE0, VELPWR

	REAL PE, SL, SH

C----------------------------------------------------------------------C

C===== VELOCITY PROPORTIONAL STOPPING BELOW VELOCITY PEO.

	PE0 = 25.
	PE = AMAX1( PE0, E )
	
	SL = ( PCOEF(1) * PE**PCOEF(2) ) + ( PCOEF(3) * PE**PCOEF(4) )
	SH = PCOEF(5) / PE**PCOEF(6) * LOG(( PCOEF(7)/PE ) + PE*PCOEF(8))
	SE = SL*SH / (SL+SH)

C----- VELPWR IS THE POWER OF VELOCITY STOPPING BELOW PE0.

	IF ( E .LE. PE0 ) THEN 
	   VELPWR = 0.45
	   IF ( Z2 .LE. 6 ) VELPWR = 0.25
	   SE = SE * (E/PE0)**VELPWR
	END IF
	
	RETURN
	END

C======================================================================C
	SUBROUTINE HESTOP(Z1,M1,Z2,M2,E,PCOEF,SE)
C
C       HELIUM ELECTRONIC STOPPING POWERS
C
C======================================================================C

	IMPLICIT NONE

	INTEGER Z1, Z2
	REAL M1, M2, E, PCOEF(8), SE

	REAL HE0

	REAL HE, HEH, SP

	REAL A, B

C----------------------------------------------------------------------C

C===== VELOCITY PROPORTIONAL STOPPING BELOW KEV/AMU HE0.

	HE0=1.

	HE = AMAX1(HE0,E)

	B = ALOG(HE)

C===== BELOW IS FORMULA FOR THE RATION OF HE TO PROTON STOPPING POWERS.

	A = .2865 + .1266*B - .001429*B*B +
     &        .02402*B*B*B - .01135*B**4 + .001475*B**5
	HEH = 1. - EXP( - AMIN1(30.,A) )

C===== ADD Z1**3 EFFECT TO HE/H STOPPING POWER RATIO HEH

	A = ( 1. + (.007 + .00005*Z2) *
     &           EXP( - ( 7.6 - AMAX1(0.,ALOG(HE)) )**2 ) )
	HEH = HEH * A*A
	CALL PSTOP(Z1,M1,Z2,M2,HE,PCOEF,SP)
	SE = SP * HEH*Z1*Z1

C===== CALC. HE VELOCITY PROPORTIONAL STOPPING

	IF ( E .LE. HE0 ) THEN
	   SE = SE * SQRT(E/HE0)
	END IF
	
	RETURN
	END

C======================================================================C
	SUBROUTINE HISTOP(Z1,M1,Z2,M2,E,EE,VFERMI,LFCTR,PCOEF,SE)
C
C       HEAVY ION ELECTRONIC STOPPING POWERS
C
C======================================================================C

	IMPLICIT NONE

	INTEGER Z1, Z2
	REAL M1, M2, E, EE, VFERMI, LFCTR, PCOEF(8), SE

	REAL YRMIN, VRMIN

	REAL V, VR, YR, Q, L, ZETA, VMIN, EEE, POWER, SP

	REAL A, B, L0, L1, Q1
	INTEGER I

C----------------------------------------------------------------------C

C===== USE VELOCITY STOPPING FOR (YRMIN=VR/Z1**.67) .LE. 0.13
C===== OR FOR VR .LE. 1.0

	YRMIN = 0.13
	VRMIN = 1.0

	V = SQRT( E/25. ) / VFERMI
 	IF ( V .GE. 1. ) THEN 
	   VR = V * VFERMI * ( 1 + 1/(5*V*V) )
	ELSE
	   VR = ( 3 * VFERMI/4 ) 
     &	            * ( 1 + ( 2 * V*V / 3. ) - V**4 / 15. )
	END IF

C===== SET YR = MAXIMUM OF (VR/Z1**.67),(VRMIN/Z1**.67) OR YRMIN.

	YR = AMAX1(YRMIN, VR / Z1**.6667)
	YR = AMAX1(YR, VRMIN / Z1**.6667)

C===== Q = IONIZATION LEVEL OF THE ION AT VELOCITY YR.

	A = -.803 * YR**0.3 + 1.3167 * YR**0.6 + 
     &                  .38157 * YR + .008983*YR*YR
	Q = AMIN1( 1., AMAX1( 0., 1. - EXP( - AMIN1( A, 50. ) ) ) )

C===== NOW WE CONVERT IONIZATION LEVEL TO EFFECTIVE CHARGE.

	B = ( AMIN1( 0.43, AMAX1( .32, .12 + .025*Z1 ) ) ) / Z1**.3333
	L0 = ( .8 - Q * ( AMIN1( 1.2, 0.6 + Z1/30. ) ) ) / Z1**.3333

	IF ( Q .LT. 0.2 ) THEN
	   L1 = 0.
	ELSE IF ( Q .LT. ( AMAX1( 0., .9 - .025*Z1 ) ) ) THEN
	   Q1 = 0.2
	   L1 = B * ( Q - .2 ) / ABS( AMAX1( 0., .9-.025*Z1 )-.2000001 )
	ELSE IF ( Q .LT. 
     &  ( AMAX1( 0., 1. - .025 * AMIN1( 16., 1.*Z1 ) ) ) ) THEN
 	   L1 = B
	ELSE
 	   L1 = B * (1.-Q) / ( .025 * AMIN1( 16., 1.*Z1 ) )
	END IF

	L = AMAX1( L1, L0*LFCTR )

C===== ADD Z1**3 EFFECT AS SHOWN IN REF. 779.

	ZETA = Q + ( 1. / ( 2. * VFERMI**2 ) ) * 
     &             ( 1. - Q ) * 
     &             ALOG( 1. + ( 4 * L * VFERMI / 1.919 )**2 )
	A = - ( 7.6 - AMAX1( 0., ALOG(E) ) )**2
	ZETA =ZETA * ( 1. + ( 1. / Z1**2 ) * ( .18 + .0015*Z2 ) *
     &		 EXP(A) )

	IF ( YR .LE. AMAX1( YRMIN, VRMIN / Z1**.6667 ) ) THEN 

C-----  CALCULATE VELOCITY STOPPING FOR YR LESS THAN YRMIN.

 	   VRMIN = AMAX1( VRMIN, YRMIN * Z1**.6667 )
	   VMIN = .5 *
     &         ( VRMIN + 
     &         SQRT( AMAX1( 0., VRMIN**2 - 0.8 * VFERMI**2 ) ) )
	   EEE = 25 * VMIN**2
	   CALL PSTOP(Z1,M1,Z2,M2,EEE,PCOEF,SP)
	   POWER = .5
	   IF ( (Z2.EQ.6) .OR. 
     &     ( ( (Z2.EQ.14) .OR. (Z2.EQ.32) ) .AND. (Z1.LE.19) ) ) THEN
	      POWER = .375
	   END IF
	   SE = ( SP * (ZETA*Z1)**2 ) * ( E / EEE )**POWER

	ELSE 

 	   CALL PSTOP(Z1,M1,Z2,M2,E,PCOEF,SP)
	   SE = SP * ( ZETA * Z1 )**2

	END IF

	RETURN
	END

C======================================================================C
