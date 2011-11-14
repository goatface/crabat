C======================================================================C
	REAL FUNCTION E2RANGEZ(Z1,M1,NKE,Z2,RE,E)
C======================================================================C
C
C	INPUTS
C         Z1	 : Z OF A PARTICLE WHICH RUNS THROUGH MEDIA
C         M1	 : MASS (AMU) OF A PARTICLE 
C                    WHICH RUNS THROUGH MEDIA
C         NKE	 : NUMBER OF KINDS OF ELMENTS IN MEDIA
C         Z2	 : Z OF EACH ELEMENT IN MEDIA
C	  RE	 : RATIO OF EACH ELEMENT IN MEDIA
C         E	 : ENERGY (MEV/AMU) OF A PARTICLE 
C                    WHICH RUNS THROUGH MEDIA
C
C	OUTPUTS
C	  E2RANGEZ : RANGE (MG/CM**2)
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1, NKE, Z2
	REAL M1, RE, E
	DIMENSION Z2(*),RE(*)


	INTEGER TABNUM
	INTEGER POINTNUM

     	REAL EE1, EE2, EE3, R1, R2, R3, AA, BB, CC, DD
	
	INTEGER CHECKTRZ
	INTEGER MAKETRZ

	INCLUDE 'TABLE_RANGE_Z.INC'

	LOGICAL FLAG_FIRSTCALL
	DATA FLAG_FIRSTCALL /.TRUE./

C----------------------------------------------------------------------C

C===== CHECK NKE =====C
	
	IF ( NKE .GT. TRZMAXNKE ) THEN
	   WRITE(*,*) 'ERROR : TOO LARGE ''NKE'' ; NKE = ',
     &                NKE, ' > ',TRZMAXNKE
	   STOP
	END IF

C===== CHECK IF FIRST CALL =====C

	IF ( FLAG_FIRSTCALL .EQV. .TRUE. ) THEN
	   CALL INITTRZ
	   FLAG_FIRSTCALL = .FALSE.
	END IF

C===== MAKE TABLE IF NEED =====C

	TABNUM = CHECKTRZ(Z1,M1,NKE,Z2,RE)

	IF ( TABNUM .EQ. 0 ) THEN
	   TABNUM =  MAKETRZ(Z1,M1,NKE,Z2,RE)
	END IF

C===== GET RANGE FROM TABLE =====C

 	DO POINTNUM=1, TRZNPOINT

	   IF ( E .GT. TRZENERGY(POINTNUM) ) GOTO 300

	   EE1 = TRZENERGY(POINTNUM-1)
	   EE2 = TRZENERGY(POINTNUM  )
	   EE3 = TRZENERGY(POINTNUM+1)

	   R1 = TRZRANGE(TABNUM,POINTNUM-1)
	   R2 = TRZRANGE(TABNUM,POINTNUM  )
	   R3 = TRZRANGE(TABNUM,POINTNUM+1)

	   DD = EE3*EE3 * (EE2-EE1) + EE2*EE2 * (EE1-EE3) +
     &              EE1*EE1 * (EE3-EE2)
	   AA = ( R3*(EE2-EE1) + R2*(EE1-EE3) + R1*(EE3-EE2) ) / DD
	   BB = ( R3 * ( EE1*EE1 - EE2*EE2 ) +
     &            R2 * ( EE3*EE3 - EE1*EE1 ) +
     &            R1 * ( EE2*EE2 - EE3*EE3) ) / DD
	   CC = ( R3 * EE2 * EE1 * (EE2-EE1) +
     &            R2 * EE1 * EE3 * (EE1-EE3) +
     &   	  R1 * EE3 * EE2 * (EE3-EE2) ) / DD
	   E2RANGEZ = AA * E*E + BB * E + CC
	   GOTO 400

 300	   CONTINUE
	
	END DO
	     
	WRITE(6,'('' ERROR : OUT OF ENERGY RANGE'')')
	STOP

C======================================================================C

 400	CONTINUE

	RETURN
	END

C======================================================================C
	REAL FUNCTION R2ENERGYZ(Z1,M1,NKE,Z2,RE,R)
C======================================================================C
C
C	INPUTS	
C         Z1	  : Z OF A PARTICLE WHICH RUNS THROUGH MEDIA
C         M1	  : MASS (AMU) OF A PARTICLE 
C                     WHICH RUNS THROUGH MEDIA
C         NKE	  : NUMBER OF KINDS OF ELMENTS IN MEDIA
C         Z2      : Z OF EACH ELEMENT IN MEDIA
C         RE	  : RATIO OF EACH ELEMENT IN MEDIA
C         R	  : RANGE (MG/CM**2)
C
C	OUTPUTS
C         R2ENERGYZ : ENERGY (MEV/AMU) OF A PARTICLE 
C                        WHICH RUNS THROUGH MEDIA
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1, NKE, Z2
	REAL M1, RE, R
	DIMENSION Z2(*),RE(*)

	REAL E,RR

	INTEGER TABNUM
	INTEGER POINTNUM

	REAL E2RANGEZ
	INTEGER CHECKTRZ

	INCLUDE 'TABLE_RANGE_Z.INC'

	REAL EE1, EE2, EE3, R1, R2, R3, AA, BB, CC, DD
	
C----------------------------------------------------------------------C
	
C===== MAKE TABLE.

	E = 18.
	RR = E2RANGEZ(Z1,M1,NKE,Z2,RE,E)

C===== GET ENERGY FROM TABLE.

	DO POINTNUM=1,TRZNPOINT
	   
	   TABNUM = CHECKTRZ(Z1,M1,NKE,Z2,RE)

	   IF ( TABNUM .EQ. 0 ) THEN
	      WRITE(*,*) 'ERROR : NO TABLE'
	      STOP
	   END IF
	   
	   IF ( R .GT. TRZRANGE(TABNUM,POINTNUM) ) GOTO 30

	   EE1 = TRZENERGY(POINTNUM-1)
	   EE2 = TRZENERGY(POINTNUM  )
	   EE3 = TRZENERGY(POINTNUM+1)
	   
	   R1 = TRZRANGE(TABNUM,POINTNUM-1)
	   R2 = TRZRANGE(TABNUM,POINTNUM  )
	   R3 = TRZRANGE(TABNUM,POINTNUM+1)

	   DD = EE3*EE3 * (EE2-EE1) + 
     &	        EE2*EE2 * (EE1-EE3) + 
     &   	EE1*EE1 * (EE3-EE2)
	   AA = ( R3 * (EE2-EE1) + 
     &            R2 * (EE1-EE3) +
     &            R1 * (EE3-EE2) ) / DD
	   BB = ( R3 * ( EE1*EE1 - EE2*EE2 ) + 
     &            R2 * ( EE3*EE3 - EE1*EE1 ) +
     &            R1 * ( EE2*EE2 - EE3*EE3 ) ) / DD
	   CC = ( R3 * EE2 * EE1 * (EE2-EE1) +
     &            R2 * EE1 * EE3 * (EE1-EE3) +
     &	          R1 * EE3 * EE2 * (EE3-EE2) ) / DD

	   R2ENERGYZ = ( -BB + SQRT( BB*BB - 4. *AA*(CC-R)) ) / 
     &                                   ( 2.*AA )

	   GOTO 40

 30	   CONTINUE

	END DO

	WRITE(6,'('' ERROR : OUT OF RANGE RANGE'')')
	STOP

 40	CONTINUE

C***** RETURN *****C

	RETURN
	END

C======================================================================C
