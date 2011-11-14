C======================================================================C
	REAL FUNCTION ENERGYNEWZ(Z1,M1,NKE,Z2,RE,E,THICK)
C======================================================================C
C
C	INPUTS
C	  Z1	       : Z OF A PARTICLE WHICH RUNS THROUGH MEDIA
C	  M1	       : MASS (AMU) OF A PARTICLE 
C		           WHICH RUNS THROUGH MEDIA
C	  NKE	       : NUMBER OF KINDS OF ELMENTS IN MEDIA
C	  Z2	       : Z OF EACH ELEMENT IN MEDIA
C	  RE	       : RATIO OF EACH ELEMENT IN MEDIA
C	  E  	       : ENERGY (MEV/AMU) OF A PARTICLE 
C		           WHICH RUNS THROUGH MEDIA
C	  THICK	       : THICKNESS (MG/CM**2) OF MEDIA
C
C	OUTPUTS
C	  ENERGYNEWZ : ENERGY (MEV/AMU) AT THE BACK OF THE MEDIA
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1,NKE,Z2
	REAL M1,RE,E,THICK
	DIMENSION Z2(*),RE(*)

	REAL RANGE,RANGEB

	REAL E2RANGEZ, R2ENERGYZ

C----------------------------------------------------------------------C

	RANGE = E2RANGEZ(Z1,M1,NKE,Z2,RE,E)
	RANGEB = RANGE - THICK

	IF ( RANGEB .GT. 0. ) THEN
	   ENERGYNEWZ = R2ENERGYZ(Z1,M1,NKE,Z2,RE,RANGEB)
	ELSE
	   ENERGYNEWZ = 0.
	END IF

	RETURN
	END

C======================================================================C
	REAL FUNCTION ENERGYOLDZ(Z1,M1,NKE,Z2,RE,E,THICK)
C======================================================================C
C
C	INPUTS
C	  Z1	       : Z OF A PARTICLE WHICH RUNS THROUGH MEDIA
C	  M1	       : MASS (AMU) OF A PARTICLE 
C		           WHICH RUNS THROUGH MEDIA
C	  NKE	       : NUMBER OF KINDS OF ELMENTS IN MEDIA
C	  Z2	       : Z OF EACH ELEMENT IN MEDIA
C	  RE	       : RATIO OF EACH ELEMENT IN MEDIA
C	  E  	       : ENERGY (MEV/AMU) OF A PARTICLE 
C		           WHICH RAN THROUGH MEDIA
C	  THICK	       : THICKNESS (MG/CM**2) OF MEDIA
C
C	OUTPUTS
C	  ENERGYNEWZ : ENERGY (MEV/AMU) IN FRONT OF THE MEDIA
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1,NKE,Z2
	REAL M1,RE,E,THICK
	DIMENSION Z2(*),RE(*)

	REAL RANGE,RANGEB

	REAL E2RANGEZ, R2ENERGYZ

C----------------------------------------------------------------------C

	RANGE = E2RANGEZ(Z1,M1,NKE,Z2,RE,E)

	RANGEB = RANGE + THICK

	IF ( RANGEB .GT. 0. ) THEN
	   ENERGYOLDZ = R2ENERGYZ(Z1,M1,NKE,Z2,RE,RANGEB)
	ELSE
	   ENERGYOLDZ = 0.
	END IF
	
	RETURN
	END

C======================================================================C


