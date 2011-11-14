C======================================================================C
C	PARAMETER
C======================================================================C

	INTEGER TRZMAXNTAB
	PARAMETER(TRZMAXNTAB = 10 )

	INTEGER TRZMAXNKE
	PARAMETER(TRZMAXNKE = 20)

	INTEGER TRZNPOINT
	PARAMETER(TRZNPOINT = 200)

C======================================================================C
C	RANGE TABLE
C======================================================================C

	INTEGER TRZNTAB

	INTEGER   TRZZ1
	REAL    TRZM1
	DIMENSION TRZZ1(TRZMAXNTAB), TRZM1(TRZMAXNTAB)

	INTEGER  TRZNKE
	INTEGER  TRZZ2
	REAL   TRZRE
	DIMENSION TRZNKE(TRZMAXNTAB), 
     &	     		TRZZ2(TRZMAXNTAB, TRZMAXNKE),
     &			TRZRE(TRZMAXNTAB, TRZMAXNKE)

	REAL  TRZENERGY
	REAL  TRZRANGE
	DIMENSION TRZENERGY(0:TRZNPOINT),
     &			TRZRANGE(TRZMAXNTAB,0:TRZNPOINT)

	COMMON /TABLERANGEZ/ TRZNTAB,
     &				TRZZ1, TRZM1,
     &                          TRZNKE, TRZZ2, TRZRE,
     &                          TRZENERGY, TRZRANGE

C======================================================================C
