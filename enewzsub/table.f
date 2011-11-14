C======================================================================C
	SUBROUTINE INITTRZ
C======================================================================C

	IMPLICIT NONE

	INTEGER POINTNUM

	INCLUDE 'TABLE_RANGE_Z.INC'

C----------------------------------------------------------------------C

	TRZNTAB = 0

	TRZENERGY(0) = 0.
	TRZENERGY(1) = 0.02
	DO POINTNUM=2,TRZNPOINT
	   TRZENERGY(POINTNUM) = TRZENERGY(POINTNUM-1)*1.0437
	END DO
	
	RETURN
	END

C======================================================================C
	INTEGER FUNCTION CHECKTRZ(Z1,M1,NKE,Z2,RE)
C======================================================================C
C
C	INPUTS
C         Z1	 : Z OF A PARTICLE WHICH RUNS THROUGH MEDIA
C         M1	 : MASS (AMU) OF A PARTICLE 
C                    WHICH RUNS THROUGH MEDIA
C         NKE	 : NUMBER OF KINDS OF ELMENTS IN MEDIA
C         Z2	 : Z OF EACH ELEMENT IN MEDIA
C	  RE	 : RATIO OF EACH ELEMENT IN MEDIA
C
C	OUTPUTS
C	  CHECKTRZ : 0      = NOT EXIST
C                      OTHERS = TABLE NUMBER
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1, NKE, Z2
	REAL M1, RE
	DIMENSION Z2(*),RE(*)
	
	INTEGER TABNUM
	INTEGER ELNUM

	INCLUDE 'TABLE_RANGE_Z.INC'
	
C----------------------------------------------------------------------C

C***** INIT *****C

	CHECKTRZ = 0

C***** SEARCH TABLE *****C
	
	DO TABNUM = 1, TRZNTAB
	   IF ( ( Z1 .EQ. TRZZ1(TABNUM) )
     &	             .AND. ( M1 .EQ. TRZM1(TABNUM) )
     &	             .AND. ( NKE .EQ. TRZNKE(TABNUM)) ) THEN
	      DO ELNUM=1,NKE
		 IF ( Z2(ELNUM) .NE. TRZZ2(TABNUM,ELNUM) ) GOTO 500
		 IF ( RE(ELNUM) .NE. TRZRE(TABNUM,ELNUM) ) GOTO 500
	      END DO
	      CHECKTRZ = TABNUM
	      GOTO 200
	   END IF
 500	   CONTINUE
	END DO
	   
C***** RETURN *****C

 200	RETURN
	END
      
C======================================================================C
	INTEGER FUNCTION MAKETRZ(Z1,M1,NKE,Z2,RE)
C======================================================================C
C
C	INPUTS
C         Z1	 : Z OF A PARTICLE WHICH RUNS THROUGH MEDIA
C         M1	 : MASS (AMU) OF A PARTICLE 
C                    WHICH RUNS THROUGH MEDIA
C         NKE	 : NUMBER OF KINDS OF ELMENTS IN MEDIA
C         Z2	 : Z OF EACH ELEMENT IN MEDIA
C	  RE	 : RATIO OF EACH ELEMENT IN MEDIA
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1, NKE, Z2
	REAL M1, RE
	DIMENSION Z2(*),RE(*)

      
	INTEGER TABNUM
	INTEGER ELNUM
	INTEGER POINTNUM


	REAL DE2, DEDX1, DEDX2, DEDX3, DEDX4, DR

	REAL DEDXZ

	INCLUDE 'TABLE_RANGE_Z.INC'

C----------------------------------------------------------------------C

C***** INCREMENT TRZNTAB *****C

	TRZNTAB = TRZNTAB + 1
	TABNUM = TRZNTAB

C***** CHECK TABLE NUMBER *****C

	IF ( TRZNTAB .GT. TRZMAXNTAB ) THEN
	   TRZNTAB = TRZMAXNTAB
	   TABNUM = 1
	END IF

C***** CREATE TABLE *****C
      	
	TRZZ1(TABNUM) = Z1
	TRZM1(TABNUM) = M1

	TRZNKE(TABNUM) = NKE
	DO ELNUM=1,NKE
	   TRZZ2(TABNUM,ELNUM) = Z2(ELNUM)
	   TRZRE(TABNUM,ELNUM) = RE(ELNUM)
	END DO
      
	TRZRANGE(TABNUM,0)=0.
	DO POINTNUM = 1, TRZNPOINT
	   DE2 = ( TRZENERGY(POINTNUM) - TRZENERGY(POINTNUM-1) ) / 2.
	   DEDX1 = DEDXZ(Z1, M1, NKE, Z2, RE,
     &        ( TRZENERGY(POINTNUM-1) + 1.33998104*DE2 ) )
	   DEDX2 = DEDXZ(Z1, M1, NKE, Z2, RE,
     &        ( TRZENERGY(POINTNUM-1) + 1.86113631*DE2 ) )
	   DEDX3 = DEDXZ(Z1, M1, NKE, Z2, RE,
     &        ( TRZENERGY(POINTNUM-1) + 0.13886369*DE2 ) )
	   DEDX4 = DEDXZ(Z1, M1, NKE, Z2, RE,
     &        ( TRZENERGY(POINTNUM-1) + 0.66001869*DE2 ) )
	   DR = DE2 * ( 0.65214515/DEDX1 + 0.34785485/DEDX2
     &        + 0.34785485/DEDX3 + 0.65214515/DEDX4 )
	   TRZRANGE(TABNUM,POINTNUM) = TRZRANGE(TABNUM,POINTNUM-1) + DR
	END DO

	MAKETRZ = TABNUM

C***** RETURN *****C

	RETURN
	END

C======================================================================C
	REAL FUNCTION DEDXZ(Z1,M1,NKE,Z2,RE,E)
C======================================================================C
C
C	INPUTS	
C         Z1	: Z OF THE PARTICLE RUNNING THROUGH MATTER
C	  M1	: MASS (AMU) OF THE PARTICLE 
C  		              RUNNING THROUGH MATTER
C         NKE  	: NUMBER OF KINDS OF ELMENTS IN MATTER
C	  Z2	: Z OF EACH ELEMENT IN MATTER
C	  RE	: RATIO OF EACH ELEMENT IN MATTER
C	  E	: ENERGY (MEV/AMU) OF THE PARTICLE 
C                             RUNNING THROUGH MATTER
C	OUTPUTS
C	  DEDXZ : DE/DX ((MEV/AMU)/(MG/CM**2))
C
C----------------------------------------------------------------------C

	IMPLICIT NONE

	INTEGER Z1, NKE, Z2
	REAL M1, RE, E
	DIMENSION Z2(*),RE(*)

	INTEGER ELNUM

	REAL SE,SN

	REAL AW2,AWMATTER

	REAL DUMMY, DUMMY8
	DIMENSION DUMMY8(8)

C----------------------------------------------------------------------C

C***** INIT *****C

	AWMATTER = 0.
	DEDXZ = 0.

C***** CALULATE DEDXZ *****C

	DO ELNUM=1,NKE
	   CALL STOP(Z1, M1, Z2(ELNUM), E*1000.*M1, SE, SN, 2)
	   CALL SCOEF(Z2(ELNUM), DUMMY, DUMMY, AW2,
     &                 DUMMY, DUMMY, DUMMY, DUMMY, DUMMY8)
	   AWMATTER = AWMATTER + AW2 * RE(ELNUM)
	   DEDXZ = DEDXZ + (SE+SN)*AW2*RE(ELNUM)
	END DO
	
	DEDXZ = DEDXZ / ( AWMATTER * M1 )

C***** RETURN *****C

 	RETURN
	END

C======================================================================C
