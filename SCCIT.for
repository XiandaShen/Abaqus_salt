C     SCCIT  
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C
      PARAMETER(IPHASEN=11,LSEED=200,PI=3.1415926D0)
C
      COMMON/PARAMETERS/AMUO,ANU,CELLR,VOIDM,VOIDSTD,DS,INTERVAL
C	  
C      PARAMETER(IPHASEN=11,INTERVAL=10,ITTIME=2000,DS=150.0D0,
C     1    AMUO=8846.0D0,ANU=0.3D0,CELLR=0.1375D0,VOIDM=0.103125D0,
C     2    VOIDSTD=0.0005D0,LSEED=200,PI=3.1415926D0)    
      DIMENSION TSTRAIN(3,3),DTSTRAIN(3,3),TSTRESS(3,3),
     1 PC(IPHASEN,3,3,3,3),FI(IPHASEN),DT(IPHASEN,IPHASEN,3,3,3,3),
     2 PSTRAIN(IPHASEN,3,3),PSTRESS(IPHASEN,3,3),PAO(IPHASEN,3,3,3,3),
     3 SAO(3,3,3,3),PA(IPHASEN,3,3,3,3),HC(3,3,3,3),
     4 EGSTRAIN(IPHASEN,3,3),PR(IPHASEN),STENSOR(3,3,3,3),
     5 PTENSOR(3,3,3,3),PSI(IPHASEN),THETA(IPHASEN),PHI(IPHASEN),
     6 Q(IPHASEN,3,3),CELLV(IPHASEN),SYMIF(3,3,3,3),AAP(3,3,3,3),
     7 ACAP(3,3,3,3),DTENSOR(IPHASEN,IPHASEN,3,3,3,3),
     8 ZEROF(3,3,3,3),DEGSTRAIN(IPHASEN,3,3),
     9 HCN(3,3,3,3),DHC(3,3,3,3)
C      
      DIMENSION TEMP21(3,3),TEMP22(3,3),TEMP23(3,3),
     1 TEMP25(3,3),TEMP26(3,3),TEMP27(3,3),TEMP28(3,3),TEMP24(3,3),
     2 TEMP41(3,3,3,3),TEMP42(3,3,3,3),TEMP43(3,3,3,3),TEMP44(3,3,3,3),
     3 TEMP45(3,3,3,3),TEMP46(3,3,3,3),TEMP47(3,3,3,3),TEMP48(3,3,3,3),
     4 TEMP49(3,3,3,3),TEMP410(3,3,3,3),TEMP411(3,3,3,3),
     5 TEMP412(3,3,3,3),TEMP413(3,3,3,3),TEMP414(3,3,3,3),
     6 TEMP415(3,3,3,3),TEMP416(3,3,3,3),TEMP417(3,3,3,3),
     7 TEMP418(3,3,3,3),TEMP419(3,3,3,3)
C	  REAL TVOID,TEMP11
C
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
C
C---------------------------------------------------------------------------
C 	                      MATERIAL PARAMETERS
C----------------------------------------------------------------------------
C
		AMUO=PROPS(1)
		ANU=PROPS(2)
		CELLR=PROPS(3)
		VOIDM=PROPS(4)
		VOIDSTD=PROPS(5)
		DS=PROPS(6)
C
C-----------------------------------------------------------------------
C                          INPUT INITIATION
C-----------------------------------------------------------------------
C    INPUT FOR EACH TIME STEP
		INTERVAL=DTIME
		CALL MAT1_MAT2(NTENS,STRAN,TSTRAIN,HALF)
		CALL MAT1_MAT2(NTENS,DSTRAN,DTSTRAIN,HALF)
		CALL MAT1_MAT2(NTENS,STRESS,TSTRESS,ONE)
C    FIRST TIME STEP EGSTRAIN INITIATION
      IF (TIME(2).EQ.0) THEN
		STATEV(:)=0.0D0
	  END IF
C
		DO I=1,IPHASEN
			EGSTRAIN(I,1,1)=STATEV(I*7-6)
			EGSTRAIN(I,2,2)=STATEV(I*7-5)	 
	 		EGSTRAIN(I,3,3)=STATEV(I*7-4)
			EGSTRAIN(I,1,2)=STATEV(I*7-3)	 
			EGSTRAIN(I,1,3)=STATEV(I*7-2)	 
			EGSTRAIN(I,2,3)=STATEV(I*7-1)	 
			CELLV(I)=STATEV(I*7)
		END DO
C
C    VOLUME FRACTION
      TEMP11=0.0D0
      DO I=1,IPHASEN
          PR(I)=1.0D0/IPHASEN
      END DO
      WRITE(6,*) 'PR=',PR
C
C	 VOID RADIUS INITATION
      IF (TIME(2).EQ.0) THEN  
        CELLV(1)=0.112D0
        CELLV(2)=0.110D0
        CELLV(3)=0.108D0
        CELLV(4)=0.106D0
        CELLV(5)=0.104D0
        CELLV(6)=0.102D0
        CELLV(7)=0.10D0
        CELLV(8)=0.098D0
        CELLV(9)=0.096D0
        CELLV(10)=0.094D0
		CELLV(11)=0.092D0
	  END IF
C	 	 
C    INCLUSION ORIENTATION	 
      DO I=1,IPHASEN
          PSI(I)=0.0D0
          THETA(I)=0.0D0
          PHI(I)=0.0D0
          CALL TRANSFORMATION(PSI(I),THETA(I),PHI(I),Q(I,:,:))
      END DO
C	 
      CALL SYMIDENDITYF(SYMIF)
      ZEROF(:,:,:,:)=0.0D0
C	  
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                      CALCULATION
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------      
C 	 
C  UPDATE TOTAL STRAIN
            CALL Aij_PLUS_Bij(TSTRAIN,DTSTRAIN,TSTRAIN)
C STIFFNESS FOR EACH INCLUSIONS
          DO I=1,IPHASEN
              FI(I)=(CELLV(I)/CELLR)**3
              CALL STIFFNESS(AMUO,FI(I),ANU,PC(I,:,:,:,:))
          END DO
C-----------------------------------------------------------------------
C HC CALCULATION
      ITERN=0
	  CALL STIFFNESS(AMUO,FI(I),ANU,HC)
      DO WHILE (0.0D0.LE.1.0D0)  
C	PTENSOR 
      CALL INVERSEFOURTH(HC,TEMP41)
      CALL ESHELBY(ANU,STENSOR)
      CALL Aijmn_Bmnkl(STENSOR,TEMP41,PTENSOR)
C
C CONSENTRATION TENSOR
          DO I=1,IPHASEN
              CALL Aijkl_MINUS_Bijkl(PC(I,:,:,:,:),PC(IPHASEN,:,:,:,:),
     1 TEMP41)
              CALL Aijmn_Bmnkl(PTENSOR,TEMP41,TEMP42)
              CALL Aijkl_PLUS_Bijkl(SYMIF,TEMP42,TEMP43)
              CALL INVERSEFOURTH(TEMP43,PAO(I,:,:,:,:))
          END DO
          SAO(:,:,:,:)=0.0D0
          DO I=1,IPHASEN
              CALL A_Bijkl(PR(I),PAO(I,:,:,:,:),3,TEMP41)
              CALL Aijkl_PLUS_Bijkl(SAO,TEMP41,SAO)
          END DO
          CALL INVERSEFOURTH(SAO,TEMP41)
          DO I=1,IPHASEN
              CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),TEMP41,PA(I,:,:,:,:))
          END DO
C HOMOGENIZED STIFFNESS
          HCN(:,:,:,:)=0.0D0
          DO I=1,IPHASEN
              CALL Aijmn_Bmnkl(PC(I,:,:,:,:),PA(I,:,:,:,:),TEMP42)
              CALL A_Bijkl(PR(I),TEMP42,3,TEMP43)
              CALL Aijkl_PLUS_Bijkl(HCN,TEMP43,HCN)
          END DO
          CALL Aijkl_MINUS_Bijkl(HCN,HC,DHC)
          DHCN=0.0D0
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      DO L=1,3
                  DHCN=DHCN+(DHC(I,J,K,L)**2)
                      END DO
                  END DO
              END DO
          END DO
          WRITE(111,*) DHCN
          IF (DHCN.LE.0.001D0) THEN
              GOTO 110
          END IF
          ITERN=ITERN+1
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      DO L=1,3
                  HC(I,J,K,L)=HCN(I,J,K,L)
                      END DO
                  END DO
              END DO
          END DO          
      END DO
110   WRITE(6,*) ITERN
C----------------------------------------------------------------
C TOTAL STRESS
          TEMP21(:,:)=0.0D0
          DO I=1,IPHASEN
              CALL Aijkl_Bkl(PC(I,:,:,:,:),EGSTRAIN(I,:,:),TEMP22)
              CALL Aij_Bijkl(TEMP22,PA(I,:,:,:,:),TEMP23)
              CALL A_Bij(PR(I),TEMP23,3,TEMP24)
              CALL Aij_MINUS_BijN(TEMP21,TEMP24,3,TEMP21)
          END DO
          IF (TIME(2).EQ.0) THEN 
              CALL Aijkl_Bkl(HC,TSTRAIN,TSTRESS)
          END IF
          CALL Aijkl_Bkl(HC,TSTRAIN,TEMP25)
          CALL Aij_PLUS_Bij(TEMP25,TEMP21,TSTRESS)
C          CALL Aij_PLUS_Bij(TSTRESS,TEMP26,TSTRESS)
C D TENSOR
          DO I=1,IPHASEN
              DO J=1,IPHASEN
                  IF (I.EQ.J) THEN
                      CALL A_Bijkl(PR(I),PA(I,:,:,:,:),3,TEMP41)
                      CALL Aijkl_MINUS_Bijkl(SYMIF,TEMP41,TEMP42)
                      CALL Aijmn_Bmnkl(TEMP42,PAO(I,:,:,:,:),TEMP43)
                      CALL Aijmn_Bmnkl(TEMP43,PTENSOR,TEMP44)
                      CALL Aijmn_Bmnkl(TEMP44,PC(I,:,:,:,:),
     1 DTENSOR(I,J,:,:,:,:))					
                  END IF 
                  IF (I.NE.J) THEN
                      CALL A_Bijkl(PR(J),PA(I,:,:,:,:),3,TEMP41)
                      CALL Aijkl_MINUS_Bijkl(ZEROF,TEMP41,TEMP42)
                      CALL Aijmn_Bmnkl(TEMP42,PAO(J,:,:,:,:),TEMP43)
                      CALL Aijmn_Bmnkl(TEMP43,PTENSOR,TEMP44)
                      CALL Aijmn_Bmnkl(TEMP44,PC(J,:,:,:,:),
     1 DTENSOR(I,J,:,:,:,:))					  
                  END IF 
              END DO
          END DO
C          
C WITH DIFFERENT P
C          AAP(:,:,:,:)=0.0D0
C          DO I=1,IPHASEN
C              CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),PTENSOR,TEMP41)
C              CALL A_Bijkl(PR(I),TEMP41,3,TEMP42)
C              CALL Aijkl_PLUS_Bijkl(AAP,TEMP42,AAP)
C          END DO
C          ACAP(:,:,:,:)=0.0D0
C          DO I=1,IPHASEN
C              CALL Aijkl_MINUS_Bijkl(HC,PC(I,:,:,:,:),TEMP41)
C              CALL Aijmn_Bmnkl(TEMP41,PAO(I,:,:,:,:),TEMP42)
C              CALL Aijmn_Bmnkl(TEMP42,PTENSOR,TEMP43)
C              CALL A_Bijkl(PR(I),TEMP43,3,TEMP44)
C              CALL Aijkl_PLUS_Bijkl(ACAP,TEMP44,ACAP)
C          END DO          
C          DO I=1,IPHASEN
C              DO J=1,IPHASEN
C                  IF (I.EQ.J) THEN
C                      CALL A_Bijkl(PR(I),PA(I,:,:,:,:),3,TEMP41)
C                      CALL Aijkl_MINUS_Bijkl(SYMIF,TEMP41,TEMP42)
C                      CALL Aijmn_Bmnkl(TEMP42,PAO(I,:,:,:,:),TEMP43)
C                      CALL Aijmn_Bmnkl(TEMP43,PTENSOR,TEMP44)
CC
C                      CALL Aijmn_Bmnkl(PA(I,:,:,:,:),AAP,TEMP45)
C                      CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),PTENSOR,TEMP46)
C                      CALL Aijkl_MINUS_Bijkl(TEMP45,TEMP46,TEMP47)
C                      CALL INVERSEFOURTH(ACAP,TEMP48)
C                      CALL Aijmn_Bmnkl(TEMP47,TEMP48,TEMP49)
C                      CALL Aijkl_MINUS_Bijkl(SYMIF,PA(I,:,:,:,:),
C     1 TEMP410)
C                      CALL TRANSPOSEF(TEMP410,TEMP411)
C                      CALL Aijkl_MINUS_Bijkl(HC,PC(I,:,:,:,:),TEMP412)
C                      CALL Aijmn_Bmnkl(TEMP412,PAO(I,:,:,:,:),TEMP413)
C                      CALL Aijmn_Bmnkl(TEMP413,PTENSOR,TEMP414)
C                      CALL Aijkl_PLUS_Bijkl(TEMP411,TEMP414,TEMP415)
C                      CALL A_Bijkl(PR(I),TEMP415,3,TEMP416)
C                      CALL Aijmn_Bmnkl(TEMP49,TEMP416,TEMP417)
C                      CALL Aijkl_PLUS_Bijkl(TEMP44,TEMP417,TEMP418)
C                      CALL Aijmn_Bmnkl(TEMP418,PC(I,:,:,:,:),
C     1 DTENSOR(I,J,:,:,:,:))
C                  END IF
C                  IF (I.NE.J) THEN
C                      CALL A_Bijkl(PR(J),PA(I,:,:,:,:),3,TEMP41)
C                      CALL Aijkl_MINUS_Bijkl(ZEROF,TEMP41,TEMP42)
C                      CALL Aijmn_Bmnkl(TEMP42,PAO(J,:,:,:,:),TEMP43)
C                      CALL Aijmn_Bmnkl(TEMP43,PTENSOR,TEMP44)
CC
C                      CALL Aijmn_Bmnkl(PA(I,:,:,:,:),AAP,TEMP45)
C                      CALL Aijmn_Bmnkl(PAO(I,:,:,:,:),PTENSOR,TEMP46)
C                      CALL Aijkl_MINUS_Bijkl(TEMP45,TEMP46,TEMP47)
C                      CALL INVERSEFOURTH(ACAP,TEMP48)
C                      CALL Aijmn_Bmnkl(TEMP47,TEMP48,TEMP49)
C                      CALL Aijkl_MINUS_Bijkl(SYMIF,PA(J,:,:,:,:),
C     1 TEMP410)
CC                      OPEN(111,FILE='DEBUG.OUT')
CC                     WRITE(111,*) TEMP410
C                      CALL TRANSPOSEF(TEMP410,TEMP411)
C                      CALL Aijkl_MINUS_Bijkl(HC,PC(J,:,:,:,:),TEMP412)
C                      CALL Aijmn_Bmnkl(TEMP412,PAO(J,:,:,:,:),TEMP413)
C                      CALL Aijmn_Bmnkl(TEMP413,PTENSOR,TEMP414)
C                      CALL Aijkl_PLUS_Bijkl(TEMP411,TEMP414,TEMP415)
C                      CALL A_Bijkl(PR(J),TEMP415,3,TEMP416)
C                      CALL Aijmn_Bmnkl(TEMP49,TEMP416,TEMP417)
C                      CALL Aijkl_PLUS_Bijkl(TEMP44,TEMP417,TEMP418)
C                      CALL Aijmn_Bmnkl(TEMP418,PC(J,:,:,:,:),
C     1 DTENSOR(I,J,:,:,:,:))
C                  END IF
C              END DO
C          END DO
C PHASE STRESS AND STRAIN
          DO J=1,IPHASEN
              TEMP21(:,:)=0.0D0
              DO I=1,IPHASEN
                  CALL Aijkl_Bkl(DTENSOR(J,I,:,:,:,:),EGSTRAIN(I,:,:),
     1 TEMP22)
                  CALL Aij_PLUS_Bij(TEMP21,TEMP22,TEMP21)
              END DO
              CALL Aijkl_Bkl(PA(J,:,:,:,:),TSTRAIN,TEMP22)
              CALL Aij_PLUS_Bij(TEMP21,TEMP22,PSTRAIN(J,:,:))
              CALL Aij_MINUS_BijN(PSTRAIN(J,:,:),EGSTRAIN(J,:,:),
     1 3,TEMP23)
              CALL Aijkl_Bkl(PC(J,:,:,:,:),TEMP23,PSTRESS(J,:,:))
          END DO 
          DO I=1,IPHASEN
              CALL CHEMOCELL(CELLV(I),DEGSTRAIN(I,:,:),PSTRESS(I,:,:),
     1 CELLR,Q(I,:,:),INTERVAL,DS)
	          CALL Aij_PLUS_Bij(EGSTRAIN(I,:,:),DEGSTRAIN(I,:,:),
     1 EGSTRAIN(I,:,:))
          END DO	 
C          OPEN(222,FILE='RVSTRAIN.OUT')
C          OPEN(333,FILE='RVSTRESS.OUT')
C          OPEN(444,FILE='RTVOID.OUT')
C          WRITE(6,*) TSTRAIN(3,3)
C          WRITE(6,*) TSTRESS(3,3) 
          TVOID=0.0D0
          DO I=1,IPHASEN
              TVOID=TVOID+4.0D0/3.0D0*PI*(CELLV(I)**3)
          END DO
		  WRITE(6,*) 'TIME=',TIME(2)		  
          WRITE(6,*) 'TVOID=',TVOID
		  WRITE(6,*) 'EGSTRAIN=',EGSTRAIN(1,:,:)
		  WRITE(6,*) 'PSTRESS1=',PSTRESS(1,:,:)
		  WRITE(6,*) 'PSTRAIN1=',PSTRAIN(1,:,:)
		  WRITE(6,*) 'CELLV1=',CELLV(1)
		  WRITE(6,*) 'PSTRESS11=',PSTRESS(11,:,:)
		  WRITE(6,*) 'PSTRAIN11=',PSTRAIN(11,:,:)
		  WRITE(6,*) 'CELLV11=',CELLV(11)
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                               OUTPUT
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C		  
      CALL MAT2_MAT1(NTENS,TSTRESS,STRESS,ONE)	 
      CALL MAT4_MAT2(HC,DDSDDE,2)
		DO I=1,IPHASEN
			STATEV(I*7-6)=EGSTRAIN(I,1,1)
			STATEV(I*7-5)=EGSTRAIN(I,2,2)	 
	 		STATEV(I*7-4)=EGSTRAIN(I,3,3)
			STATEV(I*7-3)=EGSTRAIN(I,1,2)	 
			STATEV(I*7-2)=EGSTRAIN(I,1,3)	 
			STATEV(I*7-1)=EGSTRAIN(I,2,3)	 
			STATEV(I*7)=CELLV(I)
		END DO
	 RETURN
	 END
C
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------	 
C
      SUBROUTINE MAT2_MAT1(NTENS,TENSOR,VECTOR,FACT)
C
C ====================================================================
C
C =================== MAT2_MAT1: TENSOR TO VECTOR=====================
C
C ====================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,NTENS
          VECTOR(I)=0.0D0
      END DO
C
      VECTOR( 1 ) = TENSOR(1 , 1)
      VECTOR( 2 ) = TENSOR(2 , 2)
      VECTOR( 3 ) = TENSOR(3 , 3)
      VECTOR( 4 ) = TENSOR(1 , 2)*FACT
      VECTOR( 5 ) = TENSOR(1 , 3)*FACT
      VECTOR( 6 ) = TENSOR(2 , 3)*FACT
C
      RETURN
      END
C	 
C-----------------------------------------------------------------------	 
      SUBROUTINE MAT1_MAT2(NTENS,VECTOR,TENSOR,FACT)
C ====================================================================
C 
C                MAT1_MAT2 : VECTOR TO TENOR  
C
C ====================================================================
       INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION VECTOR(NTENS),TENSOR(3,3)
C
      DO I=1,3
        DO J=1,3
          TENSOR(I,J)=0.0D0
        END DO
      END DO
C
      TENSOR(1 , 1) = VECTOR( 1 )
      TENSOR(2 , 2) = VECTOR( 2 )
        TENSOR(3 , 3) = VECTOR( 3 )
        TENSOR(1 , 2) = VECTOR( 4 )*FACT
        TENSOR(2 , 1) = TENSOR(1 , 2)
        TENSOR(1 , 3) = VECTOR( 5 )*FACT
        TENSOR(3 , 1) = TENSOR(1 , 3)
        TENSOR(2 , 3) = VECTOR( 6 )*FACT
        TENSOR(3 , 2) = TENSOR(2 , 3)
C
      RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C                STIFFNESS 
C
C ====================================================================
      SUBROUTINE STIFFNESS(AMUO,FI,ANU,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
        AMU=AMUO*(1.0D0-FI)
        ALAMDA=2.0D0*AMU*ANU/(1-2*ANU)
        DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
                    R(I,J,K,L)=0.0D0
                END DO
            END DO
        END DO
        END DO
        R(1,1,1,1)=2.0*amu+alamda
        R(2,2,2,2)=2.0*amu+alamda
        R(3,3,3,3)=2.0*amu+alamda
        R(1,1,2,2)=alamda
        R(1,1,3,3)=alamda
        R(2,2,1,1)=alamda
        R(2,2,3,3)=alamda
        R(3,3,1,1)=alamda
        R(3,3,2,2)=alamda
        R(2,3,2,3)=amu
        R(2,3,3,2)=amu
        R(3,2,3,2)=amu
        R(3,2,2,3)=amu
        R(2,1,2,1)=amu
        R(2,1,1,2)=amu
        R(1,2,2,1)=amu
        R(1,2,1,2)=amu
        R(1,3,1,3)=amu
        R(1,3,3,1)=amu
        R(3,1,1,3)=amu
        R(3,1,3,1)=amu
		RETURN
      END
C-----------------------------------------------------------------------
C ====================================================================
C 
C                INVERSEFOURTH 
C
C ====================================================================  
      SUBROUTINE INVERSEFOURTH(F,FINV)
            INCLUDE 'ABA_PARAM.INC' 
      DIMENSION F(3,3,3,3),FMATRIX(3,3),W(3,3),FINV(3,3,3,3) 
      FMATRIX(1,1)=F(1,1,1,1)
      FMATRIX(1,2)=F(1,1,2,2)
      FMATRIX(1,3)=F(1,1,3,3)
      FMATRIX(2,1)=F(2,2,1,1)
      FMATRIX(2,2)=F(2,2,2,2)
      FMATRIX(2,3)=F(2,2,3,3)
      FMATRIX(3,1)=F(3,3,1,1)
      FMATRIX(3,2)=F(3,3,2,2)
      FMATRIX(3,3)=F(3,3,3,3)
      CALL INVERSE(FMATRIX,3,3,W)
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      DO L=1,3
                          FINV(I,J,K,L)=0.0D0
                      END DO
                  END DO
              END DO
          END DO
          DO I=1,3
              DO J=1,3
              FINV(I,I,J,J)=W(I,J)
              END DO
          END DO
        FINV(1,2,1,2)=0.25/F(1,2,1,2)
        FINV(2,3,2,3)=0.25/F(2,3,2,3)
        FINV(1,3,1,3)=0.25/F(1,3,1,3)
        FINV(1,2,2,1)=FINV(1,2,1,2)
        FINV(2,1,2,1)=FINV(1,2,1,2)
        FINV(2,1,1,2)=FINV(1,2,1,2)
        FINV(3,2,3,2)=FINV(2,3,2,3)
        FINV(3,2,2,3)=FINV(2,3,2,3)
        FINV(2,3,3,2)=FINV(2,3,2,3)
        FINV(1,3,3,1)=FINV(1,3,1,3)
        FINV(3,1,1,3)=FINV(1,3,1,3)
        FINV(3,1,3,1)=FINV(1,3,1,3)
        return
      END
C-----------------------------------------------------------------------	  
C ====================================================================
C 
C                ESHELBY 
C
C ====================================================================  
      SUBROUTINE ESHELBY(ANU,R)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION R(3,3,3,3)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      DIJ=0.0D0
                      DKL=0.0D0
                      DIK=0.0D0
                      DJL=0.0D0
                      DIL=0.0D0
                      DJK=0.0D0
              IF(I.EQ.J) THEN
              DIJ=1.0D0
              END IF
              IF(K.EQ.L) THEN
              DKL=1.0D0
              END IF              
              IF(I.EQ.K) THEN
              DIK=1.0D0
              END IF              
              IF(J.EQ.L) THEN
              DJL=1.0D0
              END IF
              IF(I.EQ.L) THEN
              DIL=1.0D0
              END IF
              IF(J.EQ.K) THEN
              DJK=1.0D0
              END IF
              R(I,J,K,L)=(5.0D0*ANU-1.0D0)/15.0D0/(1.0D0-ANU)*DIJ*DKL
     1+(-5.0D0*ANU+4.0D0)/15.0D0/(1.0D0-ANU)*(DIK*DJL+DIL*DJK)
                  END DO
              END DO 
          END DO
      END DO
	  RETURN
      END
C-----------------------------------------------------------------------	  
        SUBROUTINE TRANSFORMATION(PSI,THETA,PHI,QQ)
CC====================================================================================================
C
C                          TRANSFORMATION MATRIX Q	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION PSI(1),THETA(1),PHI(1),QQ(3,3)
C
C
      QQ(1,1)=COS(PSI(1))*COS(THETA(1))*COS(PHI(1))
     1-SIN(PSI(1))*SIN(PHI(1))
      QQ(1,2)=SIN(PSI(1))*COS(THETA(1))*COS(PHI(1))
     1+COS(PSI(1))*SIN(PHI(1))
      QQ(1,3)=-SIN(THETA(1))*COS(PHI(1))
      QQ(2,1)=-COS(PSI(1))*COS(THETA(1))*SIN(PHI(1))
     1-SIN(PSI(1))*COS(PHI(1))
      QQ(2,2)=-SIN(PSI(1))*COS(THETA(1))*SIN(PHI(1))
     1+COS(PSI(1))*COS(PHI(1))
      QQ(2,3)=SIN(THETA(1))*SIN(PHI(1))
      QQ(3,1)=COS(PSI(1))*SIN(THETA(1))
      QQ(3,2)=SIN(PSI(1))*SIN(THETA(1))
      QQ(3,3)=COS(THETA(1))
C
      RETURN
      END
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          SYMIDENDITYF	  
C
C=====================================================================================================
C
      SUBROUTINE  SYMIDENDITYF(R)
            INCLUDE 'ABA_PARAM.INC'
         DIMENSION R(3,3,3,3)
         R(:,:,:,:)=0.0D0
        DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
        R(1,1,1,1)=1.0D0
        R(2,2,2,2)=1.0D0
        R(3,3,3,3)=1.0D0
        R(2,3,2,3)=0.5D0
        R(2,3,3,2)=0.5D0
        R(3,2,3,2)=0.5D0
        R(3,2,2,3)=0.5D0
        R(2,1,2,1)=0.5D0
        R(2,1,1,2)=0.5D0
        R(1,2,2,1)=0.5D0
        R(1,2,1,2)=0.5D0
        R(1,3,1,3)=0.5D0
        R(1,3,3,1)=0.5D0
        R(3,1,1,3)=0.5D0
        R(3,1,3,1)=0.5D0
                  END DO
              END DO
          END DO
        END DO
      RETURN
      end 
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          TRANSPOSEF	  
C
C=====================================================================================================
C	
      SUBROUTINE TRANSPOSEF(A,B)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION A(3,3,3,3),B(3,3,3,3)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
                      B(I,J,K,L)=A(K,L,I,J)
                  END DO
              END DO
          END DO
      END DO
      RETURN
      END 	
C-----------------------------------------------------------------------
CC====================================================================================================
C
C                          CHEMOCELL	  
C
C=====================================================================================================
C	  
      SUBROUTINE CHEMOCELL(CELLV,GVCE,GSIGMA,AAAA,Q,INTERVAL,DSK)
      INCLUDE 'ABA_PARAM.INC'
C  DIMENSION OF INPUT AND OUT PUT 
C TIME UNIT S, LENGTH UNIT MM, COMPRESSIVE STRESS AND WATER PRESSURE IS NEGATIVE
      DIMENSION  GVCE(3,3),GSIGMA(3,3),Q(3,3),
     1 CLVCE(3,3),CLSIGMA(3,3)
C DIMENSION OF VARIABLES IN THIS SUBROUTINE
      DIMENSION QGSIGMA(3,3),CVCE(3,3),CSIGMA(3,3),QLVCE(3,3)
     1 ,A(INTERVAL+1),CINA(INTERVAL+1)
      INTEGER INTERP,LNUM 
C INITIAL ZERO VECTOR 
C      INTERP=INTERVAL+1
      DIMENSION ZEROV(INTERVAL+1),
     1 TVC(INTERVAL+1),TVA(INTERVAL+1),TVB(INTERVAL+1),DR(INTERVAL)
C CONSTANTS
      DATA ZERO,HALF,ONE,TWO,FOUR,PI,DSO,IVCTEST
     1 /0.0D0,0.5D0,1.0D0,2.0D0,4.0D0,3.1415926D0,1.0D-9,0/
C WATER PRESSURE IS 1 ATM, 0.101MPA
      DATA POREPRESSURE /-0.101D0/
C CONCENTRATION mol/mm^3 R mJ/mol/K TEMPERATURE K
      DATA CPORE /6.48D-6/
      DATA R /8.37D3/
      DATA T /293.0D0/
C % mm^3/mol
      DATA OMEGA /2.7D4/  
C  INITIALIZATION
      CLSIGMA(:,:)=0.0D0
	  CLVCE(:,:)=0.0D0
C  LOCAL STRESS    
      CALL Aik_Bkj(Q,GSIGMA,QGSIGMA)            
      CALL Aik_Bkj(QGSIGMA,TRANSPOSE(Q),CLSIGMA)
C AXIAL EFFECTIVE STRESS OF CRACKS (MPA)
      SIGMA1=-CLSIGMA(1,1)*(AAAA**2)/(AAAA**2-CELLV**2)
     1-POREPRESSURE
      SIGMA2=-CLSIGMA(2,2)*(AAAA**2)/(AAAA**2-CELLV**2)
     1-POREPRESSURE
      SIGMA3=-CLSIGMA(3,3)*(AAAA**2)/(AAAA**2-CELLV**2)
     1-POREPRESSURE
C make DS to be e-18
      DS=DSK*DSO
C 
      CINA(1)=CELLV
      DO I=1,INTERVAL+1
          ZEROV(I)=0.0D0
          TVB(I)=0.0D0
          TVA(I)=0.0D0
          TVC(I)=0.0D0
      END DO
C
      DO I=1,INTERVAL
          A(I)=AAAA
          VI=-TWO*((A(I)**2)-(CINA(I)**2))*(OMEGA**2)*DS/R/T/FOUR*CPORE/
C     1    ((A(I)**4)*(CINA(I)-A(I))-(A(I)**2)*((CINA(I)**2)
     1    ((A(I)**4)*(LOG(CINA(I))-LOG(A(I)))-(A(I)**2)*((CINA(I)**2)
     2    -(A(I)**2))+((CINA(I)**4)-(A(I)**4))/FOUR)
C
          VVI=(A(I)**2)*PI-(CINA(I)**2)*PI
C DISOLUTION ON EACH CRACK
      IF (SIGMA3.GT.0.) THEN    
          VC=SIGMA3*VI
          VVC=VC*VVI
          TVC(I+1)=VC+TVC(I)
      ELSE
          VC=0.0D0
          VVC=0.0D0
          TVC(I+1)=VC+TVC(I)
      END IF
      IF (SIGMA1.GT.0.) THEN    
          VA=SIGMA1*VI
          VVA=VA*VVI
          TVA(I+1)=VA+TVA(I)
      ELSE
          VA=0.0D0
          VVA=0.0D0
          TVA(I+1)=VA+TVA(I)
      END IF
      IF (SIGMA2.GT.0.) THEN    
          VB=SIGMA2*VI
          VVB=VB*VVI
          TVB(I+1)=VB+TVB(I)
      ELSE
          VB=0.0D0
          VVB=0.0D0
          TVB(I+1)=VB+TVB(I)
      END IF
C 
      VVT=VVA+VVB+VVC
      CSURFACE=FOUR*PI*(CINA(I)**2)
      DR(I)=VVT/CSURFACE
C
      IF ((CINA(I)-DR(I)).LT.0) THEN
          IVCTEST=1
          CINA(I)=0.0000001D0
          LNUM=I
          EXIT
      END IF 
      A(I+1)=A(I)
      CINA(I+1)=CINA(I)-DR(I)
      LNUM=I+1
C      WRITE(111,*) vvt,SIGMA1,GSIGMA(1,1),QGSIGMA(1,1)
      END DO
C
        CRR=A(LNUM)-CINA(LNUM)
        CELLV=CINA(LNUM)
        CLVCE(1,1)=-TVA(LNUM)/AAAA
        CLVCE(2,2)=-TVB(LNUM)/AAAA
        CLVCE(3,3)=-TVC(LNUM)/AAAA
C      
      CALL Aik_Bkj(TRANSPOSE(Q),CLVCE,QLVCE)            
      CALL Aik_Bkj(QLVCE,Q,GVCE)
      DFI=(((AAAA-CRR)**3)-(CINA(I)**3))/(AAAA**3)
C 
      RETURN
      END

      SUBROUTINE MAT4_MAT2(TENSOR,DMATRIX,ICOE)
C
C ====================================================================
C                        MAT4_MAT2                                                                  I
C I        THIS PROGRAM TRANSFORMS THE FOURTH ORDER COMPLIANCE       I
C I        TENSOR TO A SECOND ORDER MATRIX                           I
C I                                                                  I
C ====================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION TENSOR(3,3,3,3),DMATRIX(6,6)
C
      DATA ZERO,TWO /0.0D0,2.0D0/
C
C     D2 = THE SECOND ORDER STIFFNESS MATRIX
C
      DO I=1,6
        DO J=1,6
          DMATRIX(I,J)=0.0D0
        END DO
      END DO

      IF (ICOE.EQ.1) THEN
         COE1=1.
         COE2=1.
      ELSEIF(ICOE.EQ.2) THEN
         COE1=2.
         COE2=4.
      END IF
C
      DMATRIX(1,1)=TENSOR(1,1,1,1)
      DMATRIX(1,2)=TENSOR(1,1,2,2)
      DMATRIX(1,3)=TENSOR(1,1,3,3)
      DMATRIX(1,4)=TENSOR(1,1,1,2)*COE1
      DMATRIX(1,5)=TENSOR(1,1,2,3)*COE1
      DMATRIX(1,6)=TENSOR(1,1,1,3)*COE1
C
      DMATRIX(2,1)=TENSOR(2,2,1,1)
      DMATRIX(2,2)=TENSOR(2,2,2,2)
      DMATRIX(2,3)=TENSOR(2,2,3,3)
      DMATRIX(2,4)=TENSOR(2,2,1,2)*COE1
      DMATRIX(2,5)=TENSOR(2,2,2,3)*COE1
      DMATRIX(2,6)=TENSOR(2,2,1,3)*COE1
C
      DMATRIX(3,1)=TENSOR(3,3,1,1)
      DMATRIX(3,2)=TENSOR(3,3,2,2)
      DMATRIX(3,3)=TENSOR(3,3,3,3)
      DMATRIX(3,4)=TENSOR(3,3,1,2)*COE1
      DMATRIX(3,5)=TENSOR(3,3,2,3)*COE1
      DMATRIX(3,6)=TENSOR(3,3,1,3)*COE1
C
      DMATRIX(4,1)=TENSOR(1,2,1,1)*COE1
      DMATRIX(4,2)=TENSOR(1,2,2,2)*COE1
      DMATRIX(4,3)=TENSOR(1,2,3,3)*COE1
      DMATRIX(4,4)=TENSOR(1,2,1,2)*COE2
      DMATRIX(4,5)=TENSOR(1,2,2,3)*COE2
      DMATRIX(4,6)=TENSOR(1,2,1,3)*COE2
C
      DMATRIX(5,1)=TENSOR(2,3,1,1)*COE1
      DMATRIX(5,2)=TENSOR(2,3,2,2)*COE1
      DMATRIX(5,3)=TENSOR(2,3,3,3)*COE1
      DMATRIX(5,4)=TENSOR(2,3,1,2)*COE2
      DMATRIX(5,5)=TENSOR(2,3,2,3)*COE2
      DMATRIX(5,6)=TENSOR(2,3,1,3)*COE2
C
      DMATRIX(6,1)=TENSOR(1,3,1,1)*COE1
      DMATRIX(6,2)=TENSOR(1,3,2,2)*COE1
      DMATRIX(6,3)=TENSOR(1,3,3,3)*COE1
      DMATRIX(6,4)=TENSOR(1,3,1,2)*COE2
      DMATRIX(6,5)=TENSOR(1,3,2,3)*COE2
      DMATRIX(6,6)=TENSOR(1,3,1,3)*COE2
C
C
      RETURN
      END
C
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------	  
      SUBROUTINE A_Bij(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION B(3,3),C(3,3)
      DO I=1,3
          DO J=1,3
          C(I,J)=A*B(I,J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE A_Bijkl(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION B(N,N,N,N),C(N,N,N,N)
      DO I=1,3
          DO J=1,3
              DO K=1,3
                  DO L=1,3
          C(I,J,K,L)=A*B(I,J,K,L)
                  END DO
              END DO
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bi(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N)
      C=0.0D0
      DO I=1,3
          C=C+A(I)*B(I)
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_Bj(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N),C(N,N)
      C(:,:)=0.0D0
      DO I=1,3
          DO J=1,3
          C(I,J)=A(I)*B(J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Ai_PLUS_Bi(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N),B(N),C(N)
      DO I=1,3
          C(I)=A(I)+B(I)
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_Bj(A,B,N,C)
      INCLUDE 'ABA_PARAM.INC' 
C      INTEGER N
      DIMENSION A(N,N)
      DIMENSION B(N),C(N)
      DO I=1,3
          C(I)=0.0D0
          DO J=1,3
          C(I)=C(I)+A(I,J)*B(J)
          END DO
      END DO 
      RETURN 
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_MINUS_BijN(A,B,N,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_BijN                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
C      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)-B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_PLUS_Bij(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aij_PLUS_BijN(A,B,N,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_BijN                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
C      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
      DO I=1,3
        DO J=1,3
          C(I,J)=A(I,J)+B(I,J)
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_Bkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aijkl_Bkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO I=1,3
        DO J=1,3
          C(I,J)=ZERO
          DO K=1,3
            DO L=1,3
              C(I,J)=C(I,J)+A(I,J,K,L)*B(K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_MINUS_Bijkl(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_MINUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
          C(I,J,K,L)=A(I,J,K,L)-B(I,J,K,L)
                END DO
            END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijkl_PLUS_Bijkl(A,B,C)
C====================================================================================
C                                                                      *
C                                  Aij_PLUS_Bij                        *
C                                                                      *
C====================================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
      DO I=1,3
        DO J=1,3
            DO K=1,3
                DO L=1,3
          C(I,J,K,L)=A(I,J,K,L)+B(I,J,K,L)
                END DO
            END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aijmn_Bmnkl(A,B,C)
CC====================================================================================================
C
C                          Aijmn_Bmnkl=Cijkl	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
C
C
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C(I,J,K,L)=0.0D0
              DO M=1,3
                DO N=1,3
            C(I,J,K,L)=C(I,J,K,L)+A(I,J,M,N)*B(M,N,K,L)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aik_Bkj(A,B,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION A(3,3),B(3,3),C(3,3)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
		RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE Aik_BkjN(A,B,N,C)
CC====================================================================================================
C
C                          Aik_Bkj=Cij	  
C
C=====================================================================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER N
      DIMENSION A(N,N),B(N,N),C(N,N)
C
C
      DO I=1,3
        DO J=1,3
          C(I,J)=0.0D0
          DO K=1,3
            C(I,J)=C(I,J)+A(I,K)*B(K,J)
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE Aij_Bijkl(A,B,C)
C========================================================================
C                                                                       =
C                              Aij_Bijkl                                =
C                                                                       =
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C
      DOUBLE PRECISION A(3,3),B(3,3,3,3),C(3,3)
      DATA ZERO /0.0D0/
C
      DO K=1,3
        DO L=1,3
          C(K,L)=0.0D0
          DO I=1,3
            DO J=1,3
              C(K,L)=C(K,L)+A(I,J)*B(I,J,K,L)
            END DO
          END DO
        END DO
      END DO
C
      RETURN
      END
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      SUBROUTINE INVERSE(A,N,NP,AINV)
C========================================================================
C
C    CALCULATE THE SECOND ORDER TENSOR A'S INVERSE, AINV
C    A^{-1} = AINV    
C    this subroutine inverses a (n x n) A matrix
C	 following a Gauss-Jordan elimination process
C
C========================================================================
      INCLUDE 'ABA_PARAM.INC'
C  
      DIMENSION A(NP,NP),IPIV(NP),INDXR(NP),INDXC(NP),
     1 A0(NP,NP),AINV(NP,NP)
C
      DO J=1,N
        IPIV(J)=0
      END DO
C
C
C     storage of the original A matrix
C
      DO I=1,N
        DO J=1,N
          A0(I,J)=A(I,J)
        END DO
      END DO
C
C	find a pivot among the rows of A that have not already been reduced
C
      DO I=1,N
        BIG=0.0D0
        DO J=1,N
          IF(IPIV(J).NE.1)THEN
            DO K=1,N
                IF(IPIV(K).EQ.0)THEN
                  IF(ABS(A(J,K)).GE.BIG)THEN
                   BIG=ABS(A(J,K))
                   IROW=J
                   ICOL=K
                   PIV=A(J,K)
                  END IF
                ELSEIF(IPIV(K).GT.1)THEN
                  write (7,*) 'Singular Matrix'
              END IF
            END DO
          END IF
        END DO
C
      IPIV(ICOL)=IPIV(ICOL)+1
      INDXR(I)=IROW
      INDXC(I)=ICOL
C	  
C     interchange the rows to put the pivot on the diagonal
C
      IF(irow.ne.icol)THEN
        DO L=1,N
          DUM=A(IROW,L)
          A(IROW,L)=A(ICOL,L)
          A(ICOL,L)=DUM
        END DO
      END IF
C     reduction of the row of the pivot
C       
      IF(PIV.EQ.0) write (7,*) 'Singular Matrix2'
C       
      PIVIN=1./PIV          ! numerical stabilization
C
      A(ICOL,ICOL)=1.       ! numerical stabilization
        DO L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVIN
        END DO
C
C     reduction of the column of the pivot
C
        DO LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.    ! numerical stabilization
            DO L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
            END DO
          END IF
        END DO
      END DO
C
C
C     unscramble the columns to get A-1
C		
      DO J=N,1,-1   ! reverse DO loop
        DO K=1,N
          DUM=A(K,INDXR(J))
          A(K,INDXR(J))=A(K,INDXC(J))
          A(K,INDXC(J))=DUM
        END DO
      END DO
C
C	restitution process of A and Ainv
C
      DO I=1,N
        DO J=1,N
          AINV(I,J)=A(I,J)
        END DO
      END DO
C
      DO I=1,N
        DO J=1,N
          A(I,J)=A0(I,J)
        END DO
      END DO
C     
      RETURN
      END
	  