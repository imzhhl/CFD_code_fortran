!==============================================================================!
!              Optimal discrete sources based on global conductance            !
!              Direct Problem + R-Direct differentiation problem    	       !
!  SIMPLE algorithm + Staggered Grid System + Difference Schemes + Field Plot  !
!                                 Fu-Yun Zhao  +  Dong-Dong Zhang  2016/03/23  !
!==============================================================================!
   PROGRAM InverseConvection

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'
									  
    COMMON/DVARIABLES/UD(NX,NY),VD(NX,NY),TD(NX,NY)
    
	COMMON/ARDVARIABLES/ARUD(NX,NY),ARVD(NX,NY),ARTD(NX,NY)
    COMMON/BRDVARIABLES/BRUD(NX,NY),BRVD(NX,NY),BRTD(NX,NY)
    COMMON/CRDVARIABLES/CRUD(NX,NY),CRVD(NX,NY),CRTD(NX,NY)
!Anything is possible! Orientation is the best starting point!
                TD = 0.0;    ARTD = 0.0;    BRTD = 0.0 
			  TMAX = 0.0;  ARTMAX = 0.0;  BRTMAX = 0.0  
		       CDT = 0.0;   ARCDT = 0.0;   BRCDT = 0.0 
			  ACSD = 0.0;   ACSDR = 0.0;    AGCF = 0.0;  AGCFR = 0.0
			  BCSD = 0.0;   BCSDR = 0.0;    BGCF = 0.0;  BGCFR = 0.0

!Governing parameters for natural convection problem in enclosures [Open]

               Ra = 5.0E+4;  Pr = 0.7; theta = 90.0  

!Iteration process of Conjugte Gradient Method

             KMAX = 1000                             !KKK ranges from 0 to KMAX, and the maximum iteration loops for CGM are (55+1)

!Convergent criterion for CGM

        ConverCGM = 1.0E-6

!input Grid Data
              OPEN(1,FILE='D:\INVDSCM\Grid10060.dat',STATUS='UNKNOWN')
              READ(1,*) NI
              READ(1,*) NJ
              READ(1,101) (X(I),I=1,NI)
              READ(1,101) (Y(J),J=1,NJ)
             CLOSE(1)
101         FORMAT(5(F9.6,1X)) 

              NIM = NI-1;  NJM = NJ-1

            XC(1) = X(1)
           XC(NI) = X(NI)
             DO I = 2,NIM
            XC(I) = 0.5*(X(I)+X(I-1))
            DX(I) = X(I)-X(I-1)
             ENDDO
             DO I = 1,NIM
           HDX(I) = XC(I+1)-XC(I)
             ENDDO

            YC(1) = Y(1)
           YC(NJ) = Y(NJ)
             DO J = 2,NJM
            YC(J) = 0.5*(Y(J)+Y(J-1))
            DY(J) = Y(J)-Y(J-1)
             ENDDO
             DO J = 1,NJM
           HDY(J) = YC(J+1)-YC(J)
             ENDDO

            FX(1) = 0.0
             DO I = 2,NIM
            FX(I) = (X(I)-XC(I))/(XC(I+1)-XC(I))
             ENDDO
            FY(1) = 0.0
             DO J = 2,NJM
            FY(J) = (Y(J)-YC(J))/(YC(J+1)-YC(J))
             ENDDO

!           @@@@@@@@@@@  Conjugate Gradient Method  @@@@@@@@       
!

! 
!STEP A: Assume the unkonwn height of fin H
!

		 ! NIYLOWA = 2;     NIYHIGHA = 52   
		 ! NIYLOWB = 50;    NIYHIGHB = 101   

		       SLA = 0.00; SLB = 0.30; SLC = 0.75

DO KCG = 0, KMAX
	  	   
		    IF(KCG==0) THEN
	         NIYDA = 2
	         NIYDB = 31
!			 NIYDC = 76
		    NIYRDA = NIYDA + 1				! +/-   '-' from the top to the bottom  
	        NIYRDB = NIYDB + 1				! +/-   '+' from the bottom to the top 
!		    NIYRDC = NIYDC - 1				! +/-
			 ELSE											
		     NIYDA = NINT(SLA*100 + 1.00)
		     NIYDB = NINT(SLB*100 + 1.00)
!			 NIYDC = NINT(SLC*100 + 1.00)
       	    NIYRDA = NIYDA + 1				! +/-
		    NIYRDB = NIYDB + 1				! +/-
!			NIYRDC = NIYDC - 1				! +/-
		     ENDIF

	   
	!  IF((NIYDA.LE.NIYLOWA).OR.(NIYDA.GE.NIYHIGHA)) EXIT	
       	   	   
	   !Report the Iteration Loops of CGM
	         
			 WRITE(*,911) KCG
911         FORMAT(1X, 'CGM Iteration Loop', 1X, I6, 1X/)

	   !Solve the direct problem for temperature Tmax
	         
			 WRITE(*,*) ' Solve the Direct Problem'

!			 CALL DIRECT(KCG, Ra, Pr, NIYDA)                                   ! one heat source
			 CALL DIRECT(KCG, Ra, Pr, NIYDA, NIYDB)							   ! two heat sources
!			 CALL DIRECT(KCG, Ra, Pr, theta, NIYDA, NIYDB, NIYDC)	                   ! three heat sources
       
	   !Calculate the global conductance
	          TMAX = -1.E30
              DO I = 1,NI 
              DO J = 1,NJ 
              TMAX = MAX(TMAX,TD(I,J))
              ENDDO
              ENDDO
!	           CDT = 0.1/TMAX                                                  ! one heat source
			   CDT = (0.1*1+0.1*1)/TMAX										   ! two heat sources
!			   CDT = (0.1*1+0.1*1+0.1*1)/TMAX								   ! three heat sources
		    CDTINV = 1/CDT
	        WRITE(*,913) CDT
913         FORMAT(1X, 'DIRECT the Gloable Conductance', 1X, E15.6, 1X/)
	      	 
		 !Record the Hight and NIY
!		 CALL AUTOSAVE(KCG, NIYDA, SLA, CDT)				                   ! one heat source
		 CALL AUTOSAVE(KCG, NIYDA, NIYDB, SLA, SLB, CDT)					   ! two heat sources
!	     CALL AUTOSAVE(KCG, NIYDA, NIYDB, NIYDC, SLA, SLB, SLC, CDT)		   ! three heat sources

!
!STEP B-1: Calculate the search direction OF the SOURCE A
!

	   !Sovle the A-RDirect problem temperature RTmax
			 WRITE(*,*) ' Solve the ARDirect differentiation problem'  !THE +1 OF THE NIY1+1 IS THE DATA(Y)
             
!			 CALL ARDIRECT(KCG, Ra, Pr, NIYRDA)                                ! one heat source
			 CALL ARDIRECT(KCG, Ra, Pr, NIYRDA, NIYDB)						   ! two heat sources
!			 CALL ARDIRECT(KCG, Ra, Pr, theta, NIYRDA, NIYDB, NIYDC)            	   ! three heat sources
!	
	   !Calculate the global conductance + data DRT
	        ARTMAX = -1.E30
             DO I = 1,NI 
             DO J = 1,NJ 
            ARTMAX = MAX(ARTMAX,ARTD(I,J))
             ENDDO
             ENDDO
!	         ARCDT = 0.1/ARTMAX												   ! one heat source
			 ARCDT = (0.1*1+0.1*1)/ARTMAX									   ! two heat sources
!			 ARCDT = (0.1*1+0.1*1+0.1*1)/ARTMAX								   ! three heat sources
		  ARCDTINV = 1/ARCDT
	        WRITE(*,914) ARCDT
914         FORMAT(1X, 'AR-DIRECT the Gloable Conductance', 1X, E15.6, 1X/)


		   !define the Gradient of Cost Functional: GCF								   
		   	  AGCF = (ARCDTINV - CDTINV)/0.01
		   !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
            IF(KCG==0) THEN
             AGAMA = 0.0
             AGCFR = AGCF    !--For next CGM iteration
	       ELSE        
		     AGAMA = (AGCF/AGCFR)**2                       !((RRT-RT)/(RRTR-RTR))**2
             AGCFR = AGCF    !--For next CGM iteration
           ENDIF 
       	   	   
	      !define the Conjugate Search Direction CSD(KCG)
	       IF(KCG==0) THEN
	          ACSD = - AGCF
	         ACSDR =   ACSD  !--For next CGM iteration
	       ELSE
	          ACSD = - AGCF + AGAMA*ACSDR
	         ACSDR =   ACSD  !--For next CGM iteration
	       ENDIF

!STEP B-2: Calculate the search direction of the SOURCE B
!

	   !Sovle the A-RDirect problem temperature RTmax
			 WRITE(*,*) ' Solve the BRDirect differentiation problem' !THE +1 OF THE NIY1+1 IS THE DATA(Y)
!
             CALL BRDIRECT(KCG, Ra, Pr, NIYDA, NIYRDB)						   ! two heat sources
!			 CALL BRDIRECT(KCG, Ra, Pr, theta, NIYDA, NIYRDB, NIYDC)                  ! three heat sources                  
	
	   !Calculate the global conductance + data DRT
	        BRTMAX = -1.E30
              DO I = 1,NI 
              DO J = 1,NJ 
            BRTMAX = MAX(BRTMAX,BRTD(I,J))
              ENDDO
              ENDDO
			 BRCDT = (0.1*1+0.1*1)/BRTMAX									   ! two heat sources
!	         BRCDT = (0.1*1+0.1*1+0.1*1)/BRTMAX								   ! three heat sources
		  BRCDTINV = 1/BRCDT
	        WRITE(*,916) BRCDT
916         FORMAT(1X, 'BR-DIRECT the Gloable Conductance', 1X, E15.6, 1X/)


		   !define the Gradient of Cost Functional: GCF								   
		   	  BGCF = (BRCDTINV - CDTINV)/0.01
		   !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
           IF(KCG==0) THEN
             BGAMA = 0.0
             BGCFR = BGCF    !--For next CGM iteration
	       ELSE        
		     BGAMA = (BGCF/BGCFR)**2                       !((RRT-RT)/(RRTR-RTR))**2
             BGCFR = BGCF    !--For next CGM iteration
           ENDIF 

       	   	   
	      !define the Conjugate Search Direction CSD(KCG)
	       IF(KCG==0) THEN
	          BCSD = - BGCF
	         BCSDR =   BCSD  !--For next CGM iteration
	       ELSE
	          BCSD = - BGCF + BGAMA*BCSDR
	         BCSDR =   BCSD  !--For next CGM iteration
	       ENDIF

!!STEP B-3: Calculate the search direction of the SOURCE C
!!

!	   !Sovle the A-RDirect problem temperature RTmax
!			 WRITE(*,*) ' Solve the CRDirect differentiation problem'
             
!			 CALL CRDIRECT(KCG, Ra, Pr, theta, NIYDA, NIYDB, NIYRDC)            !THE +1 OF THE NIY1+1 IS THE DATA(Y)
	
!	   !Calculate the global conductance + data DRT
!	        CRTMAX = -1.E30
!              DO I = 1,NI 
!              DO J = 1,NJ 
!            CRTMAX = MAX(CRTMAX,CRTD(I,J))
!              ENDDO
!              ENDDO
!	         CRCDT = (0.1*1+0.1*1+0.1*1)/CRTMAX
!		  CRCDTINV = 1/CRCDT
!	        WRITE(*,918) CRCDT
!918         FORMAT(1X, 'CR-DIRECT the Gloable Conductance', 1X, E15.6, 1X/)


!		   !define the Gradient of Cost Functional: GCF								   
!		   	  CGCF = (CRCDTINV - CDTINV)/0.01
!		   !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
!           IF(KCG==0) THEN
!             CGAMA = 0.0
!             CGCFR = CGCF    !--For next CGM iteration
!	       ELSE        
!		     CGAMA = (CGCF/CGCFR)**2                       !((RRT-RT)/(RRTR-RTR))**2
!             CGCFR = CGCF    !--For next CGM iteration
!           ENDIF 

       	   	   
!	      !define the Conjugate Search Direction CSD(KCG)
!	       IF(KCG==0) THEN
!	          CCSD = - CGCF
!	         CCSDR =   CCSD  !--For next CGM iteration
!	       ELSE
!	          CCSD = - CGCF + CGAMA*CCSDR
!	         CCSDR =   CCSD  !--For next CGM iteration
!	       ENDIF
!
! 
!STEP C: Calculate the step size: ALFA 
!
			 
			IF(CDT.LT.ARCDT) THEN
              ALFA = 0.01
			ELSE
			  ALFA = 0.00000001
			ENDIF
			 
			IF(CDT.LT.BRCDT) THEN 
			  ALFB = 0.01
			ELSE
			  ALFB = 0.00000001
			ENDIF

!			IF(CDT.LT.CRCDT) THEN 
!			  ALFC = 0.001
!			ELSE
!			  ALFC = 0.00000001
!			ENDIF

! 
!STEP D: Update the location of the discrete heat source
!
		       SLA = SLA + ALFA*ACSD   ! +/-
			   SLB = SLB + ALFB*BCSD   ! +/-
!			   SLC = SLC - ALFC*CCSD   ! +/-
			
	   WRITE(*,912) SLA
912    FORMAT(1X, '(11)The location of Source-A could be', 1X, F11.6, 1X/)

	   WRITE(*,915) SLB
915    FORMAT(1X, '(12)The location of Source-B could be', 1X, F11.6, 1X/)        
		 
!	   WRITE(*,917) SLC
!917    FORMAT(1X, '(12)The location of Source-C could be', 1X, F11.6, 1X/)        

		 IF((CDT.GT.ARCDT).AND.(CDT.GT.BRCDT)) EXIT								 ! two heat sources
!		 IF((CDT.GT.ARCDT).AND.(CDT.GT.BRCDT).AND.(CDT.GT.CRCDT)) EXIT			 ! three heat sources

ENDDO
END

!* * * * * * * * * * * * * * * *  the End of the Inverse Convection Program * * * *	* * * * * * * *!

!
!Record the optimal location of the discrete heat source
!
!SUBROUTINE AUTOSAVE(LSAVE, NIYDA, SLA, CDT)							         ! one heat source
 SUBROUTINE AUTOSAVE(LSAVE, NIYDA, NIYDB, SLA, SLB, CDT)						 ! two heat sources
!SUBROUTINE AUTOSAVE(LSAVE, NIYDA, NIYDB, NIYDC, SLA, SLB, SLC, CDT)             ! three heat sources

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

! CHARACTER FIFlux*6, IC*1, JC*1, KC*1, LC*1

              NIM = NI-1;  NJM = NJ-1
              
              OPEN(5, FILE='D:\INVDSCM\Ra5E4.dat', STATUS='OLD', POSITION='APPEND')
!			 WRITE(5,501) LSAVE, NIYDA, SLA, CDT									! one heat source
			 WRITE(5,501) LSAVE, NIYDA, NIYDB, SLA, SLB, CDT						! two heat sources
!			 WRITE(5,501) LSAVE, NIYDA, NIYDB, NIYDC, SLA, SLB, SLC, CDT			! three heat sources
             CLOSE(5)
!501         FORMAT(1X, I6, 1X, I3, 1X, F11.4, 1X, F11.4)							! one heat source
 501         FORMAT(1X, I6, 1X, I3, 1X, I3, 1X, F11.4, 1X, F11.4, 1X, F11.4)		! two heat sources
!501         FORMAT(1X, I6, 1X, I3, 1X, I3, 1X, I3, 1X, F11.4, 1X, F11.4, 1X, F11.4, 1X, F11.4)	 ! three heat sources

RETURN
END    
