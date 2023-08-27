!==============================================================================!
!        Optimal contaminant concentration based on average concentration      !
!              Direct Problem + R-Direct differentiation problem    	       !
!  SIMPLE algorithm + Staggered Grid System + Difference Schemes + Field Plot  !
!                                Dong-Dong Zhang   Fu-Yun Zhao     2016/09/17  !
!==============================================================================!
   PROGRAM InverseConvection

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'
									  
    COMMON/DVARIABLES/UD(NX,NY), VD(NX,NY), TD(NX,NY), CD(NX,NY)
    
	COMMON/ARDVARIABLES/ARUD(NX,NY),ARVD(NX,NY),ARTD(NX,NY),ARCD(NX,NY)
    COMMON/BRDVARIABLES/BRUD(NX,NY),BRVD(NX,NY),BRTD(NX,NY),BRCD(NX,NY)
    COMMON/CRDVARIABLES/CRUD(NX,NY),CRVD(NX,NY),CRTD(NX,NY),CRCD(NX,NY)

!Orientation is the best starting point!
                CD = 0.0;    ARCD = 0.0;    BRCD = 0.0 

			  ACSD = 0.0;   ACSDR = 0.0;    AGCF = 0.0;  AGCFR = 0.0
			  BCSD = 0.0;   BCSDR = 0.0;    BGCF = 0.0;  BGCFR = 0.0

!Governing parameters for mixed convection problem in enclosures [Open]

                Re = 200;      Pr = 0.71;     Gr = 1.0E+05 

!Iteration process of Conjugte Gradient Method

              KMAX = 1000                             !KKK ranges from 0 to KMAX, and the maximum iteration loops for CGM are (55+1)

!Convergent criterion for CGM

        ConverCGM = 1.0E-6

!input Grid Data
              OPEN(1,FILE='D:\Contaminant\Grid10x5.dat',STATUS='UNKNOWN')
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
!STEP A: Assume the positions of ventilation ports
!
			   ISX = 51;   ISY = 16;    SM = 1.0

		       SLA = 0.8622; SLB = 0.90   !; SLC = 1.60

DO KCG = 0, KMAX
	  	   
		    IF(KCG==0) THEN
	         NIYDA = 44
	         NIYDB = 46
		    NIYRDA = NIYDA + 1				! +/-   '-' from the top to the bottom  
	        NIYRDB = NIYDB - 1				! +/-   '+' from the bottom to the top 
			 ELSE											
		     NIYDA = NINT(SLA*50 + 1.00)
		     NIYDB = NINT(SLB*50 + 1.00)
       	    NIYRDA = NIYDA + 1				! +/-
		    NIYRDB = NIYDB - 1				! +/-
		     ENDIF

	   
	!  IF((NIYDA.LE.NIYLOWA).OR.(NIYDA.GE.NIYHIGHA)) EXIT	
       	   	   
	   !Report the Iteration Loops of CGM
	         
			 WRITE(*,911) KCG
911         FORMAT(1X, 'CGM Iteration Loop', 1X, I6, 1X/)

	   !Solve the direct problem for average global concentration
	         
			 WRITE(*,*) 'Solve the Direct Problem'

			 CALL DIRECT(KCG, Re, Pr, Gr, ISX, ISY, SM, NIYDA, NIYDB)							   ! two free ports
       
	   !Calculate the average global concentration
			    CAV = 0.0; CAVT = 0.0; SVE1 = 0.0
			   DO I = 2,NIM
			   DO J = 2,NJM
			    CAV = CD(I,J)*DX(I)*DY(J)
		       CAVT = CAVT + CAV
			   ENDDO
			   ENDDO
    		   SVE1 = (CAVT*0.1)/( SM*DX(ISX)*DY(ISY) * (X(NI)-X(1)) * (Y(NJ)-Y(1)) )      !1/CAVT

!The center of gravity for the contaminant distribution - XG, YG 
			     CX = 0.0; CY = 0.0; GX = 0.0; GY = 0.0; XG = 0.0; YG = 0.0
			   DO I = 2,NIM
			   DO J = 2,NJM
			     CX = CD(I,J)*XC(I)
				 CY = CD(I,J)*YC(J)
		         GX = GX + CX*DX(I)*DY(J)
				 GY = GY + CY*DY(J)*DX(I)
			   ENDDO
			   ENDDO
				 XG = GX/CAVT
				 YG = GY/CAVT 
!Mean radius of contaminant diffusion - SVE2
			   CXY = 0.0; SVE = 0.0; SVE2 = 0.0
			   DO J = 2,NJM
			   DO I = 2,NIM
				CXY = CD(I,J)*((XC(I)-XG)**2+(YC(J)-YG)**2)
				SVE = SVE + CXY*DY(J)*DX(I)
			   ENDDO
			   ENDDO
			   SVE2 = SQRT(SVE/CAVT) 	
!object function
          	  FOBJ = 0.0; SVE1OPT = 0.7794; SVE2OPT = 0.6094 

			  FOBJ = 0.5*SVE1/SVE1OPT+0.5*SVE2/SVE2OPT

	        WRITE(*,913) SVE1
913         FORMAT(1X, 'OBJECT fucntion of DIRECT', 1X, F11.6, 1X/)
	      	 
		 !Record the positions
!		 CALL AUTOSAVE(KCG, NIYDA, SLA, CDT)				                  
		 CALL AUTOSAVE(KCG, NIYDA, NIYDB, SLA, SLB, SVE1, SVE2, FOBJ)					  
!	     CALL AUTOSAVE(KCG, NIYDA, NIYDB, NIYDC, SLA, SLB, SLC, CAVT)		  

!
!STEP B-1: Calculate the search direction OF the ventration port A
!

	   !Sovle the A-RDirect problem Re-average global concentration
			 WRITE(*,*) ' Solve the ARDirect differentiation problem'  !THE +1 OF THE NIY1+1 IS THE DATA(Y)
             
!			 CALL ARDIRECT(KCG, Ra, Pr, NIYRDA)                                
			 CALL ARDIRECT(KCG, Re, Pr, Gr, ISX, ISY, SM, NIYRDA, NIYDB)						   
!			 CALL ARDIRECT(KCG, Re, Pr, Gr, NIYRDA, NIYDB, NIYDC)            
	
	   !Calculate the average global concentration + data DRT
			  ARCAV = 0.0; ARCAVT = 0.0; ARSVE1 = 0.0
			   DO I = 2,NIM
			   DO J = 2,NJM
			  ARCAV = ARCD(I,J)*DX(I)*DY(J)
		     ARCAVT = ARCAVT + ARCAV
			   ENDDO
			   ENDDO
		    ARSVE1 = (ARCAVT*0.1)/( SM*DX(ISX)*DY(ISY) * (X(NI)-X(1)) * (Y(NJ)-Y(1)) )     !1/ARCAVT	        

!The center of gravity for the contaminant distribution - ARXG, ARYG 
			     ARCX = 0.0; ARCY = 0.0; ARGX = 0.0; ARGY = 0.0; ARXG = 0.0; ARYG = 0.0
			   DO I = 2,NIM
			   DO J = 2,NJM
			     ARCX = ARCD(I,J)*XC(I)
				 ARCY = ARCD(I,J)*YC(J)
		         ARGX = ARGX + ARCX*DX(I)*DY(J)
				 ARGY = ARGY + ARCY*DY(J)*DX(I)
			   ENDDO
			   ENDDO
				 ARXG = ARGX/ARCAVT
				 ARYG = ARGY/ARCAVT 

!Mean radius of contaminant diffusion - ARSVE2
			   ARCXY = 0.0; ARSVE = 0.0; ARSVE2 = 0.0
			   DO J = 2,NJM
			   DO I = 2,NIM
				ARCXY = ARCD(I,J)*((XC(I)-XG)**2+(YC(J)-YG)**2)
				ARSVE = ARSVE + ARCXY*DY(J)*DX(I)
			   ENDDO
			   ENDDO
			   ARSVE2 = SQRT(ARSVE/ARCAVT) 		      	 

			WRITE(*,914) ARSVE1
914         FORMAT(1X, 'OBJECT fucntion of AR-DIRECT', 1X, F11.6, 1X/)

!object function
          	  ARFOBJ = 0.0; ARSVE1OPT = 0.7794; ARSVE2OPT = 0.6094 

			  ARFOBJ = 0.5*ARSVE1/ARSVE1OPT+0.5*ARSVE2/ARSVE2OPT
		     
		   !define the Gradient of Cost Functional: GCF								   
		   	  AGCF = (ARFOBJ - FOBJ)/0.02
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
!
!STEP B-2: Calculate the search direction of the ventilation port B
!

	   !Sovle the A-RDirect problem temperature RTmax
			 WRITE(*,*) ' Solve the BRDirect differentiation problem' !THE +1 OF THE NIY1+1 IS THE DATA(Y)

!            CALL BRDIRECT(KCG, Ra, Pr, NIYDA, NIYRDB)						          ! two heat sources
			 CALL BRDIRECT(KCG, Re, Pr, Gr, ISX, ISY, SM, NIYDA, NIYRDB)                  ! three heat sources                  
	
	   !Calculate the global conductance + data DRT
			  BRCAV = 0.0; BRCAVT = 0.0; BRSVE1 = 0.0
			   DO I = 2,NIM
			   DO J = 2,NJM
			  BRCAV = BRCD(I,J)*DX(I)*DY(J)
		     BRCAVT = BRCAVT + BRCAV
			   ENDDO
			   ENDDO
		    BRSVE1 = (BRCAVT*0.1)/( SM*DX(ISX)*DY(ISY) * (X(NI)-X(1)) * (Y(NJ)-Y(1)) )     !1/ARCAVT	        

!The center of gravity for the contaminant distribution - ARXG, ARYG 
			     BRCX = 0.0; BRCY = 0.0; BRGX = 0.0; BRGY = 0.0; BRXG = 0.0; BRYG = 0.0
			   DO I = 2,NIM
			   DO J = 2,NJM
			     BRCX = BRCD(I,J)*XC(I)
				 BRCY = BRCD(I,J)*YC(J)
		         BRGX = BRGX + BRCX*DX(I)*DY(J)
				 BRGY = BRGY + BRCY*DY(J)*DX(I)
			   ENDDO
			   ENDDO
				 BRXG = BRGX/BRCAVT
				 BRYG = BRGY/BRCAVT 

!Mean radius of contaminant diffusion - BRSVE2
			   BRCXY = 0.0; BRSVE = 0.0; BRSVE2 = 0.0
			   DO J = 2,NJM
			   DO I = 2,NIM
				BRCXY = BRCD(I,J)*((XC(I)-XG)**2+(YC(J)-YG)**2)
				BRSVE = BRSVE + BRCXY*DY(J)*DX(I)
			   ENDDO
			   ENDDO
			   BRSVE2 = SQRT(BRSVE/BRCAVT) 		      	 

	        WRITE(*,916) BRSVE1
916         FORMAT(1X, 'OBJECT fucntion of BR-DIRECT', 1X, F11.6, 1X/)

!object function
          	  BRFOBJ = 0.0; BRSVE1OPT = 0.7794; BRSVE2OPT = 0.6094 

			  BRFOBJ = 0.5*BRSVE1/BRSVE1OPT+0.5*BRSVE2/BRSVE2OPT
		     
		   !define the Gradient of Cost Functional: GCF								   
		   	  BGCF = (BRFOBJ - FOBJ)/0.02
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

! 
!STEP C: Calculate the step size: ALFA 
!
			 
			IF(FOBJ.GT.ARFOBJ) THEN
              ALFA = 0.005
			ELSE
			  ALFA = 0.00000001
			ENDIF
			 
			IF(FOBJ.GT.BRFOBJ) THEN 
			  ALFB = 0.005
			ELSE
			  ALFB = 0.00000001
			ENDIF

! 
!STEP D: Update the location of the discrete heat source
!
		       SLA = SLA + ALFA*ACSD   ! +/-
			   SLB = SLB - ALFB*BCSD   ! +/-
			
	   WRITE(*,912) SLA
912    FORMAT(1X, '(11)The location of Source-A could be', 1X, F11.6, 1X/)

	   WRITE(*,915) SLB
915    FORMAT(1X, '(12)The location of Source-B could be', 1X, F11.6, 1X/)        
		 

		 IF((FOBJ.LT.ARFOBJ).AND.(FOBJ.LT.BRFOBJ)) EXIT			 ! three ventilation ports

ENDDO
END

!* * * * * * * * * * * * * * * *  the End of the Inverse Convection Program * * * *	* * * * * * * *!

!
!Record the optimal location of the ventilation ports
!
!SUBROUTINE AUTOSAVE(LSAVE, NIYDA, SLA, CDT)							          ! one port
 SUBROUTINE AUTOSAVE(LSAVE, NIYDA, NIYDB, SLA, SLB, SVE1, SVE2, FOBJ)						  ! two ports
!SUBROUTINE AUTOSAVE(LSAVE, NIYDA, NIYDB, NIYDC, SLA, SLB, SLC, CAVT)             ! three ports

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

! CHARACTER FIFlux*6, IC*1, JC*1, KC*1, LC*1

              NIM = NI-1;  NJM = NJ-1
              
              OPEN(5, FILE='D:\INVCT_two\K11.dat', STATUS='OLD', POSITION='APPEND')
!			 WRITE(5,501) LSAVE, NIYDA, SLA, CDT									! one port
			 WRITE(5,501) LSAVE, NIYDA, NIYDB, SLA, SLB, SVE1, SVE2, FOBJ						! two ports
!			 WRITE(5,501) LSAVE, NIYDA, NIYDB, NIYDC, SLA, SLB, SLC, CAVT			! three ports
             CLOSE(5)
!501         FORMAT(1X, I6, 1X, I3, 1X, F11.4, 1X, F11.4)							! one port
501         FORMAT(1X, I6, 1X, I3, 1X, I3, 1X, 5(F11.4, 1X))		    ! two ports
!501         FORMAT(1X, I6, 1X, I3, 1X, I3, 1X, I3, 1X, F11.4, 1X, F11.4, 1X, F11.4, 1X, F11.4)	 ! three ports

RETURN
END    
