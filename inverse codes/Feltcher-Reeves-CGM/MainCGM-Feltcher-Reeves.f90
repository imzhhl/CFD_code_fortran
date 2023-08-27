!==============================================================================!
!                       Inverse Natural Convection Problem                     !
!              Direct Problem + Sensitivity Problem + Adjoint Problem    	   !
!  SIMPLE algorithm + Staggered Grid System + Difference Schemes + Field Plot  !
!                         Fu-Yun Zhao  +  Di Liu 							   !
!                   Hunan University + Changsha + Aug-13-2007                  !        
!==============================================================================!
! Powell-Beale's Conjugate Gradient Method

   PROGRAM InverseConvection

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/DVARIABLES/UD(NX,NY),VD(NX,NY),TD(NX,NY)
    
	COMMON/AVARIABLES1/UA1(NX,NY),VA1(NX,NY),TA1(NX,NY)
    COMMON/AVARIABLES2/UA2(NX,NY),VA2(NX,NY),TA2(NX,NY)
	
	COMMON/SVARIABLES1/US1(NX,NY),VS1(NX,NY),TS1(NX,NY)
    COMMON/SVARIABLES2/US2(NX,NY),VS2(NX,NY),TS2(NX,NY)
	
	COMMON/CGMCoefficients/CSD1(KMAX,NY), GCF1(KMAX,NY), CSD2(KMAX,NY), GCF2(KMAX,NY)
! CSD CSDR: Conjugate Search Direction
! GCF GCFR: Gradient of the Cost Functional
 
! DIMENSION Res1(NX,NY), Tm1(NX,NY), Res2(NX,NY), Tm2(NX,NY)
! Res: Residual of measurement values (Tm) and calculated ones (T)
 DIMENSION  Tm1(NX,NY), Tm2(NX,NY)
! Tm:  Temperatures at measurement points or Obervations or Design points
 DIMENSION Tm(NX,NY), Res(NX,NY),FIQ(NX,NY)
! Res: Residual of measurement values (Tm) and calculated ones (T)
! Tm  Temperatures at ALL OF TEH measurement points
! Inner Heat Sources' for adjoint problem
 DIMENSION VCF(KMAX)
! Value of cost functional

! DIMENSION FIQ1(NX,NY),FIQ2(NX,NY)
! 'Inner Heat Sources' for adjoint problem

 DIMENSION Q1(NY), Qexact1(NY),Q2(NY),Qexact2(NY)
! Q: Unknown function
! Qexact: Exact function

 DIMENSION GAMA1(KMAX), GAMA2(KMAX), CHI1(KMAX), CHI2(KMAX) !CHI IS A GREEK ALPHABET SIMBLE LOOKS LIKE "X". GAMA MEANS GAMMA
 DIMENSION Qerror1(300), Qerror2(300)!, VCF1(300), VCF2(300), 
! VCF: Value of Cost Functional
! Qerror: estimation error between the exact and the estimated functions
 INTEGER Idelta
! Idelta denotes the difference between the two unknown boundaries' sensors along the direction 
!    which is perpendicular with the dirction along the doundary 

!Anything is possible! Orientation is the best starting point!
               TD = 0.0;     TA1 = 0.0;     TA2 = 0.0;    TS1 = 0.0;	TS1 = 0.0

!              Res1 = 0.0;    Tm1 = 0.0;    Res2 = 0.0;    Tm2 = 0.0;   FIQ1 = 0.0; FIQ2 = 0.0

		      Tm1 = 0.0;     Tm2 = 0.0
              
		  Qexact1 = 0.0;      Q1 = 0.0; Qexact2 = 0.0;     Q2 = 0.0
		     ! CSD1 = 0.0;   CSDR1 = 0.0;    GCF1 = 0.0;  GCFR1 = 0.0
		     ! CSD2 = 0.0;   CSDR2 = 0.0;    GCF2 = 0.0;  GCFR2 = 0.0
			 CSD1 = 0.0;    GCF1 = 0.0
		     CSD2 = 0.0;    GCF2 = 0.0
          Qerror1 = 0.0; Qerror2 = 0.0;    VCF1 = 0.0;   VCF2 = 0.0 
		    GAMA1 = 0.0;   GAMA2 = 0.0;    CHI1 = 0.0;   CHI2 = 0.0
			
			Res  = 0.0;    Tm  = 0.0;    FIQ  = 0.0;    VCF = 0.0

!Governing parameters for natural convection problem in enclosures [Open]
               Ra = 1.0E+5;  Pr = 0.72 

!Iteration process of Conjugte Gradient Method
!             KMAX = 299    !KKK ranges from 0 to KMAX, and the maximum iteration loops for CGM are (55+1)

!Convergent criterion for CGM
        ConverCGM = 1.0E-9  

!input Grid Data
              OPEN(1,FILE='grid80.dat',STATUS='UNKNOWN')
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


!input the exact function
              OPEN(2,FILE='Qexact1.dat',STATUS='UNKNOWN')
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (1)
              READ(2,102) YC(J), Qexact1(J)
             ENDDO
             CLOSE(2)
102         FORMAT(1X, F9.6, 1X, F11.4)

!input the exact function
             OPEN(21,FILE='Qexact2.dat',STATUS='UNKNOWN')
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (2)
              READ(21,202) YC(J), Qexact2(J)
             ENDDO
             CLOSE(21)
202         FORMAT(1X, F9.6, 1X, F11.4)



!input the measurement values
!NOTE[1]: In the present work, they are obtained from direct convection problem.
!NOTE[2]: The sensors should be located in the interior of fluid field.



              Imw1 = 80; Ime1 = 80
              Jms =  2; Jmn = NJM
              OPEN(31,FILE='T-measured80.dat',STATUS='UNKNOWN')
            DO Im = Imw1,Ime1 
            DO Jm = Jms-1,Jmn+1  
              READ(31,301) YC(Jm), Tm1(Im,Jm)
			   Tm(Im,Jm) = Tm1(Im,Jm)
            ENDDO
  		    ENDDO
             CLOSE(31)
301         FORMAT(1X, F9.6, 1X, F11.4)

              Imw2 = 2; Ime2 = 2
              Jms =  2; Jmn = NJM
              OPEN(32,FILE='T-measured2.dat',STATUS='UNKNOWN')
           DO Im = Imw2,Ime2      
           DO Jm = Jms-1,Jmn+1  
             READ(32,302) YC(Jm), Tm2(Im,Jm)
			 Tm(Im,Jm) = Tm2(Im,Jm)
           ENDDO
		    ENDDO
            CLOSE(32)
302         FORMAT(1X, F9.6, 1X, F11.4)

           Idelta =  Imw1 - Ime2


!---Exact or expected or average temperature      Tex (From above Tm)
!---Standard deviation of the measurement errors  Delta
!---Random number of Gaussian distribution        OMIGA
!          Delta = 0.0
            Delta = 0.05
!            Delta = 0.10
         ConverCGM = 1.20*(Delta**2.0)/2.0  !Number of sensors 160.0
      
	    DO Im = Imw2,Imw1, Idelta   
            DO Jm = Jms,Jmn
               TOMIGA = 0.00
 	           DO NR = 1, 192
                  CALL RANDOM(rand)
				  
    	          TOMIGA = TOMIGA + rand
    	       ENDDO
    	       OMIGA = 2.576*(TOMIGA - 96.0)/96.0
               WRITE(*,*) OMIGA
		       Tm(Im,Jm) = Tm(Im,Jm) + OMIGA*Delta
            ENDDO
 		ENDDO
!             OPEN(33,FILE='F:\INCP\randTm_Case3_Ra1E+5_L090_Delta001.dat',STATUS='UNKNOWN')
        OPEN(33,FILE='Delta00.dat',STATUS='UNKNOWN')
        DO Im = Imw2,Imw1, Idelta     
            DO Jm = Jms-1,Jmn+1  
               WRITE(33,1033) YC(Jm), Tm(Im,Jm)
            ENDDO
  		ENDDO
        CLOSE(33)
1033        FORMAT(1X, F9.6, 1X, F11.4)
!
!
!           @@@@@@@@@@@  Conjugate Gradient Method  @@@@@@@@       
! 
!STEP A: Assume the unkonwn heat flux Q
!
                Q1 = 0.0 ; Q2 = 0.0
!Record the iterative history of CGM
!
   OPEN(4,FILE='Ra1e5ERROR.dat',STATUS='UNKNOWN') 
!********************************************************************!
   KZ  = 1
DO KCG = 1, KMAX
   
       !Report the Iteration Loops of CGM
	         WRITE(*,911) KCG
911         FORMAT(1X, 'CGM Iteration Loop', 1X, I6, 1X/)

!Record the estimated heat flux
	         CALL AUTOSAVE(KCG, Q1, Q2)

!Error of exact function and estimated function
          ERRORFZ1 = 0.0; ERRORFM1 = 0.0
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (1)
          ERRORFZ1 = ERRORFZ1 + (Q1(J)-Qexact1(J))*(Q1(J)-Qexact1(J))
          ERRORFM1 = ERRORFM1 + Qexact1(J)*Qexact1(J)
             ENDDO
		  
      Qerror1(KCG) = ERRORFZ1/ERRORFM1 		  !---Preparing for the history plot{1}
    
!Error of exact function and estimated function
          ERRORFZ2 = 0.0; ERRORFM2 = 0.0
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (2)
          ERRORFZ2 = ERRORFZ2 + (Q2(J)-Qexact2(J))*(Q2(J)-Qexact2(J))
          ERRORFM2 = ERRORFM2 + Qexact2(J)*Qexact2(J)
             ENDDO
      Qerror2(KCG) = ERRORFZ2/ERRORFM2 		  !---Preparing for the history plot{2}

!***********************************************************************************************************!
!***********************************************************************************************************!
	   
!Solve the direct problem for temperature TD
	         WRITE(*,*) ' Solve the Direct Problem'
             CALL DIRECT(KCG, Ra, Pr, Q1, Q2)

!Residual T - Tm

	    DO Im = Imw2,Imw1, Idelta 
            DO Jm = Jms,Jmn
!             DO Jm = Jms,Jmn, 2
!             DO Jm = Jms,Jmn, 5
!             DO Jm = Jms,Jmn, 10
              Res(Im,Jm) = TD(Im,Jm)- Tm(Im,Jm)
		    ENDDO
	    ENDDO
!Convergent criterion (Cost Functional/Performance Function/Error)		
		STemp = 0.0
        DO Im = Imw2,Imw1, Idelta
            DO Jm = Jms,Jmn
!             DO Jm = Jms,Jmn, 2
!             DO Jm = Jms,Jmn, 5
!             DO Jm = Jms,Jmn, 10
		      STemp = STemp + Res(Im,Jm)*Res(Im,Jm)
		    ENDDO
		ENDDO
        STemp = STemp/2.0
        VCF(KCG) = STemp                     !---Preparing for the history plot{2}

	    WRITE(*,*) 'Residual of cost functional', STemp

        IF(STemp.LE.ConverCGM) EXIT
		
!--- Define the inner Heat Sources for Adjoint problem
        FIQ = 0.0
        DO Im = Imw2,Imw1, Idelta
            DO Jm = Jms,Jmn
!              DO Jm = Jms,Jmn, 2
!              DO Jm = Jms,Jmn, 5
!              DO Jm = Jms,Jmn, 10
 	           FIQ(Im,Jm) = Res(Im,Jm)*DX(Im)*DY(Jm)
		    ENDDO
		ENDDO

!
!STEP B1: Calculate the Conjugate Search Direction: CSD CSDR
!
       !Sovle the adjoint problem backward in time for adjoint temperature TA
		     WRITE(*,*) ' Solve the Adjoint1 Problem'
             CALL ADJOINT1(KCG, Ra, Pr, FIQ)
			 
!!!**********************************************************************************!!!!!!!!!!!!!!  	      
!!!! Fletcher-Reeves CGM
! !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
! !define the Gradient of Cost Functional: GCF
           DO J = 2, NJM			          !--------------------------Depending on the unknown function (3)
             GCF1(KCG,J) = TA1(NI,J)
		   ENDDO
		   IF(KCG==1) THEN
             GAMA1(KCG) = 0.0
			 DO J = 2, NJM
		       CSD1(KCG,J) = - GCF1(KCG,J)
		     ENDDO   
	       ELSE
             GAMAFZ1 = 0.0;  GAMAFM1 = 0.0
             DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
!              GAMAFZ1 = GAMAFZ1 + GCF1(KCG,J)*(GCF1(KCG,J)-GCF1(KCG-1,J))   !PR(P) CGM
!              GAMAFM1 = GAMAFM1 + GCF1(KCG-1,J)*GCF1(KCG-1,J)			     !PR(P) CGM
               GAMAFZ1 = GAMAFZ1 + GCF1(KCG,J)* GCF1(KCG,J)                  ! Fletcher-Reeves CGM
               GAMAFM1 = GAMAFM1 + GCF1(KCG-1,J)*GCF1(KCG-1,J)		         ! Fletcher-Reeves CGM
	         ENDDO
             GAMA1(KCG) = GAMAFZ1/GAMAFM1 
 !define the Conjugate Search Direction CSD(KCG)
             DO J = 2, NJM
	           CSD1(KCG,J) = - GCF1(KCG,J) + GAMA1(KCG)*CSD1(KCG-1,J) !+ CHI1(KZ)*CSD1(KZ,J)
	         ENDDO
			 
	       ENDIF
		   
!!!**********************************************************************************!!!!!!!!!!!!!!	
! !!!!!!!!! THE POWELL-BEALE'S CGM      
! !define the Gradient of Cost Functional: GCF	 
	 ! DO J = 2, NJM			       !--------------------------Depending on the unknown function (3)
        ! GCF1(KCG,J) = TA1(NI,J)
     ! ENDDO
! !Factor GAMA(KCG)/.or.CHI(KCG)/ (Factor for Search Direction) IN THE POWELL-BEALE'S CGM	 
	 ! IF(KCG==1) THEN
	     ! GAMA1(KCG) = 0.0
		 ! DO J = 2, NJM
		  ! CSD1(KCG,J) = - GCF1(KCG,J)
		 ! ENDDO   
	 ! ELSE
! ! TEST IF THE GRADIENTS AT SUCCESSIVE ITERATIONS TEND TO BE NON-ORTHOGONAL	   
	   ! VARIABLE1 = 0.0; VARIABLE2 = 0.0
	   ! DO J=2, NJM
	     ! VARIABLE1 = VARIABLE1 + GCF1(KCG-1,J)*GCF1(KCG,J)
	     ! VARIABLE2 = VARIABLE2 + GCF1(KCG  ,J)*GCF1(KCG,J)
	   ! ENDDO
	   
       ! IF((ABS(VARIABLE1)).GE.(0.2*VARIABLE2)) THEN
         ! KZ = KCG-1
       ! ENDIF
	   
       ! GAMAFZ1 = 0.0;  GAMAFM1 = 0.0
       ! DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
! !         GAMAFZ1 = GAMAFZ1 + GCF1(J)*(GCF1(J)-GCFR1(J))  !PR(P) CGM
! !         GAMAFM1 = GAMAFM1 +GCFR1(J)*GCFR1(J)			  !PR(P) CGM
         ! GAMAFZ1 = GAMAFZ1 + (GCF1(KCG,J) - GCF1(KCG-1,J))* GCF1(KCG,  J)       ! Powell-Beale's CGM
         ! GAMAFM1 = GAMAFM1 + (GCF1(KCG,J) - GCF1(KCG-1,J))* CSD1(KCG-1,J)		! Powell-Beale's CGM
	   ! ENDDO
		   
       ! GAMA1(KCG) = GAMAFZ1/GAMAFM1
		
	   ! IF((KZ.EQ.0).OR.(KCG.EQ.(KZ+1)))THEN
	      ! CHI1(KZ) = 0.0
	     ! ELSE
	       ! CHIFZ1 = 0.0;  CHIFM1 = 0.0
           ! DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
             ! CHIFZ1 = CHIFZ1 + (GCF1(KZ+1,J) - GCF1(KZ,J))* GCF1(KCG,J)      ! Powell-Beale's CGM
             ! CHIFM1 = CHIFM1 + (GCF1(KZ+1,J) - GCF1(KZ,J))* CSD1(KZ,J)	     ! Powell-Beale's CGM
	       ! ENDDO
           ! CHI1(KCG) = CHIFZ1/CHIFM1
	   ! ENDIF
	
	   ! DO J = 2, NJM
	     ! CSD1(KCG,J) = - GCF1(KCG,J) + GAMA1(KCG)*CSD1(KCG-1,J) + CHI1(KZ)*CSD1(KZ,J)
	   ! ENDDO
	   	   
! !!! TEST IF SUFFICIENTLY DOWNHILL
        
	   ! IF(KCG.GT.(KZ-1)) THEN	
          ! VARIABLE3 = 0.0 ; VARIABLE4 = 0.0; VARIABLE5 = 0.0
		  ! DO J = 2, NJM
		    ! VARIABLE3 = VARIABLE3 + GCF1(KCG,J)**2
		    ! VARIABLE4 = VARIABLE4 + CSD1(KCG,J)*GCF1(KCG,J)
   ! !		VARIABLE5 = VARIABLE5 + GCF1(KCG,J)**2
          ! ENDDO
		  ! IF(((-1.2*VARIABLE3).GE.VARIABLE4).AND.(VARIABLE4.GE.(-0.8*VARIABLE3))) THEN
		    ! KZ = KCG - 1
		    ! CHI1(KZ) = 0.0
			
			! DO J = 2, NJM
	          ! CSD1(KCG,J) = - GCF1(KCG,J) + GAMA1(KCG)*CSD1(KCG-1,J) + CHI1(KZ)*CSD1(KZ,J)
		    ! ENDDO
			
		  ! ENDIF
		! ENDIF
		
	  ! ENDIF
		


!
!STEP B2: Calculate the Conjugate Search Direction: CSD2(KCG)
!!!!!!!!!!!!!!!!!!!**********************************************************!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Fletcher-Reeves CGM
! !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
           DO J = 2, NJM			          !--------------------------Depending on the unknown function (3)
             GCF2(KCG,J) = TA1(1,J)
		   ENDDO
		   IF(KCG==1) THEN
             GAMA2(KCG) = 0.0
			 DO J = 2, NJM
		       CSD2(KCG,J) = - GCF2(KCG,J)
		     ENDDO   
	       ELSE
             GAMAFZ2 = 0.0;  GAMAFM2 = 0.0
             DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
!              GAMAFZ2 = GAMAFZ2 + GCF2(KCG,J)*(GCF2(KCG,J)-GCF2(KCG-1,J))      !  PR(P) CGM
!              GAMAFM2 = GAMAFM2 + GCF2(KCG-1,J)*GCF2(KCG-1,J)			        !  PR(P) CGM
               GAMAFZ2 = GAMAFZ2 + GCF2(KCG,J)* GCF2(KCG,J)                     ! Fletcher-Reeves CGM
               GAMAFM2 = GAMAFM2 + GCF2(KCG-1,J)*GCF2(KCG-1,J)		            ! Fletcher-Reeves CGM
	         ENDDO
             GAMA2(KCG) = GAMAFZ2/GAMAFM2 
 !define the Conjugate Search Direction CSD(KCG)
             DO J = 2, NJM
	           CSD2(KCG,J) = - GCF2(KCG,J) + GAMA2(KCG)*CSD2(KCG-1,J) !+ CHI2(KZ)*CSD2(KZ,J)
	         ENDDO
			 
	       ENDIF
!!!**********************************************************************************!!!!!!!!!!!!!!
! !!!!!!!!! THE POWELL-BEALE'S CGM
! !define the Gradient of Cost Functional: GCF
   	! DO J = 2, NJM			          !--------------------------Depending on the unknown function (3)
       ! GCF2(KCG,J) = TA1(1,J)
    ! ENDDO
! !Factor GAMA(KCG)/.or.CHI(KCG)/ (Factor for Search Direction) IN THE POWELL-BEALE'S CGM
    ! IF(KCG==1) THEN
	     ! GAMA2(KCG) = 0.0
		 ! DO J = 2, NJM
		  ! CSD2(KCG,J) = - GCF2(KCG,J)
		 ! ENDDO
	 ! ELSE

! !TEST IF THE GRADIENTS AT SUCCESSIVE ITERATIONS TEND TO BE NON-ORTHOGONAL
       
	   ! VARIABLE1 = 0.0; VARIABLE2 = 0.0
	   
	   ! DO J=2, NJM
	     ! VARIABLE1 = VARIABLE1 + GCF1(KCG-1,J)*GCF1(KCG,J)
	     ! VARIABLE2 = VARIABLE2 + GCF1(KCG  ,J)*GCF1(KCG,J)
	   ! ENDDO
	   
       ! IF((ABS(VARIABLE1)).GE.(0.2*VARIABLE2)) THEN
         ! KZ = KCG-1
       ! ENDIF
	   
       ! GAMAFZ2 = 0.0;  GAMAFM2 = 0.0
	   
       ! DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
         ! GAMAFZ2 = GAMAFZ2 + (GCF2(KCG,J) - GCF2(KCG-1,J))* GCF2(KCG,  J)       ! Powell-Beale's CGM
         ! GAMAFM2 = GAMAFM2 + (GCF2(KCG,J) - GCF2(KCG-1,J))* CSD2(KCG-1,J)		! Powell-Beale's CGM
	   ! ENDDO
	   
       ! GAMA2(KCG) = GAMAFZ2/GAMAFM2
	   
	   ! IF(KCG==(KZ+1))THEN
	      ! CHI2(KZ) = 0.0
	     ! ELSE
	       ! CHIFZ2 = 0.0;  CHIFM2 = 0.0
           ! DO J = 2, NJM			       !--------------------------Depending on the unknown function (4)
             ! CHIFZ2 = CHIFZ2 + (GCF2(KZ+1,J) - GCF2(KZ,J))* GCF2(KCG,J)        ! Powell-Beale's CGM
             ! CHIFM2 = CHIFM2 + (GCF2(KZ+1,J) - GCF2(KZ,J))* CSD2(KZ,J)         ! Powell-Beale's CGM
	       ! ENDDO
           ! CHI2(KCG) = CHIFZ2/CHIFM2
	   ! ENDIF
	
	   ! DO J = 2, NJM
	     ! CSD2(KCG,J) = - GCF2(KCG,J) + GAMA2(KCG)*CSD2(KCG-1,J) + CHI2(KZ)*CSD2(KZ,J)
	   ! ENDDO
			 
! !!! TEST IF SUFFICIENTLY DOWNHILL
        
	    ! IF(KCG.GT.(KZ-1)) THEN	
          ! VARIABLE3 = 0.0 ; VARIABLE4 = 0.0; VARIABLE5 = 0.0
		  ! DO J = 2, NJM
		    ! VARIABLE3 = VARIABLE3 + GCF2(KCG,J)**2
		    ! VARIABLE4 = VARIABLE4 + CSD2(KCG,J)*GCF2(KCG,J)
   ! !		VARIABLE5 = VARIABLE5 + GCF1(KCG,J)**2
          ! ENDDO
		  ! IF(((-1.2*VARIABLE3).GE.VARIABLE4).AND.(VARIABLE4.GE.(-0.8*VARIABLE3))) THEN
		    ! KZ = KCG - 1
		    ! CHI2(KZ) = 0.0
			
			! DO J = 2, NJM
	          ! CSD2(KCG,J) = - GCF2(KCG,J) + GAMA2(KCG)*CSD2(KCG-1,J) + CHI2(KZ)*CSD2(KZ,J)
		    ! ENDDO
			
		  ! ENDIF
		! ENDIF
    ! ENDIF
! 
!STEP C1: Calculate the step size: ALFA1 
!

       !Solve the sensitivity problem for sensitivity temperature TS
		     WRITE(*,*) ' Solve the Sensitivity1 Problem'
             CALL SENSITIVITY1(KCG, Ra, Pr)

       !define the Step Size ALFA (Prudhomme method)
           ALFAFZ1 = 0.0; ALFAFM1 = 0.0

            DO Im = Imw2,Imw1,Idelta
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
           ALFAFZ1 = ALFAFZ1 - Res(Im,Jm)*TS1(Im,Jm)
           ALFAFM1 = ALFAFM1 +  TS1(Im,Jm)*TS1(Im,Jm)
			ENDDO
 	        ENDDO

             ALFA1 = ALFAFZ1/ALFAFM1

!STEP C2: Calculate the step size: ALFA2 
!

       !Solve the sensitivity problem for sensitivity temperature TS
		     WRITE(*,*) ' Solve the Sensitivity2 Problem'
             CALL SENSITIVITY2(KCG, Ra, Pr)

       !define the Step Size ALFA (Prudhomme method)
           ALFAFZ2 = 0.0; ALFAFM2 = 0.0

            DO Im = Imw2,Imw1,Idelta
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
           ALFAFZ2 = ALFAFZ2 - Res(Im,Jm)*TS2(Im,Jm)
           ALFAFM2 = ALFAFM2 +  TS2(Im,Jm)*TS2(Im,Jm)
			ENDDO
 	        ENDDO

             ALFA2 = ALFAFZ2/ALFAFM2
			 
!  CORRECTED CALCULATION OF THE STEP SIZE. ALPHA1 & ALPHA2 
!  By Marcelo J.C and Helcio R.B. Orlande 2004 Int.J. Heat&Mass Transfer 

! !Solve the sensitivity problem for sensitivity temperature TS
        ! WRITE(*,*) ' Solve the Sensitivity1 Problem'
        ! CALL SENSITIVITY1(KCG, Ra, Pr)
! !Solve the sensitivity problem for sensitivity temperature TS
        ! WRITE(*,*) ' Solve the Sensitivity2 Problem'
        ! CALL SENSITIVITY2(KCG, Ra, Pr)

! ! BELOWS ARE THE COEFFICIENTS NEEDED FOR THE CALCULATION OF THE STEP SIZE
        ! V1 = 0.0 ;  V2 = 0.0 ;  V3 = 0.0 ;  V4 = 0.0 ;  V5 = 0.0
		
		! DO Im = Imw2,Imw1, Idelta
          ! DO Jm = Jms,Jmn
            ! V1 = V1 + Res(Im,Jm)*TS2(Im,Jm)
            ! V2 = V2 + TS1(Im,Jm)**2
            ! V3 = V3 + TS1(Im,Jm)*TS2(Im,Jm)
			! V4 = V4 + Res(Im,Jm)*TS1(Im,Jm)
			! V5 = V5 + TS2(Im,Jm)**2
		  ! ENDDO
		! ENDDO
		
		! ALPHA1 = (V1*V2-V3*V4)/(V5*V2-V3**2)
		! ALPHA2 = (V4*V5-V3*V1)/(V5*V2-V3**2)
!!!!!!!!!!!!THERE IS SOMETHING WRONG WITH THE RESULTS CALCULATED BY V1,V2,V3,V4,V5.		

! 
!STEP D: Update the unknown heat flux
!
         DO J = 2 , NJM
			
		    Q1(J) = Q1(J) + ALFA1*CSD1(KCG,J)
		    Q2(J) = Q2(J) + ALFA2*CSD2(KCG,J)
			
	     ENDDO
	! WRITE(*, 128) v1,v2,v3,v4,v5
! 128 FORMAT(5(E11.6E2,2X))	
    WRITE(*, 129) ALPHA1, ALPHA2
129 FORMAT(2(E15.6E3,2X))
    DO J = 2 , NJM
	  WRITE(*,130) CSD1(KCG,J),CSD2(KCG,J)
130 FORMAT(2(E15.6E3,2X))		
    ENDDO
!***************************************************************************************!
!Record the iterative history of CGM in File=4
             WRITE(4,104) KCG, Qerror1(KCG), Qerror2(KCG), VCF(KCG)

104         FORMAT(1X,I6,1X,5(E15.6E3,2X))

   KZ = KZ + 1
ENDDO
!***********************************************************!
! close the file for recording the iterative history of CGM
    CLOSE(4)
!****the end***************!
!
END
!* * * * * * * * * * * * * * * * * *  the End of the Inverse Convection Program     * *

!
!Record the unknown heat fluxes: Q
!
SUBROUTINE AUTOSAVE(LSAVE, Q1, Q2) 

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

 DIMENSION Q1(NY),Q2(NY)

 CHARACTER FIFlux*10, IC*1, JC*1, KC*1, LC*1


              NIM = NI-1;  NJM = NJ-1

!View loop numbers as output file name. ! The maximum loop is no more than 10000.
              III =  LSAVE/1000
	          JJJ = (LSAVE-1000*III)/100
	          KKK = (LSAVE-1000*III-100*JJJ)/10
	          LLL =  LSAVE-1000*III-100*JJJ-10*KKK
	           IC = CHAR(48+III)
	           JC = CHAR(48+JJJ)
	           KC = CHAR(48+KKK)
	           LC = CHAR(48+LLL)
           FIFlux = 'UQ'//IC//JC//KC//LC//'.dat'	! Plot the inverse-calculated unknown heat flux


!---Dynamic record the unknown heat flux
             Q1(1) = 0;  Q1(NJ) = 0          !Fit the curve at both corners 
		     Q2(1) = 0;  Q2(NJ) = 0							   !Present work would not concern the "end condition of adjoint problem"
              
			  OPEN(5, FILE=FIFlux, STATUS='UNKNOWN')
             DO J = 1, NJ
             WRITE(5,501) YC(J), Q1(J), Q2(J)
             ENDDO
             CLOSE(5)
501         FORMAT(1X, F9.6, 1X, F11.4,1X, F11.4)

RETURN
END    
