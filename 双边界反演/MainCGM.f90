!==============================================================================!
!                       Inverse Natural Convection Problem                     !
!              Direct Problem + Sensitivity Problem + Adjoint Problem    	   !
!  SIMPLE algorithm + Staggered Grid System + Difference Schemes + Field Plot  !
!                         Fu-Yun Zhao  +  Di Liu 							   !
!                   Hunan University + Changsha + Aug-13-2007                  !        
!==============================================================================!
   PROGRAM InverseConvection

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/DVARIABLES/UD(NX,NY),VD(NX,NY),TD(NX,NY)
    
	COMMON/AVARIABLES1/UA1(NX,NY),VA1(NX,NY),TA1(NX,NY)
    COMMON/AVARIABLES2/UA2(NX,NY),VA2(NX,NY),TA2(NX,NY)
	
	COMMON/SVARIABLES1/US1(NX,NY),VS1(NX,NY),TS1(NX,NY)
    COMMON/SVARIABLES2/US2(NX,NY),VS2(NX,NY),TS2(NX,NY)
 
 DIMENSION Res1(NX,NY), Tm1(NX,NY), Res2(NX,NY), Tm2(NX,NY)
! Res: Residual of measurement values (Tm) and calculated ones (T)
! Tm:  Temperatures at measurement points or Obervations or Design points

 DIMENSION FIQ1(NX,NY),FIQ2(NX,NY)
! 'Inner Heat Sources' for adjoint problem

 DIMENSION Q1(NY), Qexact1(NY),Q2(NY),Qexact2(NY)
! Q: Unknown function
! Qexact: Exact function

 DIMENSION CSD1(NY), CSDR1(NY), GCF1(NY), GCFR1(NY), CSD2(NY), CSDR2(NY), GCF2(NY), GCFR2(NY)
! CSD CSDR: Conjugate Search Direction
! GCF GCFR: Gradient of the Cost Functional

 DIMENSION VCF1(106), VCF2(106), Qerror1(106), Qerror2(106)
! VCF: Value of Cost Functional
! Qerror: estimation error between the exact and the estimated functions


!Anything is possible! Orientation is the best starting point!
               TD = 0.0;     TA1 = 0.0;     TA2 = 0.0;    TS1 = 0.0;	TS1 = 0.0

              Res1 = 0.0;    Tm1 = 0.0;    Res2 = 0.0;    Tm2 = 0.0;   FIQ1 = 0.0; FIQ2 = 0.0

		  Qexact1 = 0.0;      Q1 = 0.0; Qexact2 = 0.0;     Q2 = 0.0
		     CSD1 = 0.0;   CSDR1 = 0.0;    GCF1 = 0.0;  GCFR1 = 0.0
		     CSD2 = 0.0;   CSDR2 = 0.0;    GCF2 = 0.0;  GCFR2 = 0.0
          Qerror1 = 0.0; Qerror2 = 0.0;    VCF1 = 0.0;   VCF2 = 0.0 


!Governing parameters for natural convection problem in enclosures [Open]
               Ra = 1.0E+6;  Pr = 0.72 

!Iteration process of Conjugte Gradient Method
             KMAX = 31    !KKK ranges from 0 to KMAX, and the maximum iteration loops for CGM are (55+1)

!Convergent criterion for CGM
        ConverCGM = 1.0E-9  

!input Grid Data
              OPEN(1,FILE='D:\CGM2\Grid50.dat',STATUS='UNKNOWN')
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
              OPEN(2,FILE='D:\DTP\YCQ1.dat',STATUS='UNKNOWN')
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (1)
              READ(2,102) YC(J), Qexact1(J)
             ENDDO
             CLOSE(2)
102         FORMAT(1X, F9.6, 1X, F11.4)

!input the exact function
             OPEN(21,FILE='D:\DTP\YCQ6.dat',STATUS='UNKNOWN')
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (2)
              READ(21,202) YC(J), Qexact2(J)
             ENDDO
             CLOSE(21)
202         FORMAT(1X, F9.6, 1X, F11.4)



!input the measurement values
!NOTE[1]: In the present work, they are obtained from direct convection problem.
!NOTE[2]: The sensors should be located in the interior of fluid field.



              Imw1 = 49; Ime1 = 49
              Jms =  2; Jmn = NJM
              OPEN(3,FILE='D:\TBDT\YCTML6R1E6-49.dat',STATUS='UNKNOWN')
            DO Im = Imw1,Ime1      
            DO Jm = Jms-1,Jmn+1  
              READ(3,103) YC(Jm), Tm1(Im,Jm)
            ENDDO
  		    ENDDO
             CLOSE(3)

              Imw2 = 4; Ime2 = 4
              Jms =  2; Jmn = NJM
              OPEN(3,FILE='D:\TBDT\YCTML6R1E6-4.dat',STATUS='UNKNOWN')
           DO Im = Imw2,Ime2      
           DO Jm = Jms-1,Jmn+1  
             READ(3,103) YC(Jm), Tm2(Im,Jm)
           ENDDO
		    ENDDO
            CLOSE(3)



103         FORMAT(1X, F9.6, 1X, F11.4)


!---Exact or expected or average temperature      Tex (From above Tm)
!---Standard deviation of the measurement errors  Delta
!---Random number of Gaussian distribution        OMIGA
           Delta = 0.0
!           Delta = 0.05
!            Delta = 0.10

!       ConverCGM = 50.0*(Delta**2.0)/2.0  !Number of sensors 50.0  

            DO Im = Imw1,Ime1      
            DO Jm = Jms,Jmn

           TOMIGA = 0.00

 	        DO NR = 1, 192
             CALL RANDOM(rand)
    	   TOMIGA = TOMIGA + rand
    	    ENDDO

    	    OMIGA = 2.576*(TOMIGA - 96.0)/96.0

         WRITE(*,*) OMIGA

		Tm1(Im,Jm) = Tm1(Im,Jm) + OMIGA*Delta

            ENDDO
 		    ENDDO


            DO Im = Imw2,Ime2      
            DO Jm = Jms,Jmn

           TOMIGA = 0.00

 	        DO NR = 1, 192
             CALL RANDOM(rand)
    	   TOMIGA = TOMIGA + rand
    	    ENDDO

    	    OMIGA = 2.576*(TOMIGA - 96.0)/96.0

         WRITE(*,*) OMIGA

		Tm2(Im,Jm) = Tm2(Im,Jm) + OMIGA*Delta

            ENDDO
 		    ENDDO

!             OPEN(33,FILE='F:\INCP\randTm_Case3_Ra1E+5_L090_Delta001.dat',STATUS='UNKNOWN')
!             OPEN(33,FILE='F:\INCP\randTm_Case3_Ra1E+5_L090_Delta005.dat',STATUS='UNKNOWN')
              OPEN(33,FILE='D:\DTP\Delta010.dat',STATUS='UNKNOWN')
            DO Im = Imw1,Ime1      
            DO Jm = Jms-1,Jmn+1  
             WRITE(33,1033) YC(Jm), Tm1(Im,Jm)
            ENDDO
  		    ENDDO
             CLOSE(33)
1033        FORMAT(1X, F9.6, 1X, F11.4)

              OPEN(34,FILE='D:\DTP\Delta011.dat',STATUS='UNKNOWN')
            DO Im = Imw2,Ime2      
            DO Jm = Jms-1,Jmn+1  
             WRITE(34,1034) YC(Jm), Tm2(Im,Jm)
            ENDDO
  		    ENDDO
             CLOSE(34)
1034        FORMAT(1X, F9.6, 1X, F11.4)

!
!           @@@@@@@@@@@  Conjugate Gradient Method  @@@@@@@@       
!

! 
!STEP A: Assume the unkonwn heat flux Q
!
                Q1 = 0.0 ; Q2 = 0.0


DO KCG = 0, KMAX

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



!************************prepare for the L1F************************************************************************************!
       
!	   !Error of exact function and estimated function
!          ERRORFZ2 = 0.0; ERRORFM2 = 0.0
!              DO J = 1,NJ				      !--------------------------Depending on the unknown function (2)
!          ERRORFZ2 = ERRORFZ2 + (Q2(J)+Qexact2(J))*(Q2(J)-Qexact2(J))
!          ERRORFM2 = ERRORFM2 + Qexact2(J)*Qexact2(J)
!             ENDDO
!      Qerror2(KCG) = ERRORFZ2/ERRORFM2 		  !---Preparing for the history plot{2}

!***********************************************************************************************************!
       
	   
	   !Solve the direct problem for temperature TD
	         WRITE(*,*) ' Solve the Direct Problem'
             CALL DIRECT(KCG, Ra, Pr, Q1, Q2)


       !Residual T - Tm
            DO Im = Imw1,Ime1
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
       Res1(Im,Jm) = TD(Im,Jm)- Tm1(Im,Jm)
		    ENDDO
			ENDDO  

       !Convergent criterion (Cost Functional/Performance Function/Error)
            STemp1 = 0.0
            DO Im = Imw1,Ime1
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!            DO Jm = Jms,Jmn, 10
		    STemp1 = STemp1 + Res1(Im,Jm)*Res1(Im,Jm)
		    ENDDO
			ENDDO
            STemp1 = STemp1/2.0
         VCF1(KCG) = STemp1                     !---Preparing for the history plot{2}

	         WRITE(*,*) 'Residual of cost functional1', STemp1

             IF(STemp1.LE.ConverCGM) EXIT	

            DO Im = Imw2,Ime2
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
       Res2(Im,Jm) = TD(Im,Jm)- Tm2(Im,Jm)
		    ENDDO
			ENDDO  

       !Convergent criterion (Cost Functional/Performance Function/Error)
            STemp2 = 0.0
            DO Im = Imw2,Ime2
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
		    STemp2 = STemp2 + Res2(Im,Jm)*Res2(Im,Jm)
		    ENDDO
			ENDDO
            STemp2 = STemp2/2.0
         VCF2(KCG) = STemp2                     !---Preparing for the history plot{2}

	         WRITE(*,*) 'Residual of cost functional2', STemp2

             IF(STemp2.LE.ConverCGM) EXIT

!--- Define the inner Heat Sources for Adjoint problem
             FIQ1 = 0.0
            DO Im = Imw1,Ime1
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
 	  FIQ1(Im,Jm) = Res1(Im,Jm)*DX(Im)*DY(Jm)
		    ENDDO
			ENDDO

             FIQ2 = 0.0
            DO Im = Imw2,Ime2
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
 	  FIQ2(Im,Jm) = Res2(Im,Jm)*DX(Im)*DY(Jm)
		    ENDDO
			ENDDO

!
!STEP B1: Calculate the Conjugate Search Direction: CSD CSDR
!
       !Sovle the adjoint problem backward in time for adjoint temperature TA
		     WRITE(*,*) ' Solve the Adjoint1 Problem'
             CALL ADJOINT1(KCG, Ra, Pr, FIQ1)

	   !define the Gradient of Cost Functional: GCF								   
   	         DO J = 2, NJM			          !--------------------------Depending on the unknown function (3)
          GCF1(J) = TA1(NI,J)
		     ENDDO

       !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
           IF(KCG==0) THEN
             GAMA1 = 0.0
             GCFR1 = GCF1    !--For next CGM iteration
	       ELSE
           GAMAFZ1 = 0.0;  GAMAFM1 = 0.0
              DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
!          GAMAFZ1 = GAMAFZ1 + GCF1(J)*(GCF1(J)-GCFR1(J))  !PR(P) CGM
!          GAMAFM1 = GAMAFM1 +GCFR1(J)*GCFR1(J)			  !PR(P) CGM
           GAMAFZ1 = GAMAFZ1 + GCF1(J)* GCF1(J)           ! Fletcher-Reeves CGM
           GAMAFM1 = GAMAFM1 + GCFR1(J)*GCFR1(J)		  ! Fletcher-Reeves CGM
	         ENDDO
             GAMA1 = GAMAFZ1/GAMAFM1
             GCFR1 = GCF1    !--For next CGM iteration
           ENDIF 

       !define the Conjugate Search Direction CSD(KCG)
	       IF(KCG==0) THEN
	          CSD1 = - GCF1
	         CSDR1 =   CSD1  !--For next CGM iteration
	       ELSE
	          CSD1 = - GCF1 + GAMA1*CSDR1
	         CSDR1 =   CSD1  !--For next CGM iteration
	       ENDIF


!
!STEP B2: Calculate the Conjugate Search Direction: CSD CSDR
!
       !Sovle the adjoint problem backward in time for adjoint temperature TA
		     WRITE(*,*) ' Solve the Adjoint2 Problem'
             CALL ADJOINT2(KCG, Ra, Pr, FIQ2)

	   !define the Gradient of Cost Functional: GCF
   	         DO J = 2, NJM			          !--------------------------Depending on the unknown function (3)
           GCF2(J) = TA2(NI,J)
		     ENDDO

       !Factor GAMA(KCG)/.or.BETA(KCG)/ (Factor for Search Direction)
           IF(KCG==0) THEN
             GAMA2 = 0.0
             GCFR2 = GCF2    !--For next CGM iteration
	       ELSE
           GAMAFZ2 = 0.0;  GAMAFM2 = 0.0
              DO J = 2, NJM			          !--------------------------Depending on the unknown function (4)
!          GAMAFZ2 = GAMAFZ2 + GCF2(J)*(GCF2(J)-GCFR2(J))  !PR(P) CGM
!          GAMAFM2 = GAMAFM2 +GCFR2(J)*GCFR2(J)			  !PR(P) CGM
           GAMAFZ2 = GAMAFZ2 + GCF2(J)* GCF2(J)           ! Fletcher-Reeves CGM
           GAMAFM2 = GAMAFM2 +GCFR2(J)*GCFR2(J)			  ! Fletcher-Reeves CGM
	         ENDDO
             GAMA2 = GAMAFZ2/GAMAFM2
             GCFR2 = GCF2    !--For next CGM iteration
           ENDIF 

       !define the Conjugate Search Direction CSD(KCG)
	       IF(KCG==0) THEN
	          CSD2 = - GCF2
	         CSDR2 =   CSD2  !--For next CGM iteration
	       ELSE
	          CSD2 = - GCF2 + GAMA2*CSDR2
	         CSDR2 =   CSD2  !--For next CGM iteration
	       ENDIF


! 
!STEP C1: Calculate the step size: ALFA1 
!

       !Solve the sensitivity problem for sensitivity temperature TS
		     WRITE(*,*) ' Solve the Sensitivity1 Problem'
             CALL SENSITIVITY1(KCG, Ra, Pr, CSD1)

       !define the Step Size ALFA (Prudhomme method)
           ALFAFZ1 = 0.0; ALFAFM1 = 0.0

            DO Im = Imw1,Ime1
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
           ALFAFZ1 = ALFAFZ1 - Res1(Im,Jm)*TS1(Im,Jm)
           ALFAFM1 = ALFAFM1 +  TS1(Im,Jm)*TS1(Im,Jm)
			ENDDO
 	        ENDDO

             ALFA1 = ALFAFZ1/ALFAFM1

!STEP C2: Calculate the step size: ALFA2 
!

       !Solve the sensitivity problem for sensitivity temperature TS
		     WRITE(*,*) ' Solve the Sensitivity2 Problem'
             CALL SENSITIVITY2(KCG, Ra, Pr, CSD2)

       !define the Step Size ALFA (Prudhomme method)
           ALFAFZ2 = 0.0; ALFAFM2 = 0.0

            DO Im = Imw2,Ime2
            DO Jm = Jms,Jmn
!           DO Jm = Jms,Jmn, 2
!           DO Jm = Jms,Jmn, 5
!           DO Jm = Jms,Jmn, 10
           ALFAFZ2 = ALFAFZ2 - Res2(Im,Jm)*TS2(Im,Jm)
           ALFAFM2 = ALFAFM2 +  TS2(Im,Jm)*TS2(Im,Jm)
			ENDDO
 	        ENDDO

             ALFA2 = ALFAFZ2/ALFAFM2


! 
!STEP D: Update the unknown heat flux
!
		        Q1 = Q1 + ALFA1*CSD1
		        Q2 = Q2 + ALFA2*CSD2

ENDDO



!
!Record the iterative history of CGM
!

              OPEN(4,FILE='E:\TBCGM2\L6R1E6ERROR.dat',STATUS='UNKNOWN')

             DO ICG = 0, KCG  !NOT KCG-1
             WRITE(4,104) ICG, Qerror1(ICG), Qerror2(ICG), VCF1(ICG), VCF2(ICG)
             ENDDO
             CLOSE(4)

104         FORMAT(1X,I6,1X,4(E15.6E3,2X))


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

 CHARACTER FIFlux*6, IC*1, JC*1, KC*1, LC*1


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
           FIFlux = 'UQ'//IC//JC//KC//LC	! Plot the inverse-calculated unknown heat flux


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
