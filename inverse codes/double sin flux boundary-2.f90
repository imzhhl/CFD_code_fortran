! Program for Mixed convection in enclosures
!-----------                       Fu-Yun Zhao, G-B Yang
!---------------      Wuhan China  May 05 2019
! Flux boundary     
PROGRAM Laminar

 PARAMETER(NX=905,NY=905)
    COMMON/GRIDN/NI,NJ
    COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)
    COMMON/VARIABLES/ U(NX,NY), V(NX,NY), T(NX,NY)!,Qu(NX,NY)
    COMMON/PRESSURES/ P(NX,NY)
    DIMENSION Qexact1(NY), Qexact2(NY) !, TM(NX,NY)
!    DIMENSION HFP1(NY)  ! DEFINE THE HEAT FLUX FUNCTION  
	INTEGER Lin,Luq
	REAL  HFP1(NY), HFP2(NY), HFP3(NY), HFP4(NY), Qu1(NY), Qu2(NY) 
!	COMMON/HEAT_FLUX_PROFILE/ HFP1(NY), HFP2(NY), HFP3(NY), HFP4(NY), Qu1(NY), Qu2(NY)
	REAL Re, Gr, Pr
	CHARACTER*50 FILE1
!   serial number for saving the exact temperature measured
    CHARACTER*3 SNUM
    INTEGER I1,I2,I3
    CHARACTER I1C,I2C,I3C
!   arguments to determine the doundary heat flux profile
	INTEGER Q1, Q2
    CHARACTER FLAG
	
	 ! CHARACTER*20 FILE1, FILE2 !, FILE3, FILE4
     ! PRINT*,'ENTER GRID FILE NAME FOR THE CLACULATE CODES'
     ! READ(*,'(A20)') FILE1
     ! PRINT*,'ENTER OUTPUT FILE NAME FOR FIELD RESULTS'
     ! READ(*,'(A20)') FILE2
	 ! PRINT*,'ENTER OUTPUT FILE NAME FOR EXACT HEAT FLUX VALUE'
     ! READ(*,'(A20)') FILE3
	 ! PRINT*,'ENTER OUTPUT FILE NAME FOR MASUREED POINTS'
     ! READ(*,'(A20)') FILE4
	 ! PRINT*,'ENTER OUTPUT FILE NAME FOR EXACT HEAT FLUX VALUE'
     ! READ(*,'(A20)') FILE5
     ! OPEN(9,FILE=FILE2,STATUS='UNKNOWN')
	 WRITE(*,*) "Enter the Control Parameter Reynolds Number: Re "
	 Read(*,*) Re
	 WRITE(*,*) "Enter the Control Parameter Grashof Number: Gr "
	 READ(*,*) Gr
	 WRITE(*,*) "Specify the boundary heat flux, Q1 for the right wall,"
	 WRITE(*,*) "Q2 for the left wall : Q1, Q2"
	 READ(*,*) Q1, Q2
	 WRITE(*,*) "Enter the File Name for Fields Results: "
	 READ(*,*) FILE1
     WRITE(*,*) "Save the measured temperature or not, Y/N "
	 READ(*,*) FLAG
     
	
	 OPEN(9,FILE=TRIM(FILE1)//'Field.dat',STATUS='UNKNOWN')
	 
	 !OPEN(9,FILE='HFP1-HFP2.dat',STATUS='UNKNOWN')

!--- Input of grid data
              OPEN(0,FILE='grid.dat',STATUS='UNKNOWN')
!             OPEN(0,FILE=FILE1,STATUS='UNKNOWN')  
				READ(0,*) NI
                READ(0,*) NJ
                READ(0,100) (X(I),I=1,NI)
                READ(0,100) (Y(J),J=1,NJ)
               CLOSE(0)
100           FORMAT(5(F9.6,1X)) 

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

!--- Governing parameters [OPEN]
!                 Ra = 1.0E+6;     Pr = 0.71
!                  Ra = 1.0E+4;     Pr = 0.72
!				  XLe = 2.00;      DRT = 0.00000000010;      DRC = 0.00000000010
               
                Lin = 11; Luq = 91
                 ! Re = 200;  Gr = 1.0e+5
                 Pr = 0.72
				  
              OMIGU = 0.25;    OMIGV = 0.25;   OMIGPC = 0.95;    OMIGT = 0.85
!			  
! Res1， and  Res2 are guradances for convergence ?
!
               Res1 = 0.000001;   Res2 = 0.001 

               WRITE(*,101)
                READ(*,*) Inlooplimit
101           FORMAT(/1X,'Input the maximum inner loops for inital solutions (20-20000)')

!!!!!DEFINE THE FUNCTION FORM OF THE HEAT FLUX
        DO J = 1, NJ
		    HFP1(J) = - SIN(3.14159265*YC(J))
		ENDDO
		
		DO J = 1,NJ
		   IF( ( YC(J).GE.0.20 ) .AND. ( YC(J) .LE.0.8 ) ) THEN
		          HFP2(J) = -1.0
		   ELSE
		          HFP2(J) = 0.0
		   ENDIF
		ENDDO
		
		DO J = 1, NJ
		   IF( ( YC(J) .LE. 0.50 ) .AND. (YC(J) .GE. 0.0) ) THEN
		           HFP3(J) = -2*YC(J)
		   ELSE
		           HFP3(J) = 2*YC(J) - 2
		   ENDIF
		ENDDO
		
		DO J = 1,NJ
		   HFP4(J) = - SIN( 2*3.14159265*YC(J) ) 
        ENDDO
        
        Qu1 = 0.0
        SELECT CASE(Q1)
		  CASE(1)
		    Qu1 = HFP1
		  CASE(2)
		    Qu1 = HFP2
		  CASE(3)
		    Qu1 = HFP3
		  CASE(4)
		    Qu1 = HFP4
		  CASE DEFAULT
		    WRITE(*,*) "Bad Input for Q1"
          END SELECT
          
          Qu2 = 0.0
          
    !    SELECT CASE(Q2)
		  !CASE(1)
		  !  Qu2 = -1.0*HFP1
		  !CASE(2)
		  !  Qu2 = -1.0*HFP2
		  !CASE(3)
		  !  Qu2 = -1.0*HFP3
		  !CASE(4)
		  !  Qu2 = -1.0*HFP4
		  !CASE DEFAULT
		  !  WRITE(*,*) "Bad Input for Q2"
    !    END SELECT
          
          SELECT CASE(Q2)
		  CASE(1)
		    Qu2 = HFP1
		  CASE(2)
		    Qu2 = HFP2
		  CASE(3)
		    Qu2 = HFP3
		  CASE(4)
		    Qu2 = HFP4
		  CASE DEFAULT
		    WRITE(*,*) "Bad Input for Q2"
          END SELECT
          
!--- Difference schemes for convective terms           
          NSCHEMEUV = 6  !1 - HYBRID; 2 - POWERLAW; 3 - UPWIND(1st-order)
           NSCHEMET = 5  !1 - HYBRID; 2 - POWERLAW; 3 - UPWIND(1st-order)

                  U = 0.0;     V = 0.0;     P = 0.0;     T = 0.0

!---------------------------------------------------------------------------------------------------------
!--- Inner iteration by SIMPLE algorithms!
!--
                        DO INLOOP = 1, Inlooplimit

!$$$$$$$$$$$$$$$$$$$$$$User Defined Functions of the boundary conditions[OPEN]$$$$$$$$$$$$$$$$$$$$$$$$$$
!TOP - the known temperature wall
                DO I = 2,NIM
                    U(I, NJ) = +0.0
                    V(I,NJM) = +0.0
                    T(I, NJ) =  +0.0
				ENDDO

!BOTTOM -the known temperature floor
                DO I = 2,NIM
                    U(I,1) = +0.0
                    V(I,1) = +0.0
                    T(I,1) = +0.0
				ENDDO

!LEFT - the known Qk wall
                DO J = 2,Lin
					U(1,J) = 1.0
					V(1,J) = 0.0
					T(1,J) = 0.0
				ENDDO
							 
				DO J = Lin+1 , Luq
                    U(1,J) = 0.0
                    V(1,J) = 0.0
                    T(1,J) = T(2,J) + HDX(1) * Qu2(J)
                ENDDO
                             
                DO J = Luq, NJM
                    U(1,J) = U(2,J)
                    ! V(1,J) =  0.0  ! ???right?
                    !V(1,J) =  V(2,J)  ! ?? 
                    T(1,J) = T(2,J)
                    !T(1,J+1) = T(2,J+1)
                    ! Attention, stagged grid
                    P(2,J+1) = 0.0
                ENDDO
							 
				! The Free Outlet
				!DO J =  Luq , NJM
                    ! V(1,J) = V(2,J)
					! U(1,J) = U(2,J) - (V(2,J) - V(2,J-1)) * (DX(2) / DY(J))
					!IF(U(1,J)>0) THEN                           
					!   T(1,J) = 0.0
					!ELSE
					!    T(1,J) = T(2,J)
					!ENDIF
				!ENDDO

!RIGHT - the unknown Qu wall
                DO J = 2,NJM
                    U(NIM,J) = 0.0
                    V(NI, J) = 0.0
                    T(NI, J) = T(NIM,J) - HDX(NIM) * Qu1(J)
				ENDDO
				T(NI, 1) = T(NI,  2)
				T(NI,NJ) = T(NI,NJM)
						 
! this is the file for SAVING exact value for the unknown Q
! firstly, above all, calculate TEH Qexact,
		 ! DO J=1,NJ
		   ! Qexact1(J)=  HFP1(J)
		 ! ENDDO
		 Qexact1 = Qu1
		 ! Is it necessary to calculate the heat flux by  - ( T(NI,J) - T(NIM,J) )/ HDX(NIM) ???
! SAVE IT
         OPEN(13,FILE=TRIM(FILE1)//'Qexact1.dat',STATUS='UNKNOWN')
             DO J = 1,NJ				      !--------------------------Depending on the unknown function (1)
              WRITE(13,113) YC(J), Qexact1(J)
             ENDDO
         CLOSE(13)
113      FORMAT(1X, F9.6, 1X, F11.4)
! firstly, above all, calculate TEH Qexact,
		 ! DO J=Lin,Luq
		   ! Qexact2(J)=  HFP2(J)
		 ! ENDDO
		 Qexact2 = Qu2
! SAVE IT
         OPEN(14,FILE=TRIM(FILE1)//'Qexact2.dat',STATUS='UNKNOWN')
             DO J = Lin,Luq				      !--------------------------Depending on the unknown function (1)
              WRITE(14,114) YC(J), Qexact2(J)
             ENDDO
         CLOSE(14)
114      FORMAT(1X, F9.6, 1X, F11.4)


!$$$$$$$$$$$$$$$$$$$$$$User Defined Functions of the boundary conditions$$$$$$$$$$$$$$$$$$$$$$$$$$
			  
			   CALL   SUBU(ERROU, OMIGU, NSCHEMEUV, Re, Gr)

               CALL   SUBV(ERROV, OMIGV, NSCHEMEUV, Re, Gr)

               CALL   SUBPC(ERROPC, OMIGPC)

               CALL   SUBT(ERROT, OMIGT, NSCHEMET, Re, Pr)

               AMAX0 = AMAX1(ERROU, ERROV, ERROPC, ERROT)


!---Check the normal flow rates on the set planes 
!        !Global fluid flow flux on the Left side  
              FlowL = 0.00  
               DO J = 2,NJM
              FlowL = FlowL + U(1,J)*DY(J)
               ENDDO
!        !Global fluid flow flux on the Right side
              FlowR = 0.00  
               DO J = 2,NJM
              FlowR = FlowR + U(NIM,J)*DY(J)
               ENDDO
!!        !Global fluid flow flux on the Bottom side
              FlowB = 0.00  
               DO I = 2,NIM
              FlowB = FlowB + V(I,1)*DX(I)
               ENDDO
!        !Global fluid flow flux on the Top side
              FlowT = 0.00  
               DO I = 2,NIM
              FlowT = FlowT + V(I,NJM)*DX(I)
               ENDDO
            ResFLOW = FlowR  -  FlowL  +  FlowT  -  FlowB !flow flux Je - Jw + Jn - Js

!---Check the overall heat transfer rates
         !Calculate Global heat Flux along the Left side  
               TNuL = 0.0  
                  I = 1
               DO J = 2,NJM
               HFXX = (Re*Pr)*U(I,J)*T(I,J) - (T(I+1,J)-T(I,J))/HDX(I)
			 !  HFXX = U(I,J)*T(I,J) - (T(I+1,J)-T(I,J))/HDX(I)
               TNuL = TNuL + HFXX*DY(J)
               ENDDO
         !Calculate Global Heat Flux on the Right side
               TNuR = 0.0  
                  I = NIM
               DO J = 2,NJM
               HFXX = (Re*Pr)*U(I,J)*T(I+1,J) - (T(I+1,J)-T(I,J))/HDX(I)  
			  ! HFXX =U(I,J)*T(I+1,J) - (T(I+1,J)-T(I,J))/HDX(I) 
               TNuR = TNuR + HFXX*DY(J)
               ENDDO
         !Calculate Global Heat Flux on the Bottom side
               TNuB = 0.0
                  J = 1
               DO I = 2,NIM			
               HFYY = (Re*Pr)*V(I,J)*T(I,J) - (T(I,J+1)-T(I,J))/HDY(J)
			 ! HFYY = V(I,J)*T(I,J) - (T(I,J+1)-T(I,J))/HDY(J)
               TNuB = TNuB + HFYY*DX(I)
               ENDDO
         !Calculate Global Heat Flux on the Top side
               TNuT = 0.0  
                  J = NJM
               DO I = 2,NIM
               HFYY = (Re*Pr)*V(I,J)*T(I,J+1) - (T(I,J+1)-T(I,J))/HDY(J)
			  ! HFYY = V(I,J)*T(I,J+1) - (T(I,J+1)-T(I,J))/HDY(J)
               TNuT = TNuT + HFYY*DX(I)
               ENDDO
            ResHeat = TNuR  -  TNuL  +  TNuT  -  TNuB    !thermal flux Je - Jw + Jn - Js

               IF((AMAX0.LE.Res1).AND.(ABS(ResFLOW).LE.Res2).AND.(ABS(ResHeat).LE.Res2)) EXIT	

               IF(INLOOP==1.OR.MOD(INLOOP,1000)==0) THEN
               WRITE(*,102) INLOOP, ERROU, ERROV, ERROPC, ERROT, ABS(ResFLOW), ABS(ResHeat)
               ENDIF
102           FORMAT(1X, I7, 4X, 4(E8.3E1,2X), 4X, 2(E8.3E1,2X))

          ENDDO
!--
!--- Inner iteration by SIMPLE algorithms!
!---------------------------------------------------------------------------------------------------------
              OPEN(11,FILE = TRIM(FILE1)//'Statistics.dat',STATUS = 'UNKNOWN')
!!Record the emission rates along the heat sources - constant T type [OPEN] 
               ANuE = 0.00
               DO J = 1,NJ 				                                             ![OPEN]
                ! ANuS = ANuS + DY(J)*( T(NI,J) - T(NIM,J))/HDX(NIM) 
                ANuE = ANuE + DY(J)* T(NI,J)
               ENDDO
               TotalNuE = 1.00 / ANuE
           	   WRITE(*, *) 'Total heat transfer rate on the right wall: ', TotalNuE
               WRITE(11,*) 'Total heat transfer rate on the right wall: ', TotalNuE
               
               ANuW = 0.00
               DO J = Lin+1,Luq				                                             ![OPEN]
               
                ANuW = ANuW + DY(J)* T(1,J)
               ENDDO
               TotalNuW = 1.00 / ANuW
           	   WRITE(*, *) 'Total heat transfer rate on the left wall: ', TotalNuW
               WRITE(11,*) 'Total heat transfer rate on the left wall: ', TotalNuW

!Record the Flow rates across the free ventilated ports (Top side) [OPEN]
             FLowTP = 0.00;  FLowTN = 0.00  
               DO J =  Luq , NJM
             FLowTP = FLowTP + AMAX1(U(2,J),0.00)*DY(J)
             FLowTN = FLowTN - AMIN1(U(2,J),0.00)*DY(J)
               ENDDO
           FlowTopP = FLowTP
		   FlowTopN = FLowTN
           	   WRITE(*, *) 'Overall mass flowing out rate', FlowTopP
               WRITE(*, *) 'Overall mass flowing inn rate', FlowTopN
               
               WRITE(11,*) 'Overall mass flowing out rate', FlowTopP
               WRITE(11,*) 'Overall mass flowing inn rate', FlowTopN
        
        CLOSE(11)


!Plot the fluid flow and thermal/pollutant dispersions
               CALL  VISULIZATION(Re, Gr, Pr) 
			   
! this the file for the SAVING value for measured points
    IF((FLAG .EQ. 'Y') .OR. (FLAG .EQ. 'y')) THEN  
        DO Im = 1, NI
           I1 = Im / 100
           I2 = (Im - 100*I1) / 10
           I3 = Im - 100*I1 - 10*I2
           I1C = CHAR(48 + I1)
           I2C = CHAR(48 + I2)
           I3C = CHAR(48 + I3)
           SNUM = I1C//I2C//I3C
           OPEN(301,FILE = TRIM(FILE1)//SNUM//'.dat',STATUS = 'UNKNOWN')
           IF(XC(I).LE.(X(NI)/2.0)) THEN
             DO Jm = 1,NJ
               WRITE(301,115) YC(Jm), T(Im,Jm)
             ENDDO
           ELSE
             DO Jm = 1,NJ
               WRITE(301,115) YC(Jm), T(Im,Jm)
             ENDDO
           ENDIF
115        FORMAT(1X, F9.6, 1X, F11.4)
! 116      FORMAT(1X, F9.6, 1X, F9.6, 1X, F11.4)
           CLOSE(301)
        ENDDO
    ENDIF
            
!             Imw = 1; Ime = 1
!              Jms =  2; Jmn = NJM
!         OPEN(301,FILE='T-measured1.dat')
!            
!			DO Im = Imw,Ime     
!            DO Jm = Jms-1,Jmn+1  
!              WRITE(301,115) YC(Jm), T(Im,Jm)
!            ENDDO
!  		    ENDDO
!             CLOSE(301) 
!
!115         FORMAT(1X, F9.6, 1X, F11.4)
!! 116         FORMAT(1X, F9.6, 1X, F9.6, 1X, F11.4)

 READ(*,*)

END
!--* * * * * * * * * THE END OF THE PROGRAM MAIN     * * * * * * * * * *





SUBROUTINE SUBU(ERROU, OMIGU, NSCHEMEUV, Re, Gr)


     PARAMETER(NX=905,NY=905)
     COMMON/GRIDN/NI,NJ
     COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

     COMMON/VARIABLES/U(NX,NY), V(NX,NY), T(NX,NY) !,Qu(NX,NY), C(NX,NY)

     COMMON/PRESSURES/P(NX,NY)

     COMMON/OBSTACLESUV/AAAU(NX,NY),AAAV(NX,NY)

     COMMON/PRECOEFFC/BU(NX,NY), BV(NX,NY)

  DIMENSION AE(NX,NY), AW(NX,NY), AN(NX,NY), AS(NX,NY), SU(NX,NY), SOU(NX,NY)

!             COEFU = Pr             ! Form A
!              COEFU = SQRT(Pr/Ra)    ! Form B
              COEFU = 1.00/Re        ! Form C
			 
                NIM = NI-2;   NJM = NJ-1

DO I = 2,NIM
DO J = 2,NJM

           	     FE = (U(I,J)+U(I+1,J))*DY(J)/2.0
       	         DE = COEFU*DY(J)/DX(I+1)
            AE(I,J) =    DE*SCHEME(FE/DE,NSCHEMEUV) + AMAX1(-FE,0.0)

  	             FW = (U(I,J)+U(I-1,J))*DY(J)/2.0
  	             DW = COEFU*DY(J)/DX(I)
            AW(I,J) =    DW*SCHEME(FW/DW,NSCHEMEUV) + AMAX1( FW,0.0)

   	             FN = (DX(I+1)*V(I,J)+DX(I)*V(I+1,J))/2.0   
      	         DN = COEFU*HDX(I)/HDY(J)
            AN(I,J) =    DN*SCHEME(FN/DN,NSCHEMEUV) + AMAX1(-FN,0.0)

	             FS = (DX(I+1)*V(I,J-1)+DX(I)*V(I+1,J-1))/2.0 
  	             DS = COEFU*HDX(I)/HDY(J-1)
            AS(I,J) =    DS*SCHEME(FS/DS,NSCHEMEUV) + AMAX1( FS,0.0)


                SAD = 0.0

!---CDS(Linear Interpolation)
        IF(NSCHEMEUV==4) THEN
 	          FECDS =             FE*0.5*U(I,J) + FE*0.5*U(I+1,J)  
 	          FEUDS =      AMAX1(FE,0.0)*U(I,J) + AMIN1(FE,0.0)*U(I+1,J)
	          FWCDS =           FW*0.5*U(I-1,J) + FW*0.5*U(I,J)
	          FWUDS =    AMAX1(FW,0.0)*U(I-1,J) + AMIN1(FW,0.0)*U(I,J)
	          FNCDS =     FN*(1.0-FY(J))*U(I,J) + FN*FY(J)*U(I,J+1)  
	          FNUDS =      AMAX1(FN,0.0)*U(I,J) + AMIN1(FN,0.0)*U(I,J+1)
	          FSCDS = FS*(1.0-FY(J-1))*U(I,J-1) + FS*FY(J-1)*U(I,J)
	          FSUDS =    AMAX1(FS,0.0)*U(I,J-1) + AMIN1(FS,0.0)*U(I,J)
                SAD =-FECDS + FEUDS + FWCDS - FWUDS - FNCDS + FNUDS + FSCDS - FSUDS
        ENDIF

!---Second-order Upwind Scheme
        IF(NSCHEMEUV==5)THEN
              SADEP = -(U(I,J)  -U(I-1,J))*AMAX1(+FE,0.0) 
        IF(I==NIM) THEN
              SADEN = 0.00
        ELSE 
              SADEN = +(U(I+1,J)-U(I+2,J))*AMAX1(-FE,0.0)
        ENDIF  	  	
        IF(I==2) THEN
              SADWP = 0.00
        ELSE
              SADWP = +(U(I-1,J)-U(I-2,J))*AMAX1(+FW,0.0)
        ENDIF
              SADWN = -(U(I,J)  -U(I+1,J))*AMAX1(-FW,0.0)
              SADNP = -(U(I,J)  -U(I,J-1))*AMAX1(+FN,0.0)
        IF(J==NJM) THEN
              SADNN = 0.00
        ELSE 
              SADNN = +(U(I,J+1)-U(I,J+2))*AMAX1(-FN,0.0)
        ENDIF 	  	
        IF(J==2) THEN
              SADSP = 0.00
        ELSE
              SADSP = +(U(I,J-1)-U(I,J-2))*AMAX1(+FS,0.0)
        ENDIF
              SADSN = -(U(I,J)  -U(I,J+1))*AMAX1(-FS,0.0)
	            SAD = (1.0/2.0)*(SADEP + SADEN + SADWP + SADWN + SADNP + SADNN + SADSP + SADSN)
        ENDIF

!---Deferrd correction QUICK Scheme
        IF(NSCHEMEUV==6)THEN
              SADEP = -(-   U(I-1,J)-2.0*U(I,J)  +3.0*U(I+1,J))*AMAX1(+FE,0.0)
        IF(I==NIM)THEN
              SADEN = 0.00
        ELSE
              SADEN = +(3.0*U(I,J)  -2.0*U(I+1,J)-    U(I+2,J))*AMAX1(-FE,0.0)
        ENDIF
        IF(I==2)THEN
              SADWP = 0.00
        ELSE
              SADWP = +(-   U(I-2,J)-2.0*U(I-1,J)+3.0*U(I,J)  )*AMAX1(+FW,0.0)
        ENDIF   	  	
              SADWN = -(3.0*U(I-1,J)-2.0*U(I,J)  -    U(I+1,J))*AMAX1(-FW,0.0)   	  	
              SADNP = -(-   U(I,J-1)-2.0*U(I,J)  +3.0*U(I,J+1))*AMAX1(+FN,0.0)  
        IF(J==NJM)THEN
              SADNN = 0.00
        ELSE
              SADNN = +(3.0*U(I,J)  -2.0*U(I,J+1)-    U(I,J+2))*AMAX1(-FN,0.0) 
        ENDIF 	  	
        IF(J==2)THEN
              SADSP = 0.00
        ELSE
              SADSP = +(-   U(I,J-2)-2.0*U(I,J-1)+3.0*U(I,J)  )*AMAX1(+FS,0.0) 
        ENDIF 	  	
              SADSN = -(3.0*U(I,J-1)-2.0*U(I,J)  -    U(I,J+1))*AMAX1(-FS,0.0)
 	            SAD = (1.0/8.0)*(SADEP + SADEN + SADWP + SADWN + SADNP + SADNN + SADSP + SADSN)
        ENDIF


           SOU(I,J) = DY(J)*(P(I,J)-P(I+1,J)) + SAD 

            SU(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J) + AAAU(I,J)

            BU(I,J) = DY(J)/SU(I,J)


ENDDO
ENDDO

        CALL  SLORXY(U, 1, ERROU, NIM, NJM, AE, AW, AN, AS, SU, SOU, OMIGU) 

RETURN
END




 SUBROUTINE SUBV(ERROV, OMIGV, NSCHEMEUV, Re, Gr)

  PARAMETER(NX=905,NY=905)
     COMMON/GRIDN/NI,NJ
     COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

     COMMON/VARIABLES/U(NX,NY), V(NX,NY), T(NX,NY)!, C(NX,NY), ,Qu(NX,NY)

     COMMON/PRESSURES/P(NX,NY)

     COMMON/OBSTACLESUV/AAAU(NX,NY),AAAV(NX,NY)

     COMMON/PRECOEFFC/BU(NX,NY), BV(NX,NY)

  DIMENSION AE(NX,NY), AW(NX,NY), AN(NX,NY), AS(NX,NY), SV(NX,NY), SOV(NX,NY)

!             COEFV = Pr              ! Form A
              ! COEFV = SQRT(Pr/Ra)     ! Form B
            COEFV = 1.00/Re         ! Form C

                NIM = NI-1;  NJM = NJ-2
	   
DO  I = 2,NIM
DO  J = 2,NJM

                 FE = (DY(J+1)*U(I,J)+DY(J)*U(I,J+1))/2.0	    
                 DE = COEFV*HDY(J)/HDX(I)
            AE(I,J) =    DE*SCHEME(FE/DE,NSCHEMEUV) + AMAX1(-FE,0.0)

  	             FW = (DY(J+1)*U(I-1,J)+DY(J)*U(I-1,J+1))/2.0	
	             DW = COEFV*HDY(J)/HDX(I-1)
            AW(I,J) =    DW*SCHEME(FW/DW,NSCHEMEUV) + AMAX1( FW,0.0)

                 FN = (V(I,J)+V(I,J+1))*DX(I)/2.0
                 DN = COEFV*DX(I)/DY(J+1)
            AN(I,J) =    DN*SCHEME(FN/DN,NSCHEMEUV) + AMAX1(-FN,0.0)

                 FS = (V(I,J)+V(I,J-1))*DX(I)/2.0
	             DS = COEFV*DX(I)/DY(J)
            AS(I,J) =    DS*SCHEME(FS/DS,NSCHEMEUV) + AMAX1( FS,0.0)


                SAD = 0.0
!---CDS(Linear Interpolation)
        IF(NSCHEMEUV==4)THEN
              FECDS =     FE*(1.0-FX(I))*V(I,J) + FE*FX(I)*V(I+1,J)  
              FEUDS =      AMAX1(FE,0.0)*V(I,J) + AMIN1(FE,0.0)*V(I+1,J)
              FWCDS = FW*(1.0-FX(I-1))*V(I-1,J) + FW*FX(I-1)*V(I,J)
	          FWUDS =    AMAX1(FW,0.0)*V(I-1,J) + AMIN1(FW,0.0)*V(I,J)
	          FNCDS =             FN*0.5*V(I,J) + FN*0.5*V(I,J+1)  
	          FNUDS =      AMAX1(FN,0.0)*V(I,J) + AMIN1(FN,0.0)*V(I,J+1)
	          FSCDS =           FS*0.5*V(I,J-1) + FS*0.5*V(I,J)
	          FSUDS =    AMAX1(FS,0.0)*V(I,J-1) + AMIN1(FS,0.0)*V(I,J)
                SAD = -FECDS + FEUDS + FWCDS - FWUDS - FNCDS + FNUDS + FSCDS - FSUDS
        ENDIF

!---Second-order Upwind Scheme
        IF(NSCHEMEUV==5)THEN
              SADEP = -(V(I,J)  -V(I-1,J))*AMAX1(+FE,0.0) 
        IF(I==NIM) THEN
              SADEN = 0.00
        ELSE 
              SADEN = +(V(I+1,J)-V(I+2,J))*AMAX1(-FE,0.0)
        ENDIF  	  	
        IF(I==2) THEN
              SADWP = 0.00
        ELSE
              SADWP = +(V(I-1,J)-V(I-2,J))*AMAX1(+FW,0.0)
        ENDIF
              SADWN = -(V(I,J)  -V(I+1,J))*AMAX1(-FW,0.0)
              SADNP = -(V(I,J)  -V(I,J-1))*AMAX1(+FN,0.0)
        IF(J==NJM) THEN
              SADNN = 0.00
        ELSE 
              SADNN = +(V(I,J+1)-V(I,J+2))*AMAX1(-FN,0.0)
        ENDIF 	  	
        IF(J==2) THEN
              SADSP = 0.00
        ELSE
              SADSP = +(V(I,J-1)-V(I,J-2))*AMAX1(+FS,0.0)
        ENDIF
              SADSN = -(V(I,J)  -V(I,J+1))*AMAX1(-FS,0.0)
	            SAD = (1.0/2.0)*(SADEP + SADEN + SADWP + SADWN + SADNP + SADNN + SADSP + SADSN)
        ENDIF

!---Deferrd correction QUICK Scheme
        IF(NSCHEMEUV==6)THEN
              SADEP = -(-   V(I-1,J)-2.0*V(I,J)  +3.0*V(I+1,J))*AMAX1(+FE,0.0)
        IF(I==NIM)THEN
              SADEN = 0.00
        ELSE
              SADEN = +(3.0*V(I,J)  -2.0*V(I+1,J)-    V(I+2,J))*AMAX1(-FE,0.0)
        ENDIF
        IF(I==2)THEN
              SADWP = 0.00
        ELSE
              SADWP = +(-   V(I-2,J)-2.0*V(I-1,J)+3.0*V(I,J)  )*AMAX1(+FW,0.0)
        ENDIF   	  	
              SADWN = -(3.0*V(I-1,J)-2.0*V(I,J)  -    V(I+1,J))*AMAX1(-FW,0.0)   	  	
              SADNP = -(-   V(I,J-1)-2.0*V(I,J)  +3.0*V(I,J+1))*AMAX1(+FN,0.0)  
        IF(J==NJM)THEN
              SADNN = 0.00
        ELSE
              SADNN = +(3.0*V(I,J)  -2.0*V(I,J+1)-    V(I,J+2))*AMAX1(-FN,0.0) 
        ENDIF 	  	
        IF(J==2)THEN
              SADSP = 0.00
        ELSE
              SADSP = +(-   V(I,J-2)-2.0*V(I,J-1)+3.0*V(I,J)  )*AMAX1(+FS,0.0) 
        ENDIF 	  	
              SADSN = -(3.0*V(I,J-1)-2.0*V(I,J)  -    V(I,J+1))*AMAX1(-FS,0.0)
	            SAD = (1.0/8.0)*(SADEP + SADEN + SADWP + SADWN + SADNP + SADNN + SADSP + SADSN)
        ENDIF


           SOV(I,J) = DX(I)*(P(I,J)-P(I,J+1)) + SAD                                              &  
  	   		        ! +          DX(I)*((T(I,J)+0.00     )*DY(J+1)+(T(I,J+1)+0.00       )*DY(J))/2.0       ! & ! For the form B
!   			        +                DX(I)*((T(I,J)+Br*C(I,J))*DY(J+1)+(T(I,J+1)+Br*C(I,J+1))*DY(J))/2.0       ! & ! For the form B
!                   + (Gr/(Re**2.0))*DX(I)*((T(I,J)+Br*C(I,J))*DY(J+1)+(T(I,J+1)+Br*C(I,J+1))*DY(J))/2.0 ! & ! For the form C
                  + (Gr/(Re**2.0))*DX(I)*((T(I,J)+0.00     )*DY(J+1)+(T(I,J+1)+0.00       )*DY(J))/2.0 ! & ! For the form C
!                   +             Ar*DX(I)*((T(I,J)+0.00     )*DY(J+1)+(T(I,J+1)+0.00       )*DY(J))/2.0 ! & ! For the form C

            SV(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J) + AAAV(NX,NY)

            BV(I,J) = DX(I)/SV(I,J)

ENDDO
ENDDO

       
         CALL SLORXY(V, 1, ERROV, NIM, NJM, AE, AW, AN, AS, SV, SOV, OMIGV) 

RETURN
END





SUBROUTINE SUBPC(ERROPC,OMIGPC)

 PARAMETER(NX=905,NY=905)
    COMMON/GRIDN/NI,NJ
    COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

    COMMON/VARIABLES/U(NX,NY), V(NX,NY), T(NX,NY)!, C(NX,NY), Qu(NX,NY)

    COMMON/PRESSURES/P(NX,NY)

    COMMON/PRECOEFFC/BU(NX,NY), BV(NX,NY)

 DIMENSION  AE(NX,NY), AW(NX,NY), AN(NX,NY), AS(NX,NY), SP(NX,NY), SOP(NX,NY), PC(NX,NY)

                 PC = 0.0

                NIM = NI-1;  NJM = NJ-1

               DO I = 2,NIM
               DO J = 2,NJM

      IF(I.NE.NIM) AE(I,J) = BU(I,J)*DY(J)                    
      IF(I.NE.2)   AW(I,J) = BU(I-1,J)*DY(J)
      IF(J.NE.NJM) AN(I,J) = BV(I,J)*DX(I)
      IF(J.NE.2)   AS(I,J) = BV(I,J-1)*DX(I)
        
                  SOP(I,J) = DY(J)*(U(I-1,J)-U(I,J))+DX(I)*(V(I,J-1)-V(I,J))

     	           SP(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J)

               ENDDO
               ENDDO

!              CALL SLORXY(PC,  1, ERROPC, NIM, NJM, AE, AW, AN, AS, SP, SOP, OMIGPC)
!              CALL SLORXY(PC,  2, ERROPC, NIM, NJM, AE, AW, AN, AS, SP, SOP, OMIGPC)
!              CALL SLORXY(PC,  4, ERROPC, NIM, NJM, AE, AW, AN, AS, SP, SOP, OMIGPC)
!              CALL SLORXY(PC, 16, ERROPC, NIM, NJM, AE, AW, AN, AS, SP, SOP, OMIGPC)
               CALL SLORXY(PC, 64, ERROPC, NIM, NJM, AE, AW, AN, AS, SP, SOP, OMIGPC)
!              CALL SLORXY(PC,128, ERROPC, NIM, NJM, AE, AW, AN, AS, SP, SOP, OMIGPC)

               DO I = 2,NIM
               DO J = 2,NJM
               IF(I.NE.NIM) U(I,J) = U(I,J) + BU(I,J)*(PC(I,J)-PC(I+1,J))  
               IF(J.NE.NJM) V(I,J) = V(I,J) + BV(I,J)*(PC(I,J)-PC(I,J+1))
                            P(I,J) = P(I,J) + PC(I,J)
               ENDDO
               ENDDO


RETURN
END



SUBROUTINE SUBT(ERROT,OMIGT,NSCHEMET,Re, Pr)


 PARAMETER(NX=905,NY=905)
     COMMON/GRIDN/NI,NJ
     COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

     COMMON/VARIABLES/U(NX,NY), V(NX,NY), T(NX,NY)!, C(NX,NY) ,Qu(NX,NY)

!     COMMON/OBSTACLESPTCA/INWA,INEA,JNSA,JNNA
!     COMMON/OBSTACLESPTCB/INWB,INEB,JNSB,JNNB
!     COMMON/OBSTACLESPTCC/INWC,INEC,JNSC,JNNC
!     COMMON/OBSTACLESPTCD/INWD,INED,JNSD,JNND

!     COMMON/COEFENERGY/CoeffT(NX,NY)

  DIMENSION AE(NX,NY), AW(NX,NY), AN(NX,NY), AS(NX,NY), ST(NX,NY), SOT(NX,NY)

                NIM = NI-1;  NJM = NJ-1

!            CoeffT = 1.00                ! Form A
!             CoeffT = 1.0/SQRT(Pr*Ra)     ! Form B
             CoeffT = 1.0/(Re*Pr)         ! Form C


! Conjugate heat convection and diffusion  [OPEN Solid]
! Heat conduction in the blocks
!               DO I = INWA,INEA
!               DO J = JNSA,JNNA
!        ROCPr(I,J) = 1.000                ! (pCp)s/(pCp)f  It has been implictly set as unity. For steady state, it has no effects.
!        CoeffT(I,J) = DRT/SQRT(Pr*Ra) 	   ! DRT = 竹(Solid)/竹(Fluid)   !Qinter/Qrefer = 1.00
!               ENDDO
!               ENDDO
!               DO I = INWB,INEB
!               DO J = JNSB,JNNB
!        ROCPr(I,J) = 1.000                ! (pCp)s/(pCp)f  It has been implictly set as unity. For steady state, it has no effects.
!        CoeffT(I,J) = DRT/SQRT(Pr*Ra)      ! DRT = 竹(Solid)/竹(Fluid)   !Qinter/Qrefer = 1.00
!               ENDDO
!               ENDDO
!               DO I = INWC,INEC
!               DO J = JNSC,JNNC
!        ROCPr(I,J) = 1.000                ! (pCp)s/(pCp)f  It has been implictly set as unity. For steady state, it has no effects.
!        CoeffT(I,J) = DRT/SQRT(Pr*Ra)      ! DRT = 竹(Solid)/竹(Fluid)   !Qinter/Qrefer = 1.00
!               ENDDO
!               ENDDO
!               DO I = INWD,INED
!               DO J = JNSD,JNND
!        ROCPr(I,J) = 1.000                ! (pCp)s/(pCp)f  It has been implictly set as unity. For steady state, it has no effects.
!        CoeffT(I,J) = DRT/SQRT(Pr*Ra)      ! DRT = 竹(Solid)/竹(Fluid)   !Qinter/Qrefer = 1.00
!               ENDDO
!               ENDDO


DO I = 2,NIM
DO J = 2,NJM

!              effee = (FX(I)/CoeffT(I,J)+(1.0-FX(I))/CoeffT(I+1,J))**(-1)     ! Harmonic Mean of P and E
!                 De =   effee*DY(J)/HDX(I)
                 De =  CoeffT*DY(J)/HDX(I)
                 Fe =  U(I,J)*DY(J)
            AE(I,J) =      De*SCHEME(Fe/De,NSCHEMET) + AMAX1(-Fe,0.0)

!              effww = (FX(I-1)/CoeffT(I-1,J)+(1.0-FX(I-1))/CoeffT(I,J))**(-1) ! Harmonic Mean of W and P
!                 Dw =   effww*DY(J)/HDX(I-1)
                 Dw =  CoeffT*DY(J)/HDX(I-1)
                 Fw =U(I-1,J)*DY(J)
            AW(I,J) =      Dw*SCHEME(Fw/Dw,NSCHEMET) + AMAX1( Fw,0.0)

!              effnn = (FY(J)/CoeffT(I,J)+(1.0-FY(J))/CoeffT(I,J+1))**(-1)     ! Harmonic Mean of P and N
!                 Dn =   effnn*DX(I)/HDY(J)
                 Dn =  CoeffT*DX(I)/HDY(J)
                 Fn =  V(I,J)*DX(I)
            AN(I,J) =      Dn*SCHEME(Fn/Dn,NSCHEMET) + AMAX1(-Fn,0.0)

!              effss = (FY(J-1)/CoeffT(I,J-1)+(1.0-FY(J-1))/CoeffT(I,J))**(-1) ! Harmonic Mean of S and P
!                 Ds =   effss*DX(I)/HDY(J-1)
                 Ds =  CoeffT*DX(I)/HDY(J-1)
                 Fs =V(I,J-1)*DX(I)
            AS(I,J) =      Ds*SCHEME(Fs/Ds,NSCHEMET) + AMAX1( Fs,0.0)


                SAD = 0.0
!---CDS(Linear Interpolation)
         IF(NSCHEMET==4) THEN
              FeCDS =   Fe*(1.0-FX(I))*T(I,J)  +     Fe*FX(I)*T(I+1,J)  
              FeUDS =    AMAX1(Fe,0.0)*T(I,J)  +AMIN1(Fe,0.0)*T(I+1,J)

              FwCDS = Fw*(1.0-FX(I-1))*T(I-1,J)+   Fw*FX(I-1)*T(I,J)
   	          FwUDS =    AMAX1(Fw,0.0)*T(I-1,J)+AMIN1(Fw,0.0)*T(I,J)

              FnCDS =   Fn*(1.0-FY(J))*T(I,J)  +     Fn*FY(J)*T(I,J+1)  
              FnUDS =    AMAX1(Fn,0.0)*T(I,J)  +AMIN1(Fn,0.0)*T(I,J+1)

              FsCDS = Fs*(1.0-FY(J-1))*T(I,J-1)+   Fs*FY(J-1)*T(I,J)
              FsUDS =    AMAX1(Fs,0.0)*T(I,J-1)+AMIN1(Fs,0.0)*T(I,J)
                SAD = -FeCDS + FeUDS + FwCDS - FwUDS - FnCDS + FnUDS + FsCDS - FsUDS
         ENDIF

!---Second-order Upwind Scheme
         IF(NSCHEMET==5)THEN
              SADEP = -(T(I,J)  -T(I-1,J))*AMAX1(+Fe,0.0) 
         IF(I==NIM) THEN
              SADEN = 0.00
         ELSE 
              SADEN = +(T(I+1,J)-T(I+2,J))*AMAX1(-Fe,0.0)
         ENDIF  	  	
         IF(I==2) THEN
              SADWP = 0.00
         ELSE
              SADWP = +(T(I-1,J)-T(I-2,J))*AMAX1(+Fw,0.0)
         ENDIF
              SADWN = -(T(I,J)  -T(I+1,J))*AMAX1(-Fw,0.0)
              SADNP = -(T(I,J)  -T(I,J-1))*AMAX1(+Fn,0.0)
         IF(J==NJM) THEN
              SADNN = 0.00
         ELSE 
              SADNN = +(T(I,J+1)-T(I,J+2))*AMAX1(-Fn,0.0)
         ENDIF 	  	
         IF(J==2) THEN
              SADSP = 0.00
         ELSE
              SADSP = +(T(I,J-1)-T(I,J-2))*AMAX1(+Fs,0.0)
         ENDIF
              SADSN = -(T(I,J)  -T(I,J+1))*AMAX1(-Fs,0.0)
	            SAD = (1.0/2.0)*(SADEP + SADEN + SADWP + SADWN + SADNP + SADNN + SADSP + SADSN)
         ENDIF


           SOT(I,J) = 0.0 + SAD	                    
   
            ST(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J)

ENDDO
ENDDO

         CALL SLORXY(T, 1, ERROT, NIM, NJM, AE, AW, AN, AS, ST, SOT, OMIGT) 

RETURN
END

FUNCTION SCHEME(PE, ISCHEME)

    IF(ISCHEME==1)THEN
       SCHEME = AMAX1(0.0, 1.0 - ABS(PE)*0.5)                   ! Hybrid Scheme
ELSEIF(ISCHEME==2)THEN
	   SCHEME = AMAX1(0.0, 1.0 - ABS(PE)*0.1)				    ! Power-law Scheme
	   SCHEME = SCHEME**5.0								   
ELSEIF(ISCHEME==3.OR.ISCHEME==4.OR.ISCHEME==5.OR.ISCHEME==6)THEN
	   SCHEME = 1.0										        ! Upwind Scheme
    ENDIF

RETURN
END

  !  FUNCTION HEATFLUX(ARGU)
  !  PARAMETER(NX=905,NY=905)
		!COMMON/HEAT_FLUX_PROFILE/ HFP1(NY), HFP2(NY), HFP3(NY), HFP4(NY), QF(NY), Qu1(NY), Qu2(NY)
  !      INTEGER ARGU
		!REAL QF(NY)
		!QF = 0.0
		!SELECT CASE(ARGU)
		!  CASE(1)
		!    QF = HFP1
		!  CASE(2)
		!    QF = HFP2
		!  CASE(3)
		!    QF = HFP3
		!  CASE(4)
		!    QF = HFP4
		!  CASE DEFAULT
		!    WRITE(*,*) "Bad Input"
  !      END SELECT
		!RETURN
  !  END

SUBROUTINE SLORXY(VARI,NITER,ERROR,NIM,NJM,AE,AW,AN,AS,CAS,SO,OMGA)

 PARAMETER(NX=905,NY=905)

 DIMENSION VARI(NX,NY)
 DIMENSION BAT(NX),DAT(NX)
 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),CAS(NX,NY),SO(NX,NY)

DO L=1, NITER

  DO I=2,NIM
    DO J=2,NJM
        SOUCE =SO(I,J)+AE(I,J)*VARI(I+1,J)+AW(I,J)*VARI(I-1,J)
       IF(J == 2)THEN 
        BAT(J)=AN(I,J)/CAS(I,J)
        DAT(J)=(AS(I,J)*VARI(I,J-1)+SOUCE)/CAS(I,J)
       ELSE
        BAT(J)=AN(I,J)/(CAS(I,J)-AS(I,J)*BAT(J-1))
        DAT(J)=(AS(I,J)*DAT(J-1)+SOUCE)/(CAS(I,J)-AS(I,J)*BAT(J-1))
       ENDIF
    ENDDO
    DO J=NJM,2,-1
        VARI(I,J)=VARI(I,J)+OMGA*(BAT(J)*VARI(I,J+1)+DAT(J)-VARI(I,J))
    ENDDO
  ENDDO

  DO J=2,NJM
    DO I=2,NIM
      SOUCE=SO(I,J)+AN(I,J)*VARI(I,J+1)+AS(I,J)*VARI(I,J-1)
     IF(I == 2) THEN
      BAT(I)=AE(I,J)/CAS(I,J)													 
      DAT(I)=(AW(I,J)*VARI(I-1,J)+SOUCE)/CAS(I,J)
     ELSE
      BAT(I)=AE(I,J)/(CAS(I,J)-AW(I,J)*BAT(I-1))
      DAT(I)=(AW(I,J)*DAT(I-1)+SOUCE)/(CAS(I,J)-AW(I,J)*BAT(I-1))
     ENDIF
    ENDDO
    DO I=NIM,2,-1
        VARI(I,J)=VARI(I,J)+OMGA*(BAT(I)*VARI(I+1,J)+DAT(I)-VARI(I,J))
    ENDDO
  ENDDO

  ERROR = 0.00
  DO J=2,NJM
  DO I=2,NIM
    VARIMAX=SO(I,J)+AN(I,J)*VARI(I,J+1)+AS(I,J)*VARI(I,J-1)            &
         +AE(I,J)*VARI(I+1,J)+AW(I,J)*VARI(I-1,J)-VARI(I,J)*CAS(I,J)
    VARIMAX=ABS(VARIMAX)
    IF(VARIMAX.GE.ERROR)   ERROR=VARIMAX
  ENDDO
  ENDDO

ENDDO

RETURN
END





SUBROUTINE VISULIZATION(Re, Gr, Pr)

 PARAMETER(NX=905,NY=905)
    COMMON/GRIDN/NI,NJ
    COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

    COMMON/VARIABLES/ U(NX,NY), V(NX,NY), T(NX,NY)!C(NX,NY),Qu(NX,NY)

    COMMON/PRESSURES/ P(NX,NY)

 DIMENSION U1(NX,NY), V1(NX,NY),PSI(NX,NY),HF(NX,NY),HFX(NX,NY),HFY(NX,NY)
!    DIMENSION Qexact(NY)!, TM(NX,NY)
                NIM = NI-1;  NJM = NJ-1

!---STAGGERED GRIDS ARE TRANSFORMED INTO NON-STGGERED ONES!
                 U1 = U;  V1 = V

               DO J = 1,NJ
            U1(1,J) = U(1,J)
           U1(NI,J) = U(NIM,J)
               ENDDO
               DO I = 2,NIM
            U1(I,1) = 0.5*(U(I,1)+U(I-1,1))
           U1(I,NJ) = 0.5*(U(I,NJ)+U(I-1,NJ))
               ENDDO

               DO I = 1,NI
            V1(I,1) = V(I,1)
           V1(I,NJ) = V(I,NJM)
               ENDDO
               DO J = 2,NJM
            V1(1,J) = 0.5*(V(1,J)+V(1,J-1))
           V1(NI,J) = 0.5*(V(NI,J)+V(NI,J-1))
               ENDDO

               DO I = 2,NIM
               DO J = 2,NJM
            U1(I,J) = 0.5*(U(I,J)+U(I-1,J))
            V1(I,J) = 0.5*(V(I,J)+V(I,J-1))
               ENDDO
               ENDDO
!---Solving the streamfunction
               CALL  Streamfunction(U1, V1, PSI, PSImax, PSImin)
			   CALL Heatfunction(U1, V1, T, HF, HFX, HFY, HFMAX, HFMIN, Re, Pr) 


!                OPEN(9,FILE=FILE2,STATUS='UNKNOWN')
               WRITE(9,901) 
901	          FORMAT(1X,'TITLE="Field.DAT"')
               WRITE(9,902) 
902	          FORMAT(1X,'VARIABLES="X","Y","T","U","V","P","PSI","Hf","Hfx","Hfy"')
               WRITE(9,903) NI,NJ
903	          FORMAT(1X,'ZONE I=',I3,3X,'J=',I3,3X,'F=POINT')
               DO J = 1,NJ
               DO I = 1,NI
               WRITE(9,904) XC(I), YC(J), T(I,J), U1(I,J), V1(I,J), P(I,J),PSI(I,J)	,Hf(I,J),Hfx(I,J),Hfy(I,J)
               ENDDO
               ENDDO
904           FORMAT(1X, 2(F9.6,1X), 1X, 8(E15.6,1X))
               CLOSE(9)


RETURN
END




SUBROUTINE Streamfunction(U1, V1, PSI, PSImax, PSImin) 

 PARAMETER(NX=905,NY=905)
    COMMON/GRIDN/NI,NJ
    COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

 DIMENSION   U1(NX,NY),V1(NX,NY),PSI(NX,NY)
 DIMENSION   AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),AP(NX,NY),SO(NX,NY)

         NIM = NI-1;  NJM = NJ-1

!---Initiation for Streamfunction!
         PSI = 0.0

!---BOUNDARY OF PSI Values! Integralation Start from Orientation  
    PSI(1,1) = 0.0	 

           J = 1
    PSI(2,J) = PSI(1,J)  -(V1(1,J)+V1(2,J))*HDX(1)/2.0
        DO I = 3,NIM
    PSI(I,J) = PSI(I-1,J)-(V1(I,J)*DX(I-1)+V1(I-1,J)*DX(I))/2.0
        ENDDO
    PSI(NI,J)= PSI(NIM,J)-(V1(NIM,J)+V1(NI,J))*HDX(NIM)/2.0

           I = 1
    PSI(I,2) = PSI(I,1)  +(U1(I,1)+U1(I,2))*HDY(1)/2.0
        DO J = 3,NJM
    PSI(I,J) = PSI(I,J-1)+(U1(I,J)*DY(J-1)+U1(I,J-1)*DY(J))/2.0
        ENDDO 
    PSI(I,NJ)= PSI(I,NJM)+(U1(I,NJM)+U1(I,NJ))*HDY(NJM)/2.0

           I = NI
    PSI(I,2) = PSI(I,1)  +(U1(I,1)+U1(I,2))*HDY(1)/2.0
        DO J = 3,NJM
    PSI(I,J) = PSI(I,J-1)+(U1(I,J)*DY(J-1)+U1(I,J-1)*DY(J))/2.0
        ENDDO
    PSI(I,NJ)= PSI(I,NJM)+(U1(I,NJM)+U1(I,NJ))*HDY(NJM)/2.0
  
           J = NJ
    PSI(2,J) = PSI(1,J)  -(V1(1,J)+V1(2,J))*HDX(1)/2.0
        DO I = 3,NIM
    PSI(I,J) = PSI(I-1,J)-(V1(I,J)*DX(I-1)+V1(I-1,J)*DX(I))/2.0
        ENDDO
    PSI(NI,J)= PSI(NIM,J)-(V1(NIM,J)+V1(NI,J))*HDX(NIM)/2.0

!..Poisson Method for solving the	streamfunction
        DO I = 2,NIM
        DO J = 2,NJM
     AE(I,J) = DY(J)/HDX(I)
     AN(I,J) = DX(I)/HDY(J)
     AW(I,J) = DY(J)/HDX(I-1)
     AS(I,J) = DX(I)/HDY(J-1)
     AP(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)
          Ve = V1(I,J)*(1-FX(I))     + V1(I+1,J)*FX(I)	
          Vw = V1(I-1,J)*(1-FX(I-1)) + V1(I,J)*FX(I-1)	
          Un = U1(I,J)*(1-FY(J))     + U1(I,J+1)*FY(J)	
          Us = U1(I,J-1)*(1-FY(J-1)) + U1(I,J)*FY(J-1)	
     SO(I,J) = (Ve-Vw)*DY(J)-(Un-Us)*DX(I)  
        ENDDO
        ENDDO

             CALL SLORXY(PSI, 80000, ERROPSI, NIM, NJM, AE, AW, AN, AS, AP, SO, 0.50)  

!---Set orientational point
!		DO I = 1, NI
!	    DO J = 1, NJ
!   PSI(I,J) = PSI(I,J) - (PSI(41,1)+PSI(42,1))/2.0  ! Depending on your stipulation for origin of visualization function
!	    ENDDO
!		ENDDO

!------Determining the minimum and maximum values of streamfunction in interior flow field
      PSIMAX = -1.E30
      PSIMIN = +1.E30
        DO I = 1,NI 
        DO J = 1,NJ 
      PSIMAX = MAX(PSIMAX,PSI(I,J))
      PSIMIN = MIN(PSIMIN,PSI(I,J))
        ENDDO
        ENDDO
 RETURN
END

SUBROUTINE Heatfunction(U1, V1, T, HF, HFX, HFY, HFMAX, HFMIN, Re, Pr) 

  PARAMETER(NX=905,NY=905)
     COMMON/GRIDN/NI,NJ
     COMMON/GRID/ X(NX),Y(NY),XC(NX),YC(NY),DX(NX),DY(NY),FX(NX),FY(NY),HDX(NX),HDY(NY)

  DIMENSION  AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),AP(NX,NY),SO(NX,NY)
  DIMENSION  U1(NX,NY),V1(NX,NY),T(NX,NY),HF(NX,NY),HFX(NX,NY),HFY(NX,NY)
  DIMENSION  COEFT(NX,NY), COEFHF(NX,NY)  

         NIM = NI-1; NJM = NJ-1
       ! COEFT = 1.00/SQRT(Pr*Ra)
      ! COEFHF =      SQRT(Pr*Ra)    ! Function coefficient is inverse to the corresponding variable coefficient
	  
	    COEFT = 1.00/(Pr*Re)
        COEFHF =     Pr*Re

!---Temperature Interpolation on four corners
     T(1,1)  = T(2,1)    + T(1,2)    - T(2,2)
     T(NI,1) = T(NI,2)   + T(NIM,1)  - T(NIM,2)
     T(1,NJ) = T(1,NJM)  + T(2,NJ)   - T(2,NJM)
    T(NI,NJ) = T(NI,NJM) + T(NIM,NJ) - T(NIM,NJM)

!---Initiation for Heatfunction!
          HF = 0.0;  HFX = 0.0;  HFY = 0.0 


    DO I = 2,NIM
    DO J = 2,NJM
	    Te   =    T(I,J)*(1-FX(I))   +  T(I+1,J)*FX(I)	
	    Ve   =   V1(I,J)*(1-FX(I))   + V1(I+1,J)*FX(I)	

	    Tw   =  T(I-1,J)*(1-FX(I-1)) +  T(I,J)*FX(I-1)	
	    Vw   = V1(I-1,J)*(1-FX(I-1)) + V1(I,J)*FX(I-1)	

	    Tn   =    T(I,J)*(1-FY(J))   +  T(I,J+1)*FY(J)	
	    Un   =   U1(I,J)*(1-FY(J))   + U1(I,J+1)*FY(J)	

	    Ts   =  T(I,J-1)*(1-FY(J-1)) +  T(I,J)*FY(J-1)	
	    Us   = U1(I,J-1)*(1-FY(J-1)) + U1(I,J)*FY(J-1)	

    HFX(I,J) = U1(I,J)*T(I,J)-COEFT(I,J)*(Te-Tw)/DX(I)
    HFY(I,J) = V1(I,J)*T(I,J)-COEFT(I,J)*(Tn-Ts)/DY(J)

        effe = (FX(I)/COEFHF(I,J)+(1.0-FX(I))/COEFHF(I+1,J))**(-1)     ! Harmonic Mean of P and E
        effw = (FX(I-1)/COEFHF(I-1,J)+(1.0-FX(I-1))/COEFHF(I,J))**(-1) ! Harmonic Mean of W and P
        effn = (FY(J)/COEFHF(I,J)+(1.0-FY(J))/COEFHF(I,J+1))**(-1)     ! Harmonic Mean of P and N
        effs = (FY(J-1)/COEFHF(I,J-1)+(1.0-FY(J-1))/COEFHF(I,J))**(-1) ! Harmonic Mean of S and P

     AE(I,J) = effe*DY(J)/HDX(I)
     AN(I,J) = effn*DX(I)/HDY(J)
     AW(I,J) = effw*DY(J)/HDX(I-1)
     AS(I,J) = effs*DX(I)/HDY(J-1)
     AP(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)

     SO(I,J) = (effe*Te*Ve-effw*Tw*Vw)*DY(J)-(effn*Tn*Un-effs*Ts*Us)*DX(I)

    ENDDO
    ENDDO



           I = 1
        DO J = 1,NJ 
    HFX(I,J) = U1(I,J)*T(I,J) - COEFT(I,J)*(T(I+1,J) -T(I,J))/HDX(1)
        ENDDO

           I = NI
        DO J = 1,NJ 
    HFX(I,J) = U1(I,J)*T(I,J) - COEFT(I,J)*(T(I,J) -T(I-1,J))/HDX(NIM)
        ENDDO

           J = 1
        DO I = 1,NI 
    HFY(I,J) = V1(I,J)*T(I,J) - COEFT(I,J)*(T(I,J+1) -T(I,J))/HDY(1)
        ENDDO

           J = NJ
        DO I = 1,NI 
    HFY(I,J) = V1(I,J)*T(I,J) - COEFT(I,J)*(T(I,J) -T(I,J-1))/HDY(NJM)
        ENDDO

!---Boundary conditons of Heatfunctions  
     HF(1,1) = 0.0   !Integralation Start from the orientation

           J = 1
     HF(2,J) = HF(1,J)  -(HFY(1,J)+HFY(2,J))*HDX(1)/2.0
        DO I = 3, NIM
     HF(I,J) = HF(I-1,J)-(HFY(I,J)*DX(I-1)+HFY(I-1,J)*DX(I))/2.0
        ENDDO
    HF(NI,J) = HF(NIM,J)-(HFY(NIM,J)+HFY(NI,J))*HDX(NIM)/2.0

           I = 1
     HF(I,2) = HF(I,1)  +(HFX(I,1)+HFX(I,2))*HDY(1)/2.0
        DO J = 3, NJM
     HF(I,J) = HF(I,J-1)+(HFX(I,J)*DY(J-1)+HFX(I,J-1)*DY(J))/2.0
        ENDDO 
    HF(I,NJ) = HF(I,NJM)+(HFX(I,NJM)+HFX(I,NJ))*HDY(NJM)/2.0

           I = NI
     HF(I,2) = HF(I,1)  +(HFX(I,1)+HFX(I,2))*HDY(1)/2.0
        DO J = 3, NJM
     HF(I,J) = HF(I,J-1)+(HFX(I,J)*DY(J-1)+HFX(I,J-1)*DY(J))/2.0
        ENDDO
    HF(I,NJ) = HF(I,NJM)+(HFX(I,NJM)+HFX(I,NJ))*HDY(NJM)/2.0

           J = NJ
     HF(2,J) = HF(1,J)  -(HFY(1,J)+HFY(2,J))*HDX(1)/2.0
        DO I = 3, NIM
     HF(I,J) = HF(I-1,J)-(HFY(I,J)*DX(I-1)+HFY(I-1,J)*DX(I))/2.0
        ENDDO
    HF(NI,J) = HF(NIM,J)-(HFY(NIM,J)+HFY(NI,J))*HDX(NIM)/2.0

!---Poisson Method for solving the	heatfunctions
             CALL SLORXY(HF, 80000, ERROHL, NIM, NJM, AE, AW, AN, AS, AP, SO, 0.50)  	!!!For higher DRT or DRC, iteration times should be required.

!---Ensure that values of heatlines equaling the Nusselt numbers!
          HF = HF*Pr*Re ; HFX = HFX*Pr*Re ; HFY = HFY*Pr*Re 

!---Set orientational point
!		DO I = 1, NI
!		DO J = 1, NJ
!    HF(I,J) = HF(I,J) - (HF(41,1)+HF(42,1))/2.0  ! Depending on your stipulation for origin of visualization function
!       ENDDO
!		ENDDO

!---Determining the Minimum and Maximum Heatfunctions in interior field.
       HFMAX = -1.E30
       HFMIN = +1.E30
        DO I = 1,NI 
        DO J = 1,NJ 
       HFMAX = MAX(HFMAX,HF(I,J))
       HFMIN = MIN(HFMIN,HF(I,J))
        ENDDO
        ENDDO

RETURN
END




