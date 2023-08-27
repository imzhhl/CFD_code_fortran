 SUBROUTINE ARDIRECT(LOOP, Re, Pr, Gr, ISX, ISY, SM, NIYRDA, NIYDB)  ! three heat sources

   INCLUDE 'PARAMETER.inc'	 
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)
    COMMON/ARDCOEFENERGY/CoeffT(NX,NY)
    COMMON/ARDCOEFPOLLUT/CoeffC(NX,NY) 


              NIM = NI-1;  NJM = NJ-1 

!--- Constant and Demo Parameters  
			  Scc = 0.7;       Br = 1.0   !;        SX = 101;	    SY = 51;      SM = 10.0

            OMIGU = 0.35;   OMIGV = 0.35;   OMIGPC = 0.75;   OMIGT = 0.85;   OMIGC = 0.85 
 
             Res1 = 1.0E-10;  Res2 = 1.0E-6

           INLOOP = 40000								


!--- Difference schemes for convective terms           
        NSCHEMEUV = 6     !6 - QUICK(3rd-order)
         NSCHEMET = 5     !5 - SUD(2nd-order)
         NSCHEMEC = 5     !5 - SUD(2nd-order) 


!--- Initializing the direct flow fields

                U = 0.0;     V = 0.0;     P = 0.0;     T = 0.0;     C = 0.0
			    
!--- Inner Iteration!
      DO ITER = 1, INLOOP
		 												   
		 CALL   SUBTCNeuumanDRA(NIYRDA,NIYDB)		   

         CALL   SUBUDRA(ERROU,OMIGU,NSCHEMEUV,Re)  

 	     CALL   SUBVDRA(ERROV,OMIGV,NSCHEMEUV,Re,Gr,Br)  

         CALL   SUBPCDRA(ERROPC,OMIGPC)

         CALL   SUBTDRA(ERROT,OMIGT,NSCHEMET,Re,Pr)          

		 CALL   SUBCDRA(ERROC, OMIGC, NSCHEMEC, Re, Scc, SM, ISX, ISY)

            AMAX0 = AMAX1(ERROU,ERROV,ERROPC,ERROT,ERROC)

!---Check the balance of mass flow rates across the domain 
!        !Global fluid flow flux on the West side  
           FlowWest = 0.00  
               DO J = 2,NJM
           FlowWest = FlowWest + U(1,J)*DY(J)
               ENDDO
!        !Global fluid flow flux on the East side
           FlowEast = 0.00  
               DO J = 2,NJM
           FlowEast = FlowEast + U(NIM,J)*DY(J)
               ENDDO
!!        !Global fluid flow flux on the South side
          FlowSouth = 0.00  
               DO I = 2,NIM
          FlowSouth = FlowSouth + V(I,1)*DX(I)
               ENDDO
!        !Global fluid flow flux on the North side
          FlowNorth = 0.00  
               DO I = 2,NIM
          FlowNorth = FlowNorth + V(I,NJM)*DX(I)
               ENDDO
            ResFlow = FlowEast  -  FlowWest  +  FlowNorth  -  FlowSouth     ! flow flux Je - Jw + Jn - Js

!---Check the balance of thermal flux across the domain
         !Calculate Global heat Flux along the West side  
		 !Please pay attention on the CoeffT(I,J) within the fluid side, while not on the boundary !
            TNuWest = 0.0  
                  I = 1
               DO J = 2,NJM
!              HFXX = SQRT(Ra*Pr)*U(I,J)*T(I,J) - SQRT(Ra*Pr)*CoeffT(I+1,J)*(T(I+1,J)-T(I,J))/HDX(I)	
               HFXX =             U(I,J)*T(I,J) -             CoeffT(I+1,J)*(T(I+1,J)-T(I,J))/HDX(I)	
            TNuWest = TNuWest + HFXX*DY(J)
               ENDDO
         !Calculate Global Heat Flux on the East side
            TNuEast = 0.0  
                  I = NIM
               DO J = 2,NJM
!              HFXX = SQRT(Ra*Pr)*U(I,J)*T(I+1,J) - SQRT(Ra*Pr)*CoeffT(I,J)*(T(I+1,J)-T(I,J))/HDX(I)   
               HFXX =             U(I,J)*T(I+1,J) -             CoeffT(I,J)*(T(I+1,J)-T(I,J))/HDX(I)   
            TNuEast = TNuEast + HFXX*DY(J)
               ENDDO
         !Calculate Global Heat Flux on the South side
           TNuSouth = 0.0
                  J = 1
               DO I = 2,NIM			
!              HFYY = SQRT(Ra*Pr)*V(I,J)*T(I,J) - SQRT(Ra*Pr)*CoeffT(I,J+1)*(T(I,J+1)-T(I,J))/HDY(J)
               HFYY =             V(I,J)*T(I,J) -             CoeffT(I,J+1)*(T(I,J+1)-T(I,J))/HDY(J)
           TNuSouth = TNuSouth + HFYY*DX(I)
               ENDDO
         !Calculate Global Heat Flux on the North side
           TNuNorth = 0.0  
                  J = NJM
               DO I = 2,NIM
!              HFYY = SQRT(Ra*Pr)*V(I,J)*T(I,J+1) - SQRT(Ra*Pr)*CoeffT(I,J)*(T(I,J+1)-T(I,J))/HDY(J)
               HFYY =             V(I,J)*T(I,J+1) -             CoeffT(I,J)*(T(I,J+1)-T(I,J))/HDY(J)
           TNuNorth = TNuNorth + HFYY*DX(I)
               ENDDO
            ResHEAT = TNuEast  -  TNuWest  +  TNuNorth  -  TNuSouth    !thermal flux Je - Jw + Jn - Js

!---Check the balance of species fluxes across the domain
         !Calculate Global Species Flux along the West side  
            TShWest = 0.0  
                  I = 1
               DO J = 2,NJM
!              CFXX = XLe*SQRT(Pr*Ra)*U(I,J)*C(I,J) - XLe*SQRT(Pr*Ra)*CoeffC(I+1,J)*(C(I+1,J)-C(I,J))/HDX(I)
               CFXX =                 U(I,J)*C(I,J) -                 CoeffC(I+1,J)*(C(I+1,J)-C(I,J))/HDX(I)
            TShWest = TShWest + CFXX*DY(J)
               ENDDO
         !Calculate Global Species Flux on the East side
            TShEast = 0.0  
                  I = NIM
               DO J = 2,NJM
!              CFXX = XLe*SQRT(Pr*Ra)*U(I,J)*C(I+1,J) - XLe*SQRT(Pr*Ra)*CoeffC(I,J)*(C(I+1,J)-C(I,J))/HDX(I)   
               CFXX =                 U(I,J)*C(I+1,J) -                 CoeffC(I,J)*(C(I+1,J)-C(I,J))/HDX(I)   
            TShEast = TShEast + CFXX*DY(J)
               ENDDO
         !Calculate Global Species Flux on the South side
           TShSouth = 0.0  
                  J = 1
               DO I = 2,NIM
!              CFYY = XLe*SQRT(Pr*Ra)*V(I,J)*C(I,J) - XLe*SQRT(Pr*Ra)*CoeffC(I,J+1)*(C(I,J+1)-C(I,J))/HDY(J)
               CFYY =                 V(I,J)*C(I,J) -                 CoeffC(I,J+1)*(C(I,J+1)-C(I,J))/HDY(J)
           TShSouth = TShSouth + CFYY*DX(I)
               ENDDO
         !Calculate Global Species Flux on the North side
           TShNorth = 0.0  
                  J = NJM
               DO I = 2,NIM
!              CFYY = XLe*SQRT(Pr*Ra)*V(I,J)*C(I,J+1) - XLe*SQRT(Pr*Ra)*CoeffC(I,J)*(C(I,J+1)-C(I,J))/HDY(J)
               CFYY =                 V(I,J)*C(I,J+1) -                 CoeffC(I,J)*(C(I,J+1)-C(I,J))/HDY(J)
           TShNorth = TShNorth + CFYY*DX(I)
               ENDDO
           ResPollu = TShEast  -  TShWest  +  TShNorth  -  TShSouth - SM*DX(ISX)*DY(ISY)   ! mass flux Je - Jw + Jn - Js - Ms

               IF((AMAX0.LE.Res1).AND.(ABS(ResFLOW).LE.Res2).AND.(ABS(ResHEAT).LE.Res2).AND.(ABS(ResPollu).LE.Res2)) EXIT	

               IF(ITER==1.OR.MOD(ITER,1000)==0) THEN
               WRITE(*,105) ITER, ERROU, ERROV, ERROPC, ERROT, ERROC, ABS(ResFlow), ABS(ResHeat), ABS(ResPollu)
               ENDIF
105           FORMAT(1X,I5,1X,5(E8.3E1,1X),1X,3(E8.3E1,1X))

      ENDDO
!End of inner iteration!

!--- Plot the direct flow fields
CALL VISULIZATIONDRA(LOOP)  

RETURN
END
!--* * * * * * * * * THE END OF THE PROGRAM MAIN     * * * * * * * * * *




 SUBROUTINE SUBTCNeuumanDRA(NIYRDA,NIYDB)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)


          NIM = NI-1;   NJM = NJ-1

!TOP - free port
               DO I = 2,NIYRDA-1
!           U(I,NJ) = 0.00
!          V(I,NJM) = 0.00
            T(I,NJ) = T(I,NJM)      
            C(I,NJ) = C(I,NJM)      
               ENDDO

              DO I = NIYRDA,NIYRDA+5
          V(I,NJM) = V(I,NJM-1) - DY(NJM)*(U(I,NJM)        - U(I-1,NJM)          )/DX(I)
           U(I,NJ) = U(I,NJM)
       IF(V(I,NJM).LT.0)  T(I,NJ)=0.0	 ! ambient air is of lower temperature   
       IF(V(I,NJM).LT.0)  C(I,NJ)=0.0	 ! ambient air is of less polluted   
       IF(V(I,NJM).GE.0)  T(I,NJ)= T(I,NJM) 
       IF(V(I,NJM).GE.0)  C(I,NJ)= C(I,NJM) 
              ENDDO 

               DO I = NIYRDA+6,NIYDB-1
!           U(I,NJ) = 0.00
!          V(I,NJM) = 0.00
            T(I,NJ) = T(I,NJM)      
            C(I,NJ) = C(I,NJM)      
               ENDDO

              DO I = NIYDB,NIYDB+5
          V(I,NJM) = V(I,NJM-1) - DY(NJM)*(U(I,NJM)        - U(I-1,NJM)          )/DX(I)
           U(I,NJ) = U(I,NJM)
       IF(V(I,NJM).LT.0)  T(I,NJ)=0.0	 ! ambient air is of lower temperature   
       IF(V(I,NJM).LT.0)  C(I,NJ)=0.0    ! ambient air is of less polluted   
       IF(V(I,NJM).GE.0)  T(I,NJ)= T(I,NJM) 
       IF(V(I,NJM).GE.0)  C(I,NJ)= C(I,NJM) 
              ENDDO 

               DO I = NIYDB+6,NIM
!           U(I,NJ) = 0.00
!          V(I,NJM) = 0.00
            T(I,NJ) = T(I,NJM)      
            C(I,NJ) = C(I,NJM)      
               ENDDO

!BOTTOM - adiabatic floor
               DO I = 2,45                  
!            U(I,1) = 0.00
!            V(I,1) = 0.00
             T(I,1) = T(I,2)
             C(I,1) = C(I,2)
               ENDDO 
!BOTTOM - adiabatic floor but heated
               DO I = 46,55                  
!            U(I,1) = 0.00
!            V(I,1) = 0.00
             T(I,1) = +1.0      !T(I,2) + 1*HDY(1)      !+1.0
             C(I,1) = C(I,2)          !+0.50	 !High level of soluting source
               ENDDO 
!BOTTOM - adiabatic floor
               DO I = 56,NIM                  
!            U(I,1) = 0.00
!            V(I,1) = 0.00
             T(I,1) = T(I,2)
             C(I,1) = C(I,2)
               ENDDO 

!!LEFT - free port AND adiabatic and impearmeable wall	
	
               DO J = 2,6
             U(1,J) = 1.0      
             V(1,J) = 0.0
             T(1,J) = 0.0 
             C(1,J) = 0.0
               ENDDO 

               DO J = 7,NJM		     
!            U(1,J) = 0.0	                                               
!            V(1,J) = 0.0
			 T(1,J) = T(2,J)
			 C(1,J) = C(2,J)      
               ENDDO 


!RIGHT - adiabatic and impearmeable wall
               DO J = 2,21
!           U(NIM,J)= 0.00   
!           V(NI,J) = 0.00   
            T(NI,J) = T(NIM,J) 
            C(NI,J) = C(NIM,J) 
               ENDDO 

               DO J = 22,31
!           U(NIM,J)= 0.00   
!           V(NI,J) = 0.00   
            T(NI,J) = +1.0     !T(NIM,J) + 1*HDX(NIM)            !+1.0 
            C(NI,J) = C(NIM,J) 
               ENDDO 
			   
			   DO J = 32,NJM
!           U(NIM,J)= 0.00   
!           V(NI,J) = 0.00   
            T(NI,J) = T(NIM,J) 
            C(NI,J) = C(NIM,J) 
               ENDDO
RETURN
END





SUBROUTINE SUBUDRA(ERROU,OMIGU,NSCHEMEUV, Re) 

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)
    COMMON/ARPREUV/ BU(NX,NY),BV(NX,NY)

 
 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SU(NX,NY),SOU(NX,NY)


!             COEFU = Pr             ! Form A
!             COEFU = SQRT(Pr/Ra)    ! Form B
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

            SU(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J) 

            BU(I,J) = DY(J)/SU(I,J)


ENDDO
ENDDO

        CALL  SLORXY(U, 1, ERROU, NIM, NJM, AE, AW, AN, AS, SU, SOU, OMIGU) 

RETURN
END





SUBROUTINE SUBVDRA(ERROV,OMIGV,NSCHEMEUV,Re,Gr,Br)     !,JNNA,JNNB)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)
    COMMON/ARPREUV/ BU(NX,NY),BV(NX,NY)


 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SV(NX,NY),SOV(NX,NY)

!             COEFV = Pr              ! Form A
!             COEFV = SQRT(Pr/Ra)     ! Form B
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
!  	   		        +          DX(I)*((T(I,J)+0.00     )*DY(J+1)+(T(I,J+1)+0.00       )*DY(J))/2.0       ! & ! For the form B
!  			        +          DX(I)*((T(I,J)+Br*C(I,J))*DY(J+1)+(T(I,J+1)+Br*C(I,J+1))*DY(J))/2.0       ! & ! For the form B
                    + (Gr/(Re**2.0))*DX(I)*((T(I,J)+Br*C(I,J))*DY(J+1)+(T(I,J+1)+Br*C(I,J+1))*DY(J))/2.0 ! & ! For the form C
!                   + (Gr/(Re**2.0))*DX(I)*((T(I,J)+0.00     )*DY(J+1)+(T(I,J+1)+0.00       )*DY(J))/2.0 ! & ! For the form C
!                   +             Ar*DX(I)*((T(I,J)+0.00     )*DY(J+1)+(T(I,J+1)+0.00       )*DY(J))/2.0 ! & ! For the form C

            SV(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J) 

            BV(I,J) = DX(I)/SV(I,J)

ENDDO
ENDDO

       
         CALL SLORXY(V, 1, ERROV, NIM, NJM, AE, AW, AN, AS, SV, SOV, OMIGV) 

RETURN
END





SUBROUTINE SUBPCDRA(ERROPC,OMIGPC)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)
    COMMON/ARPREUV/ BU(NX,NY),BV(NX,NY)

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SP(NX,NY),SOP(NX,NY),PC(NX,NY)

       PC = 0.0

      NIM = NI-1;  NJM = NJ-1

     DO I = 2,NIM
     DO J = 2,NJM

     IF(I.NE.NIM) AE(I,J) = BU(I,J)*DY(J)                    
     IF(I.NE.2)   AW(I,J) = BU(I-1,J)*DY(J)
     IF(J.NE.NJM) AN(I,J) = BV(I,J)*DX(I)
     IF(J.NE.2)   AS(I,J) = BV(I,J-1)*DX(I)
        
                 SOP(I,J) = DY(J)*(U(I-1,J)-U(I,J))+DX(I)*(V(I,J-1)-V(I,J))
    	          SP(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)

     ENDDO
     ENDDO

     CALL SLORXY(PC, 4,ERROPC,NIM,NJM,AE,AW,AN,AS,SP,SOP,OMIGPC) 
	
          DO I = 2,NIM
          DO J = 2,NJM
     IF(I.NE.NIM) U(I,J) = U(I,J) + BU(I,J)*(PC(I,J)-PC(I+1,J))  
     IF(J.NE.NJM) V(I,J) = V(I,J) + BV(I,J)*(PC(I,J)-PC(I,J+1))
                  P(I,J) = P(I,J) + PC(I,J)
          ENDDO
          ENDDO

RETURN
END



SUBROUTINE SUBTDRA(ERROT,OMIGT,NSCHEMET,Re,Pr) 

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
	COMMON/ARDCOEFENERGY/CoeffT(NX,NY)

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),ST(NX,NY),SOT(NX,NY)
 
! DIMENSION COEFT(NX,NY),COKR(NX,NY)

         NIM = NI-1;  NJM = NJ-1

!            CoeffT = 1.00                ! Form A
!            CoeffT = 1.0/SQRT(Pr*Ra)     ! Form B
             CoeffT = 1.0/(Re*Pr)         ! Form C


DO I = 2,NIM
DO J = 2,NJM

              effee = (FX(I)/CoeffT(I,J)+(1.0-FX(I))/CoeffT(I+1,J))**(-1)     ! Harmonic Mean of P and E
                 De =   effee*DY(J)/HDX(I)
!                De =  CoeffT*DY(J)/HDX(I)
                 Fe =  U(I,J)*DY(J)
            AE(I,J) =      De*SCHEME(Fe/De,NSCHEMET) + AMAX1(-Fe,0.0)

              effww = (FX(I-1)/CoeffT(I-1,J)+(1.0-FX(I-1))/CoeffT(I,J))**(-1) ! Harmonic Mean of W and P
                 Dw =   effww*DY(J)/HDX(I-1)
!                Dw =  CoeffT*DY(J)/HDX(I-1)
                 Fw =U(I-1,J)*DY(J)
            AW(I,J) =      Dw*SCHEME(Fw/Dw,NSCHEMET) + AMAX1( Fw,0.0)

              effnn = (FY(J)/CoeffT(I,J)+(1.0-FY(J))/CoeffT(I,J+1))**(-1)     ! Harmonic Mean of P and N
                 Dn =   effnn*DX(I)/HDY(J)
!                Dn =  CoeffT*DX(I)/HDY(J)
                 Fn =  V(I,J)*DX(I)
            AN(I,J) =      Dn*SCHEME(Fn/Dn,NSCHEMET) + AMAX1(-Fn,0.0)

              effss = (FY(J-1)/CoeffT(I,J-1)+(1.0-FY(J-1))/CoeffT(I,J))**(-1) ! Harmonic Mean of S and P
                 Ds =   effss*DX(I)/HDY(J-1)
!                Ds =  CoeffT*DX(I)/HDY(J-1)
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


SUBROUTINE SUBCDRA(ERROC,OMIGC,NSCHEMEC,Re,Scc,SM,ISX,ISY) 

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
	COMMON/ARDCOEFPOLLUT/CoeffC(NX,NY)

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SC(NX,NY),SOC(NX,NY)
 

         NIM = NI-1;  NJM = NJ-1

!            CoeffC = 1.00/XLe                 ! Form A
!            CoeffC = 1.00/(XLe*SQRT(Pr*Ra))   ! Form B
             CoeffC = 1.00/(Re*Scc)            ! Form C

DO I = 2,NIM
DO J = 2,NJM

              effee = (FX(I)/CoeffC(I,J)+(1.0-FX(I))/CoeffC(I+1,J))**(-1)     ! Harmonic Mean of P and E
                 De =   effee*DY(J)/HDX(I)
!                De =  CoeffC*DY(J)/HDX(I)
                 Fe =  U(I,J)*DY(J)
            AE(I,J) =      De*SCHEME(Fe/De,NSCHEMEC) + AMAX1(-Fe,0.0)

              effww = (FX(I-1)/CoeffC(I-1,J)+(1.0-FX(I-1))/CoeffC(I,J))**(-1) ! Harmonic Mean of W and P
                 Dw =   effww*DY(J)/HDX(I-1)
!                Dw =  CoeffC*DY(J)/HDX(I-1)
                 Fw =U(I-1,J)*DY(J)
            AW(I,J) =      Dw*SCHEME(Fw/Dw,NSCHEMEC) + AMAX1( Fw,0.0)

              effnn = (FY(J)/CoeffC(I,J)+(1.0-FY(J))/CoeffC(I,J+1))**(-1)  	  ! Harmonic Mean of P and N
                 Dn =   effnn*DX(I)/HDY(J)
!                Dn =  CoeffC*DX(I)/HDY(J)
                 Fn =  V(I,J)*DX(I)
            AN(I,J) =      Dn*SCHEME(Fn/Dn,NSCHEMEC) + AMAX1(-Fn,0.0)

              effss = (FY(J-1)/CoeffC(I,J-1)+(1.0-FY(J-1))/CoeffC(I,J))**(-1) ! Harmonic Mean of S and P
                 Ds =   effss*DX(I)/HDY(J-1)
!                Ds =  CoeffC*DX(I)/HDY(J-1)
                 Fs =V(I,J-1)*DX(I)
            AS(I,J) =      Ds*SCHEME(Fs/Ds,NSCHEMEC) + AMAX1( Fs,0.0)


                SAD = 0.0
!---CDS(Linear Interpolation)
         IF(NSCHEMEC==4) THEN
              FeCDS =   Fe*(1.0-FX(I))*C(I,J)  +     Fe*FX(I)*C(I+1,J)  
              FeUDS =    AMAX1(Fe,0.0)*C(I,J)  +AMIN1(Fe,0.0)*C(I+1,J)

              FwCDS = Fw*(1.0-FX(I-1))*C(I-1,J)+   Fw*FX(I-1)*C(I,J)
	          FwUDS =    AMAX1(Fw,0.0)*C(I-1,J)+AMIN1(Fw,0.0)*C(I,J)

              FnCDS =   Fn*(1.0-FY(J))*C(I,J)  +     Fn*FY(J)*C(I,J+1)  
              FnUDS =    AMAX1(Fn,0.0)*C(I,J)  +AMIN1(Fn,0.0)*C(I,J+1)

              FsCDS = Fs*(1.0-FY(J-1))*C(I,J-1)+   Fs*FY(J-1)*C(I,J)
              FsUDS =    AMAX1(Fs,0.0)*C(I,J-1)+AMIN1(Fs,0.0)*C(I,J)
                SAD = -FeCDS + FeUDS + FwCDS - FwUDS - FnCDS + FnUDS + FsCDS - FsUDS
         ENDIF

!---Second-order Upwind Scheme
         IF(NSCHEMEC==5)THEN
              SADEP = -(C(I,J)  -C(I-1,J))*AMAX1(+Fe,0.0) 
         IF(I==NIM) THEN
              SADEN = 0.00
         ELSE 
              SADEN = +(C(I+1,J)-C(I+2,J))*AMAX1(-Fe,0.0)
         ENDIF  	  	
         IF(I==2) THEN
              SADWP = 0.00
         ELSE
              SADWP = +(C(I-1,J)-C(I-2,J))*AMAX1(+Fw,0.0)
         ENDIF
              SADWN = -(C(I,J)  -C(I+1,J))*AMAX1(-Fw,0.0)
              SADNP = -(C(I,J)  -C(I,J-1))*AMAX1(+Fn,0.0)
         IF(J==NJM) THEN
              SADNN = 0.00
         ELSE 
              SADNN = +(C(I,J+1)-C(I,J+2))*AMAX1(-Fn,0.0)
         ENDIF 	  	
         IF(J==2) THEN
              SADSP = 0.00
         ELSE
              SADSP = +(C(I,J-1)-C(I,J-2))*AMAX1(+Fs,0.0)
         ENDIF
              SADSN = -(C(I,J)  -C(I,J+1))*AMAX1(-Fs,0.0)
	            SAD = (1.0/2.0)*(SADEP + SADEN + SADWP + SADWN + SADNP + SADNN + SADSP + SADSN)
         ENDIF

		 IF(I == ISX.AND.J == ISY) THEN
		   SOC(I,J) = 0.0 + SAD + SM*DX(I)*DY(J)
		 ELSE
		   SOC(I,J) = 0.0 + SAD
		 ENDIF                    
   
			SC(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J)

ENDDO
ENDDO

         CALL SLORXY(C, 1, ERROC, NIM, NJM, AE, AW, AN, AS, SC, SOC, OMIGC) 

RETURN
END




SUBROUTINE VISULIZATIONDRA(LSAVE)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY),C(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)

 DIMENSION U1(NX,NY),V1(NX,NY)

 CHARACTER FIField*6, IC*1, JC*1, KC*1, LC*1

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

!---Pressure interpolation on four corners
        DO J = 2,NJM
      P(1,J) = P(2,J)    
     P(NI,J) = P(NIM,J)
        ENDDO 
        DO I = 2,NIM
      P(I,1) = P(I,2)   
     P(I,NJ) = P(I,NJM)
        ENDDO     
      P(1,1) = P(2,1)   +P(1,2)   -P(2,2)
     P(NI,1) = P(NI,2)  +P(NIM,1) -P(NIM,2)
     P(1,NJ) = P(1,NJM) +P(2,NJ)  -P(2,NJM)
    P(NI,NJ) = P(NI,NJM)+P(NIM,NJ)-P(NIM,NJM)


!---Save the loop numbers as output file name. !The maximum loop is no more than 10000.
         III =  LSAVE/1000
	     JJJ = (LSAVE-1000*III)/100
	     KKK = (LSAVE-1000*III-100*JJJ)/10
	     LLL =  LSAVE-1000*III-100*JJJ-10*KKK
	      IC = CHAR(48+III)
	      JC = CHAR(48+JJJ)
	      KC = CHAR(48+KKK)
	      LC = CHAR(48+LLL)
     FIField = 'ARDF'//IC//JC//KC//LC	   !Plot the direct convection field

!---Dynamic save the flow fields
         OPEN(1, FILE=FIField, STATUS='UNKNOWN')
        WRITE(1,101) 
101	   FORMAT(1X,'TITLE="Field.DAT"')
        WRITE(1,102) 
102	   FORMAT(1X,'VARIABLES="X","Y","T","C","U","V"')
        WRITE(1,103) NI,NJ
103	   FORMAT(1X,'ZONE I=',I3,3X,'J=',I3,3X,'F=POINT')
        DO J = 1,NJ
        DO I = 1,NI
        WRITE(1,104) XC(I),YC(J),T(I,J),C(I,J),U1(I,J),V1(I,J)
        ENDDO
        ENDDO
        CLOSE(1)
104    FORMAT(1X,2(F9.6,1X),1X,4(E15.6,1X))


RETURN
END