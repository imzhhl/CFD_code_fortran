!SUBROUTINE ARDIRECT(LOOP, Ra, Pr, NIYRDA)				       ! one heat source
 SUBROUTINE ARDIRECT(LOOP, Ra, Pr, NIYRDA, NIYDB)		       ! two heat sources
!SUBROUTINE ARDIRECT(LOOP, Ra, Pr, theta, NIYRDA, NIYDB, NIYDC)       ! three heat sources
   INCLUDE 'PARAMETER.inc'	 
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)

!    COMMON/ARDOBSTACLESPTCA/INWAR,INEAR,JNSAR !,JNNA
!    COMMON/ARDOBSTACLESPTCB/INWBR,INEBR,JNSBR !,JNNB

! DIMENSION AU(NX,NY), AV(NX,NY)
              NIM = NI-1;  NJM = NJ-1    !;  IRDNIY = 72 - NIY2

!--- Partitions on the bottom and top side
!			 INWAR = 31;   INEAR = 32;    JNSAR = 2;       JNNAR = NIY2
!			 INWBR = 71;   INEBR = 72;    JNSBR = NIY2+11; JNNBR = NJM
!--- Constant and Demo Parameters  

            OMIGU = 0.35;    OMIGV = 0.35;   OMIGPC = 0.75;    OMIGT = 0.85     !;  AU = 0.0;  AV = 0.0

             Res1 = 1.0E-10;  Res2 = 1.0E-6

           INLOOP = 25000								

!--- Difference schemes for convective terms           
        NSCHEMEUV = 6     !1 - HYBRID;	2 - POWERLAW; 3 - UPWIND(1st-order); 4 - CDS(Linear Interpolation); 5 - SUD(2nd-order); 6 - QUICK(3rd-order)
         NSCHEMET = 6     !1 - HYBRID;	2 - POWERLAW; 3 - UPWIND(1st-order); 4 - CDS(Linear Interpolation); 5 - SUD(2nd-order); 6 - QUICK(3rd-order)
           ISTATE = 1     !1 - Steady;  2 - Transient


!--- Initializing the direct flow fields

                U = 0.0;     V = 0.0;     P = 0.0;     T = 0.0

!--- Inner Iteration!
      DO ITER = 1, INLOOP
		 
! 		 CALL   SUBTCNeuumanDRA(NIYRDA)                         ! one heat source
		 CALL   SUBTCNeuumanDRA(NIYRDA,NIYDB)					! two heat sources
!        CALL   SUBTCNeuumanDRA(NIYRDA,NIYDB,NIYDC)				! three heat sources

         CALL   SUBUDRA(ERROU,OMIGU,ISTATE,NSCHEMEUV,Pr,Ra,theta)     !,JNNAR,JNNBR)

 	     CALL   SUBVDRA(ERROV,OMIGV,ISTATE,NSCHEMEUV,Pr,Ra,theta)     !,JNNAR,JNNBR)

         CALL   SUBPCDRA(ERROPC,OMIGPC)

         CALL   SUBTDRA(ERROT,OMIGT,NSCHEMET,Pr,Ra)              !,JNNAR,JNNBR)


            AMAX0 = AMAX1(ERROU,ERROV,ERROPC,ERROT)

         !Calculate Global Heat Fluxes on three sides to check the thermal balance
              TNuB = 0.0;  TNuT = 0.0;  TNuB = 0.0;  TNuR = 0.0  
          
!---Check the overall heat transfer rates
         !Calculate Global heat Flux along the Left side  
               TNuL = 0.0  
                  I = 1
               DO J = 2,NJM
               HFXX = SQRT(Ra*Pr)*U(I,J)*T(I,J) - (T(I+1,J)-T(I,J))/HDX(I)
			 !  HFXX = U(I,J)*T(I,J) - (T(I+1,J)-T(I,J))/HDX(I)
               TNuL = TNuL + HFXX*DY(J)
               ENDDO
         !Calculate Global Heat Flux on the Right side
               TNuR = 0.0  
                  I = NIM
               DO J = 2,NJM
               HFXX = SQRT(Ra*Pr)*U(I,J)*T(I+1,J) - (T(I+1,J)-T(I,J))/HDX(I)  
			  ! HFXX =U(I,J)*T(I+1,J) - (T(I+1,J)-T(I,J))/HDX(I) 
               TNuR = TNuR + HFXX*DY(J)
               ENDDO
         !Calculate Global Heat Flux on the Bottom side
               TNuB = 0.0
                  J = 1
               DO I = 2,NIM			
               HFYY = SQRT(Ra*Pr)*V(I,J)*T(I,J) - (T(I,J+1)-T(I,J))/HDY(J)
			 ! HFYY = V(I,J)*T(I,J) - (T(I,J+1)-T(I,J))/HDY(J)
               TNuB = TNuB + HFYY*DX(I)
               ENDDO
         !Calculate Global Heat Flux on the Top side
               TNuT = 0.0  
                  J = NJM
               DO I = 2,NIM
               HFYY = SQRT(Ra*Pr)*V(I,J)*T(I,J+1) - (T(I,J+1)-T(I,J))/HDY(J)
			  ! HFYY = V(I,J)*T(I,J+1) - (T(I,J+1)-T(I,J))/HDY(J)
               TNuT = TNuT + HFYY*DX(I)
               ENDDO
            ResHeat = TNuR  -  TNuL  +  TNuT  -  TNuB    !thermal flux Je - Jw + Jn - Js



             IF((AMAX0.LE.Res1).AND.(Resheat.LE.Res2)) EXIT	

             IF(ITER==1.OR.MOD(ITER,1000)==0) THEN
             WRITE(*,110)  ITER, ERROU,ERROV,ERROPC,ERROT, ResHeat
             ENDIF
110         FORMAT(2X,I6,4X,5(E8.4E1,3X))

      ENDDO
!End of inner iteration!
      WRITE(*,*) 'Thermal Balance for RDirect problem', ResHeat

!--- Plot the direct flow fields
CALL VISULIZATIONDRA(LOOP)  

RETURN
END
!--* * * * * * * * * THE END OF THE PROGRAM MAIN     * * * * * * * * * *




!SUBROUTINE SUBTCNeuumanDRA(NIYRDA)										! one heat source				 
 SUBROUTINE SUBTCNeuumanDRA(NIYRDA,NIYDB)								! two heat sources
!SUBROUTINE SUBTCNeuumanDRA(NIYRDA,NIYDB,NIYDC)							! three heat sources

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)


          NIM = NI-1;   NJM = NJ-1
!TOP - the adiabatic wall
                             DO I = 2,NIM
!                         U(I,NJ) = +0.0
!                        V(I,NJM) = +0.0
                          T(I,NJ) = T(I,NJM)
						     ENDDO

!BOTTOM -the adiabatic floor
                             DO I = 2,NIM
!                          U(I,1) = +0.0
!                          V(I,1) = +0.0
                           T(I,1) = T(I,2)
						     ENDDO

!LEFT - the unknown Qu2 wall
							DO J = 2,NIYRDA-1
						  T(1,J) = T(2,J)
					         ENDDO							

							DO J = NIYRDA,NIYRDA+9
						  T(1,J) = T(2,J) + HDX(1)
					         ENDDO							

							DO J = NIYRDA+10,NJM                       !NIYDB-1       !NJM
						  T(1,J) = T(2,J)
					         ENDDO

							DO J = NIYDB,NIYDB+9
						  T(1,J) = T(2,J) + HDX(1)
					         ENDDO							

							DO J = NIYDB+10,NJM              !NIYDC-1
						  T(1,J) = T(2,J)
					         ENDDO																					

!							DO J = NIYDC,NIYDC+10
!						  T(1,J) = T(2,J) + HDX(1)
!					         ENDDO							

!							DO J = NIYDC+11,NJM
!						  T(1,J) = T(2,J)
!							ENDDO
!RIGHT - the unknown Qu1 wall
                            DO J = 2,NJM
     				     T(NI,J) = 0.00
					         ENDDO
RETURN
END





SUBROUTINE SUBUDRA(ERROU,OMIGU,ISTATE,NSCHEMEUV,Pr,Ra,theta)  !,JNNAR,JNNBR)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)
    COMMON/ARPREUV/ BU(NX,NY),BV(NX,NY)

!    COMMON/RDOBSTACLESPTCA/INWAR,INEAR,JNSAR !,JNNA
!    COMMON/RDOBSTACLESPTCB/INWBR,INEBR,JNSBR !,JNNB

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SU(NX,NY),SOU(NX,NY)
! DIMENSION AU(NX,NY)


       COEFU = SQRT(Pr/Ra)    ! Form B

         NIM = NI-2;   NJM = NJ-1
		  
DO I = 2,NIM
DO J = 2,NJM

 	      FE = (U(I,J)+U(I+1,J))*DY(J)/2.0
  	      DE = COEFU*DY(J)/DX(I+1)
     AE(I,J) = DE*SCHEME(FE/DE,NSCHEMEUV)+AMAX1(-FE,0.0)

	      FW = (U(I,J)+U(I-1,J))*DY(J)/2.0
	      DW = COEFU*DY(J)/DX(I)
     AW(I,J) = DW*SCHEME(FW/DW,NSCHEMEUV)+AMAX1( FW,0.0)

	      FN = (DX(I+1)*V(I,J)+DX(I)*V(I+1,J))/2.0   
 	      DN = COEFU*HDX(I)/HDY(J)
     AN(I,J) = DN*SCHEME(FN/DN,NSCHEMEUV)+AMAX1(-FN,0.0)

	      FS = (DX(I+1)*V(I,J-1)+DX(I)*V(I+1,J-1))/2.0 
	      DS = COEFU*HDX(I)/HDY(J-1)
     AS(I,J) = DS*SCHEME(FS/DS,NSCHEMEUV)+AMAX1( FS,0.0)

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
         SAD =-FECDS+FEUDS+FWCDS-FWUDS-FNCDS+FNUDS+FSCDS-FSUDS
     ELSE
	     SAD = 0.0
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
	     SAD = (1.0/2.0)*(SADEP+SADEN+SADWP+SADWN+SADNP+SADNN+SADSP+SADSN)
     ELSE
	     SAD = 0.0
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
	     SAD = (1.0/8.0)*(SADEP+SADEN+SADWP+SADWN+SADNP+SADNN+SADSP+SADSN)
     ELSE
	     SAD = 0.0
     ENDIF

!---Relaxation through inertia 
     IF(ISTATE==1)THEN
	   STIME = 0.0
	   CTIME = 0.0
     ELSE
	   STIME = 0.1*U(I,J)
	   CTIME = 0.1
     ENDIF

    SOU(I,J) = DY(J)*(P(I,J)-P(I+1,J))+STIME+SAD-COSD(theta)*DX(I)*(T(I,J)*DY(J+1) + T(I,J+1)*DY(J))/2.0 

     SU(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)+CTIME 

     BU(I,J) = DY(J)/SU(I,J)

ENDDO
ENDDO

CALL SLORXY(U, 1, ERROU, NIM, NJM, AE, AW, AN, AS, SU, SOU, OMIGU) 

RETURN
END





SUBROUTINE SUBVDRA(ERROV,OMIGV,ISTATE,NSCHEMEUV,Pr,Ra,theta) ! ,JNNAR,JNNBR)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/ARDPRESSURE/ P(NX,NY)
    COMMON/ARPREUV/ BU(NX,NY),BV(NX,NY)
	
!	COMMON/RDOBSTACLESPTCA/INWAR,INEAR,JNSAR !,JNNA
!    COMMON/RDOBSTACLESPTCB/INWBR,INEBR,JNSBR !,JNNB

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SV(NX,NY),SOV(NX,NY)
! DIMENSION AV(NX,NY)

!      COEFV = Pr              ! Form A
       COEFV = SQRT(Pr/Ra)     ! Form B

         NIM = NI-1;  NJM = NJ-2
	   
DO  I = 2,NIM
DO  J = 2,NJM

          FE = (DY(J+1)*U(I,J)+DY(J)*U(I,J+1))/2.0	    
  	      DE = COEFV*HDY(J)/HDX(I)
     AE(I,J) = DE*SCHEME(FE/DE,NSCHEMEUV)+AMAX1(-FE,0.0)

  	      FW = (DY(J+1)*U(I-1,J)+DY(J)*U(I-1,J+1))/2.0	
	      DW = COEFV*HDY(J)/HDX(I-1)
     AW(I,J) = DW*SCHEME(FW/DW,NSCHEMEUV)+AMAX1( FW,0.0)

          FN = (V(I,J)+V(I,J+1))*DX(I)/2.0
          DN = COEFV*DX(I)/DY(J+1)
     AN(I,J) = DN*SCHEME(FN/DN,NSCHEMEUV)+AMAX1(-FN,0.0)

          FS = (V(I,J)+V(I,J-1))*DX(I)/2.0
	      DS = COEFV*DX(I)/DY(J)
     AS(I,J) = DS*SCHEME(FS/DS,NSCHEMEUV)+AMAX1( FS,0.0)

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
         SAD = -FECDS+FEUDS+FWCDS-FWUDS-FNCDS+FNUDS+FSCDS-FSUDS
     ELSE
         SAD = 0.0
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
	     SAD = (1.0/2.0)*(SADEP+SADEN+SADWP+SADWN+SADNP+SADNN+SADSP+SADSN)
     ELSE
	     SAD = 0.0
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
	     SAD = (1.0/8.0)*(SADEP+SADEN+SADWP+SADWN+SADNP+SADNN+SADSP+SADSN)
     ELSE
	     SAD = 0.0
     ENDIF

!---Relaxation through inertia 
     IF(ISTATE==1)THEN
	   STIME = 0.0
	   CTIME = 0.0
     ELSE
	   STIME = 0.1*V(I,J)
	   CTIME = 0.1
     ENDIF

    SOV(I,J) = DX(I)*(P(I,J)-P(I,J+1))+STIME+SAD                              &  
   			 +                SIND(theta)*DX(I)*(T(I,J)*DY(J+1)+T(I,J+1)*DY(J))/2.0       ! For the form B

     SV(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)+CTIME    !+ AV(I,J)

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

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)
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
     IF(I.NE.NIM) U(I,J) = U(I,J)+BU(I,J)*(PC(I,J)-PC(I+1,J))  
     IF(J.NE.NJM) V(I,J) = V(I,J)+BV(I,J)*(PC(I,J)-PC(I,J+1))
                  P(I,J) = P(I,J)+0.6*PC(I,J)
          ENDDO
          ENDDO

RETURN
END



SUBROUTINE SUBTDRA(ERROT,OMIGT,NSCHEMET,Pr,Ra)       !,JNNAR,JNNBR)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)
	
!	COMMON/ARDOBSTACLESPTCA/INWAR,INEAR,JNSAR !,JNNA
!   COMMON/ARDOBSTACLESPTCB/INWBR,INEBR,JNSBR !,JNNB

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),ST(NX,NY),SOT(NX,NY)
! DIMENSION COEFT(NX,NY),COKR(NX,NY)

         NIM = NI-1;  NJM = NJ-1

!      COEFT = 1.00                ! Form A
       COEFT = 1.0/SQRT(Pr*Ra)     ! Form B

DO I = 2,NIM
DO J = 2,NJM

          De =COEFT*DY(J)/HDX(I)
          Fe = U(I,J)*DY(J)
     AE(I,J) = De * SCHEME(Fe/De,NSCHEMET)+AMAX1(-Fe,0.0)

          Dw =COEFT*DY(J)/HDX(I-1)
          Fw = U(I-1,J)*DY(J)
     AW(I,J) = Dw * SCHEME(Fw/Dw,NSCHEMET)+AMAX1( Fw,0.0)

          Dn =COEFT*DX(I)/HDY(J)
          Fn = V(I,J)*DX(I)
     AN(I,J) = Dn * SCHEME(Fn/Dn,NSCHEMET)+AMAX1(-Fn,0.0)

          Ds =COEFT*DX(I)/HDY(J-1)
          Fs = V(I,J-1)*DX(I)
     AS(I,J) = Ds * SCHEME(Fs/Ds,NSCHEMET)+AMAX1( Fs,0.0)

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
         SAD = -FeCDS+FeUDS+FwCDS-FwUDS-FnCDS+FnUDS+FsCDS-FsUDS
  ELSE
         SAD = 0.0
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
	     SAD = (1.0/2.0)*(SADEP+SADEN+SADWP+SADWN+SADNP+SADNN+SADSP+SADSN)
  ELSE
         SAD = 0.0
  ENDIF

!---Deferrd correction QUICK Scheme
  IF(NSCHEMET==6)THEN
       SADEP = -(-   T(I-1,J)-2.0*T(I,J)  +3.0*T(I+1,J))*AMAX1(+Fe,0.0)
    IF(I==NIM)THEN
       SADEN = 0.00
    ELSE
       SADEN = +(3.0*T(I,J)  -2.0*T(I+1,J)-    T(I+2,J))*AMAX1(-Fe,0.0)
    ENDIF
    IF(I==2)THEN
       SADWP = 0.00
    ELSE
       SADWP = +(-   T(I-2,J)-2.0*T(I-1,J)+3.0*T(I,J)  )*AMAX1(+Fw,0.0)
    ENDIF   	  	
       SADWN = -(3.0*T(I-1,J)-2.0*T(I,J)  -    T(I+1,J))*AMAX1(-Fw,0.0)   	  	
       SADNP = -(-   T(I,J-1)-2.0*T(I,J)  +3.0*T(I,J+1))*AMAX1(+Fn,0.0)  
    IF(J==NJM)THEN
       SADNN = 0.00
    ELSE
       SADNN = +(3.0*T(I,J)  -2.0*T(I,J+1)-    T(I,J+2))*AMAX1(-Fn,0.0) 
    ENDIF 	  	
    IF(J==2)THEN
       SADSP = 0.00
    ELSE
       SADSP = +(-   T(I,J-2)-2.0*T(I,J-1)+3.0*T(I,J)  )*AMAX1(+Fs,0.0) 
    ENDIF 	  	
       SADSN = -(3.0*T(I,J-1)-2.0*T(I,J)  -    T(I,J+1))*AMAX1(-Fs,0.0)
         SAD = (1.0/8.0)*(SADEP+SADEN+SADWP+SADWN+SADNP+SADNN+SADSP+SADSN)
  ELSE
         SAD = 0.0
  ENDIF

    SOT(I,J) = 0.0+SAD	                    
   
     ST(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)

ENDDO
ENDDO

        CALL SLORXY(T, 1, ERROT, NIM, NJM, AE, AW, AN, AS, ST, SOT, OMIGT) 

RETURN
END





SUBROUTINE VISULIZATIONDRA(LSAVE)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/ARDVARIABLES/U(NX,NY),V(NX,NY),T(NX,NY)
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
     FIField = 'RDF'//IC//JC//KC//LC	   !Plot the direct convection field

!---Dynamic save the flow fields
         OPEN(2, FILE=FIField, STATUS='UNKNOWN')
        WRITE(2,201) 
201	   FORMAT(1X,'TITLE="Field.DAT"')
        WRITE(2,202) 
202	   FORMAT(1X,'VARIABLES="X","Y","T","U","V"')
        WRITE(2,203) NI,NJ
203	   FORMAT(1X,'ZONE I=',I3,3X,'J=',I3,3X,'F=POINT')
        DO J = 1,NJ
        DO I = 1,NI
        WRITE(2,204) XC(I),YC(J),T(I,J),U1(I,J),V1(I,J)
        ENDDO
        ENDDO
        CLOSE(2)
204    FORMAT(1X,2(F9.6,1X),1X,3(E15.6,1X))


RETURN
END