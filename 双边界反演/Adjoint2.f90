SUBROUTINE ADJOINT2(LOOP, Ra, Pr, FIQ2)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/APRESSURE2/ P(NX,NY)

    COMMON/DVARIABLES/UD(NX,NY),VD(NX,NY),TD(NX,NY)

 DIMENSION FIQ2(NX,NY)


              NIM = NI-1;  NJM = NJ-1

!--- Constant and Demo Parameters  
            OMIGU = 0.60;    OMIGV = 0.60;   OMIGPC = 0.90;    OMIGT = 0.80

     ConvergenceA = 1.0E-10

           INLOOP = 8500								


!--- Difference schemes for convective terms driven by direct flows           
        NSCHEMEUV = 3     !1 - HYBRID;	2 - POWERLAW; 3 - UPWIND(1st-order)
         NSCHEMET = 3     !1 - HYBRID;	2 - POWERLAW; 3 - UPWIND(1st-order)
!--- Difference schemes for time terms
           ISTATE = 1     !1 - Œ»Ã¨;  2 - ∂ØÃ¨


!--- Initialization and Dirichlet boundary conditions
                U = 0.0;   V = 0.0;   P = 0.0;   T = 0.0

!--- Obtain the direct flow field
!             OPEN(1,FILE='F:\INCP\DirectFlowTemp.dat',STATUS='UNKNOWN')
!            DO J = 1,NJ
!            DO I = 1,NI
!             READ(1,911) UD(I,J),VD(I,J),TD(I,J)
!            ENDDO
!            ENDDO
!            CLOSE(1)
!911        FORMAT(1X,3(E15.6,1X))

!--- Define the Dirichlet boundary conditions!
!       DO J = 2,NJM               !---Depending both on the grid arrangement and segment location <1>
!     T(1,J) = 1.00   
! 	    ENDDO


!--- Inner Iteration!
      DO ITER = 1, INLOOP

         CALL   SUBTCNeuumanA2

         CALL   SUBUA2(ERROU,OMIGU,ISTATE,NSCHEMEUV,Pr,Ra)

 	     CALL   SUBVA2(ERROV,OMIGV,ISTATE,NSCHEMEUV,Pr,Ra)

         CALL  SUBPCA2(ERROPC,OMIGPC)

         CALL   SUBTA2(ERROT,OMIGT,NSCHEMET,Pr,Ra,FIQ2)

                AMAXA = AMAX1(ERROU,ERROV,ERROPC,ERROT)
 
             IF(AMAXA.LE.ConvergenceA) EXIT	

             IF(ITER==1.OR.MOD(ITER,1000)==0) THEN
             WRITE(*,110)  ITER, ERROU,ERROV,ERROPC,ERROT
             ENDIF
110         FORMAT(2X,I6,3X,4(E15.6E3,3X))

      ENDDO
!End of inner iteration!


!---Record the Flow, Heatflow, and Isotherm Fields 
IF(LOOP==0.OR.MOD(LOOP,10)==0)THEN
CALL VISULIZATIONA2(LOOP)
ENDIF  

RETURN
END
!--* * * * * * * * * THE END OF THE PROGRAM MAIN     * * * * * * * * * *





SUBROUTINE SUBTCNeuumanA2

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)

          NIM = NI-1;   NJM = NJ-1

         DO J = 2,NJM
       T(1,J) = T(2,J)
         ENDDO

         DO J = 2,NJM                        !---Unknown function
      T(NI,J) = T(NIM,J)
	     ENDDO

RETURN
END





SUBROUTINE SUBUA2(ERROU,OMIGU,ISTATE,NSCHEMEUV,Pr,Ra)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

   COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)
   COMMON/APRESSURE2/ P(NX,NY)
   COMMON/APREUV2/ BU(NX,NY),BV(NX,NY)

   COMMON/DVARIABLES/UD(NX,NY),VD(NX,NY),TD(NX,NY)

DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SU(NX,NY),SOU(NX,NY)

!      COEFU = Pr             ! Form A
       COEFU = SQRT(Pr/Ra)    ! Form B

         NIM = NI-2;   NJM = NJ-1

DO I = 2,NIM
DO J = 2,NJM

!---Adjoint Convective terms by Direct flow (UDS scheme) SOU2
!---Diffusion terms of component UA (CDS scheme)		 SOU5
 	     FED =-(UD(I,J)+UD(I+1,J))*DY(J)/2.0
  	      DE = COEFU*DY(J)/DX(I+1)
     AE(I,J) = DE*SCHEME(FED/DE,NSCHEMEUV)+AMAX1(-FED,0.0)

	     FWD =-(UD(I-1,J)+UD(I,J))*DY(J)/2.0
	      DW = COEFU*DY(J)/DX(I)
     AW(I,J) = DW*SCHEME(FWD/DW,NSCHEMEUV)+AMAX1( FWD,0.0)

	     FND =-(DX(I+1)*VD(I,J)+DX(I)*VD(I+1,J))/2.0   
 	      DN = COEFU*HDX(I)/HDY(J)
     AN(I,J) = DN*SCHEME(FND/DN,NSCHEMEUV)+AMAX1(-FND,0.0)

	     FSD =-(DX(I+1)*VD(I,J-1)+DX(I)*VD(I+1,J-1))/2.0 
	      DS = COEFU*HDX(I)/HDY(J-1)
     AS(I,J) = DS*SCHEME(FSD/DS,NSCHEMEUV)+AMAX1( FSD,0.0)

!---Transpose matrix of vector UD
		 VPP = (1.0-FX(I))*(V(I,J-1)+V(I,J))/2.0 + FX(I)*(V(I+1,J-1)+V(I+1,J))/2.0
		VDYP = (VD(I+1,J-1) + VD(I+1,J) - VD(I,J-1) - VD(I,J))/2.0  
	    SOU3 = U(I,J)*(UD(I+1,J)-UD(I-1,J))*DY(J)/2.0 + VPP*VDYP*DY(J)

!---Thermal gradient terms along X-direction
         TPP = (1.0-FX(I))*T(I,J) + FX(I)*T(I+1,J)
        SOU6 = TPP*(TD(I+1,J)-TD(I,J))*DY(J)

!---Relaxation through inertia 
     IF(ISTATE==1)THEN
	   STIME = 0.0
	   CTIME = 0.0
     ELSE
	   STIME = 0.1*U(I,J)
	   CTIME = 0.1
     ENDIF


    SOU(I,J) = DY(J)*(P(I,J)-P(I+1,J)) + STIME - SOU3 - SOU6 

     SU(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J) + CTIME

     BU(I,J) = DY(J)/SU(I,J)

ENDDO
ENDDO

     CALL SLORXY(U, 1, ERROU, NIM, NJM, AE, AW, AN, AS, SU, SOU, OMIGU) 

RETURN
END





SUBROUTINE SUBVA2(ERROV,OMIGV,ISTATE,NSCHEMEUV,Pr,Ra)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/APRESSURE2/ P(NX,NY)
    COMMON/APREUV2/ BU(NX,NY),BV(NX,NY)

    COMMON/DVARIABLES/UD(NX,NY),VD(NX,NY),TD(NX,NY)

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),SV(NX,NY),SOV(NX,NY)

!      COEFV = Pr              ! Form A
       COEFV = SQRT(Pr/Ra)     ! Form B

         NIM = NI-1;  NJM = NJ-2
	   
DO  I = 2,NIM
DO  J = 2,NJM

!---Adjoint Convective terms by Direct flow (UDS scheme) SOV2
!---Diffusion terms of component VA (CDS scheme)		 SOV5
         FED =-(DY(J+1)*UD(I,J)+DY(J)*UD(I,J+1))/2.0	    
  	      DE = COEFV*HDY(J)/HDX(I)
     AE(I,J) = DE*SCHEME(FED/DE,NSCHEMEUV)+AMAX1(-FED,0.0)

  	     FWD =-(DY(J+1)*UD(I-1,J)+DY(J)*UD(I-1,J+1))/2.0
	      DW = COEFV*HDY(J)/HDX(I-1)
     AW(I,J) = DW*SCHEME(FWD/DW,NSCHEMEUV)+AMAX1( FWD,0.0)

         FND =-(VD(I,J)+VD(I,J+1))*DX(I)/2.0
          DN = COEFV*DX(I)/DY(J+1)
     AN(I,J) = DN*SCHEME(FND/DN,NSCHEMEUV)+AMAX1(-FND,0.0)

         FSD =-(VD(I,J-1)+VD(I,J))*DX(I)/2.0
	      DS = COEFV*DX(I)/DY(J)
     AS(I,J) = DS*SCHEME(FSD/DS,NSCHEMEUV)+AMAX1( FSD,0.0)

!---Transpose matrix of vector VD
         UPP = (1.0-FY(J))*(U(I-1,J)+U(I,J))/2.0 + FY(J)*(U(I-1,J+1)+U(I,J+1))/2.0
        UDYP = (UD(I-1,J+1) + UD(I,J+1) - UD(I-1,J) - UD(I,J))/2.0
        SOV3 = UPP*UDYP*DX(I) + V(I,J)*(VD(I,J+1)-VD(I,J-1))*DX(I)/2.0

!---Thermal gradient terms along X-direction
         TPP = (1.0-FY(J))*T(I,J) + FY(J)*T(I,J+1)
        SOV6 = TPP*(TD(I,J+1)-TD(I,J))*DX(I)

!---Relaxation through inertia 
     IF(ISTATE==1)THEN
	   STIME = 0.0
	   CTIME = 0.0
     ELSE
	   STIME = 0.1*V(I,J)
	   CTIME = 0.1
     ENDIF


    SOV(I,J) = DX(I)*(P(I,J)-P(I,J+1)) + STIME - SOV3 - SOV6  				  
!   SOV(I,J) = DX(I)*(P(I,J)-P(I,J+1)) + STIME - SOV3  				          &
!            -          Ra*Pr*DX(I)*(T(I,J)*DY(J+1)+T(I,J+1)*DY(J))/2.0       ! For the form A 
!  			 -                DX(I)*(T(I,J)*DY(J+1)+T(I,J+1)*DY(J))/2.0       ! For the form B

     SV(I,J) = AE(I,J)+AW(I,J)+AN(I,J)+AS(I,J)+CTIME

     BV(I,J) = DX(I)/SV(I,J)

ENDDO
ENDDO
       
CALL SLORXY(V, 1, ERROV, NIM, NJM, AE, AW, AN, AS, SV, SOV, OMIGV) 

RETURN
END





SUBROUTINE SUBPCA2(ERROPC,OMIGPC)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/APRESSURE2/ P(NX,NY)
    COMMON/APREUV2/ BU(NX,NY),BV(NX,NY)

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





SUBROUTINE SUBTA2(ERROT,OMIGT,NSCHEMET,Pr,Ra,FIQ2)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/DVARIABLES2/UD(NX,NY),VD(NX,NY),TD(NX,NY)

 DIMENSION FIQ2(NX,NY)

 DIMENSION AE(NX,NY),AW(NX,NY),AN(NX,NY),AS(NX,NY),ST(NX,NY),SOT(NX,NY)


         NIM = NI-1;  NJM = NJ-1

!      COEFT = 1.00                ! Form A
       COEFT = 1.0/SQRT(Pr*Ra)     ! Form B

DO I = 2,NIM
DO J = 2,NJM

!---Adjoint Convective terms by Direct flow (UDS scheme) SOT2
!---Diffusion terms of component TA (CDS scheme)		 SOT3
          De = COEFT*DY(J)/HDX(I)
         FeD =-UD(I,J)*DY(J)
     AE(I,J) = De * SCHEME(FeD/De,NSCHEMET)+AMAX1(-FeD,0.0)

          Dw = COEFT*DY(J)/HDX(I-1)
         FwD =-UD(I-1,J)*DY(J)
     AW(I,J) = Dw * SCHEME(FwD/Dw,NSCHEMET)+AMAX1( FwD,0.0)

          Dn = COEFT*DX(I)/HDY(J)
         FnD =-VD(I,J)*DX(I)
     AN(I,J) = Dn * SCHEME(FnD/Dn,NSCHEMET)+AMAX1(-FnD,0.0)

          Ds = COEFT*DX(I)/HDY(J-1)
         FsD =-VD(I,J-1)*DX(I)
     AS(I,J) = Ds * SCHEME(FsD/Ds,NSCHEMET)+AMAX1( FsD,0.0)

!---Vertical Velocity term
        SOT4 = 0.5*(V(I,J-1) + V(I,J))*DX(I)*DY(J)

!---Measurement points Res("Inner Heat Sources")
        SOT5 = FIQ2(I,J)


    SOT(I,J) = 0.0 + SOT4 + SOT5	                    
   
     ST(I,J) = AE(I,J) + AW(I,J) + AN(I,J) + AS(I,J)

ENDDO
ENDDO

        CALL SLORXY(T, 1, ERROT, NIM, NJM, AE, AW, AN, AS, ST, SOT, OMIGT) 

RETURN
END





SUBROUTINE VISULIZATIONA2(LSAVE)

   INCLUDE 'PARAMETER.inc'
   INCLUDE 'GRID.inc'
   INCLUDE 'MESH.inc'

    COMMON/AVARIABLES2/U(NX,NY),V(NX,NY),T(NX,NY)
    COMMON/APRESSURE2/ P(NX,NY)

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


!Save the loop numbers as output file name. !The maximum loop is no more than 10000.
         III =  LSAVE/1000
	     JJJ = (LSAVE-1000*III)/100
	     KKK = (LSAVE-1000*III-100*JJJ)/10
	     LLL =  LSAVE-1000*III-100*JJJ-10*KKK
	      IC = CHAR(48+III)
	      JC = CHAR(48+JJJ)
	      KC = CHAR(48+KKK)
	      LC = CHAR(48+LLL)
     FIField = 'AF'//IC//JC//KC//LC	   !Plot the Adjoint field

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