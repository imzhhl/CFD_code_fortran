FUNCTION SCHEME(PE,ISCHEME)

    IF(ISCHEME==1)THEN
       SCHEME = AMAX1(0.,1.-ABS(PE)*0.5)                   ! Hybrid Scheme
ELSEIF(ISCHEME==2)THEN
	   SCHEME = AMAX1(0.,1.-ABS(PE)*0.1)				   ! Power-law Scheme
	   SCHEME = SCHEME**5.								   
ELSEIF(ISCHEME==3.OR.ISCHEME==4.OR.ISCHEME==5.OR.ISCHEME==6)THEN
	   SCHEME = 1.										   ! Upwind Scheme
    ENDIF

RETURN
END
