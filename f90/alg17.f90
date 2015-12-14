SUBROUTINE alg17 (istak,pltsze,itrig,title,ikdum,ifplot)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: istak
 REAL, INTENT(IN)                         :: pltsze
 INTEGER, INTENT(IN OUT)                  :: itrig
 REAL, INTENT(IN OUT)                     :: title(18)
 INTEGER, INTENT(IN OUT)                  :: ikdum
 INTEGER, INTENT(IN OUT)                  :: ifplot
 
 
 plttit=pltsze*.1
 IF (istak < 2) GO TO 10
 bal=.35*pltsze
 xlen1=.3*pltsze
 xlen2=xlen1
 ylen1=.25*pltsze
 ylen2=-1.*ylen1
 xback1=-1.9
 xback2=-6.2
 GO TO 50
 10    IF (istak == 0) GO TO 20
 xlen1=.70*pltsze
 xlen2=.15*pltsze
 xback1=-1.9-.20*pltsze
 xback2=-6.2-.20*pltsze
 IF (ikdum == 1) GO TO 30
 GO TO 40
 20    CONTINUE
 xlen1=.15*pltsze
 xlen2=.70*pltsze
 xback1=-1.9+.20*pltsze
 xback2=-6.2+.20*pltsze
 IF (ikdum == 1) GO TO 40
 30    bal=.25*pltsze
 ylen1=.50*pltsze
 ylen2=-.15*pltsze
 GO TO 50
 40    bal=.50*pltsze
 ylen1=.15*pltsze
 ylen2=-.50*pltsze
 50    CONTINUE
 yback1=-(.35+bal)
 yback2=yback1-.01*pltsze-.175
 SELECT CASE ( itrig )
   CASE (    1)
     GO TO 60
   CASE (    2)
     GO TO 70
 END SELECT
 60    CONTINUE
 GO TO 80
 70    xback1=xback1+0.35
 80    CONTINUE
 RETURN
END SUBROUTINE alg17
