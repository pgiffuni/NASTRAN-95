SUBROUTINE ofpgpw (*,FILE,out,from)
     
!     PRINT GRID POINT WEIGHT GENERATORN TABLE
!     (SOURCE PROGRAM ORIGINALLY CODED IN OFP)
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 DOUBLE PRECISION, INTENT(IN OUT)         :: out(1)
 INTEGER, INTENT(OUT)                     :: from
 INTEGER :: flag, of(5)
 
 COMMON /system/  ibuf,l,dummy(10),line
 COMMON /zzzzzz/  core(1)
 EQUIVALENCE      (l1,of(1),core(1)), (l2,of(2)), (l3,of(3)),  &
     (l4,of(4)), (l5,of(5))
 
!     FOR GRIDPOINT WEIGHT OUTPUT ONLY ONE DATA VECTOR OF 78 WORDS
!     IS EXPECTED AND IT IS THUS READ AND OUTPUT EXPLICITLY
!     (CHANGED TO D.P. BY G.CHAN/UNISYS, AND THEREFORE 156 WORDS.
!     THIS RECORD IS SENT OVER BY GPWG1B, WHICH IS NOW A D.P. ROUTINE)
 
 from = 345
 CALL READ (*2020,*60,FILE,out(1),90,0,flag)
 l1 = 0
 l2 = 0
 l3 = 202
 l4 = 0
 l5 = 0
 CALL ofp1
 line = line + 44
 WRITE  (l,350) (out(i),i=1,45)
 350 FORMAT (37X, 'MO - RIGID BODY MASS MATRIX IN BASIC COORDINATE SYSTEM',  &
     /16X,3H***,93X,3H***, /6(16X,1H*,1P,6D16.8,2H *,/),16X,  &
     3H***,93X,3H***, /40X,  &
     51HS - transformation matrix for scalar mass partition,  &
     /2(40X,3H***,5X), /3(40X,1H*,1P,3D16.8,2H *,/),2(40X,3H***,  &
     5X),  /25X,9HDIRECTION, /20X,20HMASS axis system (s),7X,  &
     4HMASS,17X,6HX-c.g.,11X,6HY-c.g.,11X,6HZ-c.g.)
 from = 355
 CALL READ (*2020,*60,FILE,out(1),66,1,flag)
 WRITE  (l,360) (out(i),i=1,12)
 360 FORMAT (28X,1HX,1P,d27.9,1P,d21.9,1P,2D17.9,/  &
     28X,1HY,1P,d27.9,1P,d21.9,1P,2D17.9,/ 28X,1HZ,1P,d27.9,1P,d21.9,1P,2D17.9)
 WRITE  (l,370) (out(i),i=13,33)
 370 FORMAT (/49X,33HI(s) - inertias relative TO c.g. , /2(38X,3H***,  &
     11X), /3(38X,1H*,1P,3D17.9,3H  *,/),2(38X,3H***,11X), /54X,  &
     25HI(q) - principal inertias, /2(38X,3H***,11X), /38X,1H*,  &
     1P,d17.9,36X,1H*, /38X,1H*,1P,d34.9,19X,1H*, /38X,1H*,1P,  &
     d51.9,3H  *, /2(38X,3H***,11X), /44X,  &
     44HQ - transformation matrix - i(q) = qt*i(s)*q, /2(38X,  &
     3H***,11X),/3(38X,1H*,1P,3D17.9,3H  *,/),2(38X,3H***,11X))
 60  RETURN
 
 2020 RETURN 1
END SUBROUTINE ofpgpw
