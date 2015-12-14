SUBROUTINE fndiam (snd1,snd2,ndstk,ndeg,lvl,lvls1,lvls2,iwk,  &
        idflt,ndlst,jwk,idim)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     FNDIAM IS THE CONTROL PROCEDURE FOR FINDING THE PSEUDO-DIAMETER
!     OF NDSTK AS WELL AS THE LEVEL STRUCTURE FROM EACH END
 
!     SND1-     ON INPUT THIS IS THE NODE NUMBER OF THE FIRST
!               ATTEMPT AT FINDING A DIAMETER.  ON OUTPUT IT
!               CONTAINS THE ACTUAL NUMBER USED.
!     SND2-     ON OUTPUT CONTAINS OTHER END OF DIAMETER
!     LVLS1-    ARRAY CONTAINING LEVEL STRUCTURE WITH SND1 AS ROOT
!     LVLS2-    ARRAY CONTAINING LEVEL STRUCTURE WITH SND2 AS ROOT
!     IDFLT-    FLAG USED IN PICKING FINAL LEVEL STRUCTURE, SET =1
!               IF WIDTH OF LVLS1 .GE. WIDTH OF LVLS2, OTHERWISE =2
!     LVL,IWK-  WORKING STORAGE
!     JWK-      WORKING STORAGE, CURRENTLY SHARING SAME SPACE WITH RENUM
!     DIMENSION OF NDLST IS THE MAX NUMBER OF NODES IN LAST LEVEL.
 
 
 INTEGER, INTENT(IN OUT)                  :: snd1
 INTEGER, INTENT(OUT)                     :: snd2
 INTEGER, INTENT(IN OUT)                  :: ndstk(1)
 INTEGER, INTENT(IN OUT)                  :: ndeg(1)
 INTEGER, INTENT(OUT)                     :: lvl(1)
 INTEGER, INTENT(OUT)                     :: lvls1(1)
 INTEGER, INTENT(OUT)                     :: lvls2(1)
 INTEGER, INTENT(IN OUT)                  :: iwk(1)
 INTEGER, INTENT(OUT)                     :: idflt
 INTEGER, INTENT(IN)                      :: ndlst(idim)
 INTEGER, INTENT(IN OUT)                  :: jwk(1)
 INTEGER, INTENT(IN OUT)                  :: idim
 INTEGER :: flag,     snd
 
 COMMON /bandb /  dum3b(3), ngrid
 COMMON /bandg /  n,        idpth
 
 flag=0
 mtw2=n
 snd=snd1
 
!     ZERO LVL TO INDICATE ALL NODES ARE AVAILABLE TO TREE
 
 20 DO  i=1,n
   lvl(i)=0
 END DO
 lvln=1
 
!     DROP A TREE FROM SND
 
 CALL tree (snd,ndstk,lvl,iwk,ndeg,lvlwth,lvlbot,lvln,maxlw,mtw2, jwk)
 IF (flag >= 1) GO TO 110
 flag=1
 70 idpth=lvln-1
 mtw1=maxlw
 
!     COPY LEVEL STRUCTURE INTO LVLS1
 
 DO  i=1,n
   lvls1(i)=lvl(i)
 END DO
 ndxn=1
 ndxl=0
 mtw2=n
 
!     SORT LAST LEVEL BY DEGREE  AND STORE IN NDLST
 
 CALL sortdg (ndlst,iwk(lvlbot),ndxl,lvlwth,ndeg)
 IF (ndxl <= idim) GO TO 100
 
!     DIMENSION EXCEEDED  . . .  STOP JOB.
 
 80 ngrid=-3
 RETURN
 
 100 CONTINUE
 snd=ndlst(1)
 GO TO 20
 110 IF (idpth >= lvln-1) GO TO 120
 
!     START AGAIN WITH NEW STARTING NODE
 
 snd1=snd
 GO TO 70
 120 IF (maxlw >= mtw2) GO TO 130
 mtw2=maxlw
 snd2=snd
 
!     STORE NARROWEST REVERSE LEVEL STRUCTURE IN LVLS2
 
 DO  i=1,n
   lvls2(i)=lvl(i)
 END DO
 130 IF (ndxn == ndxl) GO TO 140
 
!     TRY NEXT NODE IN NDLST
 
 ndxn=ndxn+1
 snd=ndlst(ndxn)
 GO TO 20
 140 idflt=1
 IF (mtw2 <= mtw1) idflt=2
 IF (idpth > idim) GO TO 80
 RETURN
END SUBROUTINE fndiam
