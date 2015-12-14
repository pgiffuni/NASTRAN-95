SUBROUTINE tree (iroot,ndstk,lvl,iwk,ndeg,lvlwth,lvlbot,lvln,  &
        maxlw,ibort,jwk)
     
!     TREE DROPS A TREE IN NDSTK FROM IROOT
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     LVL-      ARRAY INDICATING AVAILABLE NODES IN NDSTK WITH ZERO
!               ENTRIES. TREE ENTERS LEVEL NUMBERS ASSIGNED
!               DURING EXECUTION OF OF THIS PROCEDURE
!     IWK-      ON OUTPUT CONTAINS NODE NUMBERS USED IN TREE
!               ARRANGED BY LEVELS (IWK(LVLN) CONTAINS IROOT
!               AND IWK(LVLBOT+LVLWTH-1) CONTAINS LAST NODE ENTERED)
!     JWK-      ON ONTPUT CONTAINS A ROW OF UNPACKED GRID NOS.
!               CURRENTLY, JWK AND RENUM SHARE SAME CORE SPACE
!     LVLWTH-   ON OUTPUT CONTAINS WIDTH OF LAST LEVEL
!     LVLBOT-   ON OUTPUT CONTAINS INDEX INTO IWK OF FIRST
!               NODE IN LAST LEVEL
!     MAXLW-    ON OUTPUT CONTAINS THE MAXIMUM LEVEL WIDTH
!     LVLN-     ON INPUT THE FIRST AVAILABLE LOCATION IN IWK
!               USUALLY ONE BUT IF IWK IS USED TO STORE PREVIOUS
!               CONNECTED COMPONENTS, LVLN IS NEXT AVAILABLE LOCATION.
!               ON OUTPUT THE TOTAL NUMBER OF LEVELS + 1
!     IBORT-    INPUT PARAM WHICH TRIGGERS EARLY RETURN IF
!               MAXLW BECOMES .GE. IBORT
 
!     INTEGER          BUNPK
 
 
 INTEGER, INTENT(IN)                      :: iroot
 INTEGER, INTENT(IN OUT)                  :: ndstk(1)
 INTEGER, INTENT(OUT)                     :: lvl(1)
 INTEGER, INTENT(OUT)                     :: iwk(1)
 INTEGER, INTENT(IN)                      :: ndeg(1)
 INTEGER, INTENT(OUT)                     :: lvlwth
 INTEGER, INTENT(OUT)                     :: lvlbot
 INTEGER, INTENT(IN OUT)                  :: lvln
 INTEGER, INTENT(OUT)                     :: maxlw
 INTEGER, INTENT(IN)                      :: ibort
 INTEGER, INTENT(IN)                      :: jwk(1)
 
 
 maxlw =0
 itop  =lvln
 inow  =lvln
 lvlbot=lvln
 lvltop=lvln+1
 lvln  =1
 lvl(iroot)=1
 iwk(itop) =iroot
 30 lvln  =lvln+1
 35 iwknow=iwk(inow)
 ndrow =ndeg(iwknow)
 CALL bunpak(ndstk,iwknow,ndrow,jwk)
 DO  j=1,ndrow
   itest=jwk(j)
   IF (lvl(itest) /= 0) CYCLE
   lvl(itest)=lvln
   itop=itop+1
   iwk(itop)=itest
 END DO
 inow=inow+1
 IF (inow < lvltop) GO TO 35
 lvlwth=lvltop-lvlbot
 IF (maxlw < lvlwth) maxlw=lvlwth
 IF (maxlw >= ibort .OR. itop < lvltop) RETURN
 lvlbot=inow
 lvltop=itop+1
 GO TO 30
END SUBROUTINE tree
