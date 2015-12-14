SUBROUTINE rsetup (lvl,lvls1,lvls2,nacum,idim)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     SETUP COMPUTES THE REVERSE LEVELING INFO FROM LVLS2 AND STORES
!     IT INTO LVLS2.  NACUM(I) IS INITIALIZED TO NODES/ITH LEVEL FOR
!     NODES ON THE PSEUDO-DIAMETER OF THE GRAPH.  LVL IS INITIALIZED TO
!     NON-ZERO FOR NODES ON THE PSEUDO-DIAM AND NODES IN A DIFFERENT
!     COMPONENT OF THE GRAPH.
 
 
 INTEGER, INTENT(OUT)                     :: lvl(1)
 INTEGER, INTENT(IN)                      :: lvls1(1)
 INTEGER, INTENT(OUT)                     :: lvls2(1)
 INTEGER, INTENT(OUT)                     :: nacum(1)
 INTEGER, INTENT(IN OUT)                  :: idim
 
 COMMON /bandb /  dum3b(3), ngrid
 COMMON /bandg /  n,        idpth
 
!     IDIM=NUMBER OF LEVELS IN A GIVEN COMPONENT.
!     NACUM IS DIMENSIONED TO IDIM IN SIZE
 
!     DIMENSION EXCEEDED  . . .  STOP JOB.
 
 IF (idpth <= idim)  GO TO 20
 ngrid=-3
 RETURN
 
 20 DO  i=1,idpth
   nacum(i)=0
 END DO
 DO  i=1,n
   lvl(i)=1
   lvls2(i)=idpth+1-lvls2(i)
   itemp=lvls2(i)
   IF (itemp > idpth) CYCLE
   IF (itemp /= lvls1(i)) GO TO 100
   nacum(itemp)=nacum(itemp)+1
   CYCLE
   100 lvl(i)=0
 END DO
 RETURN
END SUBROUTINE rsetup
