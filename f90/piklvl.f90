SUBROUTINE piklvl (*,lvls1,lvls2,ccstor,idflt,isdir,xc,nhigh,  &
        nlow,nacum,size,stpt)
     
 
 INTEGER, INTENT(IN)                      :: lvls1(1)
 INTEGER, INTENT(IN OUT)                  :: lvls2(1)
 INTEGER, INTENT(IN)                      :: ccstor(1)
 INTEGER, INTENT(IN)                      :: idflt
 INTEGER, INTENT(OUT)                     :: isdir
 INTEGER, INTENT(IN)                      :: xc
 INTEGER, INTENT(OUT)                     :: nhigh(1)
 INTEGER, INTENT(OUT)                     :: nlow(1)
 INTEGER, INTENT(IN OUT)                  :: nacum(1)
 INTEGER, INTENT(IN OUT)                  :: size(1)
 INTEGER, INTENT(IN OUT)                  :: stpt(1)
 INTEGER :: END, temp
 
 COMMON /bandg /  idum,      idpth
 
!     THIS ROUTINE IS USED ONLY BY GIBSTK OF BANDIT MODULE
 
!     PIKLVL CHOOSES THE LEVEL STRUCTURE  USED IN NUMBERING GRAPH
 
!     LVLS1-    ON INPUT CONTAINS FORWARD LEVELING INFO
!     LVLS2-    ON INPUT CONTAINS REVERSE LEVELING INFO
!               ON OUTPUT THE FINAL LEVEL STRUCTURE CHOSEN
!     CCSTOR-   ON INPUT CONTAINS CONNECTED COMPONENT INFO
!     IDFLT-    ON INPUT =1 IF WDTH LVLS1'WDTH LVLS2, =2 OTHERWISE
!     NHIGH     KEEPS TRACK OF LEVEL WIDTHS FOR HIGH NUMBERING
!               DIMENSION OF NHIGH IS MAXIMUM ALLOWABLE NUMBER OF LEVELS
!     NLOW-     KEEPS TRACK OF LEVEL WIDTHS FOR LOW NUMBERING
!     NACUM-    KEEPS TRACK OF LEVEL WIDTHS FOR CHOSEN LEVEL STRUCTURE
!     XC-       NUMBER OF MAXIMUM ALLOWABLE CONNECTED COMPONENTS
!               (IS THE DIMENSION FOR SIZE AND STPT)
!     SIZE(I)-  SIZE OF ITH CONNECTED COMPONENT
!     STPT(I)-  INDEX INTO CCSTORE OF 1ST NODE IN ITH CON COMPT
!     ISDIR-    FLAG WHICH INDICATES WHICH WAY THE LARGEST CONNECTED
!               COMPONENT FELL.  =+1 IF LOW AND -1 IF HIGH
 
 
!     PART 1 -
!     ========
!     SORTS SIZE AND STPT HERE, IN DECENDING ORDER
!     (PREVIOUS SORT2 ROUTINE IS NOW MOVED INTO HERE.
!     THE ORIGINAL BUBBLE SORT HAS BEEN REPLACED BY THE MODIFIED SHELL
!     SORT WHICH IS MUCH FASTER   /G.CHAN,  MAY 1988)
 
 IF (xc == 0) RETURN 1
 m=xc
 10 m=m/2
 IF (m == 0) GO TO 70
 j=1
 k=xc-m
 20 i=j
 30 n=i+m
 IF (size(n)-size(i) > 0.0) THEN
   GO TO    50
 ELSE
   GO TO    60
 END IF
 50 temp   =size(i)
 size(i)=size(n)
 size(n)=temp
 temp   =stpt(i)
 stpt(i)=stpt(n)
 stpt(n)=temp
 i=i-m
 IF (i >= 1) GO TO 30
 60 j=j+1
 IF (j-k > 0) THEN
   GO TO    10
 ELSE
   GO TO    20
 END IF
 70 CONTINUE
 
 
!     PART 2 -
!     ========
!     CHOOSES THE LEVEL STRUCTURE USED IN NUMBERING GRAPH
 
 
!     FOR EACH CONNECTED COMPONENT DO
 
 DO  i=1,xc
   j  =stpt(i)
   END=size(i)+j-1
   
!     SET NHIGH AND NLOW EQUAL TO NACUM
   
   DO  k=1,idpth
     nhigh(k)=nacum(k)
     nlow(k) =nacum(k)
   END DO
   
!     UPDATE NHIGH AND NLOW FOR EACH NODE IN CONNECTED COMPONENT
   
 DO  k=j,END
 inode=ccstor(k)
 lvlnh=lvls1(inode)
 nhigh(lvlnh)=nhigh(lvlnh)+1
 lvlnl=lvls2(inode)
 nlow(lvlnl)=nlow(lvlnl)+1
END DO
MAX1=0
max2=0

!     SET MAX1=LARGEST NEW NUMBER IN NHIGH
!     SET MAX2=LARGEST NEW NUMBER IN NLOW

DO  k=1,idpth
  IF (2*nacum(k) == nlow(k)+nhigh(k)) CYCLE
  IF (nhigh(k) > MAX1) MAX1=nhigh(k)
  IF (nlow(k) > max2) max2=nlow(k)
END DO

!     SET IT= NUMBER OF LEVEL STRUCTURE TO BE USED

it=1
IF (MAX1 > max2) it=2
IF (MAX1 == max2) it=idflt
IF (it == 2) GO TO 250
IF (i == 1) isdir=-1

!     COPY LVLS1 INTO LVLS2 FOR EACH NODE IN CONNECTED COMPONENT

DO  k=j,END
inode=ccstor(k)
lvls2(inode)=lvls1(inode)
END DO

!     UPDATE NACUM TO BE THE SAME AS NHIGH

DO  k=1,idpth
  nacum(k)=nhigh(k)
END DO
CYCLE

!     UPDATE NACUM TO BE THE SAME AS NLOW

250 DO  k=1,idpth
  nacum(k)=nlow(k)
END DO
END DO
RETURN
END SUBROUTINE piklvl
