SUBROUTINE gibstk (ndstk,iold,renum,ndeg,lvl,lvls1,lvls2,ccstor,  &
        jump,icrit,nhigh,nlow,nacum,size,stpt,un,idim)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     GIBBSTOCK USES GRAPH THEORETICAL METHODS TO PRODUCE A PERMUTATION
!     OF AN INPUT ARRAY WHICH REDUCES ITS BANDWITH
 
!     THE FOLLOWING INPUT PARAMETERS ARE REQUIRED--NDSTK,N,IDEG,IOLD
 
!     THESE INTEGER ARRAYS MUST BE DIMENSIONED IN THE CALLING PROGRAM--
!     NDSTK(NR,D1),RENUM(D2+1),NDEG(D2),IOLD(D2),LVL(D2),LVLS1(D2),
!     LVLS2(D2),CCSTOR(D2)   WHERE D1 .GE. MAX DEGREE OF ANY NODE AND
!     D2 AND NR ARE .GE. THE TOTAL NUMBER OF NODES IN THE GRAPH.
 
!     EXPLANATION OF PARAMETERS--
!     NDSTK   - ADJACENCY ARRAY REPRESENTING GRAPH TO BE PROCESSED
!               NDSTK(I,J) = NODE NUMBER OF JTH CONNECTION TO NODE
!               NUMBER I.  A CONNECTION OF A NODE TO ITSELF IS NOT
!               LISTED.  EXTRA POSITIONS MUST HAVE ZERO FILL.
!     NR      - ROW DIMENSION ASSIGNED NDSTK IN CALLING PROGRAM = II1
!     IOLD(I) - RENUMBERING OF ITH NODE BEFORE GIBBSTOCK PROCESSING
!               IF NO RENUMBERING EXISTS THEN ILD(1)=1,ILD(2)=2, ETC.
!     N       - NUMBER OF NODES IN GRAPH BEING PROCESSED
!     IDEG    - MAX DEGREE OF ANY NODE IN GRAPH BEING PROCESSED
!     JUMP   IS SET TO 0 IF EITHER CRITERION IS REDUCED.
!     ICRIT   - RESEQUENCING CRITERION, SET BY BANDIT
!               1 RMS WAVEFRONT, 2 BANDWIDTH, 3 PROFILE, 4 MAX.WAVEFRONT
 
!     ON OUTPUT THESE VARIABLES CONTAIN THE FOLLOWING INFORMATION--
!     RENUM(I)- THE NEW NUMBER FOR THE ITH NODE
!     NDEG(I) - THE DEGREE OF THE ITH NODE
!     IDPTH   - NUMBER OF LEVELS IN GIBBSTOCK LEVEL STRUCTURE
!     IBW2    - THE BANDWITH AFTER RENUMBERING
!     IPF2    - THE PROFILE AFTER RENUMBERING
 
!     THE FOLLOWING ONLY HAVE MEANING IF THE GRAPH WAS ALL ONE COMPONENT
!     LVL(I)  - INDEX INTO LVLS1 TO THE FIRST NODE IN LEVEL I
!               LVL(I+1)-LVL(I)= NUMBER OF NODES IN ITH LEVEL
!     LVLS1   - LEVEL STRUCTURE CHOSEN BY GIBBSTOCK
!     LVLS2(I)- THE LEVEL ASSIGNED TO NODE I BY GIBBSTOCK
 
 
 INTEGER, INTENT(IN OUT)                  :: ndstk(1)
 INTEGER, INTENT(IN)                      :: iold(1)
 INTEGER, INTENT(OUT)                     :: renum(1)
 INTEGER, INTENT(IN)                      :: ndeg(1)
 INTEGER, INTENT(IN OUT)                  :: lvl(1)
 INTEGER, INTENT(IN OUT)                  :: lvls1(1)
 INTEGER, INTENT(IN OUT)                  :: lvls2(1)
 INTEGER, INTENT(IN OUT)                  :: ccstor(1)
 INTEGER, INTENT(OUT)                     :: jump
 INTEGER, INTENT(IN OUT)                  :: icrit
 INTEGER, INTENT(IN OUT)                  :: nhigh(1)
 INTEGER, INTENT(IN OUT)                  :: nlow(1)
 INTEGER, INTENT(IN OUT)                  :: nacum(1)
 INTEGER, INTENT(OUT)                     :: size(1)
 INTEGER, INTENT(OUT)                     :: stpt(1)
 REAL, INTENT(IN OUT)                     :: un(1)
 INTEGER, INTENT(IN)                      :: idim
 INTEGER :: stnode,   rvnode, xc,       sumwb, stnum, sbnum,  &
     obw,      op,       xcmax
 REAL :: im1,      im2
 
 COMMON /banda /  dum5a(5), method
 COMMON /bandb /  dum3b(3), ngrid
 COMMON /bandd /  obw,      nbw,      op,       np,       ncm, nzero
 COMMON /bandg /  n,        idpth,    ideg
 COMMON /bandw /  maxw0,    rms0,     maxw1,    rms1,     i77, brms0,    brms1
 COMMON /bands /  nn,       mm
 COMMON /system/  ibuf,     nout,     dum6s(6), nlpp
 
!     OLD AND NEW MAX AND RMS WAVEFRONT FOR ENTIRE PROBLEM,
!     NOT JUST GIBSTK.
!     DIMENSIONS OF NHIGH, NLOW, AND NACUM ARE IDIM EACH
!     SIZE AND STPT HAVE DIMENSION IDIM/2 AND SHOULD BE CONTIGUOUS IN
!     CORE WITH SIZE FIRST.
!     XC = NUMBER OF SUB-COMPONENTS RESULTING AFTER REMOVING DIAMETER
!     FROM ONE COMPONENT OF ORIGINAL GRAPH.
 
 xcmax = idim/2
 ncm   = 0
 n     = nn
 ibw2  = 0
 ipf2  = 0
 
!     SET RENUM(I) = 0 FOR ALL I TO INDICATE NODE I IS UNNUMBERED
!     THEN COMPUTE DEGREE OF EACH NODE AND ORIGINAL B AND P.
 
 DO  i = 1,n
   renum(i) = 0
 END DO
 CALL dgree (ndstk,ndeg,iold,ibw1,ipf1,un)
 
!     ORIGINAL ACTIVE COLUMN DATA IN MAXW1 AND RMS1, COMPUTED BY SCHEME
 
 IF (method /= 0) GO TO 35
 maxwa = maxw1
 rmsa  = rms1
 brmsa = brms1
 GO TO 38
 35 maxwa = maxw0
 rmsa  = rms0
 brmsa = brms0
 38 CONTINUE
 
!     NUMBER THE NODES OF DEGREE ZERO
!     SBNUM = LOW  END OF AVAILABLE NUMBERS FOR RENUMBERING
!     STNUM = HIGH END OF AVAILABLE NUMBERS FOR RENUMBERING
 
 sbnum = 1
 stnum = n
 DO  i = 1,n
   IF (ndeg(i) > 0) CYCLE
   renum(i) = stnum
   stnum = stnum-1
 END DO
 
!     NODES OF ZERO DEGREE APPEAR LAST IN NEW SEQUENCE.
 
 nzero = n - stnum
 ncm   = nzero
 
!     FIND AN UNNUMBERED NODE OF MIN DEGREE TO START ON
 
 50 lowdg = ideg + 1
 ncm   = ncm + 1
 nflg  = 1
 isdir = 1
 DO  i = 1,n
   IF (ndeg(i) >= lowdg .OR. renum(i) > 0) CYCLE
   lowdg  = ndeg(i)
   stnode = i
 END DO
 
!     FIND PSEUDO-DIAMETER AND ASSOCIATED LEVEL STRUCTURES.
!     STNODE AND RVNODE ARE THE ENDS OF THE DIAM AND LVLS1 AND LVLS2
!     ARE THE RESPECTIVE LEVEL STRUCTURES.
 
 CALL fndiam (stnode,rvnode,ndstk,ndeg,lvl,lvls1,lvls2,ccstor,  &
     idflt,size,un,idim)
 IF (ngrid == -3) RETURN
 IF (ndeg(stnode) <= ndeg(rvnode)) GO TO 75
 
!     NFLG INDICATES THE END TO BEGIN NUMBERING ON
 
 nflg   =-1
 stnode = rvnode
 75 CALL rsetup (lvl,lvls1,lvls2,           nacum,idim)
!                                  NHIGH,NLOW,    <===== NEW
 IF (ngrid == -3) RETURN
 
!     FIND ALL THE CONNECTED COMPONENTS  (XC COUNTS THEM)
 
 xc    = 0
 lroot = 1
 lvln  = 1
 DO  i = 1,n
   IF (lvl(i) /= 0) CYCLE
   xc = xc + 1
   IF (xc <= xcmax) GO TO 80
   
!     DIMENSION EXCEEDED.  STOP JOB.
   
   ngrid =-3
   RETURN
   
   80 stpt(xc) = lroot
   CALL tree (i,ndstk,lvl,ccstor,ndeg,lvlwth,lvlbot,lvln,maxlw,n,un)
   size(xc) = lvlbot + lvlwth - lroot
   lroot = lvlbot + lvlwth
   lvln  = lroot
 END DO
 CALL piklvl (*90,lvls1,lvls2,ccstor,idflt,isdir,xc,nhigh,nlow,  &
     nacum,size,stpt)
 
!     ON RETURN FROM PIKLVL, ISDIR INDICATES THE DIRECTION THE LARGEST
!     COMPONENT FELL.  ISDIR IS MODIFIED NOW TO INDICATE THE NUMBERING
!     DIRECTION.  NUM IS SET TO THE PROPER VALUE FOR THIS DIRECTION.
 
 90 isdir = isdir*nflg
 num   = sbnum
 IF (isdir < 0) num = stnum
 
 CALL NUMBER (stnode,num,ndstk,lvls2,ndeg,renum,lvls1,lvl,nflg,  &
     ibw2,ipf2,ccstor,isdir,nhigh,nlow,nacum,size,un,idim)
 IF (ngrid == -3) RETURN
 
!     UPDATE STNUM OR SBNUM AFTER NUMBERING
 
 IF (isdir < 0) stnum = num
 IF (isdir > 0) sbnum = num
 IF (sbnum <= stnum) GO TO 50
 
!     COMPUTE THE NEW BANDWIDTH, PROFILE, AND WAVEFRONT.
 
 CALL wavey (ndstk,renum,lvl,0,lvls2,lvls1,maxb,maxwb,averwb,  &
     sumwb,rmsb,brmsb,un)
 
 ibw2 = maxb
 ipf2 = sumwb
 IF (nlpp > 50) WRITE (nout,100) maxb,sumwb,maxwb,averwb, rmsb,brmsb
 100 FORMAT (/31X,66HAFTER resequencing by gibbs-poole-stockmeyer (gps)  &
     algorithm - - -, /40X,13HBANDWIDTH    ,i9,  /40X,13HPROFILE      ,i9,  &
     /40X,13HMAX wavefront,i9,  /40X,13HAVG wavefront,f9.3,  &
     /40X,13HRMS wavefront,f9.3,/40X,13HRMS bandwidth,f9.3)
 
!     CHECK NEW NUMBERING AGAINST OLD NUMBERING.
 
 SELECT CASE ( icrit )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 130
   CASE (    4)
     GO TO 140
 END SELECT
 110 im1   = rmsa
 im2   = ipf1
 crit1 = rmsb
 crit2 = ipf2
 GO TO 150
 120 im1   = ibw1
 im2   = ipf1
 crit1 = ibw2
 crit2 = ipf2
 GO TO 150
 130 im1   = ipf1
 im2   = ibw1
 crit1 = ipf2
 crit2 = ibw2
 GO TO 150
 140 im1   = maxwa
 im2   = rmsa
 crit1 = maxwb
 crit2 = rmsb
 
 150 IF (crit1-im1 < 0.0) THEN
   GO TO   210
 ELSE IF (crit1-im1 == 0.0) THEN
   GO TO   160
 ELSE
   GO TO   170
 END IF
 160 IF (crit2 < im2) GO TO 210
 
!     IF ORIGINAL NUMBERING IS BETTER THAN NEW ONE, SET UP TO RETURN IT
 
 170 DO  i = 1,n
   renum(i) = iold(i)
 END DO
 ibw2  = ibw1
 ipf2  = ipf1
 maxwb = maxwa
 rmsb  = rmsa
 brmsb = brmsa
 GO TO 220
 
!     EQUATE CORRESPONDING GPS AND BANDIT VARIABLES.
 
 210 jump  = 0
 220 nbw   = ibw2
 np    = ipf2
 maxw1 = maxwb
 rms1  = rmsb
 brms1 = brmsb
 RETURN
END SUBROUTINE gibstk
