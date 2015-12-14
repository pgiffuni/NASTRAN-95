SUBROUTINE NUMBER (snd,num,ndstk,lvls2,ndeg,renum,lvlst,lstpt,  &
        nflg,ibw2,ipf2,ipfa,isdir,stka,stkb,stkc,stkd,nu,idim)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     NUMBER PRODUCES THE NUMBERING OF THE GRAPH FOR MIN BANDWIDTH
 
!     SND-      ON INPUT THE NODE TO BEGIN NUMBERING ON
!     NUM-      ON INPUT AND OUTPUT, THE NEXT AVAILABLE NUMBER
!     LVLS2-    THE LEVEL STRUCTURE TO BE USED IN NUMBERING
!     RENUM-    THE ARRAY USED TO STORE THE NEW NUMBERING
!     LVLST-    ON OUTPUT CONTAINS LEVEL STRUCTURE
!     LSTPT(I)- ON OUTPUT, INDEX INTO LVLST TO FIRST NODE IN ITH LVL
!               LSTPT(I+1) - LSTPT(I) = NUMBER OF NODES IN ITH LVL
!     NFLG-     =+1 IF SND IS FORWARD END OF PSEUDO-DIAM
!               =-1 IF SND IS REVERSE END OF PSEUDO-DIAM
!     IBW2-     BANDWIDTH OF NEW NUMBERING COMPUTED BY NUMBER
!     IPF2-     PROFILE OF NEW NUMBERING COMPUTED BY NUMBER
!     IBW2 AND IPF2 HERE DO NOT INCLUDE DIAGONAL TERMS.
!     IPFA-     WORKING STORAGE USED TO COMPUTE PROFILE AND BANDWIDTH
!     ISDIR-    INDICATES STEP DIRECTION USED IN NUMBERING(+1 OR -1)
!     STACKS HAVE DIMENSION OF IDIM
!     NU-       WORK SPACE FOR BUNPAK
 
 
 INTEGER, INTENT(IN OUT)                  :: snd
 INTEGER, INTENT(IN OUT)                  :: num
 INTEGER, INTENT(IN OUT)                  :: ndstk(1)
 INTEGER, INTENT(IN)                      :: lvls2(1)
 INTEGER, INTENT(IN)                      :: ndeg(1)
 INTEGER, INTENT(OUT)                     :: renum(1)
 INTEGER, INTENT(OUT)                     :: lvlst(1)
 INTEGER, INTENT(OUT)                     :: lstpt(1)
 INTEGER, INTENT(IN)                      :: nflg
 INTEGER, INTENT(OUT)                     :: ibw2
 INTEGER, INTENT(OUT)                     :: ipf2
 INTEGER, INTENT(OUT)                     :: ipfa(1)
 INTEGER, INTENT(IN)                      :: isdir
 INTEGER, INTENT(OUT)                     :: stka(1)
 INTEGER, INTENT(OUT)                     :: stkb(1)
 INTEGER, INTENT(OUT)                     :: stkc(1)
 INTEGER, INTENT(OUT)                     :: stkd(1)
 INTEGER, INTENT(IN)                      :: nu(1)
 INTEGER, INTENT(IN OUT)                  :: idim
 INTEGER :: xa,       xb,       xc,       xd,       END, cx, test
 
 COMMON /bandb /  dum3(3),  ngrid
 COMMON /bandg /  n,        idpth,    ideg
 COMMON /bands /  dums(4),  maxgrd,   maxdeg
 COMMON /system/  ibuf,     nout
 
!     SET UP LVLST AND LSTPT FROM LVLS2
 
 DO  i=1,n
   ipfa(i)=0
 END DO
 nstpt=1
 DO  i=1,idpth
   lstpt(i)=nstpt
   DO  j=1,n
     IF (lvls2(j) /= i) CYCLE
     lvlst(nstpt)=j
     nstpt=nstpt+1
   END DO
 END DO
 lstpt(idpth+1)=nstpt
 
!     THIS ROUTINE USES FOUR STACKS, A,B,C,AND D, WITH POINTERS
!     XA,XB,XC, AND XD.  CX IS A SPECIAL POINTER INTO STKC WHICH
!     INDICATES THE PARTICULAR NODE BEING PROCESSED.
!     LVLN KEEPS TRACK OF THE LEVEL WE ARE WORKING AT.
!     INITIALLY STKC CONTAINS ONLY THE INITIAL NODE, SND.
 
 lvln=0
 IF (nflg < 0) lvln=idpth+1
 xc=1
 stkc(xc)=snd
 20 cx=1
 xd=0
 lvln=lvln+nflg
 lst=lstpt(lvln)
 lnd=lstpt(lvln+1)-1
 
!     BEGIN PROCESSING NODE STKC(CX)
 
 25 ipro=stkc(cx)
 renum(ipro)=num
 num=num+isdir
 END=ndeg(ipro)
 xa=0
 xb=0
 
!     CHECK ALL ADJACENT NODES
 
 CALL bunpak(ndstk,ipro,END,nu)
DO  i=1,END
test =nu(i)
inx=renum(test)

!     ONLY NODES NOT NUMBERED OR ALREADY ON A STACK ARE ADDED

IF (inx == 0) GO TO 30
IF (inx < 0) CYCLE

!     DO PRELIMINARY BANDWIDTH AND PROFILE CALCULATIONS

nbw=(renum(ipro)-inx)*isdir
IF (isdir > 0) inx=renum(ipro)
IF (ipfa(inx) < nbw) ipfa(inx)=nbw
CYCLE
30 renum(test)=-1

!     PUT NODES ON SAME LEVEL ON STKA, ALL OTHERS ON STKB

IF (lvls2(test) == lvls2(ipro)) GO TO 35
xb=xb+1
IF (xb > idim) GO TO 100
stkb(xb)=test
CYCLE
35 xa=xa+1
IF (xa > idim) GO TO 100
stka(xa)=test
END DO

!     SORT STKA AND STKB INTO INCREASING DEGREE AND ADD STKA TO STKC
!     AND STKB TO STKD

IF (xa == 0) GO TO 50
IF (xa == 1) GO TO 45
CALL sortdg (stkc,stka,xc,xa,ndeg)
GO TO 50
45 xc=xc+1
IF (xc > idim) GO TO 100
stkc(xc)=stka(xa)
50 IF (xb == 0) GO TO 65
IF (xb == 1) GO TO 60
CALL sortdg (stkd,stkb,xd,xb,ndeg)
GO TO 65
60 xd=xd+1
IF (xd > idim) GO TO 100
stkd(xd)=stkb(xb)

!     BE SURE TO PROCESS ALL NODES IN STKC

65 cx=cx+1
IF (xc >= cx) GO TO 25

!     WHEN STKC IS EXHAUSTED LOOK FOR MIN DEGREE NODE IN SAME LEVEL
!     WHICH HAS NOT BEEN PROCESSED

MAX=ideg+1
snd=n+1
DO  i=lst,lnd
  test=lvlst(i)
  IF (renum(test) /= 0) CYCLE
  IF (ndeg(test) >= MAX) CYCLE
  renum(snd)=0
  renum(test)=-1
  MAX=ndeg(test)
  snd=test
END DO
IF (snd == n+1) GO TO 75
xc=xc+1
IF (xc > idim) GO TO 100
stkc(xc)=snd
GO TO 25

!     IF STKD IS EMPTY WE ARE DONE, OTHERWISE COPY STKD ONTO STKC
!     AND BEGIN PROCESSING NEW STKC

75 IF (xd == 0) GO TO 90
DO  i=1,xd
  stkc(i)=stkd(i)
END DO
xc=xd
GO TO 20

!     DO FINAL BANDWIDTH AND PROFILE CALCULATIONS

90 DO  i=1,n
  IF (ipfa(i) > ibw2) ibw2=ipfa(i)
  ipf2=ipf2+ipfa(i)
END DO
RETURN

!     DIMENSION EXCEEDED  . . .  STOP JOB.

100 ngrid=-3
RETURN
END SUBROUTINE NUMBER
