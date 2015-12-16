SUBROUTINE ampc (djh1,djh2,djh,ajjl,qjh,qjho,qjhua,scr1,scr2,scr3,  &
        scr4,scr5,scr6)
     
!     THE PURPOSE OF THIS ROUTINE IS TO COMPUTE (OR RETRIEVE QJH)
 
!     IF QJH MUST BE COMPUTED
 
!     1. FORM DJH FOR THIS K (IF IDJH.EQ.0)
!        DJH = DJH1 + I*K*DJH2
!     2. FOR EACH CONSTANT THEORY
!        A. RETRIEVE AJJ PORTION = AJJTH
!        B. PERFORM THEORY FOR QJH
!           1) DOUBLET LATTICE
!              A) DECOMPOSE AJJTH
!              B) FIND PROPER DJH PORTION DJHTH
!              C) FBS FOR QJHTH
!              D) ADD TO BOTTOM OF QJHUA(CYCLE)
!           6) COMPRESSOR BLADES  (IONCE = 1).
!              A) COMPUTE QJHTH = (AJJ)*DJH.
!              B) QJHUA = QJHTH SINCE ONLY ONE BLADE AND GROUP (NGP = 1)
!           7) SWEPT TURBOPROPS   (IONCE = 1).
!              A) COMPUTE QJHTH = (AJJ)*DJH.
!              B) QJHUA = QJHTH SINCE ONLY ONE BLADE AND GROUP (NGP = 1)
 
 
 INTEGER, INTENT(IN OUT)                  :: djh1
 INTEGER, INTENT(IN OUT)                  :: djh2
 INTEGER, INTENT(IN)                      :: djh
 INTEGER, INTENT(IN)                      :: ajjl
 INTEGER, INTENT(IN OUT)                  :: qjh
 INTEGER, INTENT(IN)                      :: qjho
 INTEGER, INTENT(IN)                      :: qjhua
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN OUT)                  :: scr2
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER, INTENT(IN)                      :: scr5
 INTEGER, INTENT(IN OUT)                  :: scr6
 INTEGER :: ajjcol,qhhcol, sysbuf,FILE,NAME(2),iblock(11),mcb(7),  &
      qjhth
 REAL :: BLOCK(11)
 COMMON /ampcom/ ncolj,nsub,xm,xk,ajjcol,qhhcol,ngp,ngpd(2,30),  &
     mcbqhh(7),mcbqjh(7),noh,idjh
 COMMON /system/ sysbuf,nout,skp(52),iprec
 COMMON /zzzzzz/ iz(1)
 COMMON /unpakx/ itc,ii,jj,incr
 COMMON /packx / itc1,itc2,ii1,jj1,incr1
 COMMON /cdcmpx/ dum32(32),ib
 EQUIVALENCE     (iblock(1),BLOCK(1))
 DATA    NAME  / 4HAMPC,4H    /
 DATA    iblock(1),iblock(7),BLOCK(2),BLOCK(3),BLOCK(8) /  &
     3,        3,     1.0,      0.,      0. /
 
!     INITIALIZE
 
 ibuf1 = korsz(iz) - sysbuf + 1
 ibuf2 = ibuf1 - sysbuf
 itc   = mcbqhh(5)
 incr  = 1
 itc1  = itc
 itc2  = itc1
 incr1 = incr
 ii1   = 1
 
!     IS QJH ON SAVE FILE
 
 IF (qhhcol == 0) GO TO 100
 
!     COPY QJH FROM OLD FILE TO QJH
 
 IF (mcbqjh(1) <= 0) GO TO 10
 CALL gopen (qjho,iz(ibuf1),0)
 CALL gopen (qjh, iz(ibuf2),3)
 k = qhhcol - 1
 IF (k == 0) GO TO 20
 FILE = qjho
 DO  i = 1,k
   CALL fwdrec (*910,qjho)
 END DO
 20 CONTINUE
 CALL cyct2b (qjho,qjh,noh,iz,mcbqjh)
 CALL CLOSE (qjho,1)
 CALL CLOSE (qjh,3)
 10 CONTINUE
 RETURN
 
!     COMPUTE QJH
 
 100 CONTINUE
 
!     HAS DJH ALREADY BEEN COMPUTED
 
 IF (idjh /= 0) GO TO 110
 BLOCK(9) = xk
 CALL ssg2c (djh1,djh2,djh,1,BLOCK)
 
!     POSITION AJJL
 
 110 CALL gopen (ajjl,iz(ibuf1),0)
 k = ajjcol - 1
 IF (k == 0) GO TO 120
 FILE = ajjl
 DO  i = 1,k
   CALL fwdrec (*910,ajjl)
 END DO
 120 CONTINUE
 CALL CLOSE (ajjl,2)
 
!     SET UP TO LOOP ON CONSTANT THEORY
 
 ngps = 1
 nth  = ngpd(1,ngps)
 ncolth = 0
 135 nclold = ncolth + 1
 140 IF (ngps > ngp) GO TO 150
 IF (ngpd(1,ngps) /= nth) GO TO 150
 ncolth = ncolth + ngpd(2,ngps)
 ngps   = ngps + 1
 GO TO 140
 
!     BRANCH ON THEORY
 
 150 CONTINUE
 ionce = 0
 IF (nclold == 1 .AND. ngps > ngp) ionce = 1
 
!     COPY AJJL TO SCR1
 
 CALL gopen (ajjl,iz(ibuf1),2)
 CALL gopen (scr1,iz(ibuf2),1)
 mcb(1) = ajjl
 CALL rdtrl (mcb)
 mcb(1) = scr1
 mcb(2) = 0
 mcb(3) = ncolth
 mcb(6) = 0
 mcb(7) = 0
 ii     = nc_lold
 jj     = ncolth
 ii1    = 1
 jj1    = ncolth - nc_lold + 1
 itc    = mcb(5)
 itc1   = itc
 itc2   = itc
 incr   = 1
 incr1  = 1
 CALL ampc1 (ajjl,scr1,ncolth,iz,mcb)
 CALL CLOSE (ajjl,2)
 CALL CLOSE (scr1,1)
 CALL wrttrl (mcb)
 SELECT CASE ( nth )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
   CASE (    4)
     GO TO 4000
   CASE (    5)
     GO TO 5000
   CASE (    6)
     GO TO 6000
   CASE (    7)
     GO TO 7000
 END SELECT
 
!     DOUBLET LATTICE WITH SLENDER BODIES
 
 1000 CONTINUE
 2000 CONTINUE
 
!     TRANSPOSE MATRIX
 
 CALL tranp1 (scr1,scr4,4,scr2,scr3,scr5,scr6,0,0,0,0)
 
!     DECOMPOSE MATRIX
 
 ib = 0
 CALL cfactr (scr4,scr2,scr3,scr1,scr5,scr6,iopt)
 
!     MACH BOX
!     PISTON
 
 
!     COMPRESSOR BLADE AND SWEPT TURBOPROP THEORIES -
!     ONE BLADE ALLOWED, ONE GROUP, USE WHOLE AJJ AND DJH MATRICES.
 
 3000 CONTINUE
 4000 CONTINUE
 5000 CONTINUE
 6000 CONTINUE
 7000 CONTINUE
 
!     COPY PROPER ROWS OF DJH TO SCR4
 
 idjha = djh
 IF (ionce /= 0) GO TO 1010
 ii  = nc_lold
 jj  = ncolth
 ii1 = 1
 jj1 = ncolth-nc_lold+1
 mcb(1) = djh
 CALL rdtrl (mcb)
 itc  = mcb(5)
 itc1 = itc
 itc2 = itc
 incr = 1
 incr1  = 1
 mcb(2) = 0
 mcb(3) = jj1
 mcb(6) = 0
 mcb(7) = 0
 mcb(1) = scr4
 CALL gopen (djh,iz(ibuf1),0)
 CALL gopen (scr4,iz(ibuf2),1)
 CALL ampc1 (djh,scr4,noh,iz,mcb)
 CALL CLOSE (djh,1)
 CALL CLOSE (scr4,1)
 CALL wrttrl (mcb)
 idjha = scr4
 1010 CONTINUE
 qjhth = scr5
 IF (ionce /= 0) qjhth = qjhua
 SELECT CASE ( nth )
   CASE (    1)
     GO TO 1001
   CASE (    2)
     GO TO 2001
   CASE (    3)
     GO TO 3001
   CASE (    4)
     GO TO 4001
   CASE (    5)
     GO TO 5001
   CASE (    6)
     GO TO 6001
   CASE (    7)
     GO TO 7001
 END SELECT
 1001 CONTINUE
 2001 CONTINUE
 
!     SOLVE FOR THIS PORTION OF QJH
 
 CALL cfbsor (scr2,scr3,idjha,qjhth,iopt)
 1020 CONTINUE
 
!     COPY ACCUMULATIVELY ONTO QJHUA
 
 IF (ionce /= 0) GO TO 8000
 CALL ampc2 (scr5,qjhua,scr1)
 IF (ngps > ngp) GO TO 8000
 GO TO 135
 
!     COMPUTE THIS PORTION OF QJH  = AJJ*DJH
 
 3001 CONTINUE
 4001 CONTINUE
 5001 CONTINUE
 6001 CONTINUE
 7001 CONTINUE
 CALL ssg2b (scr1,idjha,0,qjhth,0,iprec,1,scr6)
 GO TO 1020
 
!     ALL GROUPS / THEORIES COMPLETE
 
 8000 CONTINUE
 GO TO 10
 
!     ERROR MESSAGES
 
 901 CALL mesage (ip1,FILE,NAME)
 910 ip1 = -2
 GO TO 901
END SUBROUTINE ampc
