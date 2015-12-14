SUBROUTINE read2 (maa,phia,scr1,norm,ia,uset,mi,lama,ipout,scr2,  &
        epsi,scr3)
     
!     COMPUTE MODAL MASS AND NORMALIZES VECTORS ACCORDING TO POINT,
!     MASS, OR MAX.  ALSO LOOKS FOR LARGE OFF DIAGONAL TERM
 
 
 INTEGER, INTENT(IN)                      :: maa
 INTEGER, INTENT(IN)                      :: phia
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(OUT)                     :: norm
 INTEGER, INTENT(IN)                      :: ia
 REAL, INTENT(IN OUT)                     :: uset
 INTEGER, INTENT(IN)                      :: mi
 INTEGER, INTENT(IN)                      :: lama
 INTEGER, INTENT(IN)                      :: ipout
 INTEGER, INTENT(IN OUT)                  :: scr2
 REAL, INTENT(OUT)                        :: epsi
 INTEGER, INTENT(IN OUT)                  :: scr3
 INTEGER :: point,sysbuf, ix(7),iphia(7), ihead(50), sturm,nam(2)
 REAL :: lfreq,core(13)
 DOUBLE PRECISION :: dcore(1),dxmax
 DIMENSION       im(7),ihead1(10)
 COMMON /condas/ consts(5)
 COMMON /zzzzzz/ icore(1)
 COMMON /system/ sysbuf
 COMMON /packx / ita1,itb1,ii1,jj1,incur1
 COMMON /unpakx/ itb,ii,jj,incur
 COMMON /output/ head(1)
 COMMON /sturmx/ sturm,shftpt,KEEP,ptshft,nr
 COMMON /givn  / givens,title1(100),lfreq,title2(4),nnv
 EQUIVALENCE     (consts(2),tphi), (ix(2),ncol), (ix(3),nrow),  &
     (core(1),icore(1),dcore(1)), (dxmax,xmax)
 DATA    ihead1/ 21,9,8*0  /
 DATA    ihead / 21,6,7*0,7,40*0/
 DATA    mass,   point     / 4HMASS,4HPOIN/
 DATA    MAX   / 4HMAX     /
 DATA    nam   / 4HREAD,1H2/
 
!     READ2  SHOULD NORMALIZE  PHIA  ACCORDING TO NORM +METHOD
 
 lcore = korsz(core)
 
!     DECIDE IF MI WANTED
 
 imi   = 0
 ix(1) = mi
 CALL rdtrl (ix)
 IF (ix(1) > 0) GO TO 10
 epsi = 0.0
 imi  = -1
 IF (norm == mass) norm = MAX
 10 ix(1) = phia
 CALL rdtrl (ix)
 CALL makmcb (iphia,phia,ix(3),ix(4),ix(5))
 
!     SET UP TO HANDLE IDENTITY MATRIX
 
 iden  = 0
 im(1) = maa
 CALL rdtrl (im)
 IF (im(4) == 8) iden = 1
 
!     FIND TYPE OF NORMALIZATION
 
 IF (norm == mass) GO TO 310
 ipont = 1
 IF (norm  == point) GO TO 30
 IF (ia < 1 .OR. ia > nrow) GO TO 20
 
!     TYPE IS  MAX
 
 20 ipont = 0
 
!     POINT
 
 30 ASSIGN 40 TO icopy
 GO TO 420
 
 40 CONTINUE
 
!     PROCESS PHIA - NORMALIZE - COPY TO PHIA
 
 lcore = lcore - sysbuf
 CALL gopen (scr1,core(lcore+1),0)
 lcore = lcore - sysbuf
 CALL gopen (phia,core(lcore+1),1)
 itb   = ix(5)
 jj    = nrow
 ii    = 1
 incur = 1
 ita1  = itb
 itb1  = itb
 incur1= 1
 DO  i = 1,ncol
   CALL unpack (*100,scr1,core(3))
   ii1 = ii
   jj1 = jj
   jjj = 1
   IF (itb == 2) GO TO 66
   DO  j = 1,nrow
     IF (ABS(core(j+2)) > ABS(core(jjj+2))) jjj = j
   END DO
   jjj = jjj + 2
   IF (ipont /= 1) GO TO 62
   jjj = ia + 2
   IF (ABS(core(jjj)) <= 1.0E-15) GO TO 90
   62 xmax = core(jjj)
   DO  j = 1,nrow
     core(j+2) = core(j+2)/xmax
   END DO
   GO TO 90
   66 DO  j = 1,nrow
     IF (DABS(dcore(j+1)) > DABS(dcore(jjj+1))) jjj = j
   END DO
   jjj = jjj + 1
   IF (ipont /= 1) GO TO 70
   jjj = ia + 1
   IF (DABS(dcore(jjj)) <= 1.0D-15) GO TO 90
   70 dxmax = dcore(jjj)
   DO  j = 1,nrow
     dcore(j+1) = dcore(j+1)/dxmax
   END DO
   90 CALL pack (core(3),phia,iphia)
   CYCLE
   100 ii1 = 1
   jj1 = 1
   CALL pack (core,phia,iphia)
 END DO
 CALL CLOSE (phia,1)
 CALL CLOSE (scr1,1)
 
!     COMPUTE MODAL MASS
 
 140 IF (imi  < 0) GO TO 170
 IF (iden == 0) GO TO 160
 ASSIGN 150 TO icopy
 GO TO 420
 150 CALL ssg2b (phia,scr1,0,mi,1,itb,1,scr3)
 GO TO 170
 
 160 CALL ssg2b (maa,phia,0,scr2,0,itb,1,scr3)
 CALL ssg2b (phia,scr2,0,mi,1,itb,1,scr3)
 
!     COMPUTE GENERALIZED STIFFNESS
 
 
!     COMPUTE FREQUENCY ETC
 
 170 itb  = 1
 ii   = 1
 jj   = ncol
 incur= 1
 imsg = 0
 CALL gopen (lama,core(lcore+1),0)
 CALL READ (*500,*172,lama,core(1),lcore,1,nlama)
 GO TO 520
 
!     NLAMA IS THE NUMBER OF EIGENVALUES FOUND   NCOL IS TH NUMBER OF
!     VECTORS
 
 
!     BRING IN THE ORDER FOUND
 
 172 kk = nlama + 2*ncol + 8
 
!     KK IS THE POINTER TO THE ORDER FOUND
!     L1 AND  L2 ARE COUNTERS FOR MISSING LOW FREQ. BELOW SHIFT POINTS
!     STURM AND KEEP WERE SAVED IN SDCOMP, SHFTPT AND PTSHFT IN FEER
!     AND INVPWR (REAL SYMMETRIC EIGENVALUE PROBLEM ONLY)
 
 CALL READ (*500,*171,lama,icore(kk+1),lcore,1,iflag)
 GO TO 520
 171 CALL CLOSE (lama,1)
 CALL gopen (lama,core(lcore+1),1)
 CALL WRITE (lama,ihead(1),50,0)
 CALL WRITE (lama,head(1),96,1)
 lcore = lcore + sysbuf
 core(nlama+6) = 0.0
 core(nlama+7) = 0.0
 IF (imi < 0) GO TO 180
 CALL gopen (mi,core(lcore+1),0)
 l1 = sturm
 l2 = KEEP
 shftpt = shftpt + 1.e-10
 ptshft = ptshft + 1.e-10
 180 DO  i = 1,nlama
   icore(nlama+1) = i
   l = kk + i
   icore(nlama+2) = icore(l)
   core(nlama+3)  = core(i)
   core(nlama+4)  = SQRT(ABS(core(i)))
   core(nlama+5)  = core(nlama+4)/tphi
   IF (core(i) > 1.e-10 .AND. core(i) <= shftpt) l1 = l1 - 1
   IF (core(i) > 1.e-10 .AND. core(i) <= ptshft) l2 = l2 - 1
   IF (imi <  0) GO TO 200
   IF (i > ncol) GO TO 195
   l = nlama + i + 7
   k = l - 1 + i
   CALL unpack (*195,mi,core(l))
   core(nlama+6) = core(k)
   core(nlama+7) = core(k)*core(nlama+3)
   core(l) = core(k)
   
!     ZERO OUT GENERALIZED MASS AND GENERALIZED STIFFNESS FOR THE RIGID
!     BODY MODE OF ZERO FREQUENCY
   
!     (G.C.  3/92
!     NEXT 4 NEW LINES CAUSED DEMO T03121A TO DIE. MORE STUDY IS NEEDED)
   
!     IF (CORE(I) .GE. 0.0) GO TO 200
!     CORE(NLAMA+3) = 0.0
!     CORE(NLAMA+4) = 0.0
!     CORE(NLAMA+5) = 0.0
   GO TO 200
   
!     NO MORE VECTORS
!     REPLACE STURM BY SMALLER OF L1 OR L2, IF NOT ALL LOWER MODES FOUND
!     SET STRUM TO   -1 IF THERE IS NOT ENOUGH INFORMATION,
!     SET STRUM TO -999 IF DIAG 37 IS REQUESTED (NOT TO PRINT MESSAGE).
   
   195 core(nlama+6) = 0.0
   core(nlama+7) = 0.0
   200 CALL WRITE (lama,core(nlama+1),7,0)
 END DO
 IF (l1 <  0) l1 = 0
 IF (l2 <  0) l2 = 0
 IF (l1 > l2) l1 = l2
 IF (sturm /= -1 .AND. l1 >= 0) sturm = l1
 IF (sturm > nr .AND. nr > 0) sturm = sturm - nr
 IF (KEEP <= 0 .AND. ptshft > 0.) sturm = -1
 CALL sswtch (37,j)
 IF (j == 1) sturm = -999
 CALL CLOSE (lama,1)
 IF (imi < 0) GO TO 220
 CALL CLOSE (mi,1)
 220 imsg  = 0
 xmax  = 0.
 xmax1 = 0.
 istor = 0
 jstor = 0
 
!     EPSI = 0 IMPLIES TO NOT CHECK MODAL MASS TERMS
 
 IF (epsi == 0.0) GO TO 270
 CALL gopen (mi,core(lcore+1),0)
 loop260:  DO  i = 1,ncol
   m    = nlama + i + 7
   mcol = m + ncol
   CALL unpack (*540,mi,core(mcol))
   IF (core(m) == 0) CYCLE loop260
   DO  j = 1,ncol
     IF (i == j) CYCLE loop260
     k  = mcol  + j - 1
     mm = nlama + j + 7
     IF (core(mm) == 0.0) CYCLE
     gm = ABS(core(k))/SQRT(ABS(core(m)*core(mm)))
     IF (gm > xmax1) GO TO 240
     230 CONTINUE
     IF (gm <= epsi) CYCLE
     imsg = imsg + 1
     xmax = AMAX1(xmax,gm)
     CYCLE
     240 xmax1 = gm
     istor = i
     jstor = j
     GO TO 230
   END DO
 END DO loop260
 
 CALL CLOSE (mi,1)
 IF (imsg   /=  0) CALL mesage (34,xmax,epsi)
 270 IF (givens == .0) GO TO 275
 IF (nnv    /=  0) GO TO 275
 IF (lfreq  > .0) GO TO 600
 275 CALL gopen (ipout,core(lcore+1),0)
 CALL READ (*510,*280,ipout,core(1),lcore,1,iflag)
 GO TO 520
 280 CALL CLOSE (ipout,1)
 CALL gopen (ipout,core(lcore+1),1)
 ihead1(3) = icore(1)
 CALL WRITE (ipout,ihead1,10,0)
 i0 = 0
 core (i0+ 9) = xmax1
 icore(i0+10) = istor
 icore(i0+11) = jstor
 icore(i0+12) = imsg
 icore(i0+13) = sturm
 CALL WRITE (ipout,core(2),40,0)
 CALL WRITE (ipout,head,96,1)
 IF (icore(1) /= 1) GO TO 290
 iflag = iflag - 12
 ihead1( 3) = 3
 ihead1(10) = 6
 CALL WRITE (ipout,ihead1,50,0)
 CALL WRITE (ipout,head,96,1)
 IF (iflag == 0) GO TO 290
 CALL WRITE (ipout,core(13),iflag,0)
 290 CALL CLOSE (ipout,1)
 ix(1) = ipout
 CALL wrttrl (ix)
 RETURN
 
!     COMPUTE UNNORMALIZED MODAL MASS
 
 310 ASSIGN 320 TO icopy
 GO TO 420
 320 IF (iden == 0) GO TO 330
 
!     MASS MATRIX IS IDENTITY
 
 CALL ssg2b (phia,scr1,0,mi,1,iphia(5),1,scr3)
 GO TO 340
 
 330 CALL ssg2b (maa,phia,0,scr2,0,iphia(5),1,scr3)
 CALL ssg2b (phia,scr2,0,mi,1,iphia(5),1,scr3)
 
!     BRING IN DIAGONALS
 
 340 lcore = lcore - sysbuf
 CALL gopen (mi,core(lcore+1),0)
 itb = iphia(5)
 ii  = 1
 jj  = ncol
 IF (itb /= 2) GO TO 356
 DO  j = 1,ncol
   CALL unpack (*348,mi,dcore(ncol+1))
   k = ncol + j
   dcore(j) = 1.0D0/DSQRT(DABS(dcore(k)))
   CYCLE
   348 dcore(j) = 0.0D0
 END DO
 GO TO 362
 356 DO  j = 1,ncol
   CALL unpack (*358,mi,core(ncol+1))
   k = ncol + j
   core(j) = 1.0/SQRT(ABS(core(k)))
   CYCLE
   358 core(j) = 0.0
 END DO
 362 CALL CLOSE (mi,1)
 
!     DIVIDE EACH TERM BY SQRT (MI)
 
 CALL gopen (scr1,core(lcore+1),0)
 lcore = lcore - sysbuf
 CALL gopen (phia,core(lcore+1),1)
 ii = 1
 jj = nrow
 incur = 1
 ita1  = itb
 itb1  = itb
 ncol2 = itb*ncol
 nrow2 = itb*nrow
 ii1   = 1
 jj1   = nrow
 incur1= 1
 DO  i = 1,ncol
   CALL unpack (*390,scr1,core(ncol2+1))
   IF (itb /= 2) GO TO 368
   DO  j = 1,nrow
     k = ncol + j
     dcore(k) = dcore(k)*dcore(i)
   END DO
   GO TO 380
   368 DO  j = 1,nrow
     k = ncol+j
     core(k) = core(k)*core(i)
   END DO
   380 CALL pack (core(ncol2+1),phia,iphia)
   CYCLE
   390 DO  j = 1,nrow2
     k = ncol2 + j
     core(k) = 0.0
   END DO
   GO TO 380
 END DO
 CALL CLOSE (phia,1)
 CALL CLOSE (scr1,1)
 GO TO 140
 
!     COPY ROUTINE - PHIA TO SCR1
 
 420 lcore = lcore - sysbuf
 CALL gopen (phia,core(lcore+1),0)
 lcore = lcore - sysbuf
 CALL gopen (scr1,core(lcore+1),1)
 dcore(1) = 0.0D+0
 itb   = ix(5)
 ita1  = itb
 itb1  = itb
 incur = 1
 incur1= 1
 DO  jjj = 1,ncol
   ii = 0
   CALL unpack (*435,phia,core(3))
   ii1 = ii
   jj1 = jj
   CALL pack (core(3),scr1,iphia)
   CYCLE
   435 ii1 = 1
   jj1 = 1
   CALL pack (core,scr1,iphia)
 END DO
 CALL CLOSE (phia,1)
 CALL CLOSE (scr1,1)
 lcore = lcore + 2*sysbuf
 GO TO icopy, (40,320,150)
 490 CALL mesage (-2,ip1,nam)
 500 ip1 = lama
 GO TO 490
 510 ip1 = ipout
 GO TO 490
 520 CALL mesage (-8,0,nam)
 530 CALL mesage (-3,lama,nam)
 540 CALL mesage (-5,mi,nam)
 
 
 ENTRY read5 (ipout)
!     ===================
 
!     PUT OUT EIGENVALUE SUMMARY IN CASE NO EIGENVALUES FOUND
 
 lcore = korsz(core) - sysbuf
 istor = 0
 jstor = 0
 imsg  = 0
 xmax1 = 0.
 ix(2) = 1
 DO  i = 3,7
   ix(i) = 0
 END DO
 GO TO 275
 
!     REARRANGE THE EIGENVALUE TABLE, IF NECESSARY, FOR GIVENS METHOD
 
 600 CALL gopen (lama,core(lcore+1),0)
 CALL skprec (lama,1)
 nwords = 7*nlama
 CALL READ (*500,*530,lama,core(1),nwords,1,nwrds)
 refreq = core(3)
 DO  i = 2,nlama
   j = 7*(i-1) + 3
   IF (core(j) >= refreq) CYCLE
   refreq = core(j)
   GO TO 660
 END DO
 GO TO 740
 660 CALL bckrec (lama)
 CALL CLOSE (lama,2)
 CALL gopen (lama,core(lcore+1),3)
 DO  i = 1,nlama
   IF (core(3) == refreq) EXIT
   t2 = core(2)
   t3 = core(3)
   t4 = core(4)
   t5 = core(5)
   t6 = core(6)
   t7 = core(7)
   DO  j = 2,nlama
     k = 7*(j-2)
     core(k+2) = core(k+ 9)
     core(k+3) = core(k+10)
     core(k+4) = core(k+11)
     core(k+5) = core(k+12)
     core(k+6) = core(k+13)
     core(k+7) = core(k+14)
   END DO
   k = 7*(nlama-1)
   core(k+2) = t2
   core(k+3) = t3
   core(k+4) = t4
   core(k+5) = t5
   core(k+6) = t6
   core(k+7) = t7
 END DO
 720 CALL WRITE (lama,core(1),nwords,1)
 740 CALL CLOSE (lama,1)
 GO TO 275
END SUBROUTINE read2
