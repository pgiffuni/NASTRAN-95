SUBROUTINE trlgc (tmldtb,trl,dit,itrl,fct,fco,tol,iflag)
     
!     THE PURPOSE OF THIS SUBROUTINE IS TO PRODUCE A MATRIX OF FUNCTIONS
!     OF TIME.  EACH COLUMN IS A TIME STEP (AS DEFINE BY TRL) AND EACH
!     TERM IN A COLUMN CORRESPONDS TO A UNIQUE FUNCTION OF TIME (EITHER
!     BY TABLE FROM TLOAD, TIME DELAY, OR QVECT)
 
!     INPUTS (3)
!         TMLDTB - TABLE SHOWING TIME DEPENDANT DATA
!         TRL    - TIME STEP LIST
!         DIT    - DIRECT INPUT TABLES
!         ITRL   - SELECTED TRL SET NUMBER FROM CASECC
 
!     OUTPUTS(3)
!         FCT    - TIME FUNCTIONS AT ALL TIMES
!         FCO    - TIME FUNCTIONS AT OUTPUT TIMES
!         IFLAG  - -1 IMPLIES ALL TIMES OUTPUT (I.E. FCO = FCT)
!         TOL    - TABLE OF OUTPUT TIMES
 
!     THE FORMAT OF THE  TMLDTB TABLE IS AS FOLLOWS
!         REC NO.  WORD  DESCRIPTION
!         0        1-2   TABLE NAME
!         1        1     TERM NUMBER
!                  2     TLOAD ID
!                  3     TLOAD TYPE(1,2)
!                  4     TAU ( FROM DELAY CARDS--REAL)
!                  5     TID (TABLES FROM TLOAD1 CARD)
!                  5     T1   CONSTANTS FROM TLOAD 2 CARDS
!                  6     T2
!                  7     F
!                  8     P
!                  9     C
!                  10    B
!                  11    QVECT POINTER INTO SECOND RECORD
 
!         WORDS  1 THRU 11 ARE REPEATED FOR EACH UNIQUE TIME FUNCTION
 
!         2        1    I1   QVECT TABLE ID'S
!                  2    I2
!                  3    I3
!                  4    V1   QVECT ORIENTATION VECTORS
!                  5    V2
!                  6    V3
!                  7    V4
!                  8    V5
!                  9    V6
 
!     CORE LAYOUT IS AS FOLLOWS $                                POINT
!     ========================================  ===============  =====
!     TERM DESCRIPTORS (11 WORDS PER TERM)      11*NTERM WORDS   ITERM
!     QVECT STUFF      (9  WORDS PER QVECT)     9*NQVECT WORDS   IQVECT
!     TRL   STUFF      (3 WORDS PER GROUP)      3*NGROUP WORDS+1 TGROUP
!     TABLE LIST       (1 WORD  PER UNIQUETAB)  NTAB WORDS+1     ITAB
!     TABLE DATA        PRETAB STORED           LTAB WORDS       ILTAB
!     TERM VALUES                               NTERM WORDS      IVS
 
!     3    BUFFERS      FCT                                      IBUF1
!                       FCO                                      IBUF2
!                       TOL                                      IBUF3
 
 
 INTEGER, INTENT(IN)                      :: tmldtb
 INTEGER, INTENT(IN)                      :: trl
 INTEGER, INTENT(IN OUT)                  :: dit
 INTEGER, INTENT(IN OUT)                  :: itrl
 INTEGER, INTENT(IN OUT)                  :: fct
 INTEGER, INTENT(IN OUT)                  :: fco
 INTEGER, INTENT(IN)                      :: tol
 INTEGER, INTENT(OUT)                     :: iflag
 LOGICAL :: dec
 INTEGER :: mcb(7),mcb1(7),NAME(2), sysbuf, FILE,itlist(13),iz(1)
 COMMON /BLANK / dummy,ncont
 COMMON /machin/ mach
 COMMON /zzzzzz/ z(1)
 COMMON /zblpkx/ za(4),ii1
 COMMON /system/ sysbuf
 COMMON /packx / it1,it2,ii,jj,incr
 COMMON /condas/ consts(5)
 EQUIVALENCE     (z(1),iz(1)),(consts(2),twopi),(consts(4),degra)
 DATA    NAME  / 4HTRLG,4HC   /
 DATA    itlist/ 4,1105,11,1,1205,12,2,1305,13,3,1405,14,4 /
 
 dec    = mach == 5 .OR. mach == 6 .OR. mach == 21
 noload = 0
 mcb(1) = tmldtb
 CALL rdtrl (mcb)
 IF (mcb(2) <= 0) noload = -1
 mcb(2) = 100
 igroup = 1
 iflag  =-1
 nz     = korsz(z)
 ibuf1  = nz    - sysbuf
 ibuf2  = ibuf1 - sysbuf
 ibuf3  = ibuf2 - sysbuf
 nz     = ibuf3 - 1
 IF (nz <= 0) CALL mesage (-8,0,NAME)
 
!     BRING IN  TIME DATA
 
 IF (noload /= 0) GO TO 30
 iterm = 1
 lrec  = 11
 FILE  = tmldtb
 CALL gopen (tmldtb,iz(ibuf1),0)
 CALL READ (*290,*10,tmldtb,iz(iterm),nz,0,ilen)
 CALL mesage (-8,0,NAME)
 10 nterm = ilen/lrec
 iqvec = iterm + ilen
 nz    = nz - ilen
 
!     BRING IN  QVECT DATA
 
 CALL READ (*290,*20,tmldtb,iz(iqvec),nz,0,ilen)
 CALL mesage (-8,0,NAME)
 20 nqvect = ilen/9
 igroup = iqvec + ilen
 nz     = nz - ilen
 CALL CLOSE (tmldtb,1)
 
!     FIND TRL STUFF FOR CORE
 
 30 FILE = trl
 CALL OPEN (*310,trl,iz(ibuf1),0)
 CALL fread (trl,iz(igroup),3,1)
 CALL skprec (trl,iz(igroup+2))
 40 CALL READ (*291,*50,trl,iz(igroup),nz,0,ilen)
 CALL mesage (-8,0,NAME)
 50 IF (iz(igroup) /= itrl) GO TO 40
 ngroup = (ilen-1)/3
 itab   = igroup + ilen
 nz     = nz - ilen
 CALL CLOSE (trl,1)
 IF (noload /= 0) GO TO 122
 
!     BUILD LIST OF UNIQUE TABLES
 
 ntabl = 1
 k     = itab + ntabl
 iz(k) = 0
 DO  i = 1,nterm
   k = iterm + lrec*(i-1) + 4
   IF (iz(k-2) /= 3) GO TO 60
   itid = iz(k)
   ASSIGN  60  TO iret
   GO TO 90
   60 k = iterm + lrec*(i-1) + 10
   IF (iz(k) == 0) CYCLE
   
!     LOOK AT QVECT  TABLE  ID S
   
   iq   = (iz(k)-1)*9 + iqvec
   itid = iz(iq)
   ASSIGN 70  TO iret
   GO TO 90
   70 itid = iz(iq+1)
   ASSIGN 80 TO iret
   GO TO 90
   80 itid = iz(iq+2)
   ASSIGN 120 TO iret
   
!     SEARCH TABLE LIST
   
   90 l = numtyp(itid)
   IF (dec .AND. itid > 16000 .AND. itid <= 99999999) l = 1
   IF (itid <= 0 .OR. l /= 1) GO TO 110
   DO  l = 1,ntabl
     k =  itab + l
     IF (iz(k) == itid) GO TO 110
   END DO
   
!     NEW TABLE
   
   ntabl = ntabl + 1
   k     = itab  + ntabl
   iz(k) = itid
   110 GO TO iret, (60,70,80,120)
 END DO
 iz(itab) = ntabl
 iltab = itab + ntabl + 1
 nz    = nz - ntabl - 1
 
!     BRING IN TABLE STUFF
 
 ltab = 0
 IF (ntabl == 1) GO TO 121
 CALL pretab (dit,iz(iltab),iz(iltab),iz(ibuf1),nz,ltab,iz(itab), itlist)
 121 CONTINUE
 nz  = nz - ltab
 ivs = iltab + ltab
 IF (nz < nterm) CALL mesage (-8,0,NAME)
 
!     SET UP FOR PACK
 
 it1 = 1
 it2 = 1
 ii  = 1
 jj  = nterm
 incr= 1
 CALL makmcb (mcb, fct,nterm,2,it2)
 CALL makmcb (mcb1,fco,nterm,2,it2)
 
!     OPEN OUTPUT FILES
 
 CALL gopen (fct,iz(ibuf1),1)
 122 CONTINUE
 FILE = tol
 TO   = 0.0
 IF (ncont <= 2) GO TO 123
 
!     BRING BACK LAST TIME FOR CONTINUE MODE
 
 CALL OPEN  (*310,tol,iz(ibuf2),0)
 CALL fread (tol,TO,-ncont-1,0)
 CALL fread (tol,TO,1,1)
 CALL CLOSE (tol,1)
 123 CONTINUE
 CALL OPEN  (*310,tol,iz(ibuf2),1)
 CALL fname (tol,za)
 CALL WRITE (tol,za,2,0)
 IF (noload /= 0) GO TO 150
 
!     DETERMINE IF ALL TIME STEPS OUTPUT
 
 DO  i = 1,ngroup
   k =  igroup + (i-1)*3 + 3
   IF (iz(k) /= 1) GO TO 140
 END DO
 iflag = -1
 GO TO 150
 140 iflag = 1
 CALL gopen (fco,iz(ibuf3),1)
 150 CONTINUE
 t   = TO
 ist = -1
 DO  i = 1,ngroup
   
!     PICK UP  TIME CONSTANTS
   
   k = igroup + (i-1)*3 + 1
   nstep  = iz(k)
   IF (i == ngroup) nstep = nstep + 1
   nout   = iz(k+2)
   deltat =  z(k+1)
   IF (i == 1) nstep = nstep + 1
   DO  j = 1,nstep
     IF (noload /= 0) GO TO 231
     DO   l = 1,nterm
       ip = iterm + (l-1)*lrec
       m  = iz(ip+2) - 2
       SELECT CASE ( m )
         CASE (    1)
           GO TO 160
         CASE (    2)
           GO TO 170
       END SELECT
       
!     TLOAD1  CARD
       
       160 tt = t - z(ip+3)
       CALL tab (iz(ip+4),tt,ft)
       GO TO 200
       
!     TLOAD2  CARD2
       
       170 tt   = t - z(ip+3) - z(ip+4)
       zrad = z(ip+7)*degra
       IF (tt == 0.0) GO TO 180
       IF (tt < 0.0 .OR. tt > z(ip+5)-z(ip+4)) GO TO 190
       ft = tt**z(ip+9)*EXP(z(ip+8)*tt)*COS(twopi*z(ip+6)*tt + zrad)
       GO TO 200
       
!     TT = 0.0  TRY  LIMITS OF EXPRESSION
       
       180 IF (z(ip+ 9) /= 0.0) GO TO 190
       ft = COS(zrad)
       GO TO 200
       
!     FT = 0.0
       
       190 ft = 0.0
       
!     NOW TRY FOR  QVECT  STUFF
       
       200 IF (iz(ip+10) == 0) GO TO 220
       
!     EVALUATE  QVECT FUNCTION
       
       iq = (iz(ip+10)-1)*9 + iqvec
       tt = t - z(ip+3)
       
!     CHECK FOR CONSTANT FLUX VALUE (FLOATING POINT).
!     IF TIME DEPENDENT, CALL TABLE LOOKUP.
       
       iq1 = iz(iq)
       q1  = z(iq)
       lx  = numtyp(iq1)
       IF (dec .AND. iq1 > 16000 .AND. iq1 <= 99999999) lx = 1
       IF (iq1 <= 0 .OR. lx /= 1) GO TO 202
       CALL tab (iq1,tt,q1)
       202 iq2 = iz(iq+1)
       q2  = z(iq+1)
       lx  = numtyp(iq2)
       IF (dec .AND. iq2 > 16000 .AND. iq2 <= 99999999) lx = 1
       IF (iq2 <= 0 .OR. lx /= 1) GO TO 204
       CALL tab (iq2,tt,q2)
       204 iq3 = iz(iq+2)
       q3  = z(iq+2)
       lx  = numtyp(iq3)
       IF (dec .AND. iq3 > 16000 .AND. iq3 <= 99999999) lx = 1
       IF (iq3 <= 0 .OR. lx /= 1) GO TO 206
       CALL tab (iq3,tt,q3)
       206 IF (z(iq+6) /= 0.0 .OR. z(iq+6) /= 0.0 .OR. z(iq+7) /= 0.0 .OR.  &
           z(iq+8) /= 0.0) GO TO 210
       
!     V2 = 0
       
       rt = q1*z(iq+3) + q2*z(iq+4) + q3*z(iq+5)
       IF (rt > 0.0) rt = 0.0
       ft = -rt*ft
       GO TO 220
       
!     V2   0
       
       210 ft = SQRT((q1*z(iq+3) + q2*z(iq+4) + q3*z(iq+5))**2 +  &
           (q1*z(iq+6) + q2*z(iq+7) + q3*z(iq+8))**2)*ft
       GO TO 220
       
!     PUT IN FT
       
       220 m    = ivs + l - 1
       z(m) = ft
     END DO
     
!     COLUMN BUILT
     
     CALL pack (z(ivs),fct,mcb)
     231 CONTINUE
     IF (i == ngroup .AND. j == nstep-1) GO TO 240
     IF (j == 1 .OR. j == nstep) GO TO 240
     IF (MOD(j+ist,nout) /= 0) GO TO 260
     
!     OUTPUT TIME
     
     240 CALL WRITE (tol,t,1,0)
     IF (iflag == -1) GO TO 250
     CALL pack (z(ivs),fco,mcb1)
     250 IF (j == nstep) deltat = z(k+4)
     260 t = t + deltat
   END DO
   ist = 0
 END DO
 
!     ALL OUTPUT
 
 CALL WRITE (tol,0,0,1)
 CALL CLOSE (tol,1)
 IF (noload /= 0) GO TO 281
 CALL CLOSE (fct,1)
 CALL wrttrl (mcb)
 IF (iflag == -1) GO TO 281
 CALL CLOSE (fco,1)
 CALL wrttrl (mcb1)
 281 CONTINUE
 mcb(1) = tol
 CALL wrttrl (mcb)
 RETURN
 
!     ERROR MESSAGES
 
 290 ip1 = -2
 300 CALL mesage (ip1,FILE,NAME)
 RETURN
 310 ip1 = -1
 GO TO 300
 
!     NO PROPER TSTEP CARD FOUND
 
 291 CALL mesage (-31,itrl,NAME)
 RETURN
END SUBROUTINE trlgc
