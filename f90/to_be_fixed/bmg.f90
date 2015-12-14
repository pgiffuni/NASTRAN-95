SUBROUTINE bmg
     
!     HYDROELASTIC BOUNDARY MATRIX GENERATOR
 
!     7/12/73 NO AXIAL SYMMETRY UPPER INTEGRATION LIMIT OF LAST
!             CIRCUMFERENTIAL GRID IS INCORRECT
 
 LOGICAL :: labfl    ,lkbfl    ,nstar    ,head
 INTEGER :: sysbuf   ,nabfl(2) ,nkbfl(2) ,subr(2)  ,FORM   ,  &
     buf(10)  ,z        ,scrt1    ,bdpool   ,cstm   ,  &
     eqexin   ,bgpdt    ,rd       ,rdrew    ,entrys ,  &
     wrt      ,wrtrew   ,clsrew   ,dmig(3)  ,cls    ,  &
     FILE     ,buf1     ,buf2     ,buf3     ,core   ,  &
     point    ,bndfl(2) ,flag     ,mones(3) ,eor
 REAL :: rbuf(10) ,rz(1)    ,kii
 DOUBLE PRECISION :: dz(1)    ,dtemp(3) ,term(3)  ,vi(3)    ,t0f(9) ,  &
     ti(9)    ,ain      ,dub
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm      ,uwm      ,uim      ,sfm
 COMMON /condas/  consts(5)
 COMMON /system/  sysbuf   ,iout     ,iskp(52) ,iprec
 COMMON /BLANK /  kflags(2),value(2)
 COMMON /zzzzzz/  z(1)
 COMMON /names /  rd,rdrew,wrt,wrtrew,clsrew,cls
 EQUIVALENCE      (consts(2),twopi)  ,(consts(4),degrad)  ,  &
     (z(1),rz(1),dz(1)) ,(buf(1),rbuf(1))
 DATA    subr  /  4HBMG ,4H    /,  nabfl / 4HABFL,4H    /
 DATA    bndfl /  9614  ,96    /,  mones / -1, -1, -1   /
 DATA    is    /  1            /,  dmig  / 114, 1, 120  /
 DATA    eor   ,  noeor/ 1,0   /,  nkbfl / 4HKBFL,4H    /
 DATA    matpol,  bgpdt,eqexin,cstm / 101,102,103,104   /
 DATA    bdpool/  201          /,  scrt1 / 301          /
 DATA    iz2   ,  iz6,iz7,iz8,iz9   /  2, 6, 7, 8, 9    /
 
!     DEFINE CORE AND BUFFER POINTERS
 
 core = korsz(z)
 buf1 = core - sysbuf - 2
 buf2 = buf1 - sysbuf - 2
 buf3 = buf2 - sysbuf - 2
 core = buf3 - 1
 IF (core < 100) CALL mesage (-8,0,subr)
 kflags(1) = -1
 kflags(2) = -1
 
!     OPEN MATPOOL AND LOCATE THE BNDFL RECORD AS PREPARED BY IFP4.
 
 CALL preloc (*10000,z(buf1),matpol)
 CALL locate (*10000,z(buf1),bndfl,flag)
 
!     THIS MODULE DOES NOTHING IF THE MATPOOL IS PURGED OR THE BNDFL
!     RECORD IS ABSENT.  NOW READ THE HEADER DATA OF THIS RECORD.
 
 FILE = matpol
 CALL READ (*10002,*10003,matpol,z(1),9,noeor,flag)
 in = 10
 nn = z(iz9) + 9
 IF (nn+5 > core) CALL mesage (-8,0,subr)
 
!     READ THE INDICES
 
 CALL READ (*10002,*10003,matpol,z(in),z(iz9),noeor,flag)
 value(1) = z(iz6)
 value(2) = 0.0
 IF (z(iz6) == 0) value(1) = 1.0
 
!     MODIFY LIST OF INDICES TO FIT THE FOLLOWING TABLE
 
!        M      S1    S2         N              N*
!        -      --    --         -              --
!        0                      ALL             ALL
 
!                               K M
!     .GE.2      S     S        ---             NONE
!                                2
 
!                             (2K+1)M
!     .GE.2      S     A      -------           NONE
!                                4
 
!                                              (2K+1)M
!     .GE.2      A     S        NONE           -------
!                                                 4
 
!                                               K M
!     .GE.2      A     A        NONE            ---
!                                                2
 
!     K MAY BE 0,1,2,..... IN ORDER TO CHECK INDICE FOR MATCH.
 
 IF (z(iz6) > 0.0) THEN
   GO TO    10
 ELSE
   GO TO    90
 END IF
 
!     M IS POSITIVE THUS CHECK FOR STAR OR NO-STAR INDICES PERMITTED.
!     DETERMINE THE FORM OF THE CHECK EQUATION.
 
!          Z(7) = S1
!          Z(8) = S2
!          Z(6) = M
 
 10 IF (z(iz7) == is) GO TO 12
 nstar = .true.
 GO TO 13
 12 nstar = .false.
 13 IF (z(iz7) == z(iz8)) GO TO 14
 FORM = 1
 GO TO 15
 14 FORM = -1
 
!     NOW FORM NEW LIST OF INDICES
 
 15 inn = nn + 1
 nnn = nn
 DO  i = in,nn
   n = (z(i)-1)/2
   IF (MOD(z(i),2) == 0.0) THEN
     GO TO    20
   ELSE
     GO TO    30
   END IF
   
!     NON-STAR CASE
   
   20 IF (nstar) CYCLE
   IF (FORM ) 40,70,70
   
!     STAR CASE
   
   30 IF (.NOT.nstar) CYCLE
   IF (FORM) 40,70,70
   
!                           K M
!     CHECK USING EQUATION  ---
!                            2
   
   40 n2 = n*2
   k  = n2/z(iz6)
   IF (k*z(iz6) /= n2) CYCLE
   
!     GOOD INDICE,  ADD IT TO THE LIST
   
   50 nnn = nnn + 1
   z(nnn) = z(i)
   CYCLE
   
!                            (2K+1)M
!     CHECK USING EQUATION   -------
!                               4
   
   70 n4 = n*4
   ik = n4 / z(iz6)
   ik = ik - 1
   k  = ik / 2
   IF ((2*k+1)*z(iz6) == n4) GO TO 50
 END DO
 
!     LIST IS COMPLETE
 
 in = inn
 nn = nnn
 90 labfl = .true.
 IF (nn < in) labfl = .false.
 
!     SET LKBFL AS A FLAG INDICATING WHETHER KBFL WILL BE GENERATED
!     ALONG WITH ABFL.  IF G IS NON-ZERO THEN KBFL WILL BE GENERATED.
 
 lkbfl = .true.
 IF (rz(iz2) == 0.0) lkbfl = .false.
 IF (lkbfl) kflags(1) = 0
 IF (labfl) kflags(2) = 0
 IF (.NOT.labfl .AND. .NOT.lkbfl) GO TO 10000
 
!     BGPDT IS NOW READ INTO CORE AS 5 WORD ENTRIES, RESERVING FIRST
!     WORD FOR THE EXTERNAL ID.
 
 FILE   = bgpdt
 ibgpdt = nn + 1
 nbgpdt = nn
 CALL gopen (bgpdt,z(buf2),rdrew)
 100 CALL READ (*10002,*120,bgpdt,z(nbgpdt+2),4,noeor,flag)
 nbgpdt = nbgpdt + 5
 IF (nbgpdt+5 > core) CALL mesage (-8,0,subr)
 GO TO 100
 120 CALL CLOSE (bgpdt,clsrew)
 
!     READ EQEXIN PLACING EXTERNAL ID ON RESPECTIVE BGPDT ENTRY.
 
 FILE = eqexin
 CALL gopen (eqexin,z(buf2),rdrew)
 130 CALL READ (*10002,*140,eqexin,buf,2,noeor,flag)
 n = 5*buf(2) - 5 + ibgpdt
 z(n) = buf(1)
 GO TO 130
 140 CALL CLOSE (eqexin,clsrew)
 lbgpdt = nbgpdt - ibgpdt + 1
 entrys = lbgpdt / 5
 
!     SORT THE BGPDT ON EXTERNAL ID
 
 CALL sort (0,0,5,1,z(ibgpdt),lbgpdt)
 
!  BLAST CSTM INTO CORE
 
 FILE = cstm
 CALL gopen (cstm,z(buf2),rdrew)
 icstm = nbgpdt + 1
 CALL READ (*10002,*160,cstm,z(icstm),core-icstm,noeor,flag)
 CALL mesage (-8,0,subr)
 160 ncstm = icstm + flag - 1
 lcstm = ncstm - icstm + 1
 CALL CLOSE (cstm,clsrew)
 
!     LOCATE THE T   MATRIX IN THE CSTM DATA BY USING CSID = CDF IN
!                 0F
 
!     THE HEADER DATA.      ( Z(1) )
 
 DO  i = icstm,ncstm,14
   IF (z(1) == z(i)) GO TO 240
 END DO
 WRITE  (iout,210) sfm,z(1)
 210 FORMAT (a25,' 4060, COORDINATE SYSTEM =',i9,  &
     ' CAN NOT BE FOUND IN CSTM DATA.')
 GO TO 9999
 240 n = i + 5
 DO  i = 1,9
   t0f(i) = DBLE(rz(n))
   n = n + 1
 END DO
 
!     OPEN BDPOOL FOR ABFL, AND SCRATCH1 FOR KBFL AND WRITE THE DMIG
!     HEADER INFORMATION.
 
 CALL gopen (bdpool,z(buf2),wrtrew)
 
!     WRITE DMIG RECORD ID
 
 CALL WRITE (bdpool,dmig,3,noeor)
 buf(1) = nabfl(1)
 buf(2) = nabfl(2)
 buf(3) = 0
 buf(4) = 1
 buf(5) = 1
 buf(6) = iprec
 buf(7) = 0
 buf(8) = 0
 buf(9) = 0
 IF (.NOT. labfl) GO TO 270
 CALL WRITE (bdpool,buf,9,noeor)
 270 IF (.NOT.lkbfl) GO TO 280
 FILE = scrt1
 CALL OPEN (*10001,scrt1,z(buf3),wrtrew)
 buf(1) = nkbfl(1)
 buf(2) = nkbfl(2)
 CALL WRITE (scrt1,buf,9,noeor)
 
!     READ SOME FLUID-PT DATA (IDF,R,Z,L,C,S,RHO)
 
 280 FILE = matpol
 CALL READ (*10002,*10003,matpol,idf,1,noeor,flag)
 285 idata = ncstm + 1
 ndata = ncstm + 6
 CALL READ (*10002,*10003,matpol,z(idata),6,noeor,flag)
 
!     START BUILDING TABLE OF CONNECTED GRID POINTS.
!     READ ID,PHI.  CREATE A 26 WORD ENTRY FOR EACH ID,PHI.
 
 itable = ndata + 1
 
!     INSURE THAT TABLE STARTS ON AN EVEN BOUNDARY FOR DOUBLE
!     PRECISION
 
 IF (MOD(itable,2) /= 1) itable = itable + 1
 ntable = itable - 1
 290 CALL READ (*10002,*10003,matpol,z(ntable+1),2,noeor,flag)
 IF (z(ntable+1) == -1) GO TO 300
 
!     CONVERT PHI TO RADIANS
 
 rz(ntable+2) = rz(ntable+2)*degrad
 ntable = ntable + 26
 IF (ntable+26 > core) CALL mesage (-8,0,subr)
 GO TO 290
 
!     COMPUTATION AND INSERTION OF PHI   AND PHI   FOR EACH ENTRY.
!                                     0         1
 
 300 DO  i = itable,ntable,26
   
!     SET UP PHI  IN THIRD SLOT OF ENTRY = (PHI  + PHI   )/2.0
!               0                              I      I-1
   
   IF (i /= itable) GO TO 310
   
!     SPECIAL CASE ON FIRST POINT, TEST M TO FIND PHI
!                                                    I-1
   
   IF (z(iz6) > 1) GO TO 320
   phil1 = rz(ntable-24) - twopi
   GO TO 350
   320 phil1 = rz(itable+1)
   GO TO 350
   310 phil1 = rz(i-25)
   350 rz(i+2) = (rz(i+1) + phil1) / 2.0
   
!     SET UP PHI  IN FOURTH SLOT OF ENTRY = (PHI  + PHI   )/2.0
!               1                               I      I+1
   
   IF (i /= ntable-25) GO TO 345
   
!     SPECIAL CASE ON LAST POINT, TEST M TO FIND PHI
!                                                   I+1
   
   IF (z(iz6) > 1) GO TO 340
   phip1 = rz(itable+1) + twopi
   GO TO 360
   340 phip1 = rz(ntable-24)
   GO TO 360
   345 phip1 = rz(i+27)
   360 rz(i+3) = (rz(i+1) + phip1) / 2.0
 END DO
 
!     PICK UP NEXT FLUID POINT IDF
 
 nextid = 0
 CALL READ (*10002,*400,matpol,nextid,1,noeor,flag)
 IF (nextid /= idf) GO TO 400
 
!     NEXTID IS SAME AS CURRENT IDF, THUS ADD ANOTHER ENTRY OF R,Z,L,C,
!     S,RH FIRST MOVE SINGLE ENTRY DOWN UNDER TABLE SO IT CAN GROW.
 
 z(ntable+1) = z(idata  )
 z(ntable+2) = z(idata+1)
 z(ntable+3) = z(idata+2)
 z(ntable+4) = z(idata+3)
 z(ntable+5) = z(idata+4)
 z(ntable+6) = z(idata+5)
 idata = ntable + 1
 ndata = ntable + 6
 380 IF (ndata+6 > core) CALL mesage (-8,0,subr)
 CALL READ (*10002,*10003,matpol,z(ndata+1),6,noeor,flag)
 ndata = ndata + 6
 
!     SKIP THE ID-PHI PAIRS AS THEY SHOULD BE IDENTICAL TO ONES ALREADY
!     IN THE TABLE.
 
 390 CALL READ (*10002,*10003,matpol,buf,2,noeor,flag)
 IF (buf(1) /= -1) GO TO 390
 
!     READ THE NEXTID
 
 nextid = 0
 CALL READ (*10002,*400,matpol,nextid,1,noeor,flag)
 IF (nextid == idf) GO TO 380
 
!     SORT THE TABLE ON FIELD ONE OF EACH ENTRY THE ID.
 
 400 CALL sort (0,0,26,1,z(itable),ntable-itable+1)
 
!                                  T
!     FOR EACH ENTRY GENERATE THE T T   MATRICE AND IF LKBFL = .TRUE.
!                                  I 0F
 
!     THE W  MATRICE.
!          I
 
 DO  i = itable,ntable,26
   
!     LOCATE THE TRANSFORMATION MATRIX IN DOUBLE PRECISION.
!     FIRST LOCATE BGPDT ENTRY
   
   kid = z(i)
   CALL bisloc (*10004,kid,z(ibgpdt),5,entrys,point)
   point = point + ibgpdt
   CALL bmgtns (z(icstm),lcstm,z(point),ti(1))
   
!     COMPUTE VI MATRIX.  (3X3)
   
   CALL gmmatd (ti(1),3,3,1, t0f(1),3,3,0, z(i+4))
   IF (.NOT. lkbfl) CYCLE
   j = (i+4)/2
   rz(i+22) = dz(j+3)
   rz(i+23) = dz(j+6)
   rz(i+24) = dz(j+9)
 END DO
 
!     GENERATION AND OUTPUT OF MATRIX COLUMNS TO THE ABFL MATRIX.
 
 IF (.NOT.labfl) GO TO 690
 DO  i = in,nn
   
!     COLUMN INDEX INFORMATION GJ,CJ FOR THIS HARMONIC COLUMN
   
   buf(1) = idf + z(i)*500000
   buf(2) = 0
   CALL WRITE (bdpool,buf,2,noeor)
   
!     TERMS OF THE COLUMN
   
   DO  j = itable,ntable,26
     
!     3 TERMS FOR THE J-TH ID ARE THE FOLLOWING SUMMATION
     
     term(1) = 0.0D0
     term(2) = 0.0D0
     term(3) = 0.0D0
     DO  k = idata,ndata,6
       
!                     N
!     COMPUTATION OF A
!                     I
       
       ain = rz(k)*rz(k+2)
       n   = (z(i) - 1) / 2
       fn  = n
       IF (n == 0) THEN
         GO TO   510
       ELSE
         GO TO   520
       END IF
       
!     N = 0
       
       510 ain = ain*DBLE(rz(j+3) - rz(j+2))
       GO TO 545
       
!     N IS POSITIVE, CHECK FOR STAR CASE = N*
       
       520 IF (MOD(z(i),2) == 0.0) THEN
         GO TO   530
       ELSE
         GO TO   540
       END IF
       530 dub = (SIN(rz(j+3)*fn) - SIN(rz(j+2)*fn)) / fn
       ain = ain*dub
       GO TO 545
       540 dub = (COS(rz(j+2)*fn) - COS(rz(j+3)*fn)) / fn
       ain = ain*dub
       
!     FORM VI MATRIX FOR THIS POINT
       
       545 dtemp(1) = rz(k+3)*COS(rz(j+1))
       dtemp(2) = rz(k+3)*SIN(rz(j+1))
       dtemp(3) = rz(k+4)
       CALL gmmatd (z(j+4),3,3,0, dtemp(1),3,1,0, vi(1))
       DO  l = 1,3
         term(l) = term(l) + ain*vi(l)
       END DO
     END DO
     
!     OUTPUT THESE 3 TERMS
     
     buf(1) = z(j)
     DO  k = 1,3
       buf(2) = k
       rbuf(3) = term(k)
       IF (rbuf(3) == 0.0) THEN
         GO TO   660
       END IF
       655 CALL WRITE (bdpool,buf,3,noeor)
     END DO
   END DO
   CALL WRITE (bdpool,mones,2,noeor)
 END DO
 
!     GENERATION AND OUTPUT OF COLUMNS TO THE KBFL MATRIX.
 
 690 IF (.NOT.lkbfl) GO TO 800
 DO  i = itable,ntable,26
   cosphi = COS(rz(i+1))
   sinphi = SIN(rz(i+1))
   angle   = rz(i+3) - rz(i+2)
   
!     PUT OUT 3 COLUMNS FOR EACH OF THESE CONNECTED GRIDPOINTS
   
!     SOLVE NOW FOR K   V  = 3X1  CONSTANT FOR THE 3 COLUMNS
!                    II  I
   
!     AND IS A SUMMATION
   
   term(1) = 0.0D0
   term(2) = 0.0D0
   term(3) = 0.0D0
   DO  j = idata,ndata,6
     kii = rz(j)*rz(j+2)*rz(j+5)*rz(iz2)*angle
     dtemp(1) = kii*rz(j+3)*cosphi
     dtemp(2) = kii*rz(j+3)*sinphi
     dtemp(3) = kii*rz(j+4)
     CALL gmmatd (z(i+4),3,3,0, dtemp(1),3,1,0, vi(1))
     DO  k = 1,3
       term(k) = term(k) + vi(k)
     END DO
   END DO
   
!     PUT OUT THE 3 COLUMNS
   
   DO  j = 1,3
     head = .false.
     l = i + j + 21
     dtemp(1) = DBLE(rz(l))*term(1)
     dtemp(2) = DBLE(rz(l))*term(2)
     dtemp(3) = DBLE(rz(l))*term(3)
     buf(1) = z(i)
     DO  k = 1,3
       buf(2) = k
       rbuf(3) = dtemp(k)
       
!     TERM IS NOT WRITTEN IF HAS A ZERO VALUE
       
       IF (rbuf(3) == 0.0) THEN
         GO TO   730
       END IF
       720 IF (head) GO TO 721
       buf(4) = z(i)
       buf(5) = j
       CALL WRITE (scrt1,buf(4),2,noeor)
       head = .true.
       721 CALL WRITE (scrt1,buf,3,noeor)
     END DO
     IF (head) CALL WRITE (scrt1,mones,2,noeor)
   END DO
 END DO
 
!     PROCESS THE NEXT FLUID POINT
 
 800 IF (nextid == 0) THEN
   GO TO   840
 END IF
 810 idf = nextid
 GO TO 285
 
!     ALL FLUID POINTS HAVE NOW BEEN PROCESSED.  APPEND THE KBFL, IF
!     ANY, DATA TO THE ABFL DATA AND WRAP UP.
 
 840 IF (labfl) CALL WRITE (bdpool,mones,2,noeor)
 IF (.NOT.lkbfl) GO TO 900
 CALL WRITE (scrt1,0,0,eor)
 CALL CLOSE (scrt1,clsrew)
 FILE = scrt1
 CALL OPEN (*10001,scrt1,z(buf3),rdrew)
 850 CALL READ (*10002,*860,scrt1,z(1),core,noeor,flag)
 CALL WRITE (bdpool,z(1),core,noeor)
 GO TO 850
 860 CALL WRITE (bdpool,z(1),flag,noeor)
 CALL WRITE (bdpool,mones,2,eor)
 900 CALL CLOSE (bdpool,clsrew)
 
!     PREPARE AND WRITE TRAILER
 
 buf(1) = bdpool
 
!     SET TRAILER BIT FOR DMIG CARDS
 
 buf(2) = 32768
 buf(3) = 0
 buf(4) = 0
 buf(5) = 0
 buf(6) = 0
 buf(7) = 0
 CALL wrttrl (buf)
 CALL CLOSE (scrt1,clsrew)
 
!     END OF PROCESSING
 
 10000 CALL CLOSE (matpol,clsrew)
 RETURN
 
!     ERROR CONDITIONS
 
 10001 CALL mesage (-1,FILE,subr)
 10002 CALL mesage (-2,FILE,subr)
 10003 CALL mesage (-3,FILE,subr)
 GO TO 10000
 10004 WRITE  (iout,10005) sfm,z(i)
 10005 FORMAT (a25,' 4061, CONNECTED FLUID POINT ID =',i10,  &
     ' IS MISSING BGPDT DATA.')
 
 9999 CALL mesage (-61,0,0)
 RETURN
 
END SUBROUTINE bmg
