SUBROUTINE fa1
     
!     FA1 IS THE DRIVER FOR PART ONE OF FLUTTER ANALYSIS
 
 INTEGER :: sysbuf,out,buff,buff1,ns(2),floop,tstart,  &
     khh,bhh,mhh,qhhl,casecc,flist,fsave,kxhh,mxhh,  &
     bxhh,scr1,rec0(8),flut(10),imeth(2),fmethd,smeth,  &
     trl(10),aero(2),flfact(2),fluter(2),  &
     sr,sm,sk,pr,pm,pk,sl,iblock(12),method(4)
 REAL :: BLOCK(12),REC(8),kfreq,rho
 DIMENSION       dlt(3),z(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,out
 COMMON /BLANK / floop,tstart,icead
 COMMON /output/ hdg(96)
 COMMON /packx / iti,ito,ij,nn,incr1
 COMMON /zzzzzz/ iz(1)
 COMMON /unpakx/ iout,ii,jj,incr
 EQUIVALENCE     (rec0(1),REC(1))
 EQUIVALENCE     (BLOCK(1),iblock(1))
 EQUIVALENCE     (iz(1),z(1))
 DATA  khh /101/,   bhh /102/,  mhh /103/, qhhl /104/, casecc /105/  &
     ,   flist /106/, fsave /201/, kxhh /202/, bxhh /203/,   mxhh /204/
 DATA  scr1/301/,    ns /4HFA1 ,4H    /
 DATA  imeth   /  4HS   ,4HL   /
 DATA  nmd /4  /, method/4HK   ,4HKE  ,4HPK  ,4HINV /
 DATA  trl     /  90,1006,7*0,6 /
 DATA  aero    /  3202,32 /,   flfact /4102,41/,   fluter /3902,39/
 
 DO  i = 1,12
   iblock(i) = 0
 END DO
 ncore = korsz(iz)
 buff  = ncore - sysbuf - 1
 buff1 = buff  - sysbuf
 IF (floop /= 0) GO TO 200
 
!     FIRST TIME THROUGH FIND FMETHOD ON CASECC
 
 ifile = casecc
 CALL gopen (casecc,iz(buff+1),0)
 CALL READ (*530,*10,casecc,iz,buff,1,nwr)
 10 lcc = nwr
 CALL CLOSE (casecc,1)
 
!     GET DATA FOR REC0 OF FSAVE
 
 CALL fname (fsave,rec0)
 ifile = flist
 CALL preloc (*480,iz(buff+1),flist)
 CALL locate (*470,iz(buff+1),aero,idum)
 CALL READ (*530,*530,flist,rec0(4),4,1,nwr)
 REC(6) = REC(6)*0.5
 CALL locate (*15,iz(buff+1),flfact,idum)
 CALL READ (*530,*20,flist,iz(lcc+1),buff,1,nwr)
 GO TO 400
 15 nwr = 0
 20 lfl = nwr + lcc
 CALL locate (*450,iz(buff+1),fluter,idum)
 30 CALL READ (*530,*450,flist,flut,10,0,nwr)
 i165 = 165
 IF (flut(1) /= iz(i165)) GO TO 30
 CALL CLOSE (flist,1)
 rec0(8) = flut(9)
 iep = flut(10)
 DO  i = 1,nmd
   IF (flut(2) == method(i)) GO TO 50
 END DO
 GO TO 490
 50 rec0(3) = i
 fmethd = i
 SELECT CASE ( i )
   CASE (    1)
     GO TO 60
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 61
   CASE (    4)
     GO TO 490
 END SELECT
 60 rec0(4) = 0
 IF (flut(7) == imeth(1)) rec0(4) = 1
 IF (flut(7) == imeth(2)) rec0(4) = 2
 IF (rec0(4) == 0) GO TO 430
 smeth = rec0(4)
 GO TO 65
 
!     PK METHOD HAS LINEAR SPLINE ONLY
 
 61 rec0(4) = 2
 smeth = 2
 65 CONTINUE
 
!     BUILD RECORDS 0,1,2,3 OF SAVE
 
 ifile = fsave
 CALL OPEN (*480,fsave,iz(buff+1),1)
 CALL WRITE (fsave,rec0,8,1)
 bref = REC(6)
 rref = REC(7)
 neiw = rec0(8)
 
!     BUILD M,K,RHO LIST FOR FLUTTER LOOP
 
 sr = 0
 sm = 0
 sk = 0
 i  = lcc
 IF (i == lfl) GO TO 410
 70 i  = i + 1
 IF (iz(i) == flut(4)) sr = i
 IF (iz(i) == flut(5)) sm = i
 IF (iz(i) == flut(6)) sk = i
 80 i = i + 1
 IF (i >= lfl) GO TO 90
 IF (iz(i) == -1) GO TO 70
 GO TO 80
 90 IF (sr == 0 .OR. sm == 0 .OR. sk == 0) GO TO 410
 nrho = 0
 pr = sr
 95 pr = pr + 1
 IF (iz(pr) == -1) GO TO 97
 nrho = nrho + 1
 GO TO 95
 97 nloops = 0
 IF (fmethd /= 3) GO TO 105
 
!     J.PETKAS/LOCKHEED      3/91
!     19 LINES OF OLD CODE FOR BUILDING ELEMENTS OF FSAVE FOR PK METHOD
!     WERE IN ERROR, AND ARE NOW REPLACED BY NEXT 29 NEW LINES
 
 pm = sm
 101 pm = pm + 1
 IF (iz(pm) == -1) GO TO 130
 dlt(1) = z(pm)
 
!     CENTER LOOP ON RHO
 
 pr = sr
 102 pr = pr + 1
 IF (iz(pr) == -1) GO TO 101
 dlt(3) = z(pr)
 
!     INNER LOOP ON VELOCITY
 
 pk = sk
 103 pk = pk + 1
 IF (iz(pk) == -1) GO TO 102
 dlt(2) = z(pk)
 nloops = nloops + 1
 CALL WRITE (fsave,dlt,3,0)
 GO TO 103
 
!     ALGORITHM FOR BUILDING ELEMENTS OF FSAVE FOR K AND KE METHODS
 
 105 CONTINUE
 
!     OUTER LOOP ON MACH NUMBER
 
 pm = sm
 107 pm = pm + 1
 IF (iz(pm) == -1) GO TO 130
 dlt(1) = z(pm)
 
!     CENTER LOOP ON KFREQ
 
 pk = sk
 110 pk = pk + 1
 IF (iz(pk) == -1) GO TO 107
 dlt(2) = z(pk)
 
!     INNER LOOP ON RHO
 
 pr = sr
 120 pr = pr + 1
 IF (iz(pr) == -1) GO TO 110
 dlt(3) = z(pr)
 nloops = nloops + 1
 CALL WRITE (fsave,dlt,3,0)
 GO TO 120
 130 CALL WRITE (fsave,0,0,1)
 
!     PICK UP M AND K FROM QHHL
 
 ifile = qhhl
 CALL OPEN (*480,qhhl,iz(buff1+1),0)
 CALL READ (*530,*140,qhhl,iz(lcc+1),buff1,1,nwr)
 GO TO 400
 140 lfl = nwr + lcc
 sl  = lcc + 5
 CALL CLOSE (qhhl,1)
 rec0(1) = qhhl
 CALL rdtrl (rec0)
 np  = MIN0(iz(sl-1),rec0(2)/rec0(3))
 lfl = MIN0(lfl,2*np+sl-1)
 np  = lfl - sl + 1
 CALL WRITE (fsave,iz(sl),np,1)
 np  = np/2
 
!     WRITE CASECC RECORD AND TRAILER
 
 CALL WRITE (fsave,iz(1),lcc,1)
 CALL CLOSE (fsave,1)
 rec0(1) = fsave
 rec0(2) = floop
 rec0(3) = nloops
 rec0(4) = np
 rec0(5) = lcc
 rec0(6) = 0
 rec0(7) = nrho
 CALL wrttrl (rec0)
 GO TO 210
 200 ifile = fsave
 CALL OPEN (*480,fsave,iz(buff+1),0)
 CALL READ (*530,*530,fsave,iz(1),8,1,nwr)
 CALL CLOSE (fsave,1)
 izx   = 0
 fmethd= iz(izx+3)
 smeth = iz(izx+4)
 bref  =  z(izx+6)
 rref  =  z(izx+7)
 neiw  = iz(izx+8)
 210 rec0(1) = fsave
 CALL rdtrl (rec0)
 
!     START OF LOOPING BUMP LOOP COUNTER SET TIME AND GO
 
 floop = floop + 1
 nloops= rec0(3)
 CALL klock (tstart)
 SELECT CASE ( fmethd )
   CASE (    1)
     GO TO 220
   CASE (    2)
     GO TO 230
   CASE (    3)
     GO TO 240
   CASE (    4)
     GO TO 490
 END SELECT
 
!     K METHOD BUILD PROPER QHH ON SCR1
 
 220 CALL fa1k (smeth,kfreq,rho,scr1,0)
 GO TO 300
 
!     KE METHOD DO INCORE EIGNVALUE EXTRACTION
 
 230 rec0(1) = bhh
 CALL rdtrl (rec0)
 IF (rec0(1) > 0 .AND. rec0(7) > 0) GO TO 510
 rec0(1) = khh
 CALL rdtrl (rec0)
 ico = rec0(2)*rec0(2)*4 + 4
 235 CALL fa1k (smeth,kfreq,rho,scr1,ico)
 CALL fa1ke (scr1,kfreq,bref,rho,rref,floop,nloops)
 IF (floop >= nloops) GO TO 350
 floop = floop + 1
 GO TO 235
 
!     PK METHOD  LINEAR INTERPOLATION  AND INCORE LOOP FOR
!     EIGENVALUE CONVERGENCE
 
 240 CALL fa1pki (fsave,qhhl)
 CALL fa1pke (khh,bhh,mhh,bxhh,fsave,nloops,bref,rref,neiw,iep)
 IF (floop >= nloops) GO TO 250
 floop = floop  + 1
 GO TO 240
 
!     PHID  - KXHH   CLAMAD - BXHH
 
 250 ibuf = buff1 - sysbuf
 trl(1) = scr1
 CALL rdtrl (trl)
 IF (trl(2) == 0) GO TO 350
 CALL OPEN (*350,scr1,z(ibuf),0)
 CALL READ (*290,*255,scr1,REC,6,1,nwr)
 255 CALL READ (*290,*260,scr1,z,ibuf,1,nwr)
 260 nn = nwr/2
 CALL gopen (kxhh,z(buff),1)
 CALL gopen (bxhh,z(buff1),1)
 CALL WRITE (bxhh,trl(1),50,0)
 CALL WRITE (bxhh,hdg,96,1)
 trl(1) = kxhh
 trl(2) = 0
 trl(3) = nn
 trl(4) = 2
 trl(5) = 3
 iti = 3
 ito = 3
 ij  = 1
 incr1 = 1
 265 CALL WRITE (bxhh,REC,6,0)
 CALL pack  (z,kxhh,trl)
 CALL READ  (*280,*270,scr1,REC,6,1,nwr)
 270 CALL READ  (*280,*265,scr1,z,ibuf,1,nwr)
 280 CALL WRITE (bxhh,0,0,1)
 CALL CLOSE (bxhh,1)
 CALL CLOSE (kxhh,1)
 CALL wrttrl (trl)
 trl(1) = bxhh
 trl(2) = 1006
 trl(7) = 0
 CALL wrttrl (trl)
 290 CALL CLOSE (scr1,1)
 GO TO 350
 
!     COPY KHH TO KXHH
 
 300 CALL gopen (khh,iz(buff+1),0)
 CALL gopen (kxhh,iz(buff1+1),1)
 rec0(1) = khh
 CALL rdtrl (rec0)
 rec0(1) = kxhh
 iout = rec0(5)
 incr = 1
 i = rec0(2)
 rec0(2) = 0
 rec0(6) = 0
 rec0(7) = 0
 CALL cyct2b (khh,kxhh,i,z,rec0)
 CALL CLOSE  (khh,1)
 CALL CLOSE  (kxhh,1)
 CALL wrttrl (rec0)
 
!     BUILD BXHH = (K/B)BHH
 
 rec0(1) = bhh
 CALL rdtrl (rec0)
 IF (rec0(1) <= 0) GO TO 310
 iblock(2) = 1
 BLOCK(3) = kfreq/bref
 CALL ssg2c (bhh,0,bxhh,0,BLOCK(2))
 310 CONTINUE
 
!                2  2
!     MXHH  =  (K /B ) MHH  + (RHO*RREF/2.0) QHH
 
 iblock(2) = 1
 BLOCK (3) = (kfreq*kfreq)/(bref*bref)
 iblock(8) = 1
 BLOCK (9) = rho*rref/2.0
 CALL ssg2c (mhh,scr1,mxhh,0,BLOCK(2))
 
!     THE END
 
 350 CONTINUE
 rec0(1) = fsave
 CALL rdtrl (rec0)
 rec0(2) = floop
 CALL wrttrl (rec0)
 IF (floop == nloops) floop = -1
 icead = 1
 IF (fmethd == 2) icead = -1
 IF (fmethd == 3) icead = -1
 GO TO 600
 
!     ERROR MESSAGES
 
 400 CALL mesage (-8,0,ns)
 410 WRITE  (out,420) ufm,flut(4),flut(5),flut(6)
 420 FORMAT (a23,', ONE OR MORE OF THE FOLLOWING FLFACT SETS WERE NOT',  &
     ' FOUND - ',3I9)
 GO TO 540
 430 WRITE  (out,440) ufm,flut(7)
 440 FORMAT (a23,' 2267, INTERPOLATION METHOD ',a4,' UNKNOWN')
 GO TO 540
 450 i165 = 165
 WRITE  (out,460) ufm,iz(i165)
 460 FORMAT (a23,' 2268, FMETHOD SET',i9,' NOT FOUND')
 GO TO 540
 470 CALL mesage (-7,0,ns)
 480 CALL mesage (-1,ifile,ns)
 490 WRITE  (out,500) ufm,flut(2)
 500 FORMAT (a23,' 2269, FLUTTER METHOD ',a4,' NOT IMPLEMENTED')
 GO TO 540
 510 WRITE  (out,520) ufm,flut(2)
 520 FORMAT (a23,', FLUTTER METHOD ',a4,' NOT IMPLEMENTED WITH B ', 'MATRIX')
 GO TO 540
 530 CALL mesage (-3,ifile,ns)
 540 CALL mesage (-61,0,ns)
 
 600 RETURN
END SUBROUTINE fa1
