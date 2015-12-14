SUBROUTINE cyct2
     
!     CYCT2 TRANSFORMS CYCLIC PROBLEMS BETWEEN SOLUTION VARIABLES AND
!     THE CYCLIC COMPONENTS
 
!     INPUT DATA BLOCKS - CYCD    CYCLIC COMPONENT CONSTRAINT DATA
!     INPUT DATA BLOCKS - KAA     MATRIX - STIFFNESS    MAY BE PURGED
!     INPUT DATA BLOCKS - MAA     MATRIX - MASS         MAY BE PURGED
!     INPUT DATA BLOCKS - V1I     MATRIX - LOAD OR DISP MAY BE PURGED
!     INPUT DATA BLOCKS - V2I     MATRIX - EIGENVECTORS MAY BE PURGED
!     INPUT DATA BLOCKS - LAMX    TABLE  - EIGENVALUES MUST EXIS IF V2I
 
!     OUTPUT DATA BLOCKS- KXX,MXX,V1O,V2O,LAMA
 
!     PARAMETERS - CDIR           INPUT,  BCD, (FORE OR BACK)
!     PARAMETERS - NSEG           INPUT,  INTEGER,NUMBER OF SEGS
!     PARAMETERS - KSEG           INPUT,  INTEGER,CYCLIC INDEX
!     PARAMETERS - CYCSEQ         INPUT,  INTEGER,ALTERNATE=-1
!     PARAMETERS - NLOAD          INPUT,  INTEGER,NUMBER OF LOAD COND
!     PARAMETERS - NOGO           OUTPUT, INTEGER,-1 = ERROR
 
!     SCRATCH FILES (6)
 
!     DEFINITION OF VARIABLES
!     LUA       LENGT OF A SET
!     ITYP      TYPE (0=ROT, 1=DIH)
!     IDIR      DIRECTION (0=FORE, 1=BACK)
!     IFLAG     1 IMPLIES KSEG = 0 OR 2*KSEG = NSEG
!     IPASS     1 IMPLIES SECOND PASS TRROUGH CYCD
!     IGC       1 IMPLIES FIRST MATRIX TYPE (GC FOR ROT)
!     ICS       1 IMPLIES FIRST COLUMN TYPE (COSINE FOR ROT)
 
 
 INTEGER :: cycd,kaa,v1i,v2i,lamx,v1o,v2o,cdir(2),cycseq,  &
     sysbuf,FILE,NAME(2),scr1,scr2,scr3,mcb(14),fore,  &
     iz(1),mcb1(7),mcb2(7),scr4,scr5,scr6
 DOUBLE PRECISION :: dz,arg,pi,COS,SIN,constd
 COMMON /unpakx/ itc,iik,jjk,incr1
 COMMON /zblpkx/ dz(2),iii
 COMMON /packx / ita,itb,ii,jj,incr
 COMMON /zzzzzz/ z(1)
 COMMON /system/ ksystm(65)
 COMMON /condad/ constd(5)
 COMMON /BLANK / cdir,nseg,kseg,cycseq,nload,nogo
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(55),iprec),  &
     (constd(1),pi    ),(kseg,kindex),  &
     (z(1),iz(1)),(mcb(1),mcb1(1)),(mcb(8),mcb2(1))
 DATA    cycd,kaa,maa,v1i,v2i,lamx,kxx,mxx,v1o,v2o,lama/  &
     101 ,102,103,104,105,106 ,201,202,203,204,205 /
 DATA    scr1,scr2,scr3,scr4,scr5,scr6  / 301 ,302 ,303 ,304 ,305 ,306   /
 DATA    NAME,fore /4HCYCT,4H2   ,4HFORE/
 
 
!IBMNB 6/93
 cycd = 101
 kaa  = 102
 maa  = 103
 v1i  = 104
 v2i  = 105
 lamx = 106
 kxx  = 201
 mxx  = 202
 v1o  = 203
 v2o  = 204
 lama = 205
 scr1 = 301
 scr2 = 302
 scr3 = 303
 scr4 = 304
 scr5 = 305
 scr6 = 306
!IBMNE
 nz    = korsz(iz)
 nogo  = 1
 v1i   = 104
 v1o   = 203
 scr3  = 303
 mcb(1)= cycd
 CALL rdtrl (mcb)
 lua   = mcb(3)
 ityp  = mcb(2) - 1
 idir  = 1
 IF (cdir(1) == fore) idir = 0
 nx    = nz
 ibuf1 = nz - sysbuf + 1
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 nz    = ibuf3 - 1
 IF (2*kseg > nseg .OR. kseg < 0 .OR. nseg <= 0) GO TO 640
 j = 2
 IF (mcb(5) == 4) j = 4
 IF (nz < j*lua) CALL mesage (-8,0,NAME)
 
!     PRODUCE GC AND GS MATRICES (ON SCR1 AND SCR2)
 
 arg = FLOAT(kseg)/FLOAT(nseg)
 arg = arg*pi
 IF (ityp == 0) arg = 2.0D0*arg
 
!     BRING IN CYCD
 
 CALL gopen (cycd,iz(ibuf1),0)
 CALL fread (cycd,iz(1),lua,1)
 CALL CLOSE (cycd,1)
 CALL gopen (scr1,iz(ibuf1),1)
 
!     COMPUTE COS AND SIN
 
 IF (ityp == 0) GO TO 30
 IF (kseg == 0) GO TO 10
 IF (2*kseg == nseg) GO TO 20
 GO TO 50
 10 COS = 1.0
 SIN = 0.0
 GO TO 60
 20 COS = 0.0
 SIN = 1.0
 GO TO 60
 30 IF (kseg   ==    0) GO TO 10
 IF (2*kseg == nseg) GO TO 40
 GO TO 50
 40 COS = -1.0
 SIN = 0.0
 GO TO 60
 50 CONTINUE
 COS = DCOS(arg)
 SIN = DSIN(arg)
 60 CONTINUE
 iflag = 0
 IF (kseg == 0 .OR. 2*kseg == nseg) iflag = 1
 IF (ityp /= 0 .OR. iflag == 0) CALL gopen (scr2,iz(ibuf2),1)
 ita = 2
 itb = 1
 incr= 1
 ii  = 1
 jj  = lua
 CALL makmcb (mcb1,scr1,lua,2,iprec)
 CALL makmcb (mcb2,scr2,lua,2,iprec)
 CALL wrttrl (mcb1)
 CALL wrttrl (mcb2)
 ipass = 0
 IF (ityp /= 0) GO TO 200
 
!     BUILD ROTATIONAL MATRICES
 
 70 l  = 1
 80 IF (iz(l) < 0) GO TO 190
 mm = iz(l)
 ip = 1
 
!     FIRST BUILD GC
 
 igc  = 1
 FILE = scr1
 
!     FIRST DO COSINE
 
 ics = 1
 IF (ipass /= 0) ics = 0
 
!     BUILD COLUMN
 
 90 CONTINUE
 CALL bldpk (2,iprec,FILE,0,0)
 IF (mm < 0) THEN
   GO TO   190
 ELSE IF (mm == 0) THEN
   GO TO   100
 ELSE
   GO TO   110
 END IF
 
!     INTERIOR POINT
 
 100 CONTINUE
 IF ((ics == 0 .AND. igc == 1) .OR. (ics == 1 .AND.  igc == 0)) GO TO 170
 iii   = l
 dz(1) = 1.0
 CALL zblpki
 GO TO 170
 
!     SIDE 1 POINTS
 
 110 IF (ics /= 0) GO TO 140
 
!     SINE COLUMN
 
 IF (igc /= 0) GO TO 160
 
!     MATRIX IS GS
 
 120 IF (l > mm) GO TO 130
 dz(1) = 1.0
 iii   = l
 CALL zblpki
 dz(1) = COS
 iii   = mm
 CALL zblpki
 GO TO 170
 130 iii   = mm
 dz(1) = COS
 CALL zblpki
 iii   = l
 dz(1) = 1.0
 CALL zblpki
 GO TO 170
 
!     COSINE COLUMN
 
 140 IF (igc /= 0) GO TO 150
 
!     MATRIX IS GS
 
 iii   = mm
 dz(1) =-SIN
 CALL zblpki
 GO TO 170
 
!     MATRIX IS GC
 
 150 GO TO 120
 
!     MATRIX IS GC
 
 160 iii   = mm
 dz(1) = SIN
 CALL zblpki
 GO TO 170
 170 CONTINUE
 CALL bldpkn (FILE,0,mcb(ip))
 IF (cycseq == 1) GO TO 180
 
!     NOW DO SINE COLUMN
 
 IF (ics   == 0) GO TO 180
 IF (iflag == 1) GO TO 190
 ics = 0
 GO TO 90
 
!     NOW DO GS
 
 180 IF (iflag == 1 .OR. ip == 8) GO TO 190
 ip  = 8
 igc = 0
 ics = 1
 FILE = scr2
 GO TO 90
 
!     CONSIDER NEXT CYCD VALUE
 
 190 l = l + 1
 IF (l <= lua) GO TO 80
 
!     GONE THRU CYCD ONCE. DONE IF CYCSEQ = -1
 
 IF (cycseq == -1) GO TO 400
 
!     MUST NOW DO SINE COLUMNS UNLESS IFLAG = 1
 
 IF (ipass == 1) GO TO 400
 IF (iflag == 1) GO TO 400
 ipass = 1
 GO TO 70
 
!     BUILD DIHEDRAL MATRICES
 
 200 ipass = 0
 210 l     = 1
 220 ip    = 1
 igc   = 1
 FILE  = scr1
 
!     FIRST DO S COLUMN
 
 ics = 1
 IF (ipass /= 0) ics = 0
 mm  = iz(l)
 IF (mm > 0 .AND. ipass == 1) GO TO 390
 230 CONTINUE
 CALL bldpk (2,iprec,FILE,0,0)
 IF (mm > 0) GO TO 280
 
!     INTERIOR POINT
 
 IF (ics /= 0) GO TO 260
 
!     A COLUMN
 
 IF (igc /= 0) GO TO 250
 
!     MATRIX IS GA  - COLUMN IS A
 
 240 dz(1) = 1.0
 iii = l
 CALL zblpki
 GO TO 370
 
!     MATRIX IS GS - COLUMN IS A
 
 250 GO TO 370
 
!     SCOLUMN
 
 260 IF (igc /= 0) GO TO 270
 
!     MATRIX IS GA - S COLUMN
 
 GO TO 370
 
!     MATRIX IS GS - COLUMN IS S
 
 270 GO TO 240
 
!     SIDE POINT
 
 280 IF (igc == 0) GO TO 350
 
!     MATRIX IS GS
 
 SELECT CASE ( mm )
   CASE (    1)
     GO TO 290
   CASE (    2)
     GO TO 320
   CASE (    3)
     GO TO 330
   CASE (    4)
     GO TO 370
 END SELECT
 290 iii   = l
 300 dz(1) = COS
 310 CALL zblpki
 GO TO 370
 320 iii   = l
 dz(1) =-SIN
 GO TO 310
 330 iii   = l
 340 dz(1) = 1.0
 GO TO 310
 
!     MATRIX IS GA
 
 350 iii = l
 SELECT CASE ( mm )
   CASE (    1)
     GO TO 360
   CASE (    2)
     GO TO 300
   CASE (    3)
     GO TO 370
   CASE (    4)
     GO TO 340
 END SELECT
 360 dz(1) = SIN
 GO TO 310
 370 CONTINUE
 CALL bldpkn (FILE,0,mcb(ip))
 IF (cycseq == 1 .OR. mm > 0) GO TO 380
 
!     NOW DO A COLUMN
 
 IF (ics == 0) GO TO 380
 ics = 0
 GO TO 230
 
!     NOW DO GA
 
 380 IF (ip == 8) GO TO 390
 ip  = 8
 igc = 0
 FILE= scr2
 ics = 1
 GO TO 230
 
!     CONSIDER NEXT CYCD VALUE
 
 390 l = l + 1
 IF (l <= lua) GO TO 220
 
!     GONE THRU CYCD ONCE - DONE IF CYCSEQ = -1
 
 IF (cycseq == -1) GO TO 400
 
!     NOW DO A COLUMNS
 
 IF (ipass == 1) GO TO 400
 ipass = 1
 GO TO 210
 
!     CLOSE UP SHOP
 
 400 CALL CLOSE (scr1,1)
 CALL CLOSE (scr2,1)
 CALL wrttrl (mcb1)
 IF (iflag == 0 .OR. ityp /= 0) CALL wrttrl (mcb2)
 itc = 1
 iik = 1
 jjk = lua
 incr1 = 1
 IF (idir /= 0) GO TO 490
 
!     FORWARD TRANSFORMATIONS
 
 
!     TRANSFORM MATRICES
 
 CALL cyct2a (kaa,kxx,scr1,scr2,scr3,scr4,scr5)
 CALL cyct2a (maa,mxx,scr1,scr2,scr3,scr4,scr5)
 
 mcb1(1) = kaa
 mcb2(1) = maa
 CALL rdtrl (mcb1(1))
 CALL rdtrl (mcb2(1))
 IF (mcb1(5) > 2 .OR.  mcb2(5) > 2) GO TO 405
 IF (mcb1(4) /= 6 .AND. mcb2(4) /= 6) GO TO 405
 mcb1(1) = kxx
 mcb2(1) = mxx
 CALL rdtrl (mcb1(1))
 CALL rdtrl (mcb2(1))
 mcb1(4) = 6
 mcb2(4) = 6
 IF (mcb1(1) > 0) CALL wrttrl (mcb1(1))
 IF (mcb2(1) > 0) CALL wrttrl (mcb2(1))
 
!     TRANSFORM LOADS
 
 405 mcb(1) = v1i
 CALL rdtrl (mcb(1))
 IF (mcb(1) <= 0) GO TO 460
 itc = mcb(5)
 IF (itc == 4 .AND. nz < 4*lua) CALL mesage (-8,0,NAME)
 CALL gopen (v1i,iz(ibuf1),0)
 CALL gopen (scr3 ,iz(ibuf2),1)
 CALL gopen (scr4,iz(ibuf3),1)
 
!     COMPUTE NUMBER OF RECORDS TO SKIP
 
 CALL makmcb (mcb1,scr3,lua,2,mcb(5))
 CALL makmcb (mcb2,scr4,lua,2,mcb(5))
 IF (kseg == 0) GO TO 420
 nskip = nload*kseg*(ityp+1)*2 - nload*(ityp+1)
 FILE  = v1i
 DO  i = 1,nskip
   CALL fwdrec (*620,v1i)
 END DO
 420 CONTINUE
 CALL cyct2b (v1i,scr3,nload,iz,mcb1)
 IF (ityp  == 0) GO TO 430
 IF (iflag /= 0) GO TO 430
 
!     COPY - PCA
 
 DO  j = 1,nload
   CALL fwdrec (*620,v1i)
 END DO
 CALL cyct2b (v1i,scr3,nload,iz,mcb1)
 
!     NOW COPY ONTO PS
 
 430 IF (ityp == 0 .AND. iflag /= 0) GO TO 440
 CALL cyct2b (v1i,scr4,nload,iz,mcb2)
 IF (iflag /= 0) GO TO 440
 IF (ityp  == 0) GO TO 440
 CALL REWIND (v1i)
 CALL fwdrec (*620,v1i)
 nlps = nskip + nload
 DO  j = 1,nlps
   CALL fwdrec (*620,v1i)
 END DO
 itc = -mcb(5)
 CALL cyct2b (v1i,scr4,nload,iz,mcb2)
 itc = mcb(5)
 
!     DONE WITH COPY
 
 440 CALL CLOSE (v1i,1)
 CALL CLOSE (scr3,1)
 CALL CLOSE (scr4,1)
 CALL wrttrl (mcb1)
 CALL wrttrl (mcb2)
 IF (iflag /= 0 .AND. ityp == 0) GO TO 450
 CALL ssg2b (scr1,scr3 ,0,scr5,1,iprec,1,v1o)
 CALL ssg2b (scr2,scr4,scr5,v1o,1,iprec,1,scr3)
 GO TO 460
 
!     NO GS
 
 450 CALL ssg2b (scr1,scr3,0,v1o,1,iprec,1,scr5)
 
!     TRANSFORM EIGENVECTORS FORWARD
 
 460 IF (v1o == v2o) RETURN
 mcb(1) = v2i
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) RETURN
 itc = mcb(5)
 IF (itc == 4 .AND. nz < 4*lua) CALL mesage (-8,0,NAME)
 IF (MOD(mcb(2),2) /= 2) CALL mesage (-7,0,NAME)
 IF (iflag /= 1 .OR. ityp /= 0) GO TO 470
 
!     IN = OUT
 
 v1o  = v2o
 scr3 = v1i
 GO TO 450
 470 CALL gopen (v2i,iz(ibuf1),0)
 CALL gopen (scr3,iz(ibuf2),1)
 CALL gopen (scr4,iz(ibuf3),1)
 ncopy = mcb(2)
 CALL makmcb (mcb1,scr3,lua,2,mcb(5))
 CALL makmcb (mcb2,scr4,lua,2,mcb(5))
 DO  i = 1,ncopy
   FILE = scr3
   ip = 1
   IF (MOD(i,2) == 0) ip = 8
   IF (MOD(i,2) == 0) FILE = scr4
   CALL cyct2b (v2i,FILE,1,iz,mcb(ip))
 END DO
 v1o = v2o
 v1i = v2i
 GO TO 440
 
!     DIRECTION IS BACK
 
 490 CONTINUE
 mcb(1) = v1i
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) GO TO 560
 iks = mcb(3)
 itc = mcb(5)
 IF (itc == 4 .AND. nz < 4*lua) CALL mesage (-8,0,NAME)
 
!     POSITION V1O
 
 mcb(1) = v1o
 IF (kindex == 0) GO TO 495
 CALL rdtrl (mcb)
 IF (mcb(2) > 0) GO TO 500
 495 CONTINUE
 CALL gopen  (v1o,iz(ibuf1),1)
 CALL CLOSE  (v1o,2)
 CALL makmcb (mcb,v1o,lua,2,mcb(5))
 CALL wrttrl (mcb)
 GO TO 510
 500 CONTINUE
 CALL gopen  (v1o,iz(ibuf1),0)
 CALL skpfil (v1o,+1)
 CALL skpfil (v1o,-1)
 CALL CLOSE  (v1o,2)
 510 CONTINUE
 IF (ityp == 0) GO TO 550
 
!     DISTRIBUTE UX1 AND UX2 FOR MULTIPLYS
 
 IF (iflag == 1) GO TO 550
 CALL makmcb (mcb1,scr3,iks,2,mcb(5))
 CALL makmcb (mcb2,scr4,iks,2,mcb(5))
 CALL gopen  (v1i,iz(ibuf1),0)
 CALL gopen  (scr3,iz(ibuf2),1)
 CALL gopen  (scr4,iz(ibuf3),1)
 CALL cyct2b (v1i,scr3,nload,iz(1),mcb1)
 CALL cyct2b (v1i,scr4,nload,iz(1),mcb2)
 CALL CLOSE  (scr3,1)
 CALL wrttrl (mcb1)
 CALL CLOSE  (scr4,1)
 CALL wrttrl (mcb2)
 CALL CLOSE  (v1i,1)
 
!     COMPUTE UCS
 
 520 CALL ssg2b  (scr1,scr3,0,scr5,0,iprec,1,scr6)
 CALL gopen  (v1o,iz(ibuf1),3)
 CALL gopen  (scr5,iz(ibuf2),0)
 mcb(1) = v1o
 CALL rdtrl  (mcb(1))
 CALL cyct2b (scr5,v1o,nload,iz(1),mcb)
 IF (ityp == 0 .AND. iflag /= 0) GO TO 540
 CALL CLOSE  (v1o,2)
 CALL CLOSE  (scr5,1)
 IF (ityp == 0 .OR. iflag /= 0) GO TO 530
 
!     COMPUTE UCA
 
 CALL ssg2b  (scr2,scr4,0,scr5,0,iprec, 0,scr6)
 CALL gopen  (v1o,iz(ibuf1),3)
 CALL gopen  (scr5,iz(ibuf2),0)
 CALL cyct2b (scr5,v1o,nload,iz(1),mcb)
 CALL CLOSE  (v1o,2)
 CALL CLOSE  (scr5,1)
 
!     COMPUTE USS
 
 CALL ssg2b  (scr1,scr4,0,scr5,0,iprec,1,scr6)
 CALL gopen  (v1o,iz(ibuf1),3)
 CALL gopen  (scr5,iz(ibuf2),0)
 CALL cyct2b (scr5,v1o,nload,iz(1),mcb)
 CALL CLOSE  (scr5,1)
 CALL CLOSE  (v1o,2)
 
!     COMPUTE USA
 
 530 CONTINUE
 CALL ssg2b  (scr2,scr3,0,scr5,0,iprec,1,scr6)
 CALL gopen  (v1o,iz(ibuf1),3)
 CALL gopen  (scr5,iz(ibuf2),0)
 CALL cyct2b (scr5,v1o,nload,iz(1),mcb)
 540 CONTINUE
 CALL CLOSE  (scr5,1)
 CALL CLOSE  (v1o,1)
 CALL wrttrl (mcb)
 GO TO 560
 
!     DO ROTATIONAL OR SPECIAL CASE DIH
 
 550 scr3 = v1i
 GO TO 520
 
!     SEE IF DONE
 
 560 mcb(1) = v2i
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) RETURN
 scr3 = 303
 itc  = mcb(5)
 IF (itc == 4 .AND. nz < 4*lua) CALL mesage (-8,0,NAME)
 
!     NOW DO EIGENVECTORS
 
 
!     COMPUTE NEW VECTORS
 
 CALL ssg2b (scr1,v2i,0,scr3,0,iprec,1,scr5)
 IF (ityp == 0 .AND. iflag == 1) GO TO 570
 CALL ssg2b (scr2,v2i,0,scr4,0,iprec,1,scr5)
 570 CONTINUE
 
!     POSITION FILES
 
 
!      SET LAMA FLAG
 
 mcb(1) = lamx
 CALL rdtrl (mcb)
 ilama  = 0
 IF (mcb(1) <= 0) ilama = 1
 CALL gopen (v2o,iz(ibuf1),1)
 IF (ilama /= 0) GO TO 571
 CALL gopen (lama,iz(ibuf2),1)
 FILE = lamx
 CALL gopen (lamx,iz(ibuf3),0)
 CALL READ  (*620,*630,lamx,iz(1),146,1,iflag)
 CALL WRITE (lama,iz(1),146,1)
 571 CONTINUE
 mcb(1) = v2i
 CALL rdtrl (mcb)
 nload = mcb(2)
 CALL makmcb (mcb,v2o,lua,2,mcb(5))
 ibuf4 = ibuf3 - sysbuf
 CALL gopen (scr3,iz(ibuf4),0)
 IF (ityp == 0 .AND. iflag == 1) GO TO 580
 ibuf5 = ibuf4 - sysbuf
 CALL gopen (scr4,iz(ibuf5),0)
 580 DO  i = 1,nload
   CALL cyct2b (scr3,v2o,1,iz(1),mcb)
   IF (ilama /= 0) GO TO 572
   CALL READ  (*620,*630,lamx,iz(1),7,0,iflag)
   CALL WRITE (lama,iz(1),7,0)
   572 CONTINUE
   IF (ityp == 0 .AND. iflag == 1) CYCLE
   IF (ilama == 0) CALL WRITE (lama,iz(1),7,0)
   CALL cyct2b (scr4,v2o,1,iz(1),mcb)
 END DO
 CALL wrttrl (mcb)
 CALL CLOSE  (v2o,1)
 CALL CLOSE  (scr3,1)
 CALL CLOSE  (scr4,1)
 IF (ilama /= 0) GO TO 573
 CALL CLOSE  (lama,1)
 CALL CLOSE  (lamx,1)
 mcb(1) = lama
 CALL wrttrl (mcb)
 573 CONTINUE
 
!     DONE
 
 RETURN
 
!     ERROR MESSAGES
 
! 600 IP1 = -1
 610 CALL mesage (ip1,FILE,NAME)
 GO TO 640
 620 ip1 = -2
 GO TO 610
 630 ip1 = -3
 GO TO 610
 640 CALL mesage (7,0,NAME)
 nogo = -1
 RETURN
END SUBROUTINE cyct2
