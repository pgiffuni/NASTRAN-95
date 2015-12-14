SUBROUTINE gp4
     
!     GP4  PERFORMS THE FOLLOWING FUNCTIONS--
!       1. READS CASECC AND MAKES ANALYSIS OF SUBCASE LOGIC
!       2. PROCESSES RIGID ELEMENTS AND ALL OTHER CONSTRAINT DATA (MPC,
!          SPC, OMIT, SUPORT, ASET, ETC.)
!       3. BUILDS THE USET FOR THE CURRENT SUBCASE
!       4. CALLS GP4SP TO EXAMINE GRID POINT SINGULARITIES
!       5. BUILDS THE RGT MATRIX AND YS VECTOR FOR CURRENT SUBCASE
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift ,rshift ,andf   ,orf    ,complf
 DIMENSION       buf(20),mpc(2) ,omit(2),suport(2)      ,spc(2) ,  &
     mpcadd(2)      ,spc1(2),spcadd(2)      ,mask(6),  &
     NAME(2),mcb(7) ,mcbust(7)      ,mcbys(7)       ,  &
     omitx1(2)      ,aset(2),aset1(2)       ,mak(4) , spcd(2),ctype(18)
 REAL :: rz(1)  ,bufr(2)
 CHARACTER (LEN=23) :: ufm
!WKBI 3/95 NCL94002
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
!WKBR 3/95 NCL94002      COMMON /XMSSG / UFM
 COMMON /xmssg / ufm    ,uwm    ,uim
 COMMON /machin/ mach   ,ihalf  ,jhalf
 COMMON /bitpos/ um     ,uo     ,ur     ,usg    ,usb    ,ul     ,  &
     ua     ,uf     ,us     ,un     ,ug
 COMMON /BLANK / luset  ,mpcf1  ,mpcf2  ,single ,omit1  ,react  ,  &
     nskip  ,repeat ,nosets ,nol    ,noa    ,idsub  , iautsp
 COMMON /gp4fil/ geomp  ,bgpdt  ,cstm   ,rgt    ,scr1
 COMMON /gp4prm/ buf    ,buf1   ,buf2   ,buf3   ,buf4   ,knkl1  ,  &
     mask16 ,nogo   ,gpoint ,kn
 COMMON /gp4spx/ mskum  ,mskuo  ,mskur  ,mskus  ,mskul  ,msksng ,  &
     spcset ,mpcset ,nauto  ,iogpst
 COMMON /names / rd     ,rdrew  ,wrt    ,wrtrew ,clsrew
 COMMON /packx / ita1   ,itb1   ,ii1    ,jj1    ,incr1
 COMMON /system/ ksystm(65)
 COMMON /two   / two(32)
 COMMON /unpakx/ itb    ,ii     ,jj     ,incr
 COMMON /zblpkx/ x(4)   ,ix
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm( 1),sysbuf), (ksystm( 2),outtap ),  &
     (ksystm(27),iaxic ), (ksystm(38),iaxif  ),  &
     (z(1)      ,rz(1) ), (buf(1)    ,bufr(1)),  &
     (ugset     ,usgset), (ib6       ,buf(6) )
 DATA    omit  / 5001,   50/, suport/ 5601,   56/,  &
     spc   / 5501,   55/, spc1  / 5481,   58/,  &
     spcadd/ 5491,   59/, omitx1/ 4951,   63/,  &
     aset  / 5561,   76/, aset1 / 5571,   77/,  &
     spcd  / 5110,   51/, mpc   / 4901,   49/,  &
     mpcadd/ 4891,   60/
 DATA    NAME  / 4HGP4  ,4H    /
 DATA    mset  / 4H m   /,   sg/4H sg   /, r/ 4H r  /
 DATA    ys    , uset   /202    ,203    /
 DATA    scr2           /302            /
 DATA    mpcax1, mpcax2 /101    ,102    /
 DATA    casecc, eqexin ,gpdt   /101    ,103  ,104  /
 DATA    ctype / 4HMPC , 4H    , 4HOMIT, 4H    ,  &
     4HOMIT, 4H1   , 4HSUPO, 4HRT  ,  &
     4HSPC1, 4H    , 4HSPC , 4H    ,  &
     4HSPCD, 4H    , 4HASET, 4H    ,  &
     4HASET, 4H1   /
 DATA    iz2,iz3,iz5,iz16,iz138/ 2, 3, 5, 16, 138   /
 
!     PERFORM GENERAL INITIALIZATION
 
!WKBI 3/95 NCL94002
 CALL sswtch ( 51, l51 )
 geomp  = 102
 bgpdt  = 105
 cstm   = 106
 rgt    = 201
 scr1   = 301
 nauto  = 0
 iogpst = -1
 buf1   = korsz(z) - sysbuf - 2
 buf2   = buf1 - sysbuf
 buf3   = buf2 - sysbuf
 buf4   = buf3 - sysbuf
 icrq   = luset- buf4
 insuff = 10
 IF (luset >= buf4) GO TO 2430
 mask16 = jhalf
 mask15 = jhalf/2
 n23    = 2
 mskum  = two(um )
 mskuo  = two(uo )
 mskur  = two(ur )
 mskusg = two(usg)
 mskusb = two(usb)
 mskul  = two(ul )
 mskua  = two(ua )
 mskuf  = two(uf )
 mskus  = two(us )
 mskun  = two(un )
 mskug  = two(ug )
 mskung = orf(mskun,mskug)
 mskfng = orf(mskuf,mskung)
 msksng = orf(mskus,mskung)
 mask(1)= orf(mskum,mskug)
 mask(2)= orf(mskuo,mskfng)
 mask(3)= orf(mskur,orf(mskua,mskfng))
 mask(4)= orf(mskusg,msksng)
 mask(5)= orf(mskusb,msksng)
 mask(6)= orf(mskul,orf(mskua,mskfng))
 mak(1) = orf(mskum,mskul)
 mak(2) = orf(mskus,mskul)
 mak(3) = orf(mskuo,mskul)
 mak(4) = orf(mskur,mskul)
 CALL makmcb (mcbys,ys,0,2,1)
 CALL makmcb (mcbust,uset,luset,0,0)
 multi  = -1
 usgset = -1
 single = -1
 omit1  = -1
 nosets = -1
 asetx  = -1
 react  = -1
 noys   =  0
 nogeom =  0
 nol    = -1
 noa    = +1
 nogo   =  0
 nogoof =  0
 dup    =  0
 iflag  =  0
 flag   =  0
 mskck  = complf(lshift(complf(0),20))
 rigid  =  0
 spcold = -1
 mpcold = -1
 l21    =  0
 l22    =  0
 mcb(1) = geomp
 CALL rdtrl (mcb(1))
 IF (mcb(1) < 0) GO TO 20
 
!     BIT ASSIGNMENTS FOR RIGID ELEMENTS -
!     CRIGD1 - 53       CRROD    - 65       CRBE1 - 68
!     CRIGD2 - 54       CRBAR    - 66       CRBE2 - 69
!     CRIGD3 - 83       CRTRPLT  - 67       CRBE3 - 70
!     CRIGDR - 82       CRSPLINE - 71
 
 IF (andf(mcb(5),two(21)) == two(21)) rigid = 1
 IF (andf(mcb(5),two(22)) == two(22)) rigid = 1
 IF (andf(mcb(7),two(19)) == two(19)) rigid = 1
 IF (andf(mcb(7),two(18)) == two(18)) rigid = 1
 i = mcb(6)
 DO  j = 17,23
   IF (andf(i,two(j)) == two(j)) rigid = 1
 END DO
 CALL makmcb (mcb,rgt,0,2,1)
 
!     SUBCASE LOGIC -- NSKIP IS 0 (SET BY PARAM MODULE) IF FIRST
!     SUBCASE. OTHERWISE NSKIP IS THE NO. OF RECORDS TO SKIP ON CASE
!     CONTROL DATA BLOCK TO REACH THE LAST SUBCASE. GP4 SETS THE
!     FOLLOWING PARAMETERS -
!     (1) MPCF1 = +1 (DO NOT PURGE OR EQUIV MCE DATA BLOCKS) = -1 (PURGE
!                 AND EQUIV TO TAKE).
!     (2) MPCF2 = +1 (EXECUTE MCE1 AND MCE2) = -1 (DO NOT EXECUTE)
!     (3) REPEAT= +1 (MORE SUBCASES AFTER THIS ONE) = -1 (LAST SUBCASE).
!     (4) NSKIP = NO. OF RECORDS TO SKIP ON CASE CONTROL TO REACH THE
!                 CURRENT SUBCASE (FOR MODULES IN REMAINDER OF LOOP).
 
 20 repeat= -1
 mpcf1 = -1
 mpcf2 = -1
 nskp1 =  1
 FILE  = casecc
 CALL gopen (casecc,z(buf1),0)
 IF (nskip > 1) CALL skprec (casecc,nskip-1)
 CALL fread (casecc,z,36,1)
 IF (nskip > 0) GO TO 30
 
!     FIRST SUBCASE - INITIALIZE.
 
 mpcset = z(iz2)
 spcset = z(iz3)
 nskip  = 1
 GO TO 50
 
!     SUBSEQUENT SUBCASE - POSITION CASE CONTROL AND INITIALIZE.
 
 30 mpcold = z(iz2)
 spcold = z(iz3)
 40 nskip  = nskip + 1
 CALL fread (casecc,z,36,1)
 IF (z(iz16) /= 0) GO TO 40
 IF (z(iz2) == mpcold .AND. z(iz3) == spcold) GO TO 40
 mpcset = z(iz2)
 spcset = z(iz3)
 
!     LOOK AHEAD TO END OF CURRENT SUBCASE AND SET PARAMETERS.
 
 50 CALL READ (*60,*2420,casecc,z,138,1,flag)
 
!     CHECK FOR SYMMETRY
 
 IF (z(iz16) /= 0) GO TO 50
 
!     CHECK FOR BUCKLING OR DIFFERENTIAL STIFFNESS
 
 IF (z(iz5) /= 0 .OR. z(iz138) /= 0) GO TO 60
 IF (z(iz2) == mpcset .AND. z(iz3) == spcset) GO TO 110
 repeat = 1
 
!     CHECK TO SEE IF MPC SET IS SELECTED OR IF RIGID ELEMENTS EXIST
 
 60 IF (mpcset == 0 .AND. rigid == 0) GO TO 70
 mpcf1 = 1
 mpcf2 = 1
 IF (nskip  ==      1) GO TO 70
 IF (mpcset == mpcold) mpcf2 = -1
 70 CALL CLOSE (casecc,clsrew)
 ASSIGN 120 TO ret
 
!     READ EQEXIN INTO CORE
 
 80 FILE = eqexin
 CALL gopen (eqexin,z(buf1),0)
 CALL READ  (*2410,*90,eqexin,z,buf4,1,kn)
 insuff = 80
 icrq   = buf4
 GO TO 2430
 90 CALL READ  (*2410,*2420,eqexin,z(kn+1),kn,1,flag)
 CALL CLOSE (eqexin, clsrew)
 km  = 2*kn
 kn2 = kn/2
 
!     FORM ARRAY OF SORTED SIL VALUES STARTING AT Z(KM+1)
 
 DO  i = 1, kn2
   j = 2*(i-1) + 2 + kn
   z(km+i) = z(j)/10
 END DO
 CALL sort (0,0,1,1,z(km+1),kn2)
 z(km+kn2+1) = luset + 1
 knkl1 = km + kn2 + 2
 
!     SET DIAG-S 21 AND 22 FOR DEGREE-OF-FREEDOM PRINTER LATER.
 
 CALL sswtch (21,l21)
 CALL sswtch (22,l22)
 GO TO ret, (120,1930,1660)
 
 110 nskp1 = nskp1 + 1
 GO TO 50
 
!     OPEN INPUT DATA FILE
 
 120 FILE = geomp
 CALL preloc (*130,z(buf1),geomp)
 nogeom = 1
 
!     CHECK TO SEE IF MPC SET IS SELECTED OR IF RIGID ELEMENTS EXIST
 
 IF (mpcset == 0 .AND. rigid == 0) GO TO 130
 
!     OPEN RGT FILE
 
 FILE = rgt
 CALL gopen (rgt,z(buf3),1)
 
!     IF RIGID ELEMENTS EXIST, GENERATE THEIR COEFFICIENTS
 
 nogoo = nogo
 nogo  = 0
 IF (rigid == 1) CALL criggp (n23)
 IF (nogo  /= 0) GO TO 2540
 nogo = nogoo
 
!     OPEN SCRATCH DATA FILE
 
 130 FILE = scr1
 CALL OPEN (*2400,scr1,z(buf2),wrtrew)
 
!     CHECK TO SEE IF GEOMP FILE EXISTS
 
 IF (nogeom == 0) GO TO 790
 
!     CHECK TO SEE IF MPC SET IS SELECTED OR IF RIGID ELEMENTS EXIST
 
 IF (mpcset == 0 .AND. rigid == 0) GO TO 610
 IF (mpcset /= 0) GO TO 140
 
!     NO MPC SET IS SELECTED
 
 multi = 0
 impc  = knkl1
 i = impc
 j = buf3 - 1
 GO TO 370
 
!     IF MPC SET IS SELECTED, DETERMINE IF SET IS ON MPCADD CARD.
!     IF NOT, SIMULATE AN MPCADD SET LIST WITH ONE SET = MPCSET.
 
 140 impcad = knkl1
 nmpcad = knkl1
 impc   = impcad + 2
 i      = impcad
 z(i)   = mpcset
 z(i+1) = 0
 FILE   = geomp
 CALL locate (*200,z(buf1),mpcadd,flag)
 150 CALL READ (*2410,*200,geomp,id,1,0,flag)
 IF (id == mpcset) GO TO 170
 160 CALL fread (geomp,buf,1,0)
 IF (buf(1) /= -1) GO TO 160
 GO TO 150
 170 CALL READ (*2410,*190,geomp,buf,1,0,flag)
 IF (buf(1) == -1) GO TO 180
 z(i  ) = buf(1)
 z(i+1) = 0
 i  = i + 2
 GO TO 170
 180 CALL fwdrec (*2410,geomp)
 190 impc   = i
 nmpcad = i - 2
 
!     READ MPC CARDS. FOR EACH EQUATION WHOSE SET ID MATCHES A SET ID
!     IN THE MPCADD SET LIST, CONVERT THE GRID POINT AND COMPONENT NO.
!     (OR SCALAR NO.) TO A SIL VALUE. COMPUTE THE ROW AND COLUMN NO.
!     FOR THE POINT AND SAVE THIS ALONG WITH ITS VALUE.
 
 200 CALL locate (*320,z(buf1),mpc,flag)
 j = buf3 - 1
 i = impc
 multi = 0
 ASSIGN 260  TO ret
 ASSIGN 2460 TO ret1
 ASSIGN 250  TO ret2
 ASSIGN 270  TO ret3
 210 CALL READ (*2410,*320,geomp,id,1,0,flag)
 DO  k = impcad,nmpcad,2
   IF (z(k) == id) GO TO 240
 END DO
 230 CALL fread (geomp,buf,3,0)
 IF (buf(1) /= -1) GO TO 230
 GO TO 210
 240 multi = multi + 1
 z(k+1)= 1
 ifl   = 0
 250 CALL fread (geomp,buf,3,0)
 IF (buf(1) == -1) GO TO 310
 gpoint = buf(1)
 GO TO 2100
 260 INDEX = 1
 icomp = buf(2)
 GO TO 2300
 270 IF (icomp /= 0) gpoint = gpoint + icomp - 1
 IF (ifl == 0) sild = gpoint
 IF (n23 == 3) GO TO 300
 IF (gpoint > mask15) GO TO 290
 z(i  ) = orf(lshift(gpoint,ihalf),sild)
 z(i+1) = buf(3)
 280 i = i + n23
 insuff = 236
 IF (i >= j) GO TO 2430
 ifl = 1
 GO TO 250
 
!     GPOINT IS TOO BIG TO BE PACKED INTO HALF A WORD.  ABANDON COL.
!     AND ROW PACKING LOGIC, AND DO IT OVER AGAIN WITHOUT PACKING.
 
 290 n23 = 3
 CALL REWIND (geomp)
 CALL fwdrec (*2410,geomp)
 GO TO 200
 300 z(i  ) = gpoint
 z(i+1) = sild
 z(i+2) = buf(3)
 GO TO 280
 
!     SAVE A LIST OF DEPENDENT SIL VALUES
 
 310 z(j)= sild
 j   = j - 1
 GO TO 210
 
!     DETERMINE IF ALL MPC SETS IN MPCADD SET LIST HAVE BEEN INPUT
 
 320 IF (nogo /= 0) GO TO 2540
 nogoo = nogo
 nogo  = 0
 igotch= 0
 DO  k = impcad,nmpcad,2
   IF (z(k+1) /= 0) GO TO 340
   nogo  = -1
   IF (z(k) == 200000000 .AND. iaxif /= 0) CYCLE
   IF (iaxic == 0) GO TO 330
   IF (z(k) == mpcax1 .OR. z(k) == mpcax2) CYCLE
   IF (z(k) == 200000000) CYCLE
   330 nogo  = +1
   buf(1)= z(k)
   buf(2)= 0
   CALL mesage (30,47,buf)
   CYCLE
   340 igotch= 1
 END DO
 IF (nogo == 0) GO TO 370
 IF (nogo == -1 .AND. igotch == 1) GO TO 360
 mpcset=  0
 multi = -1
 mpcf1 = -1
 mpcf2 = -1
 IF (nogo == -1 .AND. nogoo == 0) nogo = 0
 GO TO 600
 360 CONTINUE
 IF (nogo == -1 .AND. nogoo == 0) nogo = 0
 
!     CHECK TO SEE IF RIGID ELEMENTS EXIST
 
 370 IF (rigid == 0) GO TO 470
 
!     EXPAND THE DEPENDENT SET BY APPENDING RIGID ELEMENT
!     DATA TO MPC DATA
 
 CALL gopen  (rgt,z(buf3),0)
 CALL skprec (rgt,1)
 i1   = buf3 - i
 CALL READ (*2410,*380,rgt,z(i),i1,1,nrigid)
 insuff = 3020
 GO TO 2430
 380 j = j - nrigid
 multi = multi + nrigid
 CALL skprec (rgt,-2)
 CALL READ (*2410,*410,rgt,z(i),i1,1,flag)
 insuff = 3030
 i2 = i1
 390 CALL bckrec (rgt)
 CALL READ (*2410,*400,rgt,z(i),-i2,0,flag)
 CALL READ (*2410,*400,rgt,z(i), i1,0,flag)
 i2 = i2 + i1
 GO TO 390
 400 flag = i2 + flag
 GO TO 440
 
!     RE-CODE COLUMN-ROW PACKED WORD IF NECESSARY FOR DATA JUST BROUGHT
!     IN FROM RIGID ELEMENTS
!     THEN READ THE LAST RECORD FROM RGT
 
 410 IF (n23 == 3) GO TO 430
 i1 = i - 1
 i2 = i1
 i3 = i1 + flag
 420 z(i2+1) = orf(lshift(z(i1+1),ihalf),z(i1+2))
 z(i2+2) = z(i1+3)
 i1 = i1 + 3
 i2 = i2 + 2
 IF (i1 < i3) GO TO 420
 flag = i2 - i + 1
 
 430 insuff = 3050
 440 i3 = i + flag
 IF (i3 < j)  GO TO 460
 WRITE  (outtap,450) i,i3,j,flag,buf3,nrigid,n23
 450 FORMAT ('  GP4/3060 I,I3,J,FLAG,BUF3,NRIGID,N23 =',7I7)
 icrq = i - j
 GO TO 2430
 460 i  = i3
 CALL READ  (*2410,*2420,rgt,z(j+1),nrigid,1,flag)
 CALL CLOSE (rgt,clsrew)
 CALL gopen (rgt,z(buf3),1)
 
!     SORT THE LIST OF DEPENDENT SIL VALUES
!     THUS FORMING THE UM SUBSET
 
 470 ii = j + 1
 m  = buf3 - ii
 nnx= buf3 - 1
 IF (m == 1) GO TO 510
 CALL sort (0,0,1,1,z(ii),m)
 
!     CHECK FOR DEPENDENT COMPONENT ERRORS IN MPC/RIGID ELEMENT DATA
 
 jj   = nnx - 1
 nold = 0
 jxx  = 0
 DO  j = ii,jj
   IF (z(j) == nold) CYCLE
   IF (z(j) /= z(j+1)) CYCLE
   nold = z(j)
   nogo = 1
   jxx  = jxx + 1
   IF (jxx > 50) CYCLE
   CALL page2 (2)
   WRITE  (outtap,480) ufm,z(j)
   480 FORMAT (a23,' 2423, DEPENDENT COMPONENT SPECIFIED MORE THAN ONCE',  &
       ' ON MPC CARDS AND/OR IN RIGID ELEMENTS.  SIL =',i9)
 END DO
 IF (jxx > 50) WRITE (outtap,500)
 500 FORMAT (//12X,12H... AND more,/)
 510 IF (nogo /= 0) GO TO 2540
 CALL WRITE (scr1,z(ii),m,1)
 
!     SORT THE LIST OF CODED COL AND ROW NOS (OR UNCODED NOS)
!     THEN BLDPK EACH COL THUS FORMING THE RG MATRIX
 
 n   = i - impc
 nmpc= i - n23
 j   = impc
 IF (n23 == 3) CALL sort2k (0,0,3,1,z(j),n)
 IF (n23 == 2) CALL sort   (0,0,2,1,z(j),n)
 
!     CHECK FOR INDEPENDENT COMPONENT ERRORS IN MPC DATA
 
 kj   = j + n - 2*n23
 nold = 0
 nogo = 0
 DO  kk = j,kj,n23
   IF (z(kk) ==      nold) CYCLE
   IF (z(kk) /= z(kk+n23)) CYCLE
   IF (n23 == 3 .AND. z(kk+1) /= z(kk+n23+1)) CYCLE
   nold = z(kk)
   nogo = 1
   jj   = nold
   IF (n23 == 2) jj = rshift(nold,ihalf)
   CALL page2 (-2)
   WRITE  (outtap,530) ufm,jj
   530 FORMAT (a23,' 3180, INDEPENDENT COMPONENT SPECIFIED MORE THAN ',  &
       'ONCE IN AN MPC RELATIONSHIP.   SIL =',i6)
 END DO
 IF (nogo /= 0) GO TO 2540
 ncol= 1
 m   = buf3 - i
 n231= n23  - 1
 550 CALL bldpk (1,1,rgt,0,0)
 560 IF (j > nmpc) GO TO 590
 jj = z(j)
 IF (n23 == 2) jj = rshift(z(j),ihalf)
 IF (jj > ncol) GO TO 590
 ix = z(j+1)
 IF (n23 == 2) ix = andf(z(j),mask16)
 x(1) = z(j+n231)
 DO  nn1 = ii,nnx
   IF (ix == z(nn1)) GO TO 580
 END DO
 GO TO 2540
 580 ix = nn1 - ii + 1
 CALL zblpki
 j  = j + n23
 GO TO 560
 590 CALL bldpkn (rgt,0,mcb)
 ncol = ncol + 1
 IF (ncol <= luset) GO TO 550
 mcb(3) = multi
 CALL wrttrl (mcb)
 600 CALL CLOSE (rgt,clsrew)
 
!     READ OMIT CARDS (IF PRESENT).
 
 610 i = knkl1
 CALL locate (*650,z(buf1),omit,flag)
 ASSIGN 630  TO ret
 ASSIGN 2470 TO ret1
 ASSIGN 620  TO ret2
 ASSIGN 640  TO ret3
 omit1 = 1
 620 CALL READ (*2410,*650,geomp,buf,2,0,flag)
 gpoint= buf(1)
 GO TO 2100
 630 INDEX = 3
 icomp = buf(2)
 GO TO 2300
 640 IF (icomp /= 0) gpoint = gpoint + icomp - 1
 z(i)= gpoint
 i   = i + 1
 IF (i <= buf3) GO TO 620
 icrq = i - buf3
 insuff = 345
 GO TO 2430
 
!     READ OMIT1 CARDS (IF PRESENT).
 
 650 IF (nogo /= 0) GO TO 2540
 CALL locate (*720,z(buf1),omitx1,flag)
 omit1 = 1
 ASSIGN 680  TO ret
 ASSIGN 2470 TO ret1
 ASSIGN 670  TO ret2
 ASSIGN 690  TO ret3
 660 CALL READ (*2410,*720,geomp,buf,1,0,flag)
 IF (buf(1) /= 0) CALL scalex (1,buf(1),buf(8))
 670 CALL READ (*2410,*720,geomp,buf(2),1,0,flag)
 IF (buf(2) == -1)  GO TO 660
 gpoint = buf(2)
 GO TO 2100
 680 INDEX = 5
 icomp = buf(1)
 GO TO 2300
 690 IF (icomp /= 0) GO TO 700
 z(i) = gpoint
 i    = i + 1
 GO TO 670
 700 gpoint = gpoint - 1
 DO  ijk = 1,6
   IF (buf(ijk+7) == 0) GO TO 670
   z(i) = gpoint+buf(ijk+7)
   i    = i + 1
 END DO
 GO TO 670
 720 IF (omit1 /= 1) GO TO 730
 IF (nogo  /= 0) GO TO 2540
 
!     SORT OMIT AND OMIT1 DATA AND WRITE IT ON SCR1.
 
 n = i - knkl1
 i = knkl1
 CALL sort (0,0,1,1,z(i),n)
 CALL WRITE (scr1,z(i),n,1)
 
!     READ SUPORT CARDS (IF PRESENT)
 
 730 CALL locate (*780,z(buf1),suport,flag)
 react = 1
 i = knkl1
 ASSIGN 750  TO ret
 ASSIGN 2480 TO ret1
 ASSIGN 740  TO ret2
 ASSIGN 760  TO ret3
 740 CALL READ (*2410,*770,geomp,buf,2,0,flag)
 gpoint = buf(1)
 GO TO 2100
 750 INDEX = 7
 icomp = buf(2)
 GO TO 2300
 760 IF (icomp /= 0) gpoint = gpoint + icomp - 1
 z(i) = gpoint
 i    = i + 1
 IF (i < buf3) GO TO 740
 icrq   = i - buf3
 insuff = 445
 GO TO 2430
 770 IF (nogo /= 0) GO TO 2540
 n = i - knkl1
 i = knkl1
 CALL sort (0,0,1,1,z(i),n)
 CALL WRITE (scr1,z(i),n,1)
 
!     READ THE GPDT AND EXTRACT CONSTRAINED POINTS (IF ANY)
 
 780 CALL CLOSE (geomp,clsrew)
 790 FILE = gpdt
 ASSIGN 810 TO ret
 CALL gopen (gpdt,z(buf1),0)
 800 CALL READ (*2400,*820,gpdt,buf,7,0,flag)
 IF (buf(7) == 0) GO TO 800
 j = buf(1) + km
 buf(1) = z(j)
 CALL scalex (buf,buf(7),buf(8))
 GO TO 2200
 810 CALL WRITE (scr1,buf(8),n,0)
 ugset = 1
 GO TO 800
 820 IF (ugset > 0) CALL WRITE (scr1,0,0,1)
 CALL CLOSE (gpdt,clsrew)
 FILE = geomp
 IF (nogeom == 0) GO TO 830
 CALL preloc (*2400,z(buf1),geomp)
 GO TO 840
 830 IF (mpcset /= 0) CALL mesage (30,47,mpcset)
 IF (spcset /= 0) CALL mesage (30,53,spcset)
 IF (mpcset /= 0 .OR. spcset /= 0) nogo = +1
 GO TO 1280
 
!     IF SPC SET IS SELECTED, READ SPCADD CARDS (IF PRESENT).
!     DETERMINE IF SET ID IS ON SPCADD CARD.
!     IF NOT, SIMULATE AN SPCADD SET LIST WITH ONE SET = SPCSET.
 
 840 IF (spcset == 0) GO TO 1150
 ispcad = knkl1
 nspcad = knkl1
 ispc   = ispcad + 2
 i      = ispcad
 z(i  ) = spcset
 z(i+1) = 0
 CALL locate (*900,z(buf1),spcadd,flag)
 850 CALL READ (*2410,*900,geomp,id,1,0,flag)
 IF (id == spcset) GO TO 870
 860 CALL fread (geomp,id,1,0)
 IF (id /= -1) GO TO 860
 GO TO 850
 870 CALL READ (*2410,*890,geomp,buf,1,0,flag)
 IF (buf(1) == -1) GO TO 880
 z(i  ) = buf(1)
 z(i+1) = 0
 i      = i + 2
 GO TO 870
 880 CALL fwdrec (*2410,geomp)
 890 ispc   = i
 nspcad = i - 2
 
!     READ SPC1 AND SPC CARDS.
!     FOR EACH SET ID WHICH IS IN THE SPCADD SET LIST,
!     CONVERT THE GRID POINT NO. AND COMPONENT VALUE (OR SCALAR NO.)
!     TO AN SIL VALUE. SAVE A LIST IN CORE OF SIL VALUES AND
!     ENFORCED DISPLACEMENT (ON SPC1 CARDS, ENF. DISPL. = 0.)
 
 900 i = ispc
 GO TO 1010
 
!     SPC1 PROCESSING EXECUTES AFTER SPC PROCESSING
 
 910 IF (nogo /= 0) GO TO 2540
 CALL locate (*1130,z(buf1),spc1,flag)
 ASSIGN 970  TO ret
 ASSIGN 2490 TO ret1
 ASSIGN 960  TO ret2
 ASSIGN 980  TO ret3
 920 CALL READ (*2410,*1130,geomp,id,1,0,flag)
 DO  k = ispcad,nspcad,2
   IF (z(k) == id) GO TO 950
 END DO
 940 CALL fread (geomp,buf,1,0)
 IF (buf(1) /= -1) GO TO 940
 GO TO 920
 950 z(k+1) = 1
 CALL fread (geomp,buf,1,0)
 single = 1
 IF (buf(1) /= 0) CALL scalex (1,buf(1),buf(8))
 960 CALL READ (*2410,*920,geomp,buf(2),1,0,flag)
 IF (buf(2) < 0) GO TO 920
 gpoint = buf(2)
 GO TO 2100
 970 INDEX = 9
 icomp = buf(1)
 GO TO 2300
 980 IF (icomp /= 0) GO TO 990
 z(i  ) = gpoint
 z(i+1) = 0
 i      = i + 2
 GO TO 960
 990 gpoint = gpoint - 1
 DO  ijk = 1,6
   IF (buf(ijk+7) == 0) GO TO 960
   z(i  ) = gpoint+buf(ijk+7)
   z(i+1) = 0
   i      = i + 2
 END DO
 GO TO 960
 
!     PROCESSING OF SPC CARDS EXECUTES FIRST.
 
 1010 CALL locate (*910,z(buf1),spc,flag)
 ASSIGN 1050  TO ret
 ASSIGN 2530 TO ret1
 ASSIGN 1020  TO ret2
 ASSIGN 1060  TO ret3
 1020 CALL READ (*2410,*1090,geomp,buf,4,0,flag)
 DO  k = ispcad,nspcad,2
   IF (z(k) == buf(1)) GO TO 1040
 END DO
 GO TO 1020
 1040 single = 1
 z(k+1) = 1
 gpoint = buf(2)
 GO TO 2100
 1050 INDEX = 11
 icomp = buf(3)
 GO TO 2300
 1060 IF (icomp /= 0) GO TO 1070
 z(i  ) = gpoint
 z(i+1) = buf(4)
 i      = i+2
 GO TO 1020
 1070 CALL scalex (gpoint,buf(3),buf(8))
 DO  ijk = 1,6
   IF (buf(ijk+7) == 0) GO TO 1020
   z(i  ) = buf(ijk+7)
   z(i+1) = buf(4)
   i      = i + 2
 END DO
 GO TO 1020
 1090 IF (nogo /= 0) GO TO 2540
 n = i - ispc
 IF (n <= 2) GO TO 910
 
!     CHECK FOR DUPLICATELY DEFINED ENFORCED DISPLACEMENTS ON SPC CARDS
 
 CALL sort (0,0,2,1,z(ispc),n)
 n    = n - 2
 nold = 0
 DO  k = 1,n,2
   IF (z(ispc+k-1) == nold) CYCLE
   IF (z(ispc+k-1) /= z(ispc+k+1)) CYCLE
   IF (z(ispc+k) == 0 .AND. z(ispc+k+2) == 0) CYCLE
   nold = z(ispc+k-1)
   nogo = 1
   CALL page2 (3)
   WRITE  (outtap,1100) ufm,nold
   1100 FORMAT (a23,' 3147, ENFORCED DISPLACEMENT ON SPC CARDS SPECIFIED',  &
       ' MORE THAN ONCE', /5X,'FOR THE SAME COMPONENT.  SIL VALUE =' ,    i10)
 END DO
 IF (nogo /= 0) GO TO 2540
 GO TO 910
 
!     FLUID PROBLEM AND NO SPC-S AT ALL.
 
 1120 spcset = 0
 GO TO 840
 1130 nspc = i - 2
 icrq = nspc - buf3
 insuff = 740
 IF (icrq > 0) GO TO 2430
 
!     DETERMINE IF ALL SPC SETS IN SPCADD SET LIST HAVE BEEN DEFINED
 
 IF (nogo /= 0) GO TO 2540
 DO  k = ispcad,nspcad,2
   IF (z(k+1) /= 0) CYCLE
   IF (iaxif /= 0 .AND. z(k) == 200000000) GO TO 1120
   nogo   = 1
   buf(1) = z(k)
   buf(2) = 0
   CALL mesage (30,53,buf)
 END DO
 IF (nogo /= 0) GO TO 2540
 
!     SORT THE SPC LIST AND WRITE IT ON SCR1
 
 n = nspc - ispc + 2
 CALL sort  (0,0,2,1,z(ispc),n)
 CALL WRITE (scr1,z(ispc),n,1)
 
!     READ ASET CARDS (IF PRESENT)
 
 1150 i = knkl1
 CALL locate (*1190,z(buf1),aset,flag)
 ASSIGN 1170 TO ret
 ASSIGN 2470 TO ret1
 ASSIGN 1160 TO ret2
 ASSIGN 1180 TO ret3
 asetx = 1
 1160 CALL READ (*2410,*1190,geomp,buf,2,0,flag)
 gpoint = buf(1)
 GO TO 2100
 1170 INDEX = 15
 icomp = buf(2)
 GO TO 2300
 1180 IF (icomp /= 0) gpoint = gpoint + icomp - 1
 z(i) = gpoint
 i    = i + 1
 IF (i <= buf3) GO TO 1160
 icrq = i - buf3
 insuff = 1445
 GO TO 2430
 
!     READ ASET1 CARDS (IF PRESENT)
 
 1190 IF (nogo /= 0) GO TO 2540
 CALL locate (*1260,z(buf1),aset1,flag)
 asetx = 1
 ASSIGN 1220 TO ret
 ASSIGN 2470 TO ret1
 ASSIGN 1210 TO ret2
 ASSIGN 1230 TO ret3
 1200 CALL READ (*2410,*1260,geomp,buf,1,0,flag)
 IF (buf(1) /= 0) CALL scalex (1,buf(1),buf(8))
 1210 CALL READ (*2410,*1260,geomp,buf(2),1,0,flag)
 IF (buf(2) == -1) GO TO 1200
 gpoint = buf(2)
 GO TO 2100
 1220 INDEX = 17
 icomp = buf(1)
 GO TO 2300
 1230 IF (icomp /= 0) GO TO 1240
 z(i) = gpoint
 i    = i + 1
 GO TO 1210
 1240 gpoint = gpoint - 1
 DO  ijk = 1,6
   IF (buf(ijk+7) == 0) GO TO 1210
   z(i) = gpoint + buf(ijk+7)
   i    = i + 1
 END DO
 GO TO  1210
 1260 IF (asetx /= 1) GO TO 1270
 IF (nogo  /= 0) GO TO 2540
 
!     SORT ASET AND ASET1 DATA AND WRITE IT ON SCR1
 
 n = i - knkl1
 i = knkl1
 CALL sort  (0,0,1,1,z(i),n)
 CALL WRITE (scr1,z(i),n,1)
 1270 CALL CLOSE (geomp,clsrew)
 1280 CALL CLOSE (scr1,clsrew)
 
!     FORM THE BASIC USET BY READING EACH OF THE SUBSETS AND
!     TURNING ON THE APPROPRIATE BIT IN THE APPROPRIATE WORD
 
 FILE = scr1
 CALL OPEN (*2400,scr1,z(buf2),rdrew)
 DO  k = 1,luset
   z(k)   = 0
 END DO
 buf(1) = multi
 buf(2) = omit1
 buf(3) = react
 buf(4) = usgset
 buf(5) = single
 buf(6) = asetx
 icount = 0
 DO  k = 1,6
   IF (buf(k) < 0) CYCLE
   IF (k < 5) icount = icount + 1
   SELECT CASE ( k )
     CASE (    1)
       GO TO 1300
     CASE (    2)
       GO TO 1310
     CASE (    3)
       GO TO 1300
     CASE (    4)
       GO TO 1300
     CASE (    5)
       GO TO 1300
     CASE (    6)
       GO TO 1310
   END SELECT
   1300 mcbust(5) = orf(mcbust(5),mask(k))
   nosets = 1
   IF (k == 5) GO TO 1350
   1310 CALL READ (*2410,*1360,scr1,j,1,0,flag)
   IF (k == 2) GO TO 1340
   IF (k == 6) GO TO 1330
   IF (andf(z(j),mask(k)) /= mask(k)) GO TO 1340
   dup = 1
   IF (iflag /= 0) GO TO 1320
   FILE = uset
   CALL OPEN (*2400,uset,z(buf1),wrtrew)
   iflag  = 1
   FILE   = scr1
   1320 buf(1) = j
   buf(2) = k
   CALL WRITE (uset,buf(1),2,0)
   GO TO 1340
   1330 IF (andf(z(j),mskua) /= 0) GO TO 1310
   1340 z(j) = orf(z(j),mask(k))
   GO TO 1310
   1350 CALL READ (*2410,*1360,scr1,buf(7),2,0,flag)
   j    = buf(7)
   z(j) = orf(z(j),mask(k))
   GO TO 1350
 END DO
 IF (dup == 0) GO TO 1370
 CALL WRITE (uset,0,0,1)
 CALL CLOSE (uset,clsrew)
 1370 CALL CLOSE (scr1,clsrew)
 
!     THE FOLLOWING CONVENTION WILL BE USED WITH REGARD TO DEGREES OF
!     FREEDOM NOT SPECIFICALLY INCLUDED OR OMITTED-
!       1. IF ASET OR ASET1 CARDS ARE PRESENT, UNSPECIFIED DEGREES OF
!          FREEDOM WILL BE OMITTED.
!       2. IF ASET OR ASET1 CARDS ARE NOT PRESENT AND OMIT OR OMIT1
!          CARDS ARE PRESENT, UNSPECIFIED DEGREES OF FREEDOM WILL BE
!          INCLUDED IN THE ANALYSIS SET.
!       3. IF NO ASET, ASET1, OMIT, OR OMIT 1 CARDS ARE PRESENT ALL
!          UNSPECIFIED DEGREES OF FREEDOM WILL BE INCLUDED IN THE
!          ANALYSIS SET.
!       4. IF BOTH ASET OR ASET1 CARDS AND OMIT OR OMIT1 CARDS ARE
!          SUPPLIED, UNSPECIFIED DEGREES OF FREEDOM WILL BE OMITTED.
 
 mskrst = mask(2)
 IF (asetx > 0) GO TO 1380
 mskrst = mask(6)
 imsk   = 0
 1380 DO  k = 1, luset
   IF (andf(mskck,z(k)) /= 0) CYCLE
   imsk = mskrst
   z(k) = orf(z(k),mskrst)
 END DO
 IF (imsk == mask(6)) asetx = 1
 IF (imsk == mask(2)) omit1 = 1
 
!     CALL SUBROUTINE GP4SP TO EXAMINE GRID POINT SINGULARITIES
 
 CALL gp4sp (buf2,buf3,buf4)
 
!     TURN ON CERTAIN FLAGS IF THERE ARE OMIT OR ASET
!     DEGREES OF FREEDOM
 
 omit1 = -1
 DO  k = 1,luset
   IF (andf(z(k),mskuo) == 0) CYCLE
   mcbust(5) = orf(mcbust(5),mask(2))
   nosets = 1
   omit1  = 1
   EXIT
 END DO
 1410 DO  k = 1,luset
   IF (andf(z(k),mskua) == 0) CYCLE
   mcbust(5) = orf(mcbust(5),mask(6))
   nol = 1
   GO TO 1430
 END DO
 
 1430 CALL OPEN (*2400,scr1,z(buf2),rdrew)
 CALL skprec (scr1,icount)
 
!     OPEN YS FILE. WRITE SPCSET IN YS HEADER.
!     IF NO USB SET (FROM SPC AND SPC1 CARDS), WRITE NULL COLUMN
!     FOR YS VECTOR. IF USB SET IS PRESENT, BUILD THE YS VECTOR.
 
 FILE = scr1
 CALL OPEN (*1440,ys,z(buf3),wrtrew)
 noys = 1
 CALL fname (ys,buf)
 buf(3) = spcset
 CALL WRITE (ys,buf,3,1)
 1440 ix = 0
 ii = 1
 IF (single > 0) GO TO 1450
 IF (nauto > 0 .OR. usgset > 0) single = 1
 IF (noys /= 0) CALL bldpk (1,1,ys,0,0)
 GO TO 1490
 1450 IF (noys /= 0) CALL bldpk (1,1,ys,0,0)
 1460 CALL READ (*2410,*1490,scr1,buf,2,0,flag)
 j = buf(1)
 IF (buf(2) == 0) GO TO 1460
 DO  k = ii,j
   IF (andf(z(k),mskus) /= 0) ix = ix + 1
 END DO
 ii   = j + 1
 x(1) = buf(2)
 IF (noys   /= 0) GO TO 1480
 IF (nogoof /= 0) GO TO 1460
 nogo   = 1
 nogoof = 1
 CALL mesage (30,132,buf)
 GO TO 1460
 1480 CALL zblpki
 GO TO 1460
 1490 IF (noys /= 0) CALL bldpkn (ys,0,mcbys)
 IF (ii > luset) GO TO 1510
 DO  k = ii,luset
   IF (andf(z(k),mskus) /= 0) ix = ix + 1
 END DO
 1510 mcbys(3) = ix
 IF (noys == 0) GO TO 1520
 CALL wrttrl (mcbys)
 CALL CLOSE (ys,clsrew)
 1520 CALL CLOSE (scr1,clsrew)
 
 IF (l21+l22 > 0 .OR. idsub > 0) CALL gp4prt (buf1)
 IF (nauto == 0) GO TO 1540
 
!     CHANGE AUTO SPC FLAGS TO BOUNDARY SPC FLAGS
 
 j = 0
 DO  k = 1,luset
   IF (andf(z(k),mskus) == 0) CYCLE
   IF (andf(z(k),mskusg) /= 0 .OR.  andf(z(k),mskusb) /= 0) CYCLE
   z(k) = mask(5)
   j = 1
 END DO
 IF (j == 1) mcbust(5) = orf(mcbust(5),mask(5))
 
 1540 FILE = uset
 IF (dup == 0) GO TO 1570
 CALL OPEN (*2400,uset,z(buf1),rdrew)
 FILE = scr1
 CALL OPEN (*2400,scr1,z(buf2),wrtrew)
 FILE = uset
 1550 CALL READ  (*1560,*1560,uset,buf(1),2,0,flag)
 CALL WRITE (scr1,buf(1),2,0)
 GO TO 1550
 1560 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (uset,clsrew)
 1570 CALL OPEN  (*2400,uset,z(buf1),wrtrew)
 CALL fname (uset,buf)
 buf(3) = spcset
 buf(4) = mpcset
 CALL WRITE (uset,buf,4,1)
 CALL WRITE (uset,z(1),luset,1)
 IF (nol == 1) mcbust(5)= orf(mcbust(5),mask(6))
 
!     SEPARATE TRAILER WORD 4 INTO TWO PARTS
 
 mcbust(4) = rshift(mcbust(5),ihalf)
 mcbust(5) = andf(mcbust(5),complf(lshift(mcbust(4),ihalf)))
 CALL wrttrl (mcbust)
 CALL CLOSE (uset,clsrew)
 
!     PROCESS USET FOR CONSISTENCY OF DISPLACEMENT SET DEFINITIONS.
!     EACH POINT IN USET MAY BELONG TO AT MOST ONE DEPENDENT SUBSET.
 
 flag    = 0
 mask(1) = mskum
 mask(2) = mskus
 mask(3) = mskuo
 mask(4) = mskur
 mskums  = orf(mskum,mskus)
 mskuor  = orf(mskuo,mskur)
 buf( 1) = orf(mskus,mskuor)
 buf( 2) = orf(mskum,mskuor)
 buf( 3) = orf(mskur,mskums)
 buf(4)  = orf(mskuo,mskums)
 mskall  = orf(mskums,mskuor)
 mskal   = orf(mskall,mskul)
 DO  i = 1,luset
   iuset = z(i)
   idepn = andf(mskal,iuset)
   DO  ik = 1,4
     IF (andf(mak(ik),idepn) == mak(ik)) GO TO 1600
   END DO
   idepn = andf(iuset,mskall)
   IF (idepn == 0) CYCLE
   DO  j = 1,4
     msk1 = mask(j)
     msk2 = buf( j)
     IF (andf(idepn,msk1) == 0) CYCLE
     IF (andf(idepn,msk2) /= 0) GO TO 1600
   END DO
   CYCLE
   1600 IF (flag /= 0 .OR. iflag /= 0) GO TO 1610
   FILE = scr1
   CALL OPEN (*2400,scr1,z(buf1),wrtrew)
   1610 buf(5) = i
   buf(6) = idepn
   flag   = 1
   CALL WRITE (scr1,buf(5),2,0)
 END DO
 1630 IF (mpcf1 > 0 .OR. single > 0 .OR. omit1 > 0 .OR. react > 0) nosets = 1
 IF (mpcf1 == -1 .AND. single == -1 .AND. omit1 == -1) noa = -1
 IF (andf(mskua,mcbust(5)) /= 0 .OR. omit1 < 0) GO TO 1650
 CALL page2 (2)
 WRITE  (outtap,1640) ufm
 1640 FORMAT (a23,' 2403, INVALID TO HAVE AN O-SET WITH A NULL A-SET.')
 nogo = 1
 1650 CONTINUE
 IF (nogo /= 0) GO TO 2540
 IF (iflag /= 0 .OR. flag /= 0) GO TO 1920
 
!     RECOMPUTE YS MATRIX TO ACCOUNT FOR SPCD CARDS
 
 
 IF (noys == 0 .OR . nogeom == 0) GO TO 1910
!     BRING EQEXIN,SIL,AND USET BACK INTO CORE
 
 ASSIGN 1660 TO ret
 GO TO 80
 1660 CALL gopen (uset,z(buf1),0)
 FILE = uset
 CALL READ (*2410,*1670,uset,z(knkl1),buf4-knkl1,1,luset)
 icrq = buf4
 insuff = 9711
 GO TO 2430
 1670 CALL CLOSE (uset,1)
 
!     CONVERT USET POINTERS INTO SILA VALUES
 
 m  = knkl1
 n  = knkl1 + luset - 1
 ix = 0
 DO  i = m,n
   IF (andf(z(i),mskus) /= 0) GO TO 1680
   z(i) = 0
   CYCLE
   1680 ix  = ix + 1
   z(i)= ix
 END DO
 
!     POSITION CASECC
 
 FILE  = casecc
 iload = n + 1
 icrq  = n + 2*nskp1 + 1 - buf4
 insuff = 977
 IF (icrq > 0) GO TO 2430
 CALL gopen  (casecc,z(buf1),0)
 CALL skprec (casecc,nskip-1)
 DO  i = 1,nskp1
   1700 CALL fread (casecc,buf,16,1)
   IF (buf(16) /= 0) GO TO 1700
   k      = iload + 2*(i-1)
   z(k  ) = buf(4)
   z(k+1) = 0
 END DO
 CALL CLOSE (casecc,clsrew)
 
!     CONVERT SPCD CARD TO SILA + VALUE AND WRITE ON SCR2
 
 CALL gopen (scr2,z(buf2),1)
 FILE = geomp
 CALL preloc (*2400,z(buf1),geomp)
 CALL locate (*1830,z(buf1),spcd,flag)
 nn    = 2*nskp1 + iload - 2
 iold  = 0
 irecn = 0
 1720 CALL READ (*2410,*1820,geomp,buf,4,0,flag)
 DO  i = iload,nn,2
   IF (buf(1) == z(i)) GO TO 1740
 END DO
 
!     GO ON TO NEXT SET
 
 GO TO 1720
 
 1740 IF (buf(1) == iold) GO TO 1760
 IF (iold /= 0) CALL WRITE (scr2,0,0,1)
 iold  = buf(1)
 irecn = irecn + 1
 DO  i = iload,nn,2
   IF (iold == z(i)) z(i+1) = irecn
 END DO
 1760 gpoint = buf(2)
 ASSIGN 1770  TO ret
 ASSIGN 2530 TO ret1
 ASSIGN 1720  TO ret2
 ASSIGN 1780 TO ret3
 GO TO 2100
 
!     FOUND SIL
 
 1770 INDEX = 13
 icomp = buf(3)
 GO TO 2300
 1780 IF (icomp /= 0) GO TO 1790
 m = knkl1 + gpoint - 1
 IF (z(m) == 0) GO TO 1810
 mcb(1) = z(m)
 mcb(2) = buf(4)
 CALL WRITE (scr2,mcb,2,0)
 GO TO 1720
 
!     BREAK UP COMPONENTS
 
 1790 CALL scalex (gpoint,buf(3),buf(8))
 DO  i = 1,6
   IF (buf(i+7) == 0) GO TO 1720
   m = knkl1 + buf(i+7) - 1
   IF (z(m) == 0) GO TO 1810
   mcb(1) = z(m)
   mcb(2) = buf(4)
   CALL WRITE (scr2,mcb,2,0)
 END DO
 GO TO 1720
 1810 n      = 108
 buf(1) = buf(2)
 buf(2) = buf(i+7) - gpoint
 GO TO 2520
 
!     END OF SPCD-S
 
 1820 IF (nogo /= 0) GO TO 2540
 CALL WRITE (scr2,0,0,1)
 1830 CALL CLOSE (geomp,1)
 CALL CLOSE (scr2,1)
 IF (single < 0) GO TO 1910
 
!     BRING IN OLD YS
 
 n = 2*nskp1
 DO  i = 1,n
   k = iload + i - 1
   z(i) = z(k)
 END DO
 ioys = n
 inys = ioys + ix
 icrq = inys + ix - buf4
 insuff = 988
 IF (icrq > 0) GO TO 2430
 mcb(1) = ys
 CALL rdtrl (mcb)
 mcb(2) = 0
 mcb(6) = 0
 mcb(7) = 0
 CALL gopen (ys,z(buf1),0)
 itb  = mcb(5)
 ita1 = itb
 itb1 = itb
 incr = 1
 incr1= 1
 ii   = 1
 ii1  = 1
 jj   = mcb(3)
 jj1  = jj
 DO  i = 1,ix
   rz(ioys+i) = 0.0
 END DO
 CALL unpack (*1860,ys,rz(ioys+1))
 1860 CALL CLOSE (ys,clsrew)
 CALL gopen (ys,z(buf1),1)
 CALL gopen (scr2,z(buf2),0)
 FILE = scr2
 DO  i = 1,n,2
   
!     COPY OLD YS TO NEW YS
   
   DO  k = 1,ix
     rz(inys+k) = rz(ioys+k)
   END DO
   IF (z(i+1) == 0) GO TO 1890
   
!     POSITION SCR2
   
   CALL skprec (scr2,z(i+1)-1)
   1880 CALL READ (*2410,*1890,scr2,buf,2,0,flag)
   k = buf(1) + inys
   rz(k) = bufr(2)
   GO TO 1880
   
!     PUT OUT COLUMN
   
   1890 CALL pack (rz(inys+1),ys,mcb)
   CALL REWIND (scr2)
   CALL fwdrec (*2410,scr2)
 END DO
 CALL CLOSE  (ys,1)
 CALL wrttrl (mcb)
 CALL CLOSE  (scr2,1)
 1910 IF (nogo /= 0) GO TO 2540
 IF (flag /= 0) GO TO 1920
 IF (iogpst == 1) CALL mesage (17,iautsp,0)
 RETURN
 
!     INCONSISTENT DISPLACEMENT SET DEFINITIONS--
!     READ EQEXIN AND SIL INTO CORE. FOR EACH INCONSISTANT DEFINITION,
!     LOOK UP EXTERNAL NUMBER AND QUEUE MESSAGE.
 
 1920 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (scr1,clsrew)
 ASSIGN 1930 TO ret
 GO TO 80
 1930 CALL OPEN (*2400,scr1,z(buf1),rdrew)
 isil = km + 1
 neqx = kn - 1
 z(knkl1) = luset + 1
 1940 CALL READ (*2080,*2080,scr1,buf(5),2,0,iflg)
 DO  i = isil,knkl1
   IF (z(i+1) > buf(5)) EXIT
 END DO
 1960 intrnl = i - km
 komp = buf(5) - z(i) + 1
 IF (z(i+1)-z(i) == 1) komp = 0
 DO  j = 1,neqx,2
   IF (z(j+1) == intrnl) EXIT
 END DO
 1980 IF (dup == 0) GO TO 2070
 IF (iflag == 0) GO TO 2070
 CALL page2 (2)
 SELECT CASE ( ib6 )
   CASE (    1)
     GO TO 1990
   CASE (    2)
     GO TO 1940
   CASE (    3)
     GO TO 2010
   CASE (    4)
     GO TO 2030
 END SELECT
 1990 IF (komp == 0) GO TO 2000
 nogo = 1
 WRITE (outtap,2050) ufm,z(j),komp,mset
 GO TO 1940
 2000 WRITE (outtap,2060) ufm,z(j),mset
 nogo = 1
 GO TO 1940
 2010 IF (komp == 0) GO TO 2020
 WRITE (outtap,2050) ufm,z(j),komp,r
 nogo = 1
 GO TO 1940
 2020 WRITE (outtap,2060) ufm,z(j),r
 nogo = 1
 GO TO 1940
 2030 IF (komp == 0) GO TO 2040
 WRITE (outtap,2050) ufm,z(j),komp,sg
 nogo = 1
 GO TO 1940
 2040 WRITE (outtap,2060) ufm,z(j),sg
 nogo = 1
 GO TO 1940
 2050 FORMAT (a23,' 2152, GRID POINT',i9,' COMPONENT',i3,  &
     ' DUPLICATELY DEFINED IN THE ',a4,5H set.)
 2060 FORMAT (a23,' 2153, SCALAR POINT',i9,' DUPLICATELY DEFINED IN ',  &
     'THE ',a4,5H set.)
 2070 buf(7) = z(j)
 buf(8) = komp
 IF (andf(buf(6),mskum) /= 0) buf(8)= buf(8) + 10
 IF (andf(buf(6),mskus) /= 0) buf(8)= buf(8) + 100
 IF (andf(buf(6),mskuo) /= 0) buf(8)= buf(8) + 1000
 IF (andf(buf(6),mskur) /= 0) buf(8)= buf(8) + 10000
 IF (andf(buf(6),mskul) /= 0) buf(8)= buf(8) + 100000
 CALL mesage (30,101,buf(7))
 GO TO 1940
 2080 IF (dup   == 0) GO TO 2090
 IF (iflag == 0) GO TO 2090
 iflag = 0
 IF (flag /= 0) GO TO 1940
 CALL CLOSE (scr1,clsrew)
 GO TO 1630
 2090 CALL CLOSE (scr1,clsrew)
 GO TO 2540
 
 
!     INTERNAL SUBROUTINE TO PERFORM BINARY SEARCH IN EQEXIN
!     AND CONVERT THE EXTERNAL NUMBER TO A SIL VALUE AND A
!     CORRESPONDING TYPE CODE
 
 2100 klo = 0
 khi = kn2
 lastk = 0
 2110 k = (klo+khi+1)/2
 IF (lastk == k) GO TO 2150
 lastk = k
 IF (gpoint-z(2*k-1) < 0.0) THEN
   GO TO  2120
 ELSE IF (gpoint-z(2*k-1) == 0.0) THEN
   GO TO  2140
 ELSE
   GO TO  2130
 END IF
 2120 khi = k
 GO TO 2110
 2130 klo = k
 GO TO 2110
 2140 k = 2*k + kn
 ipoint = gpoint
 gpoint = z(k)/10
 icode  = z(k) - 10*gpoint
 GO TO ret, (260,630,680,750,970,1050,1770,1170,1220)
 2150 GO TO ret1, (2460,2470,2480,2490,2530)
 
 
!     INTERNAL SUBROUTINE TO SORT THE SCALAR COMPONENTS
 
 2200 DO  ii = 1,6
   IF (buf(ii+7) == 0) GO TO 2220
 END DO
 ii = 7
 2220 n  = ii - 1
 IF (n == 0) GO TO ret, (810)
 DO  ii = 1,n
   ijk = luset + 1
   DO  jj = ii,n
     IF (buf(jj+7) >= ijk) CYCLE
     ijk = buf(jj+7)
     jjx = jj
   END DO
   buf(jjx+7) = buf(ii+7)
   buf(ii +7) = ijk
 END DO
 GO TO ret, (810)
 
!     CHECK TO SEE IF GRID AND SCALAR POINTS HAVE BEEN PROPERLY USED
!     ON CONSTRAINT CARDS
 
 2300 IF (icode == 2) GO TO 2320
 
!     GRID POINTS ARE CHECKED HERE
 
 IF (icomp > 0) GO TO 2350
 nogo = 1
 CALL page2 (2)
 WRITE  (outtap,2310) ufm,ipoint,ctype(INDEX),ctype(INDEX+1)
 2310 FORMAT (a23,' 3145, COMPONENT 0 (OR BLANK) SPECIFIED FOR GRID ',  &
     'POINT',i9,4H on ,2A4,6HCARDS.)
 GO TO 2340
 
!     SCALAR POINTS ARE CHECKED HERE
 
 2320 IF (icomp <= 1) GO TO 2350
 nogo = 1
 CALL page2 (2)
 WRITE  (outtap,2330) ufm,ipoint,ctype(INDEX),ctype(INDEX+1)
 2330 FORMAT (a23,' 3146, ILLEGAL COMPONENT SPECIFIED FOR SCALAR POINT',  &
     i9,4H on ,2A4,6HCARDS.)
 2340 GO TO ret2, (250,620,670,740,960,1020,1720,1160,1210)
 2350 GO TO ret3, (270,640,690,760,980,1060,1780,1180,1230)
 
 
!     FATAL ERROR MESSAGES
 
 2400 j = -1
 GO TO 2450
 2410 j = -2
 GO TO 2450
 2420 j = -3
 GO TO 2450
 2430 j = -8
 WRITE  (outtap,2440) insuff
 2440 FORMAT (/33X,'GP4 INSUFFICIENT CORE AT ',i5)
 FILE = icrq
 2450 CALL mesage (j,FILE,NAME)
 2460 buf(1) = gpoint
 buf(2) = mpcset
 n = 48
 gpoint = 1
 GO TO 2520
 2470 buf(1) = gpoint
 gpoint = 1
 n = 49
 GO TO 2510
 2480 buf(1) = gpoint
 gpoint = 1
 n = 50
 GO TO 2510
 2490 n = 51
 2500 buf(1) = gpoint
 buf(2) = spcset
 gpoint = 1
!WKBNB 3/95 NCL94002
 IF ( l51 == 0 ) GO TO 2520
 WRITE ( outtap, 9001 ) uwm, 2051, buf(1), spcset
 9001  FORMAT( a25,i5,' UNDEFINED GRID POINT ',i6,' IN SINGLE-POINT'  &
     ,' CONSTRAINT SET ',i8)
 GO TO 2521
!WKBNE 3/95 NCL94002
 2510 buf(2) = 0
 2520 nogo   = 1
 CALL mesage (30,n,buf)
!WKBI  3/95 NCL94002
 2521 CONTINUE
 GO TO ret2, (250,620,670,740,960,1020,1720,1160,1210)
 2530 n = 52
 GO TO 2500
 2540 IF (l21+l22 > 0 .OR. idsub > 0) CALL gp4prt (-buf4)
 j = -37
 GO TO 2450
END SUBROUTINE gp4
