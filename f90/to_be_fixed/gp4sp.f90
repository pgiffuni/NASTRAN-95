SUBROUTINE gp4sp (ibuf1,ibuf2,ibuf3)
     
!     ROUTINE TO LOOK AT GPST TO ELIMINATE SINGULARITIES
 
 
 
 INTEGER, INTENT(IN OUT)                  :: ibuf1
 INTEGER, INTENT(IN)                      :: ibuf2
 INTEGER, INTENT(IN)                      :: ibuf3
 EXTERNAL        andf  ,orf   ,complf,lshift
 INTEGER :: andf  ,orf   ,complf,eqexin,gpst  ,ogpst ,scr2  ,  &
     mcb(7),omit1 ,spcset,ogpst1(10)
 DIMENSION       iponts(9)    ,jponts(9)    ,indxms(9)    ,  &
     iexcld(9)    ,isubnm(2)    ,iword(8)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm
 COMMON /BLANK / luset ,mpcf1 ,mpcf2 ,single,omit1 ,react ,nskip ,  &
     repeat,nosets,nol   ,noa   ,idsub ,iautsp
 COMMON /gp4fil/ dum(3),irgt
 COMMON /gp4spx/ mskum ,mskuo ,mskur ,mskus ,mskul ,  &
     msksng,spcset,mpcset,nauto ,iogpst
 COMMON /output/ head(1)
 COMMON /system/ isysbf,ioutpt,jdum(6),nlpp ,kdum(2),line, dd(78),ipunch
 COMMON /unpakx/ itypot,iiii  ,jjjj  ,incr
 COMMON /zzzzzz/ iz(1)
 DATA    eqexin, gpst, ogpst, scr2 /103, 107, 205, 302/
 DATA    isubnm /4HGP4S, 4HP       /
 DATA    ncard,  ierror / 2*0      /
 DATA    iscr2,  ieqexn / 2*-1     /
 
 
 INDEX = IABS (iautsp)
 IF (INDEX > 2) GO TO 730
 IF (INDEX == 2 .AND. omit1 > 0) INDEX = 3
 mskms = orf(mskum,mskus)
 IF (iautsp == 0) GO TO 6
 multi = 0
 IF (mpcf1 == -1) GO TO 6
 mskin = lshift(1,12)
 mskxx = complf(mskin)
 CALL gopen (irgt,iz(ibuf1),0)
 itypot = 1
 iiii   = 1
 jjjj   = 1
 incr   = 1
 DO  i = 1,luset
   CALL unpack (*4,irgt,idum)
   IF (andf(iz(i),mskms) /= 0) CYCLE
   multi = 1
   iz(i) = orf(iz(i),mskin)
 END DO
 CALL CLOSE (irgt,1)
 6 CONTINUE
 CALL OPEN (*610,gpst,iz(ibuf1),0)
 ifile = gpst
 CALL fwdrec (*700,gpst)
 
 10 CALL READ (*480,*480,gpst,iordr,1,0,iflag)
 initl = iordr
 ims   = 0
 ispc  = 0
 DO  i = 1, 9
   indxms(i) = 0
   iexcld(i) = 0
   iponts(i) = 0
   jponts(i) = 0
 END DO
 CALL fread (gpst,npts,1,0)
 CALL fread (gpst,iponts,npts,0)
 ibase = iponts(1)
 
!     SET VARIOUS FLAGS FOR THE SINGULARITIES
 
 DO  i = 1,npts
   ii = iponts(i)
   j  = iz(ii)
   IF (andf(j,mskms) /= 0) indxms(i) = 1
   IF (iautsp        == 0) CYCLE
   IF (andf(j,mskuo) /= 0 .AND. andf(j,mskul) /= 0) GO TO 740
   IF (andf(j,mskuo) /= 0 .AND. andf(j,mskur) /= 0) GO TO 740
   IF (andf(j,mskur) /= 0 .AND. andf(j,mskul) /= 0) GO TO 740
   IF (andf(j,mskur) /= 0) iexcld(i) = 1
   IF (multi == 0 .OR. indxms(i) /= 0) GO TO 25
   IF (andf(j,mskin) /= 0) iexcld(i) = 1
   25 IF (INDEX         < 3) CYCLE
   IF (andf(j,mskul) /= 0) iexcld(i) = 1
 END DO
 
!     DETERMINE THE ORDER OF SINGULARITY
 
 IF (iordr-2 < 0) THEN
   GO TO   230
 ELSE IF (iordr-2 == 0) THEN
   GO TO   260
 ELSE
   GO TO   410
 END IF
 
 
 40 logic = 100
 IF (ispc > initl) GO TO 750
 IF (ieqexn ==   1) GO TO 60
 ieqexn = 1
 
!     BRING IN EQEXIN
 
 CALL gopen (eqexin,iz(ibuf2),0)
 CALL skprec (eqexin,1)
 mcb(1) = eqexin
 CALL rdtrl (mcb)
 icore = luset + 2*(mcb(2)+1) - ibuf2
 IF (icore >= 0) GO TO 720
 ifile = eqexin
 CALL READ (*700,*50,eqexin,iz(luset+1),ibuf2-luset,0,neqexn)
 GO TO 720
 50 CALL CLOSE (eqexin,1)
 CALL sort (0,0,2,2,iz(luset+1),neqexn)
 iz(luset+neqexn+1) = 0
 iz(luset+neqexn+2) = 10*(luset+1)
 neqexn = neqexn + 2
 
!     LOOK UP SIL IN EQEXIN
 
 istart = 2
 60 kk = ibase
 DO  i = istart,neqexn,2
   k = luset + i
   isil = iz(k)/10
   IF (kk < isil) GO TO 80
 END DO
 logic = 110
 GO TO 750
 
!     PICK UP POINT ID AND TYPE (GRID OR SCALAR) FROM EQEXIN
 
 80 igpid = iz(k-3)
 isil  = iz(k-2)/10
 ityp  = iz(k-2) - 10*isil
 istart= i - 2
 IF (ityp == 1) GO TO 90
 
!     SCALAR POINT
 
 iordr = 0
 npts  = 0
 
 
 90 IF (ispc == 0) GO TO 140
 logic = 120
 IF (ityp == 2 .AND. ispc > 1) GO TO 750
 IF (ityp  == 2) jponts(1) = 0
 IF (iscr2 == 1) GO TO 100
 iscr2    = 1
 iword(1) = spcset
 IF (iword(1) <= 0) iword(1) = 1
 
!     INITIALIZE SCR2
 
 CALL gopen (scr2,iz(ibuf3),1)
 
!     WRITE AUTOMATICALLY GENERATED SPC1 DATA ON SCR2
 
 100 DO  i = 1,ispc
   IF (ityp      == 2) GO TO 120
   IF (jponts(i) > 0) GO TO 110
   logic = 130
   GO TO 750
   110 jponts(i) = jponts(i) - isil + 1
   120 CALL WRITE (scr2,jponts(i),1,0)
   CALL WRITE (scr2,igpid,    1,0)
 END DO
 IF (ispc+ims >= initl) GO TO 10
 
 
 140 IF (iogpst == 1) GO TO 150
 iogpst = 1
 
!     INITIALIZE OGPST
 
 CALL gopen (ogpst,iz(ibuf2),1)
 ogpst1( 1) = 0
 ogpst1( 2) = 8
 ogpst1( 3) = spcset
 ogpst1( 4) = mpcset
 ogpst1(10) = 12
 CALL WRITE (ogpst,ogpst1, 10,0)
 CALL WRITE (ogpst,iz,     40,0)
 CALL WRITE (ogpst,head(1),96,1)
 
!     PUT OUT ERROR RECORDS ON OGPST
 
 150 CALL WRITE (ogpst,igpid,1,0)
 CALL WRITE (ogpst,ityp ,1,0)
 CALL WRITE (ogpst,iordr,1,0)
 iordr = iordr + 1
 IF (iordr == 1) GO TO 180
 DO  i = 1,npts
   IF (iponts(i) > 0) GO TO 160
   logic = 140
   GO TO 750
   160 iponts(i) = iponts(i) - isil + 1
 END DO
 logic = 150
 SELECT CASE ( iordr )
   CASE (    1)
     GO TO 750
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 210
   CASE (    4)
     GO TO 220
 END SELECT
 
!     SCALAR
 
 180 DO  i = 1,9
   iponts(i) = 0
 END DO
 GO TO 220
 
!     FIRST ORDER OUTPUT
 
 200 iponts(4) = iponts(2)
 iponts(7) = iponts(3)
 iponts(2) = 0
 iponts(3) = 0
 iponts(5) = 0
 iponts(6) = 0
 iponts(8) = 0
 iponts(9) = 0
 GO TO 220
 
!     SECOND ORDER OUTPUT
 
 210 iponts(8) = iponts(6)
 iponts(7) = iponts(5)
 iponts(5) = iponts(4)
 iponts(4) = iponts(3)
 iponts(3) = 0
 iponts(6) = 0
 iponts(9) = 0
 
!     THIRD ORDER OUTPUT
 
 220 CALL WRITE (ogpst,iponts,9,0)
 GO TO 10
 
!     FIRST ORDER SINGULARITY
 
 230 DO  i = 1,npts
   IF (indxms(i) /= 0) GO TO 10
 END DO
 IF (iautsp    == 0) GO TO 40
 DO  i = 1,npts
   IF (iexcld(i) /= 0) CYCLE
   ii     = iponts(i)
   iz(ii) = msksng
   nauto  = nauto + 1
   ispc   = ispc  + 1
   jponts(ispc) = ii
   GO TO 40
 END DO
 GO TO 40
 
!     SECOND ORDER SINGULARITY
 
 260 iloop = 1
 270 DO  i = 1,npts,2
   ii = iponts(i)
   IF (ii        == 0) GO TO 310
   IF (indxms(i) /= 0) GO TO 280
   IF (iloop     == 1) GO TO 310
   IF (iexcld(i) /= 0) GO TO 310
   iz(ii) = msksng
   nauto  = nauto + 1
   ispc   = ispc  + 1
   jponts(ispc) = ii
   GO TO 290
   280 ims   = ims + 1
   290 iordr = 1
   DO  iii = 1,npts
     IF (iponts(iii) == ii) iponts(iii) = 0
   END DO
   ii = 0
   310 jj = iponts(i+1)
   IF (jj          == 0) GO TO 330
   IF (indxms(i+1) /= 0) GO TO 320
   IF (iloop       == 1) CYCLE
   IF (iexcld(i+1) /= 0) CYCLE
   iz(jj) = msksng
   nauto  = nauto + 1
   ispc   = ispc  + 1
   jponts(ispc) = jj
   GO TO 330
   320 ims = ims + 1
   330 IF (ii /= 0) GO TO 340
   logic = 160
   IF (ispc+ims < 2) GO TO 750
   IF (ispc     == 0) GO TO 10
   GO TO 380
   340 iordr = 1
   DO  iii = 1, npts
     IF (iponts(iii) == jj) iponts(iii) = 0
   END DO
 END DO
 IF (iautsp == 0) GO TO 370
 IF (iloop  == 2) GO TO 370
 iloop = 2
 GO TO 270
 370 IF (iordr == 1) GO TO 380
 IF (iordr == 2) GO TO 40
 logic = 170
 GO TO 750
 380 iok = 0
 DO  i = 1, npts
   IF (iponts(i) == 0) CYCLE
   iok = iok + 1
   iponts(iok) = iponts(i)
   IF (iok /=  i) iponts(i) = 0
   IF (i == npts) CYCLE
   ii = i + 1
   DO  j = ii,npts
     IF (iponts(j) == 0) CYCLE
     IF (iponts(j) == iponts(iok)) iponts(j) = 0
   END DO
 END DO
 npts = iok
 IF (npts == 0) GO TO 40
 
 logic = 180
 IF (npts > 2) GO TO 750
 logic = 190
 IF (iponts(1) == iponts(2)) GO TO 750
 GO TO 40
 
!     THIRD ORDER SINGULARITY
 
 410 iok = 0
 DO  i = 1,npts
   IF (indxms(i) /= 0) GO TO 430
   IF (iautsp    == 0) GO TO 420
   IF (iexcld(i) /= 0) GO TO 420
   ii     = iponts(i)
   iz(ii) = msksng
   nauto  = nauto + 1
   ispc   = ispc  + 1
   jponts(ispc) = ii
   GO TO 440
   420 iok = 1
   CYCLE
   430 ims   = ims   + 1
   440 iordr = iordr - 1
   iponts(i) = 0
 END DO
 IF (iok == 1) GO TO 460
 logic = 200
 IF (ispc+ims /= 3) GO TO 750
 IF (ispc     == 0) GO TO 10
 GO TO 40
 460 iok = 0
 DO  i = 1,npts
   IF (iponts(i) == 0) CYCLE
   iok = iok + 1
   iponts(iok) = iponts(i)
   IF (iok /= i) iponts(i) = 0
 END DO
 npts = iok
 GO TO 40
 
 480 CALL CLOSE (gpst,1)
 IF (iogpst /= 1) GO TO 490
 CALL CLOSE (ogpst,1)
 IF (ierror /= 0) GO TO 490
 CALL makmcb (ogpst1,ogpst,0,0,0)
 ogpst1(2) = 8
 CALL wrttrl (ogpst1)
 490 IF (iautsp == 0) GO TO 610
 IF (nauto  > 0) GO TO 500
 logic = 210
 IF (iscr2  == 1) GO TO 750
 IF (iogpst == 1 .AND. INDEX < 3) WRITE (ioutpt,810) uwm
 IF (iogpst == 1 .AND. INDEX == 3) WRITE (ioutpt,815) uwm
 GO TO 610
 500 logic = 220
 IF (iscr2 /= 1) GO TO 750
 CALL WRITE (scr2,0,0,1)
 CALL CLOSE (scr2,1)
 IF (ierror /= 0) GO TO 610
 IF (iogpst /= 1) WRITE (ioutpt,800) uim
 IF (iogpst == 1) WRITE (ioutpt,805) uim
 IF (iogpst == 1 .AND. INDEX < 3) WRITE (ioutpt,820) uwm
 IF (iogpst == 1 .AND. INDEX == 3) WRITE (ioutpt,825) uwm
 
!     PRINT OUT AND, IF REQUESTED, PUNCH OUT
!     AUTOMATICALLY GENERATED SPC DATA CARDS
 
 CALL gopen (scr2,iz(ibuf3),0)
 ifile = scr2
 CALL READ (*700,*510,scr2,iz(luset+1),ibuf3-luset,0,iflag)
 icore = luset + 2*nauto - ibuf3
 GO TO 720
 510 logic = 230
 IF (iflag /= 2*nauto) GO TO 750
 CALL sort (0,0,2,1,iz(luset+1),iflag)
 i    = luset + 1
 iold = -1
 ist  = i
 520 j = 0
 530 IF (i > luset+iflag) GO TO 540
 IF (iold >= 0 .AND. iz(i) /= iold) GO TO 540
 iold = iz(i)
 j = j + 2
 i = i + 2
 GO TO 530
 540 CALL sort (0,0,2,-2,iz(ist),j)
 IF (i > luset+iflag) GO TO 550
 iold = iz(i)
 ist  = i
 GO TO 520
 
 550 i = luset + 1
 iold = -1
 CALL page1
 WRITE (ioutpt,830)
 line = line + 6
 560 ii = 2
 DO  j = 1,6
   IF (i > luset+iflag) EXIT
   IF (iold >= 0 .AND. iz(i) /= iold) EXIT
   iold = iz(i)
   iword(ii+1) = iz(i+1)
   ii = ii + 1
   i  = i  + 2
 END DO
 580 iword(2) = iold
 IF (line <= nlpp) GO TO 590
 CALL page1
 WRITE (ioutpt,830)
 line = line + 6
 590 ncard = ncard + 1
 WRITE (ioutpt,840) ncard,(iword(j),j=1,ii)
 line = line + 1
 IF (iautsp < 0) WRITE (ipunch,850) (iword(j),j=1,ii)
 IF (i > luset+iflag) GO TO 600
 iold = iz(i)
 GO TO 560
 600 CALL CLOSE (scr2,1)
 610 IF (iautsp == 0 .OR. multi == 0) RETURN
 DO  i = 1,luset
   iz(i) = andf(iz(i),mskxx)
 END DO
 RETURN
 
!     ERROR MESSAGES
 
 700 num = -2
 710 CALL mesage (num,ifile,isubnm)
 720 num = -8
 ifile = icore
 GO TO 710
 730 ierror = 1
 WRITE (ioutpt,870) uwm
 GO TO 480
 740 ierror = 2
 WRITE (ioutpt,880) uwm
 GO TO 480
 750 WRITE (ioutpt,860) sfm,logic
 CALL mesage (-61,0,0)
 
 800 FORMAT (a29,' 2435, AT USER''S REQUEST, ALL POTENTIAL ',  &
     'SINGULARITIES HAVE BEEN REMOVED BY THE', /5X,  &
     'APPLICATION OF SINGLE POINT CONSTRAINTS.  REFER TO PRINT'  &
     ,       'OUT OF AUTOMATICALLY GENERATED SPC1 CARDS FOR DETAILS.')
 805 FORMAT (a29,' 2436, AT USER''S REQUEST, ONE OR MORE POTENTIAL ',  &
     'SINGULARITIES HAVE BEEN REMOVED BY THE', /5X,  &
     'APPLICATION OF SINGLE POINT CONSTRAINTS.  REFER TO PRINT'  &
     ,       'OUT OF AUTOMATICALLY GENERATED SPC1 CARDS FOR DETAILS.')
 810 FORMAT (a25,' 2437A, IN SPITE OF THE USER''S REQUEST, NONE OF ',  &
     'THE POTENTIAL SINGULARITIES HAS BEEN REMOVED', /5X,  &
     'BECAUSE OF THG PRESENCE OF SUPORT CARDS AND/OR MULTI',  &
     'POINT CONSTRAINTS OR RIGID ELEMENTS.', /5X,  &
     'REFER TO THE GRID POINT SINGULARITY TABLE FOR DETAILS.')
 815 FORMAT (a25,' 2437A, IN SPITE OF THE USER''S REQUEST, NONE OF ',  &
     'THE POTENTIAL SINGULARITIES HAS BEEN REMOVED', /5X,  &
     'BECAUSE OF THG PRESENCE OF SUPORT CARDS AND/OR MULTI',  &
     'POINT CONSTRAINTS OR RIGID ELEMENTS', /5X,'OR BECAUSE ',  &
     'THE SINGULARITIES ARE NOT PART OF THE OMIT SET (O-SET) ',  &
     'DEGREES OF FREEDOM.', /5X,  &
     'REFER TO THE GRID POINT SINGULARITY TABLE FOR DETAILS.')
 820 FORMAT (a25,' 2437, ONE OR MORE POTENTIAL SINGULARITIES HAVE NOT',  &
     ' BEEN REMOVED', /5X,'BECAUSE OF THG PRESENCE OF SUPORT ',  &
     'CARDS AND/OR MULTIPOINT CONSTRAINTS OR RIGID ELEMENTS.',  &
     /5X,'REFER TO THE GRID POINT SINGULARITY TABLE FOR DETAILS.')
 825 FORMAT (a25,' 2437, ONE OR MORE POTENTIAL SINGULARITIES HAVE NOT',  &
     ' BEEN REMOVED', /5X,'BECAUSE OF THG PRESENCE OF SUPORT ',  &
     'CARDS AND/OR MULTIPOINT CONSTRAINTS OR RIGID ELEMENTS',  &
     /5X,'OR BECAUSE THE SINGULARITIES ARE NOT PART OF THE ',  &
     'OMIT SET (O-SET) DEGREES OF FREEDOM.', /5X,  &
     'REFER TO THE GRID POINT SINGULARITY TABLE FOR DETAILS.')
 830 FORMAT (//32X, 'A U T O M A T I C A L L Y   ', 'G E N E R A T E D   ',  &
     'S P C 1   C A R D S', /, 16X, 'CARD ',8X, /,  &
     16X, 'COUNT',8X, '---1--- +++2+++ ---3--- +++4+++ ---5--- ',  &
     '+++6+++ ---7--- +++8+++ ---9--- +++10+++',/)
 840 FORMAT (15X, i5, '-', 8X, 'SPC1    ',8I8)
 850 FORMAT (                  'SPC1    ',8I8)
 860 FORMAT (a25,' 2438, LOGIC ERROR NO.',i4,  &
     ' IN SUBROUTINE GP4SP IN MODULE GP4')
 870 FORMAT (a25,' 2439, ILLEGAL VALUE INPUT FOR PARAMETER AUTOSPC - ',  &
     'SINGULARITY PROCESSING SKIPPED IN MODULE GP4')
 880 FORMAT (a25,' 2440, SINGULARITY PROCESSING SKIPPED IN MODULE GP4',  &
     ' BECAUSE OF INCONSISTENT SET DEFINITION')
 
 RETURN
END SUBROUTINE gp4sp
