SUBROUTINE grbvec
     
!     THIS SUBROUITNE IS THE MAIN DRIVER FOR THE VECGRB MODULE
!     WHICH GENERATES
 
!     (1) THE GEOMETRIC RIGID BODY VECTORS ABOUT THE INDICATED GRID
!         POINT OR ORIGIN.
!         THIS SET OF VECTORS CONSISTS OF UNIT DISPLACEMENTS IN ZERO
!         COORDINATE SYSTEM ABOUT THE SPECIFIED GRID IN GLOBAL COORD.
!         FOR EASE OF ASSEMBLY THE VECTOR IS GENERATED IN THE TRANSPOSED
!         FORM, THAT IS, WITH SIX ROW, ONE FOR EACH OF THE SIX UNIT
!         MOTIONS AND G-SET COLUMNS, ONE FOR EACH DOF'S CORRESPONDING
!         MOTION. THIS SET OF VECTORS WOULD BE EXACTLY EQUAL TO A UNIT
!         DISPLACEMENT CHECK IF ALL THE GRIDS HAD STIFFNESS BUT WERE
!         NOT GROUNDED.
 
!     (2) A g-SET SIZED CSTM FROM BASIC TO GOLBAL
 
!     DMAP SEQUENCE -
 
!     VECGRB    BGPDT,EQEXIN,CSTM/OUTVEC/P1/P2/P3   $
 
!     WHERE     P1 = 1, GENERATE CSTM FROM BASIC TO GLOBAL
!                  = 2, GENERATE PHIRBT
!               P2 = REFERENCE GRID FOR PHIRB (0=BASIC, DEFAULT)
!               P3 = CURRENTLY NOT USED
 
!     EXAMPLES -
 
!     (1) G-SET EQUILIBRIUM CHECK
!     THIS CHECK MULTIPLIES THE STIFFNESS MATRIX TIMES THE GEOMETRIC
!     RIGID BODY SHAPES PENERATED BY VECRGB. THE FORCES OBATINED FROM
!     THIS MULTIPLICATION SHOULD BE ZERO.
 
!     VECGRB   BGPDT,CSTM,EQEXIN/PHIRBT/2/0  $ CREATE TRANSPOSE OF RIGID
!     TRNSP    PHIRBT/PHIRB                  $ BODY VECTORS, THEN TRNSP
!     MPYAD    KGG,PHIRB,/KPHIG/0            $ MULTIPLY BY STIFFNESS.
!     MPYAD    PHIRBT,KPHIG,/KPHG6/0         $ SUM FORCES AND PRINT
!     MATPRN   KPHG6,,,, //                  $ 6X6 SUMMATION. PRINT ALL
!     MATGPR   GPL,USET,SIL,KPHIG//*G*/*G*// $ FORCES OVER 0.0001
!              .0001                         $
 
!     (2) COORDINATE SYSTEM TRANSFORMATION
 
!     VECRGB   BGPDT,CSTM,EQEXIN/BCSTM/1     $ TRANSFORM GLOBAL KGG TO
!     TRNSP    BCSTM/BCSTMT/                 $ BASIC
!     MPYAD    BCSTM,KGG,/BGKGG/0            $
!     MPYAD    BGKGG,BCSTMT,/BBKGG/0         $
 
!     THIS SUBROUTINE WAS ORIGINALLY CALLED CSTMX, AND WAS WRITTEN BY
!     P.KIRCHMAN/SWALES, 2/1993, WITH THE DMAP MODULE OF THE SAME NAME
 
!     THE DMAP MODLUE IS RENAMED TO GEOMETRIC RIGID BODY VECTOR, VECGRB,
!     AND THE SUBROUTINE GRBVEC. THE ORIGINAL SUBROUTINE WAS RE-CODED BY
!     G.CHAN/UNISIS, USING NASTRAN TRADITIONAL FORTRAN STYLE, AND THE
!     SUBSTITUTION OF GMMATD ROUTINE FOR DP3X3M. ALSO, THE ORDER OF 2ND
!     AND 3RD INPUT DATA BLOCKS IS INTERCHANGED.
 
!     THE ORIGINAL CSTMX IS INCLUDED IN THE 1993 RELEASE. IT IS ONLY
!     FOR BACKUP PURPOSE. CSTMX WILL BE DELETED IN NEXT NASTRAN RELEASE
!     (THE ORIGINAL CSTMX ROUTINE PRODUCED HUNDREDS OF FORTRAN ERRORS
!     ON CDC MACHINE WITH FTN5 COMPILER. 3/93)
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: rbv
 INTEGER :: trl(7),NAME(2),sub(2)
 REAL :: rx(1)
 DOUBLE PRECISION :: rbvr(3),v3(3),zero,one,xg,yg,zg,rad,xl,tu(9),  &
     t1(9),t2(9),t(9),rvec(9),rbvec(9),v(3),vout(6)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /BLANK /  p1,p2,p3
 COMMON /zzzzzz/  ix(1)
 COMMON /system/  ibuff,nout
 COMMON /packx /  tin,tou,ii,jj,incr
 EQUIVALENCE      (ix(1),rx(1))
 DATA    bgpdt ,  eqexin, cstm / 101,102,103  /,  &
     outvec/  201 /,  sub  / 4HGRBV,4HC   /
 DATA    zero  ,  one /  0.0D+0, 1.0D+0       /
 DATA    cstmx ,  eqe,xin,       bgp,dt       /  &
     4HCSTM,  4HEQEX,4HIN  , 4HBGPD,4HT   /
 
!     CHECK FOR THE PRESENCE OF OUTPUT DATA BLOCK
 
 trl(1) = outvec
 CALL rdtrl (trl)
 IF (trl(1) <= 0) GO TO 1000
 
!     INITIALIZATION
 
 lcor  = korsz(ix(1)) - ibuff
 buf1  = lcor - 1
 tu(1) = one
 tu(2) = zero
 tu(3) = zero
 tu(4) = zero
 tu(5) = one
 tu(6) = zero
 tu(7) = zero
 tu(8) = zero
 tu(9) = one
 
!     CHECK THE PRESENCE OF BGPDT AND CSTM FILES
 
 trl(1)= bgpdt
 CALL rdtrl (trl)
 IF (trl(1) <= 0) GO TO 1020
 nent  = trl(2)
 trl(1)= cstm
 CALL rdtrl (trl)
 ncst  = trl(3)
 IF (trl(1) <= 0) ncst = 0
 nent4 = nent*4
 ncst14= ncst*14
 IF (nent4+ncst14 > lcor) GO TO 1100
 rbv   = p1 == 2
 
!     CHECK IF THIRD INPUT FILE IS PRESENT, THEN OPEN EQEXIN FILE AND
!     READ THE FIRST TABLE INTO CORE IF APPROPRIEATE
 
 grdpnt = 0
 IF (p1 /= 2 .OR. p2 == 0) GO TO 60
 trl(1) = eqexin
 CALL rdtrl (trl)
 IF (trl(1) <= 0) GO TO 1020
 FILE = eqexin
 CALL fname (eqexin,NAME)
 IF (NAME(1) /= eqe .OR. NAME(2) /= xin) GO TO 1040
 CALL OPEN (*1200,eqexin,ix(buf1),0)
 CALL fwdrec (*1300,eqexin)
 CALL READ (*1300,*20,eqexin,ix(1),buf1-1,1,flag)
 CALL mesage (-8,0,sub)
 20 trl2 = trl(2)*2
 IF (flag /= trl2) GO TO 1320
 j = 1
 DO  i = 1,trl2,2
   IF (p2 /= ix(j)) GO TO 30
   grdpnt = ix(j+1)
   GO TO 50
   30 j = j + 2
 END DO
 grdpnt = 0
 WRITE  (nout,40) uwm,p2
 40 FORMAT (a25,' - ID ',i8,' IS NOT A GRID POINT.  THE ORIGIN WILL ',  &
     'BE USED.')
 50 CALL CLOSE (eqexin,1)
 
!     OPEN AND READ BGPDT TABLE INTO BEGINNING OF CORE
 
 60 FILE = bgpdt
 CALL fname (bgpdt,NAME)
 IF (NAME(1) /= bgp .AND. NAME(2) /= dt) GO TO 1040
 CALL OPEN (*1200,bgpdt,ix(buf1),0)
 CALL fwdrec (*1300,bgpdt)
 CALL READ (*1300,*70,bgpdt,ix(1),buf1-1,1,flag)
 CALL mesage (-8,0,sub)
 70 IF (flag /= nent4) GO TO 1330
 CALL CLOSE (bgpdt,1)
 
!     OPEN AND READ CSTM FIRST TABLE INOT CORE AFTER BGPDT
 
 IF (ncst == 0) GO TO 90
 FILE = cstm
 CALL fname (cstm,NAME)
 IF (NAME(1) /= cstmx) GO TO 1040
 CALL OPEN (*1200,cstm,ix(buf1),0)
 CALL fwdrec (*1300,cstm)
 CALL READ (*1300,*80,cstm,ix(nent4+1),buf1-nent4-1,1,flag)
 CALL mesage (-8,0,sub)
 80 IF (flag /= ncst14) GO TO 1340
 CALL CLOSE (cstm,1)
 
!     USE BGPDT INFO TO FIGURE OUT THE g-SET SIZE FOR OUTPUT
 
 90 size = 0
 i    = 1
 DO  j = 1,nent
   IF (ix(i) < 0) size = size + 1
   IF (ix(i) >= 0) size = size + 6
   i = i + 4
 END DO
 
!     STORE RIGID BODY REFERENCE VECTOR
 
 rbvr(1) = 0.
 rbvr(2) = 0.
 rbvr(3) = 0.
 IF (.NOT.rbv .OR. grdpnt == 0) GO TO 110
 rbvr(1) = rx(grdpnt*4-2)
 rbvr(2) = rx(grdpnt*4-1)
 rbvr(3) = rx(grdpnt*4  )
 
!     OPNE OUTPUT FILE AND FILL OUTPUT TRAILER
 
 110 CALL fname (outvec,NAME)
 CALL OPEN (*1200,outvec,ix(buf1),1)
 CALL WRITE (outvec,NAME,2,1)
 IF (     rbv) CALL makmcb (trl(1),outvec,6,2,2)
 IF (.NOT.rbv) CALL makmcb (trl(1),outvec,size,2,2)
 
!     INITIALIZE PACK COMMONS
 
 tin  = 2
 tou  = 2
 incr = 1
 col  = 1
 
!     BEGIN LOOP FOR NUMBER OF ENTRIES IN BGPDT
 
 e4  = 0
 DO  ENTRY = 1,nent
   e4  = e4 + 4
   ocid= ix(e4-3)
   xg  = rx(e4-2)
   yg  = rx(e4-1)
   zg  = rx(e4  )
   
!     RBVEC IS A VECTOR BASED ON UNIT ROTATIONS OF A VECTOR FROM THE
!     REFERENCE GRID TO THE GRID IN QUESTION. THE TRANSFORMATION IS
!     FROM BASIC ROTATIONS TO BASIC TRANSLATIONS.
   
   IF (.NOT.rbv) GO TO 130
   rvec(1) = zero
   rvec(2) = (zg-rbvr(3))
   rvec(3) =-(yg-rbvr(2))
   rvec(4) =-(zg-rbvr(3))
   rvec(5) = zero
   rvec(6) = (xg-rbvr(1))
   rvec(7) = (yg-rbvr(2))
   rvec(8) =-(xg-rbvr(1))
   rvec(9) = zero
   
!     IF THIS ENTRY IS A SCALAR AND A RIGID BODY VECTOR HAS BEEN
!     REQUESTED, STORE A ZERO COLUMN
   
   130  IF (ocid /= -1 .OR. .NOT.rbv) GO TO 140
   ii  = 1
   jj  = 1
   CALL pack (zero,outvec,trl)
   col = col + 1
   CYCLE
   
!     IF THIS ENTRY IS A SCALAR AND A CSTM HAS BEEN REQUESTED, SIMPLY
!     PLACE A ONE ON THE DIAGONAL AND CONTINUE
   
   140 IF (ocid /= -1) GO TO 150
   ii  = col
   jj  = col
   CALL pack (one,outvec,trl)
   col = col + 1
   CYCLE
   
!     IF THIS ENTRY IS ALREADY IN BASIC COORDINATES, STORE AN IDENTITY
!     IN THE APPOPRIATE SIX BY SIX
   
   150 IF (ocid /= 0 .OR. .NOT.rbv) GO TO 190
   ii  = 1
   jj  = 1
   DO  i = 1,3
     i3  = 0
     DO  j = 1,3
       vout(j  ) =   tu(j+i*3)
       vout(j+3) = rvec(j+i*3)
     END DO
     i3  = i3 + 3
     ii  = 1
     jj  = 6
     CALL pack (vout,outvec,trl)
     col = col + 1
   END DO
   DO  i = 1,3
     ii  = i + 3
     jj  = i + 3
     CALL pack (one,outvec,trl)
   END DO
   CYCLE
   
   190 IF (ocid /= 0) GO TO 210
   DO  i = 1,6
     ii  = col
     jj  = col
     CALL pack (one,outvec,trl)
     col = col + 1
   END DO
   CYCLE
   
!     CSTM MUST BE MISSING
   
   210 IF (ncst /= 0) GO TO 220
   trl(1) = cstm
   GO TO 1020
   
!     SET UP VECTORS AND MATRICES COMMON TO ALL COORDINATE SYSTEM
!     TRANSFORMATIONS
   
!     FIND COORDINATE SYSTEM
   
   220 DO  icst = 1,ncst
     IF (ix(icst*14-13+nent4) == ocid) GO TO 240
   END DO
   GO TO 1400
   
!     GET COORDINATE SYSTEM TYPE AND
!     TRANSFORMATION FROM BASIC TO COORDINATE SYSTEM ORIGIN TRIAD
   
   240 ocidt = ix(icst*14-12+nent4)
   t1(1) = rx(icst*14- 8+nent4)
   t1(4) = rx(icst*14- 7+nent4)
   t1(7) = rx(icst*14- 6+nent4)
   t1(2) = rx(icst*14- 5+nent4)
   t1(5) = rx(icst*14- 4+nent4)
   t1(8) = rx(icst*14- 3+nent4)
   t1(3) = rx(icst*14- 2+nent4)
   t1(6) = rx(icst*14- 1+nent4)
   t1(9) = rx(icst*14   +nent4)
   
   ik    = icst*14 + nent4
   v3(1) = rx(ENTRY*4-2) - rx(ik-11)
   v3(2) = rx(ENTRY*4-1) - rx(ik-10)
   v3(3) = rx(ENTRY*4  ) - rx(ik- 9)
   
   v(1)  = rx(ik-8)*v3(1) + rx(ik-5)*v3(2) + rx(ik-2)*v3(3)
   v(2)  = rx(ik-7)*v3(1) + rx(ik-4)*v3(2) + rx(ik-1)*v3(3)
   v(3)  = rx(ik-6)*v3(1) + rx(ik-3)*v3(2) + rx(ik  )*v3(3)
   
!     SPECIAL CHECKS FOR ZERO RADIUS CYLINDRICAL OR SPHERICAL COORDINATE
!     SYSTEM. IF SO TREAT AS RECTANGULAR.
   
   rad = SQRT(v(1)**2 + v(2)**2)
   IF (rad == 0.) ocidt = 1
   
!     PERFORM INDIVIDUAL COORDINATE SYSTEM TRANSFORMATION AND GENERATE
!     T2
   
   SELECT CASE ( ocidt )
     CASE (    1)
       GO TO 250
     CASE (    2)
       GO TO 330
     CASE (    3)
       GO TO 340
   END SELECT
   
!     RECTANGULAR, T = T1
   
   250 IF (.NOT.rbv) GO TO 290
   INDEX = 1
   CALL gmmatd (rvec,3,3,0, t1,3,3,0, rbvec)
   
!     ADD RIGID BODY INFORMATION TO LOWER OFF DIAGONAL 3X3 IF REQUESTED
   
   DO  i = 1,3
     i3  = 0
     DO  j = 1,3
       vout(j  ) =    t1(j+i3)
       vout(j+3) = rbvec(j+i3)
     END DO
     i3  = i3 + 3
     ii  = 1
     jj  = 6
     CALL pack (vout,outvec,trl)
     col = col + 1
   END DO
   GO TO 310
   
!     OR SIMPLY PACK THE TRANSFORMATION
   
   290 INDEX = col
   DO  i = 1,3
     ii  = INDEX
     jj  = INDEX + 2
     CALL pack (t1(i*3-2),outvec,trl)
     col = col + 1
   END DO
   
!    STORE LOWER 3X3, AND GET NEXT GRID
   
   310 DO  i = 1,3
     ii  = INDEX + 3
     jj  = INDEX + 5
     CALL pack (t1(i*3-2),outvec,trl)
     col = col + 1
   END DO
   CYCLE
   
!     CYLINDRICAL
   
   330 t2(1) = v(1)/rad
   t2(4) =-v(2)/rad
   t2(7) = zero
   t2(2) =-t2(4)
   t2(5) = t2(1)
   t2(8) = zero
   t2(3) = zero
   t2(6) = zero
   t2(9) = one
   GO TO 350
   
!     SPHERICAL
   
   340 xl    = SQRT(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
   IF (xl <= 0.0) GO TO 1060
   t2(1) = v(1)/xl
   t2(4) =(v(1)*v(3))/(rad*xl)
   t2(7) =-v(2)/rad
   t2(2) = v(2)/xl
   t2(5) =(v(2)*v(3))/(rad*xl)
   t2(8) = v(1)/rad
   t2(3) = v(3)/xl
   t2(6) =-rad/xl
   t2(9) = zero
   
   350 CALL gmmatd (t1,3,3,0, t2,3,3,0, t)
   IF (.NOT.rbv) GO TO 400
   
!     ADD RIGID BODY INFORMATION TO LOWER OFF DIAGONAL 3X3 IF REQUESTED
!     THEN PACK
   
   INDEX = 1
   CALL gmmatd (rvec,3,3,0, t,3,3,0, rbvec)
   DO  i = 1,3
     i3  = 0
     DO  j = 1,3
       vout(j  ) =     t(j+i3)
       vout(j+3) = rbvec(j+i3)
     END DO
     i3  = i3 + 3
     ii  = 1
     jj  = 6
     CALL pack (vout,outvec,trl)
     col = col + 1
   END DO
   GO TO 420
   
!     OR SIMPLY PACK THE TRANSFORMATION
   
   400 INDEX = col
   DO  i = 1,3
     ii  = INDEX
     jj  = INDEX + 2
     CALL pack (t(i*3-2),outvec,trl)
     col = col + 1
   END DO
   
!     STORE LOWER 3X3
   
   420 DO  i = 1,3
     ii  = INDEX + 3
     jj  = INDEX + 5
     CALL pack (t(i*3-2),outvec,trl)
     col = col + 1
   END DO
   
 END DO
 
 CALL CLOSE  (outvec,1)
 CALL wrttrl (trl)
 RETURN
 
!     ERRORS
 
 1000 WRITE  (nout,1010) ufm
 1010 FORMAT (a23,'.  MISSING REQUIRED OUTPUT FILE')
 GO TO  1500
 1020 WRITE  (nout,1030) ufm,trl(1)
 1030 FORMAT (a23,'.  MISSING REQUIRED INPUT FILE',i4)
 GO TO  1500
 1040 WRITE  (nout,1050) ufm,NAME
 1050 FORMAT (a23,'. INPUT FILE ',2A4,' ERROR')
 GO TO  1500
 1060 WRITE  (nout,1070) ufm
 1070 FORMAT (a23,' FROM GRBVEC. ZERO RADIAL LENGTH, ERROR AT 340')
 GO TO  1500
 1100 j = -8
 GO TO  1490
 1200 j = -1
 GO TO  1490
 1300 j = -2
 GO TO  1490
 1320 j = trl2
 GO TO  1350
 1330 j = nent4
 GO TO  1350
 1340 j = ncst14
 1350 WRITE  (nout,1360) sfm,NAME,j,flag
 1360 FORMAT (a25,'. EXPECTED RECORD LENGTH DOES NOT MATCH ACTUAL ',  &
     ' RECORD LENGTH ON INPUT FILE ',2A4, /5X,2I10)
 GO TO  1500
 1400 WRITE  (nout,1410) ufm,ocid
 1410 FORMAT (a23,'. UNABLE TO FIND COORDINATE SYSTEM ',i8)
 GO TO  1500
 
 1490 CALL mesage (j,FILE,sub)
 1500 CALL mesage (-61,0,sub)
 RETURN
END SUBROUTINE grbvec
