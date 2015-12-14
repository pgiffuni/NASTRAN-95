SUBROUTINE hbdyd
     
!     THIS IS THE BOUNDARY CONDITION (HEAT) ELEMENT ROUTINE
!     IT PRODUCES THE STIFFNESS AND OR DAMPING ELEMENT MATRICES.
 
 LOGICAL :: heat     ,nogo
 INTEGER :: ngrids(7),necpt(53),outpt    ,siltab(8)  ,  &
     set1(8)  ,set2(4)  ,sils     ,dict(13)   , elid     ,estid
 REAL :: ecpt(53)
 DOUBLE PRECISION :: c(16)    ,cc(4,4)  ,pi       ,master(8,8),  &
     mast(64) ,ke       ,me       ,itemp      ,  &
     a1(5)    ,a2(3)    ,a3(3)    ,a4(3)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm      ,uwm      ,uim      ,sfm
 COMMON /system/  ksystm(100)
 COMMON /emgest/  ecpt1    ,iflag    ,sils(8)  ,v(3)       ,  &
     ecpt14   ,matflg   ,af       ,emiss      ,  &
     absorp   ,r1       ,r2       ,csid(4,8)  , avgtmp
 COMMON /hmtout/  hx       ,cpx
 COMMON /matin /  matid    ,inflag   ,eltemp
 COMMON /emgprm/  d15(15)  ,kmbgg(3) ,iprec    ,nogo       , heat     ,icmbar
 COMMON /emgdic/  dmmm(2)  ,nlocs    ,elid     ,estid
 COMMON /condad/  pi
 EQUIVALENCE      (necpt(1), ecpt(1)), (set2(1), set1(5))  ,  &
     (ecpt1   , ecpt(1)), (dict5  , dict(5))  ,  &
     (ksystm(2),outpt  ), (cc(1,1), c(1)   )  , (master(1,1),mast(1))
 DATA    ngrids/  1, 2, 2, 3, 4, 2 ,2/
 
!     EST ENTRY FOR -CHBDY- ELEMENT
!     ======================================================
!     ECPT( 1)  = EL-ID       ELEMENT ID
!     ECPT( 2)  = IFLAG       ELEM. TYPE FLAG = (1,2,3,4,5,6,7)
!     ECPT( 3)  = SIL-1       SCALER INDICES
!     ECPT( 4)  = SIL-2
!     ECPT( 5)  = SIL-3
!     ECPT( 6)  = SIL-4
!     ECPT( 7)  = SIL-5
!     ECPT( 8)  = SIL-6
!     ECPT( 9)  = SIL-7
!     ECPT(10)  = SIL-8
!     ECPT(11)  = V1          ORIENTATION VECTOR
!     ECPT(12)  = V2
!     ECPT(13)  = V3
!     ECPT(14)  = ECPT14
!     ECPT(15)  = MATFLG      MAT ID FOR MAT4, MAT5 DATA
!     ECPT(16)  = AF          AREA FACTOR
!     ECPT(17)  = EMISS       EMISSIVITY COEFF
!     ECPT(18)  = ABSORP      ABSORPTIVITY COEFF
!     ECPT(19)  = R1          RADII OF ELIPTICAL CYLINDER
!     ECPT(20)  = R2
!     ECPT(21)  = CSID-1      COORDINATE SYSTEM ID AND
!     ECPT(22)  = X1          COORDINATE GRID POINTS
!     ECPT(23)  = Y1          (1-4 ARE ELEMENT POINTS,
!     ECPT(24)  = Z1
!     ECPT(25)  = CSID-2
!     ECPT(26)  = X2
!     ECPT(27)  = Y2
!     ECPT(28)  = Z2
!     ECPT(29)  = CSID-3
!     ECPT(30)  = X3
!     ECPT(31)  = Y3
!     ECPT(32)  = Z3
!     ECPT(33)  = CSID-4
!     ECPT(34)  = X4
!     ECPT(35)  = Y4
!     ECPT(36)  = Z4
!     ECPT(37)  = CSID-5       5-8 ARE POINTS IN THE FLUID)
!       -ETC-     -ETC-
!     ECPT(53)  = AVGTMP      AVERAGE ELEM. TEMPERATURE
 
!     GENERAL INITIALIZATION
 
 IF (.NOT. heat) RETURN
 imhere = 0
 IF (iflag < 1 .OR. iflag > 7) GO TO 470
 IF (iflag == 7) af = pi*(DBLE(r1)+DBLE(r2))
 n = ngrids(iflag)
 dict(1) = estid
 dict(2) = 1
 dict(4) = 1
 dict5   = 0.0
 
!     MASTER OUTPUT MATRIX OF SIZE UP TO 8 X 8 IS FORMED.  DUPLICATE
!     SILS ARE SUPERIMPOSED RESULTING IN A POSSIBLY SMALLER OUTPUT MATRX
 
!     FOR A GIVEN ELEMENT-ID THE MATRIX OUTPUT WILL BE OF ORDER EQUAL
!     TO THE NUMBER OF UNIQUE SILS PRESENT.
 
!     IFLAG = 1 WILL BE 1X1 OR 2X2    *
!     IFLAG = 2 WILL BE 2X2 UP TO 4X4  *
!     IFLAG = 3 WILL BE 2X2 UP TO 4X4   * (DEPENDING ON GROUDING AND
!     IFLAG = 4 WILL BE 3X3 UP TO 6X6  *   DUPLICATE SILS.)
!     IFLAG = 5 WILL BE 4X4 UP TO 8X8 *
 
!     -SET1- WILL BE A MAP OF OUTPUT POSITIONS FOR SILS 1 THRU 4
!     -SET2- WILL BE A MAP OF OUTPUT POSITIONS FOR SILS 5 THRU 8
 
 
!     FIRST FORM THE TABLE OF UNIQUE SILS.
 
 isize = 0
 loop50:  DO  i = 1,8
   IF (sils(i) <= 0) CYCLE loop50
   IF (isize   <= 0) GO TO 40
   DO  j = 1,isize
     IF (sils(i) == siltab(j)) CYCLE loop50
   END DO
   40 isize = isize + 1
   siltab(isize) = sils(i)
 END DO loop50
 CALL sort (0,0,1,1,siltab(1),isize)
 imhere = 50
 IF (isize <= 0) GO TO 470
 
!     BUILD -SET1- AND -SET2- MAPS OF WHERE OUTPUTS GO IN MASTER OUTPUT.
 
 DO  i = 1,8
   j = 8
   IF (sils(i) <= 0) GO TO 90
   DO  j = 1,isize
     IF (sils(i) == siltab(j)) GO TO 90
   END DO
   imhere = 80
   GO TO 470
   90 set1(i) = j
 END DO
 dict(3) = isize
 
!     FORM STIFFNESS -HEAT- IF REQUESTED.
 
 IF (kmbgg(1) == 0) GO TO 360
 inflag = 1
 eltemp = avgtmp
 matid  = matflg
 IF (matid == 0) GO TO 360
 CALL hmat (necpt)
 cp = cpx
 h  = hx
 IF (h == 0.0) GO TO 360
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 120
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 140
   CASE (    4)
     GO TO 210
   CASE (    5)
     GO TO 240
   CASE (    6)
     GO TO 130
   CASE (    7)
     GO TO 130
 END SELECT
 
!     IFLAG = 1, (POINT), 1 GRID-POINT.  (1 X 1)  C = H * AF
 
 120 c(1) = h
 c(2) = af
 c(1) = c(1)*c(2)
 GO TO 300
 
!     IFLAG = 2, (LINE OR ELLIPTIC CYL. )    **    **
!             2 GRID POINTS           H*AF*L * 2  1 *
!                           (2X2)  C =------ *      *
!                                       6    * 1  2 *
!                                            **    **
 
 130 c(1) = h
 c(2) = af
 c(3) = ecpt(26) - ecpt(22)
 c(4) = ecpt(27) - ecpt(23)
 c(5) = ecpt(28) - ecpt(24)
 c(1) = c(1)*c(2)*DSQRT(c(3)**2 + c(4)**2 + c(5)**2)/3.0D0
 c(2) = c(1)/2.0D0
 c(5) = c(2)
 c(6) = c(1)
 GO TO 300
 
!     IFLAG = 3, (REVOLUTION), 2 GRID-POINTS     **                **
!                                                *(3X +X )  (X + X )*
!                                        H*2PI*L *   1  2     1   2 *
!                             (2X2)  C = ------- *                  *
!                                          12    *(X + X )  (X +3X )*
!                                                *  1   2     1   2 *
!                                                **                **
 
 140 IF (ecpt(22) <= 0.0 .OR. ecpt(26) <= 0.0) GO TO 180
 IF (ecpt(23) /= 0.0 .OR. ecpt(27) /= 0.0) GO TO 180
 GO TO 200
 180 WRITE  (outpt,190) ufm,necpt(1)
 190 FORMAT (a23,' 3088, ILLEGAL GEOMETRY FOR REVOLUTION ELEMENT',i14)
 nogo = .true.
 GO TO 490
 
!     FILL CONDUCTIVITIY MATRIX
 
 200 c(1) = h
 c(2) = pi
 c(3) = ecpt(26) - ecpt(22)
 c(4) = ecpt(28) - ecpt(24)
 
!     NOTE Y2 AND Y1 ARE 0 FOR REVOLUTION ELEMENT.
 
 c(1) = c(1)*c(2)*DSQRT(c(3)**2 + c(4)**2)/6.0D0
 c(2) = c(1)*DBLE(ecpt(22) + ecpt(26))
 c(5) = c(2)
 c(6) = c(1)*DBLE(ecpt(22) + 3.0*ecpt(26))
 c(1) = c(1)*DBLE(3.0*ecpt(22) + ecpt(26))
 GO TO 300
 
!     IFLAG = 4, (TRIANGLE), 3 GRID-POINTS.       **       **
!                                                 * 2  1  1 *
!                                           H * A *         *
!                                (3X3) C =  ----- * 1  2  1 *
!                                            24   *         *
!                                                 * 1  1  2 *
!                                                 **       **
 
 
!     COMPUTE AREA -A- OF TRIANGLE   GET R2-R1 AND R3-R2
 
 210 c(1) = ecpt(26) - ecpt(22)
 c(2) = ecpt(27) - ecpt(23)
 c(3) = ecpt(28) - ecpt(24)
 c(4) = ecpt(30) - ecpt(26)
 c(5) = ecpt(31) - ecpt(27)
 c(6) = ecpt(32) - ecpt(28)
 
!     (R2-R1) X (R3-R2)  INTO  C(1),C(2),C(3)
 
 CALL daxb (c(1),c(4),c(1))
 c(7) = DSQRT(c(1)**2 + c(2)**2 + c(3)**2)
 IF (c(7) <= 0.0) GO TO 220
 c(2) = c(7)*DBLE(h)/24.0D0
 c(1) = 2.0D0*c(2)
 c(3) = c(2)
 c(5) = c(2)
 c(6) = c(1)
 c(7) = c(2)
 c(9) = c(2)
 c(10)= c(2)
 c(11)= c(1)
 GO TO 300
 220 WRITE  (outpt,230) ufm,necpt(1)
 230 FORMAT (a23,' 3089, ILLEGAL GEOMETRY FOR TRIANGLE ELEMENT',i14)
 nogo = .true.
 GO TO 490
 
!     IFLAG = 5, (QUADRILATERAL), 4 GRID-POINTS.
 
!               ***                                              ***
!               * 2(A2+A3+A4)  (A3+A4)     (A2+A4)      (A2+A3)    *
!               *                                                  *
!               *              2(A1+A3+A4) (A1+A4)      (A1+A3)    *
!   (4X4)  C  = *                                                  *
!               *                           2(A1+A2+A4) (A1+A2)    *
!               *     -SYM-                                        *
!               *                                       2(A1+A2+A3)*
!               ***                                              ***
 
!     R  =  XI, YI, ZI
!      I
 
!     A1 = MAG((R3-R2) X (R4-R3))
!     A2 = MAG((R4-R3) X (R1-R4))
!     A3 = MAG((R1-R4) X (R2-R1))
!     A4 = MAG((R2-R1) X (R3-R2))
 
 
!     R3-R2
 
 240 c( 1) = ecpt(30) - ecpt(26)
 c( 2) = ecpt(31) - ecpt(27)
 c( 3) = ecpt(32) - ecpt(28)
 
!     R4-R3
 
 c( 4) = ecpt(34) - ecpt(30)
 c( 5) = ecpt(35) - ecpt(31)
 c( 6) = ecpt(36) - ecpt(32)
 
!     R1-R4
 
 c( 7) = ecpt(22) - ecpt(34)
 c( 8) = ecpt(23) - ecpt(35)
 c( 9) = ecpt(24) - ecpt(36)
 
!     R2-R1
 
 c(10) = ecpt(26) - ecpt(22)
 c(11) = ecpt(27) - ecpt(23)
 c(12) = ecpt(28) - ecpt(24)
 
 
 CALL daxb (c( 1),c( 4),a1(1))
 CALL daxb (c( 4),c( 7),a2(1))
 CALL daxb (c( 7),c(10),a3(1))
 CALL daxb (c(10),c( 1),a4(1))
 
 c(1) = a1(1)*a2(1) + a1(2)*a2(2) + a1(3)*a2(3)
 c(2) = a1(1)*a3(1) + a1(2)*a3(2) + a1(3)*a3(3)
 c(3) = a1(1)*a4(1) + a1(2)*a4(2) + a1(3)*a4(3)
 IF (c(1)*c(2)*c(3) <= 0.0D0) GO TO 280
 a1(1) = DSQRT(a1(1)**2 + a1(2)**2 + a1(3)**2)
 a1(2) = DSQRT(a2(1)**2 + a2(2)**2 + a2(3)**2)
 a1(3) = DSQRT(a3(1)**2 + a3(2)**2 + a3(3)**2)
 a1(4) = DSQRT(a4(1)**2 + a4(2)**2 + a4(3)**2)
 a1(5) = a1(1) + a1(2) + a1(3) + a1(4)
 itemp = DBLE(h)/48.0D0
 DO  i = 1,4
   ic = 4*(i-1)
   DO  j = 1,4
     ij = ic + j
     IF (i == j) GO TO 260
     c(ij) = itemp*(a1(5) - a1(i) - a1(j))
     CYCLE
     260 c(ij) = itemp*(2.0D0*(a1(5) - a1(i)))
   END DO
 END DO
 GO TO 300
 280 WRITE  (outpt,290) ufm,necpt(1)
 290 FORMAT (a23,' 3090, ILLEGAL GEOMETRY FOR QUAD. ELEMENT',i14)
 nogo =.true.
 GO TO 490
 
!     HERE WHEN -C- MATRIX OF SIZE N X N IS READY FOR INSERTION (MAPING)
!     INTO MASTER OUTPUT MATRIX OF SIZE ISIZE X ISIZE.
 
 300 DO  i = 1,64
   mast(i) = 0.0D0
 END DO
 
 DO  i = 1,n
   i1 = set1(i)
   i2 = set2(i)
   DO  j = 1,n
     j1 = set1(j)
     j2 = set2(j)
     ke = cc(i,j)
     master(i1,j1) = master(i1,j1) + ke
     master(i1,j2) = master(i1,j2) - ke
     master(i2,j1) = master(i2,j1) - ke
     master(i2,j2) = master(i2,j2) + ke
   END DO
 END DO
 
!     CONDENSE (ISIZE X ISIZE) MATRIX IN (8 X 8) MASTER ARRAY INTO A
!     SINGLE STRAND FOR OUTPUT TO EMGOUT
 
 k = 0
 DO  jcol = 1,isize
   DO  irow = 1,isize
     k = k + 1
     mast(k) = master(irow,jcol)
   END DO
 END DO
 
!     OUTPUT VIA EMGOUT THE TRIANGLE IN GLOBAL FOR STIFFNESS MATRIX
 
 CALL emgout (mast(1),mast(1),k,1,dict,1,iprec)
 
!     FORM DAMPING -HEAT- IF REQUESTED.
 
 360 IF (kmbgg(3) == 0) GO TO 490
 inflag = 4
 eltemp = avgtmp
 matid  = matflg
 IF (matid == 0) GO TO 490
 CALL hmat (necpt)
 cp = hx
 IF (cp == 0.0) GO TO 490
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 380
   CASE (    2)
     GO TO 390
   CASE (    3)
     GO TO 400
   CASE (    4)
     GO TO 410
   CASE (    5)
     GO TO 420
   CASE (    6)
     GO TO 390
   CASE (    7)
     GO TO 390
 END SELECT
 
!     IFLAG = 1, (POINT), 1 GRID-POINT.  (1 X 1)  C = CP* AF
 
 380 c(1) = cp
 c(2) = af
 c(1) = c(1)*c(2)
 GO TO 440
 
!     IFLAG = 2, (LINE OR ELLIPTIC CYL. )
!             2 GRID POINTS           CP*AF*L*      *
!                                  C = ------*1 , 1 *
!                                        2   *      *
 
 390 c(1) = cp
 c(2) = af
 c(3) = ecpt(26) - ecpt(22)
 c(4) = ecpt(27) - ecpt(23)
 c(5) = ecpt(28) - ecpt(24)
 c(1) = c(1)*c(2)*DSQRT(c(3)**2 + c(4)**2 + c(5)**2)/2.0D0
 c(2) = c(1)
 GO TO 440
 
!     IFLAG = 3, (REVOLUTION), 2 GRID-POINTS
!                                               CP*PI*L *              *
!                                           C = ------- *2X +X , 2X +X *
!                                                  3    *  1  2    2  1*
 
 400 c(1) = cp
 c(2) = pi
 c(3) = ecpt(26) - ecpt(22)
 c(4) = ecpt(28) - ecpt(24)
 
!     NOTE Y2 AND Y1 ARE 0 FOR REVOLUTION ELEMENT.
 
 c(1) = c(1)*c(2)*DSQRT(c(3)**2 + c(4)**2)/3.0D0
 c(2) = c(1)*DBLE(ecpt(22) + 2.0*ecpt(26))
 c(1) = c(1)*DBLE(2.0*ecpt(22) + ecpt(26))
 GO TO 440
 
!     IFLAG = 4, (TRIANGLE), 3 GRID-POINTS.
!                                          CP*A *         *
!                                      C = ---- * 1, 1, 1 *
!                                           3   *         *
 
 
!     COMPUTE AREA -A- OF TRIANGLE   GET R2-R1 AND R3-R2
 
 410 c(1) = ecpt(26) - ecpt(22)
 c(2) = ecpt(27) - ecpt(23)
 c(3) = ecpt(28) - ecpt(24)
 c(4) = ecpt(30) - ecpt(26)
 c(5) = ecpt(31) - ecpt(27)
 c(6) = ecpt(32) - ecpt(28)
 
!     (R2-R1) X (R3-R2)  INTO  C(1),C(2),C(3)
 
 CALL daxb (c(1),c(4),c(1))
 c(7) = DSQRT(c(1)**2 + c(2)**2 + c(3)**2)
 c(1) = c(7)*DBLE(cp)/6.0D0
 c(2) = c(1)
 c(3) = c(1)
 GO TO 440
 
!     IFLAG = 5, (QUADRILATERAL), 4 GRID-POINTS.
 
!                                CP *                                  *
!                            C = -- * A +A +A , A +A +A , A +A +A , ETC*
!                                6  *  2  3  4   3  4  1   4  1  2     *
 
!     R  =  XI, YI, ZI
!      I
 
!     A1 = MAG((R3-R2) X (R4-R3))
!     A2 = MAG((R4-R3) X (R1-R4))
!     A3 = MAG((R1-R4) X (R2-R1))
!     A4 = MAG((R2-R1) X (R3-R2))
 
 
!     R3-R2
 
 420 c( 1) = ecpt(30) - ecpt(26)
 c( 2) = ecpt(31) - ecpt(27)
 c( 3) = ecpt(32) - ecpt(28)
 
!     R4-R3
 
 c( 4) = ecpt(34) - ecpt(30)
 c( 5) = ecpt(35) - ecpt(31)
 c( 6) = ecpt(36) - ecpt(32)
 
!     R1-R4
 
 c( 7) = ecpt(22) - ecpt(34)
 c( 8) = ecpt(23) - ecpt(35)
 c( 9) = ecpt(24) - ecpt(36)
 
!     R2-R1
 
 c(10) = ecpt(26) - ecpt(22)
 c(11) = ecpt(27) - ecpt(23)
 c(12) = ecpt(28) - ecpt(24)
 
 
 CALL daxb (c( 1),c( 4),a1(1))
 CALL daxb (c( 4),c( 7),a2(1))
 CALL daxb (c( 7),c(10),a3(1))
 CALL daxb (c(10),c( 1),a4(1))
 
 a1(1) = DSQRT(a1(1)**2 + a1(2)**2 + a1(3)**2)
 a1(2) = DSQRT(a2(1)**2 + a2(2)**2 + a2(3)**2)
 a1(3) = DSQRT(a3(1)**2 + a3(2)**2 + a3(3)**2)
 a1(4) = DSQRT(a4(1)**2 + a4(2)**2 + a4(3)**2)
 a1(5) = a1(1) + a1(2) + a1(3) + a1(4)
 itemp = DBLE(cp)/12.0D0
 DO  i = 1,4
   c(i) = itemp*(a1(5) - a1(i))
 END DO
 GO TO 440
 
!     HERE WHEN DIAGONAL C MATRIX OF SIZE 1 X N IS READY FOR INSERTION
!     (MAPING) INTO MASTER DIAGONAL OUTPUT MATRIX OF SIZE 1 X ISIZE.
 
 440 DO  i = 1,8
   mast(i) = 0.0D0
 END DO
 
 DO  i = 1,n
   i1 = set1(i)
   i2 = set2(i)
   me = c(i)
   mast(i1) = mast(i1) + me
   mast(i2) = mast(i2) + me
 END DO
 
!     OUTPUT VIA EMGOUT THE DIAGONAL MATRIX IN GLOBAL
 
 dict(2) = 2
 CALL emgout (mast(1),mast(1),isize,1,dict,3,iprec)
 GO TO 490
 
!     LOGIC ERROR
 
 470 WRITE  (outpt,480) sfm,imhere,necpt(1),sils
 480 FORMAT (a25,' 3037 FROM HBDYD.', /5X,  &
     'LOGIC ERROR,  IMHERE =',i5,'  ELEMENT ID = ',i10, /5X, 'SILS =',8I10)
 nogo = .true.
 490 RETURN
END SUBROUTINE hbdyd
