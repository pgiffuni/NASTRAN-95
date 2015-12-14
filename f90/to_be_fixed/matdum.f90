SUBROUTINE matdum (ia,iprc,npl,nout)
     
!     THIS ROUTINE IS CALLED ONLY BY MATPRN TO PRINT UP TO 5 MATRICES
 
!     IF IPRC = 0, MATRICES ARE PRINTED IN THEIR ORIG. PRECISION FORMAT
!     IF IPRC = 1, MATRICES ARE PRINTED IN SINGLE PRECISION E FORMAT
!     IF IPRC = 2, MATRICES ARE PRINTED IN DOUBLE PRECISION D FORMAT
!     IF IPRC =-1, ONLY THE DIAGONAL ELEMENTS OF THE MATRICES ARE
!                  PRINTED IN THEIR ORIG. PRECISION FORMAT
 
!     INPUT MATRIX IA(1) CAN BE IN S.P., D.P., S.P.CMPLX OR D.P.CMPLX
 
!     NPL IS THE NO. OF DATA VALUES PRINTED PER OUTPUT LINE
!     FOR S.P. REAL  DEFAULT IS 8, MAX IS 14
!     FOR D.P. REAL  DEFAULT IS 6, MAX IS 12
!     EVEN NUMBER ONLY FOR COMPLEX
 
!     P3, P4, P5 ARE PRINTOUT CONTROLS
!     P3 = m, MATRIX COLUMNS, 1 THRU m, WILL BE PRINTED.
!             DEFAULT = 0, ALL MATRIX COLUMNS WILL BE PRINTED.
!        =-m, SEE P4 = -n
!     P4 = n, LAST n MATRIX COLUMNS ARE PRINTED. DEFAULT = 0
!        =-n, AND P3 = -m, EVERY OTHER n MATRIX COLUMNS WILL BE PRINTED,
!             STARTIN FROM COLUMN m.
!     P5 = k, EACH PRINTED COLUMN WILL NOT EXCEED k LINES LONG AND THE
!             REMAINING DATA WILL BE OMITTED.
!     NOUT = P6, FORTRAN UNIT (SEE MATPRN)
 
 
 INTEGER, INTENT(IN)                      :: ia(7)
 INTEGER, INTENT(IN)                      :: iprc
 INTEGER, INTENT(IN)                      :: npl
 INTEGER, INTENT(IN OUT)                  :: nout
 LOGICAL :: jump
 INTEGER :: sysbuf,iblnk,p12,p3,p4,p5,px(5),icol(1),FILE(14)
 DOUBLE PRECISION :: dcol(1)
 DIMENSION  TYPE(10),FORM(18)
 CHARACTER (LEN=15) :: rfmt,fmtr(2,7)
 CHARACTER (LEN=35) :: cfmt,fmtc(2,4)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /BLANK /  p12(2),p3,p4,p5
!ZZ   COMMON /ZZTBPR/  COL(1)
 COMMON /zzzzzz/  col(20000)
 COMMON /unpakx/  it,k,l,incr
 COMMON /system/  sysbuf,iout,inx(6),nlpp,inx1(2),line
 COMMON /output/  head1(96),head2(96)
 EQUIVALENCE      (col(1),dcol(1),icol(1)), (iblnk,BLANK)
 DATA    TYPE  /  4HS.p.,4HREAL,4HD.p.,4HREAL,4HCOMP,4HLEX ,4HCMP ,  &
     4HD.p.,4HILL ,4HDEFN       /
 DATA    FORM  /  4HSQUA,4HRE  ,4HRECT,4HANG ,4HDIAG,4HONAL,4HLOW ,  &
     4HTRI ,4HUPP ,4HTRI ,4HSYMM,4HETRC,4HVECT,4HOR  ,  &
     4HIDEN,4HTITY,4HILL ,4HDEFN/
 DATA    BLANK ,  xmatr ,xix   ,cont  ,xinue ,ddx   /  &
     4H    ,  4HMATR,4HIX  ,4HCONT,4HINUE,4HD   /
 DATA    FILE  /  4HUT1 ,4HUT2 ,4HN/a ,4HINPT,4HINP1,4HINP2,4HINP3,  &
     4HINP4,4HINP5,4HINP6,4HINP7,4HINP8,4HINP9,4HINPT/
 DATA    fmtr  /  '(1X,1P, 8E16.8)',  '(1X,1P,6D21.12)',  &
     '(1X,1P, 9E14.6)',  '(1X,1P,7D18.10)',  &
     '(1X,1P,10E13.5)',  '(1X,1P, 8D16.8)',  &
     '(1X,1P,11E11.3)',  '(1X,1P, 9D14.6)',  &
     '(1X,1P,12E10.2)',  '(1X,1P,10D13.4)',  &
     '(1X,1P,13E10.2)',  '(1X,1P,11D11.3)', '(1X,1P,14E 9.1)',  '(1X,1P,12D10.2)'/
 DATA    fmtc  /  '(4(1X,1P,E14.7,1HR,  1P,E15.7,1HI))',  &
     '(3(1X,1P,D20.13,1HR,1P,D21.13,1HI))',  &
     '(5(1X,1P,E11.4,1HR,  1P,E12.4,1HI))',  &
     '(4(1X,1P,D14.7,1HR,  1P,D15.7,1HI))',  &
     '(6(1X,1P,E 9.2,1HR,  1P,E10.2,1HI))',  &
     '(5(1X,1P,D11.4,1HR,  1P,D12.4,1HI))',  &
     '(7(1X,1P,E 7.0,1HR,  1P,E 8.0,1HI))', '(6(1X,1P,D 9.2,1HR,  1P,D10.2,1H)))'/
 
 namea= ia(1)
 ncol = ia(2)
 nrow = ia(3)
 IF   = ia(4)
 it   = ia(5)
 IF (IF /= 7) GO TO 5
 
!     ROW VECTOR
!     A ROW OF MATRIX ELEMENTS STORED IN COLUMN FORMAT
!     INTERCHANGE ROW AND COLUMN FOR PRINTING
 
 j    = ncol
 ncol = nrow
 nrow = j
 5 IF (it <= 0 .OR. it > 4) it = 5
 IF (IF <= 0 .OR. IF > 8) IF = 9
 IF (nout /= iout) WRITE (iout,7) uim,FILE(nout-10),nout
 7 FORMAT (a29,', MATRIX PRINTOUT SAVED IN ',a4,' (FORTRAN UNIT',i4, 1H))
 IF (iprc == -1) GO TO 80
 
!     SET UP FORMAT FOR OUTPUT PRINT LINE
 
 j = iprc
 IF (it >= 3) j = iprc + 2
 SELECT CASE ( j )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 30
   CASE (    4)
     GO TO 40
 END SELECT
 10 j = npl - 7
 rfmt = fmtr(1,j)
 GO TO 50
 20 j = npl - 5
 rfmt = fmtr(2,j)
 GO TO 50
 30 j = (npl/2) - 3
 cfmt = fmtc(1,j)
 GO TO 50
 40 j = (npl/2) - 2
 cfmt = fmtc(2,j)
 50 npl1 = npl - 1
 
!     SET UP P3 AND P4 PRINTOUT OPTIONS
 
 mm = p3
 nn = ia(2)
 IF (p3 <= 0) mm = ia(2)
 IF (p4 < 0) GO TO 60
 jump = .false.
 nn = ia(2) - p4
 GO TO 70
 60 jump = .true.
 jmp4 = -p4
 jmp3 = IABS(p3)
 IF (p3 == 0) jmp3 = 1
 70 nplp5 = ia(3)
 IF (p5 /= 0) nplp5 = npl*p5
!WKBI SPR 93013
 IF ( it > 2 ) nplp5 = 2*nplp5
 
 80 DO  i = 1,96
   head2(i) = BLANK
 END DO
 head2(1) = xmatr
 head2(2) = xix
 head2(6) = cont
 head2(7) = xinue
 head2(8) = ddx
 lcol = korsz(col) - sysbuf
 incr = 1
 CALL gopen (namea,col(lcol+1),0)
 CALL page1
 CALL fname (namea,head2(3))
 WRITE (nout,90) head2(3),head2(4),namea,TYPE(2*it-1),TYPE(2*it),  &
     ncol,nrow,FORM(2*IF-1),FORM(2*IF)
 90 FORMAT (1H0,6X,7HMATRIX ,2A4,11H (gino NAME,i4,2H ),6H is a ,2A4,  &
     1X,i6,10H column x ,i6,5H row ,2A4,8H matrix.)
 IF (it == 5 .OR. ncol == 0 .OR. nrow == 0) GO TO 460
 
!     IF = 3, DIAGONAL MATRIX
!        = 7, ROW VECTOR
!        = 8, IDENTITY MATRIX
 
 IF (IF-8) 100,440,460
 100 IF (iprc == -1) GO TO 510
 IF (IF /= 3 .AND. IF /= 7) GO TO 110
 ncol = 1
 nrow = ia(3)
 110 inull= 0
 ASSIGN 150 TO ihop
 jj = 1
 120 k  = 0
 l  = 0
 CALL unpack (*330,namea,col)
 IF (jj <= mm .OR. jj >= nn) GO TO 130
 k  = nn - mm - 1
 jj = jj + k
 IF (jj > ncol) GO TO 340
 CALL skprec (namea,k)
 GO TO 340
 130 IF (.NOT.jump) GO TO 140
 IF (MOD(jj,jmp4) /= jmp3) GO TO 340
 140 IF (inull == 1) GO TO 490
 150 nrow = l - k + 1
 SELECT CASE ( IF )
   CASE (    1)
     GO TO 160
   CASE (    2)
     GO TO 160
   CASE (    3)
     GO TO 360
   CASE (    4)
     GO TO 160
   CASE (    5)
     GO TO 160
   CASE (    6)
     GO TO 160
   CASE (    7)
     GO TO 380
 END SELECT
 160 WRITE (nout,170) jj,k,l
 line = line + 3
 IF (line >= nlpp) CALL page
 170 FORMAT (8H0COLUMN ,i6,5X,6H rows ,i6,6H thru ,i6,5X,50(1H-),/,1H )
 IF (it > 2) nrow = 2*nrow
 180 k = 0
 190 j = k + 1
 IF (j >  nrow) GO TO 340
 k = j + npl1
 IF (k >  nrow) k = nrow
 IF (k > nplp5) GO TO 340
 kj = k - j
 SELECT CASE ( it )
   CASE (    1)
     GO TO 210
   CASE (    2)
     GO TO 240
   CASE (    3)
     GO TO 270
   CASE (    4)
     GO TO 300
 END SELECT
 
 200 ln = (kj+npl1)/npl
 line = line + ln
 IF (line >= nlpp) CALL page
 GO TO 190
 
!     REAL SINGLE PRECISION
 
 210 IF (iprc == 2) GO TO 220
 WRITE  (nout,rfmt) (col(i),i=j,k)
! 215 FORMAT (1X,1P,10E13.5)
!     LN = (KJ+10)/10
 GO TO 200
 220 i = k
 DO  ln = j,k
   dcol(i) = col(i)
   i = i - 1
 END DO
 
!     REAL DOUBLE PRECISION
 
 240 IF (iprc == 1) GO TO 250
 WRITE  (nout,rfmt) (dcol(i),i=j,k)
! 245 FORMAT (1X,1P,8D16.8)
!     LN = (KJ+8)/8
 GO TO 200
 250 DO  i = j,k
   col(i) = dcol(i)
 END DO
 GO TO 210
 
!     COMPLEX SINGLE
 
 270 IF (iprc == 2) GO TO 280
 WRITE  (nout,cfmt) (col(i),i=j,k)
! 275 FORMAT (1X,5(1P,E12.4,1HR,1P,E12.4,1HI))
!     LN = (KJ+10)/10
 GO TO 200
 280 i = k
 DO  ln = j,k
   dcol(i) = col(i)
   i = i - 1
 END DO
 
!     COMPLEX DOUBLE
 
 300 IF (iprc == 1) GO TO 310
 WRITE  (nout,cfmt) (dcol(i),i=j,k)
! 305 FORMAT (1X,4(1P,D15.8,1HR,1P,D15.8,1HI))
!     LN = (KJ+8)/8
 GO TO 200
 310 DO  i = j,k
   col(i) = dcol(i)
 END DO
 GO TO 270
 
 330 IF (inull == 1) GO TO 340
 inull = 1
 ibegn = jj
 340 jj    = jj + 1
 IF (jj  <= ncol) GO TO 120
 ASSIGN 350 TO ihop
 IF (inull == 1) GO TO 490
 350 CALL CLOSE (namea,1)
 GO TO 400
 360 WRITE (nout,370) k,l
 line = line + 2
 370 FORMAT ('0DIAGONAL ELEMENTS FOR COLUMNS',i6,4H TO ,i6,4H are,///)
 GO TO 180
 380 WRITE (nout,390) k,l
 line = line + 2
 390 FORMAT ('0ROW ELEMENTS FOR COLUMNS',i6,4H TO ,i6,4H are,///)
 GO TO 180
 400 WRITE  (nout,410) ia(6)
 410 FORMAT ('0THE NUMBER OF NON-ZERO WORDS IN THE LONGEST RECORD =', i8)
 ia7a = ia(7) / 100
 ia7c = ia(7) - 100*ia7a
 ia7b = ia7c / 10
 ia7c = ia7c - 10*ia7b
 WRITE  (nout,420) ia7a,ia7b,ia7c
 420 FORMAT ('0THE DENSITY OF THIS MATRIX IS ',i3,1H.,2I1,' PERCENT.')
 GO TO 750
 
 440 WRITE  (nout,450)
 450 FORMAT ('0IDENTITY MATRIX')
 460 CALL CLOSE (namea,1)
 
!     FUNNY MATRIX - SAVE MODULE PARAMETERS AND TABLE PRINT IT
 
 DO  i = 1,5
   px(i) = p12(i)
 END DO
 p12(1) = iblnk
 p12(2) = iblnk
 p3     = 3
 p4     = 3
 CALL tabprt (namea)
 DO  i = 1,5
   p12(i) = px(i)
 END DO
 GO TO 750
 490 ifin = jj - 1
 WRITE (nout,500) ibegn,ifin
 inull = 0
 line  = line + 2
 IF (line >= nlpp) CALL page
 500 FORMAT ('0COLUMNS ',i7,6H thru ,i7,' ARE NULL.')
 GO TO ihop, (150,350)
 
!     PRINT ONLY THE DIAGONAL ELEMENTS, IPRC = -1
!     TO CHECKOUT THE DIAGONALS FOR POSSIBLE MATRIX SINGULARITY
 
 510 WRITE  (nout,520)
 520 FORMAT (/23X,'(ELEMENTS ON DIAGONAL ONLY)')
 IF (ncol /= nrow) WRITE (nout,530)
 530 FORMAT (23X,'*** MATRIX IS NOT SQUARE ***')
 WRITE  (nout,540)
 540 FORMAT (1X)
 nn = MIN0(ncol,nrow)
 jj = 0
 DO  i = 1,nn
   k  = i
   l  = i
   CALL unpack (*550,namea,col(jj+1))
   GO TO 570
   550 DO  j = 1,it
     col(jj+j) = 0.0
   END DO
   570 jj = jj + it
 END DO
 CALL CLOSE (namea,1)
 SELECT CASE ( it )
   CASE (    1)
     GO TO 580
   CASE (    2)
     GO TO 600
   CASE (    3)
     GO TO 620
   CASE (    4)
     GO TO 640
 END SELECT
 580 WRITE  (nout,590) (col(j),j=1,jj)
 590 FORMAT (1X,1P,10E13.6)
 GO TO 660
 600 jj = jj/2
 WRITE  (nout,610) (dcol(j),j=1,jj)
 610 FORMAT (1X,1P,10D13.6)
 GO TO 660
 620 WRITE  (nout,630) (col(j),j=1,jj)
 630 FORMAT ((1X,5(1P,e12.5,1HR,1P,e12.5,1HI)))
 GO TO 660
 640 jj = jj/2
 WRITE  (nout,650) (dcol(j),j=1,jj)
 650 FORMAT ((1X,5(1P,d12.5,1HR,1P,d12.5,1HI)))
 660 kj = it
 IF (it >= 3) kj = it - 2
 nn = 0
 mm = 1
 loop680:  DO  j = 1,jj
   ln = mm + it - 1
   DO  i = mm,ln
     IF (col(i) /= 0.0) GO TO 680
   END DO
   nn = nn + 1
   icol(nn) = j
   mm = mm + kj
 END DO loop680
 IF (nn == 0) GO TO 710
 mm = MIN0(nn,200)
 WRITE  (nout,690) (icol(i),i=1,mm)
 690 FORMAT ('0*** ZERO DIAGONALS IN THE FOLLOWING COLUMNS -', /,(1X,20I6))
 IF (nn > 200) WRITE (nout,700)
 700 FORMAT (' ...AND MORE')
 GO TO 730
 710 WRITE  (nout,720)
 720 FORMAT ('0*** NO ZERO ON DIAGONALS')
 730 WRITE  (nout,740) ia
 740 FORMAT (/5X,'GINO FILE',i5,'   TRAILER =',6I7)
 line = line + nlpp
 
 750 RETURN
END SUBROUTINE matdum
