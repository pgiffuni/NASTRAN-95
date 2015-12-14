SUBROUTINE gpstg
     
!     THIS SUBROUTINE GENERATES THE GRID POINT SINGULARITY TABLE
!     BY EXAMINING THE TRANSLATIONAL AND ROTATIONAL 3 X 3
!     SUBMATRICES ALONG THE LEADING DIAGONAL OF THE INPUT
!     STIFFNESS MATRIX
 
 DIMENSION   iarray(8), isubnm(2)
 
 INTEGER :: gpst     , ttlwds
 
 DOUBLE PRECISION :: b(9), fl(3), d
 DOUBLE PRECISION :: m(3), r(3) , temp, fm, fr, det, const, dtol
 
 COMMON /gpstgx/ gpst , igpst, npvt, nsing , ibuf2
 COMMON /gpstgy/ d(18)
 COMMON /system/ isys(69)    , tolel
 COMMON /zzzzzz/ iz(1)
 
 EQUIVALENCE (iorder, iarray(1)), (nwds, iarray(2))
 
 DATA isubnm / 4HGPST,4HG   /
 
 dtol = tolel
 
! AT THIS POINT, BOTH TRANSLATIONAL AND ROTATIONAL DIAGONAL 3X3 S ARE
! STORED IN THE D ARRAY.  HENCE WE PROCESS THEM.
 
 ip = npvt - 1
 ASSIGN 470 TO igoto
 ASSIGN 20 TO iback
 DO  i = 1,9
   b(i) = d(i)
 END DO
 GO TO 90
 20 DO  i = 1,9
   b(i) = d(i+9)
 END DO
 
! INSURE THE SYMMETRY OF THE B MATRIX
 
 IF (b(2) /= 0.0D0 .AND. b(4) /= 0.0D0) GO TO 40
 b(2) = 0.0D0
 b(4) = 0.0D0
 GO TO 50
 40 temp = (b(2) + b(4)) / 2.0D0
 b(2) = temp
 b(4) = temp
 50 IF (b(3) /= 0.0D0 .AND. b(7) /= 0.0D0) GO TO 60
 b(3) = 0.0D0
 b(7) = 0.0D0
 GO TO 70
 60 temp = (b(3) + b(7)) / 2.0D0
 b(3) = temp
 b(7) = temp
 70 IF (b(6) /= 0.0D0 .AND. b(8) /= 0.0D0) GO TO 80
 b(6) = 0.0D0
 b(8) = 0.0D0
 GO TO 90
 80 temp = (b(6) + b(8)) / 2.0D0
 
! SCALE THE MATRIX BY DIVIDING EACH ELEMENT OF B BY THE LARGEST ELEMENT.
! IF THE LARGEST ELEMENT IS NON-POSITIVE, THE SINGULARITY IS OF ORDER 3.
 
 90 temp = b(1)
 DO  i = 2,9
   IF (b(i) > temp) temp = b(i)
 END DO
 IF (temp <= 0.0D0) GO TO 430
 DO  i = 1,9
   b(i) = b(i) / temp
 END DO
 
! FIND THE SQUARES OF THE MAGNITUDES OF THE VECTORS OF THE ROWS OF THE
! B MATRIX.
 
 iorder = 0
 j = 0
 DO  i = 1,9,3
   j = j + 1
   fl(j) = b(i)**2 + b(i+1)**2 + b(i+2)**2
   IF (fl(j) == 0.0D0) iorder = iorder + 1
 END DO
 IF (iorder == 2)  GO TO 410
 IF (iorder == 0)  GO TO 260
 
! AT THIS POINT ONE AND ONLY ONE FL(I) IS ZERO.
 
 DO  i = 1,3
   isave = i
   IF (fl(i) == 0. 0D0) THEN
      SELECT CASE ( isave )
       CASE (    1)
         GO TO 140
       CASE (    2)
         GO TO 150
       CASE (    3)
         GO TO 160
     END SELECT
   END IF
 END DO
 CALL mesage (-30,26,isubnm)
 140 fm = b(5) * b(9)  -  b(6) * b(8)
 fr = DSQRT( (b(5)**2  +  b(6)**2)  *  (b(8)**2  +  b(9)**2) )
 GO TO 170
 150 fm = b(1) * b(9)  -  b(3) * b(7)
 fr = DSQRT( (b(1)**2  +  b(3)**2)  *  (b(7)**2  +  b(9)**2) )
 GO TO 170
 160 fm = b(1) * b(5)  -  b(2) * b(4)
 fr = DSQRT( (b(1)**2  +  b(2)**2)  *  (b(4)**2  +  b(5)**2) )
 170 IF ( fm   == 0.0D0 ) GO TO 180
 IF ( fr   <= 0.0D0 ) GO TO 250
 IF ( fm/fr >= dtol ) GO TO 250
 
! HERE WE HAVE THAT THE ORDER OF THE SINGULARITY IS 2.
 
 180 iorder = 2
 nwds   = 0
 ttlwds  = 2
 SELECT CASE ( isave )
   CASE (    1)
     GO TO 190
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 210
 END SELECT
 190 k1   = 5
 k2   = 9
 inc1 = 1
 inc2 = 3
 inc3 = 2
 GO TO 220
 200 k1   = 1
 k2   = 9
 inc1 = 2
 inc2 = 3
 inc3 = 1
 GO TO 220
 210 k1   = 1
 k2   = 5
 inc1 = 3
 inc2 = 2
 inc3 = 1
 220 IF (b(k1) <= 0.0D0  .AND.  b(k2) <= 0.0D0) GO TO 430
 IF (b(k1) <= 0.0D0) GO TO 230
 nwds      = 2
 ttlwds     = 4
 iarray(3) = ip + inc1
 iarray(4) = ip + inc2
 ipoint    = 5
 GO TO 240
 230 ipoint = 3
 240 IF (b(k2) <= 0.0D0) GO TO 440
 nwds = nwds + 2
 ttlwds = ttlwds + 2
 iarray(ipoint)   = ip + inc1
 iarray(ipoint+1) = ip + inc3
 GO TO 440
 
! AT THIS POINT WE HAVE THAT ONE AND ONLY ONE FL IS ZERO BUT THAT ORDER
! OF THE SINGULARITY IS 1.
 
 250 iorder    = 1
 nwds      = 1
 ttlwds    = 3
 iarray(3) = ip + isave
 GO TO 440
 
! AT STATEMENT NO. 260, WE HAVE THAT ALL THE FL(I) ARE .GT. 0.0D0, SO
! THAT THE DETERMINANT, DET, OF B MUST BE COMPUTED.
 
 260 det = b(1) * ( b(5)*b(9) - b(6)*b(8) )  &
     - b(2) * ( b(4)*b(9) - b(6)*b(7) ) + b(3) * ( b(4)*b(8) - b(5)*b(7) )
 const = 0.05D0*dtol * fl(1) * fl(2) * fl(3)
 IF (det > const) GO TO 460
 
! COMPUTE M(I) AND R(I)
 
 m(1) = b(5) * b(9) - b(6) * b(8)
 m(2) = b(1) * b(9) - b(3) * b(7)
 m(3) = b(1) * b(5) - b(2) * b(4)
 r(1) = DSQRT ( b(5)**2 + b(6)**2 ) * DSQRT ( b(8)**2 + b(9)**2 )
 r(2) = DSQRT ( b(1)**2 + b(3)**2 ) * DSQRT ( b(7)**2 + b(9)**2 )
 r(3) = DSQRT ( b(1)**2 + b(2)**2 ) * DSQRT ( b(4)**2 + b(5)**2 )
 
! FIND I1,J1,K1 SUCH THAT M(I1)/R(I1) .GE. M(J1)/R(J1) .GE. M(K1)/R(K1)
 
 i1 = 1
 j1 = 2
 k1 = 3
 IF (m(1)*r(2) >= m(2)*r(1)) GO TO 270
 i1 = 2
 j1 = 1
 270 IF (m(i1)*r(k1) >= m(k1)*r(i1)) GO TO 280
 itemp = i1
 i1    = k1
 k1    = itemp
 280 IF (m(j1)*r(k1) >= m(k1)*r(j1)) GO TO 290
 itemp = j1
 j1    = k1
 k1    = itemp
 290 IF (m(i1) >= r(i1)*dtol) GO TO 400
 
! HERE THE SINGULARITY IS OF ORDER 2.
 
 nwds   = 0
 ttlwds = 2
 iorder = 2
 
! FIND II, JJ, KK SUCH THAT B(II) .GE. B(JJ) .GE. B(KK)
 
 ii = 1
 jj = 5
 kk = 9
 IF (b(1) >= b(5)) GO TO 300
 ii = 5
 jj = 1
 300 IF (b(ii) >= b(kk)) GO TO 310
 itemp = ii
 ii    = kk
 kk    = itemp
 310 IF (b(jj) >= b(kk)) GO TO 320
 itemp = jj
 jj    = kk
 kk    = itemp
 320 ll    = ii
 kount = 0
 ipoint= 3
 330 IF (b(ll) <= 0.0D0) GO TO 440
 nwds   = nwds + 2
 ttlwds = ttlwds + 2
 IF (ll - 5 < 0) THEN
   GO TO   340
 ELSE IF (ll - 5 == 0) THEN
   GO TO   350
 ELSE
   GO TO   360
 END IF
 340 inc1 = 2
 inc2 = 3
 GO TO 370
 350 inc1 = 1
 inc2 = 3
 GO TO 370
 360 inc1 = 1
 inc2 = 2
 370 iarray(ipoint)   = ip + inc1
 iarray(ipoint+1) = ip + inc2
 ipoint = ipoint + 2
 kount  = kount  + 1
 IF (kount - 2 < 0) THEN
   GO TO   380
 ELSE IF (kount - 2 == 0) THEN
   GO TO   390
 ELSE
   GO TO   440
 END IF
 380 ll = jj
 GO TO 330
 390 ll = kk
 GO TO 330
 
! AT THIS POINT THE SINGULARITY IS OF ORDER 1.
 
 400 iorder = 1
 nwds   = 1
 ttlwds = 3
 iarray(3) = ip + i1
 IF (m(j1) < r(j1)*dtol) GO TO 440
 nwds   = 2
 ttlwds = 4
 iarray(4) = ip + j1
 IF (m(k1) < r(k1)*dtol) GO TO 440
 nwds   = 3
 ttlwds = 5
 iarray(5) = ip + k1
 GO TO 440
 
! AT THIS POINT 2 ROWS OF THE B MATRIX ARE IDENTICALLY ZERO.
 
 410 nwds   = 2
 ttlwds = 4
 ipoint = 2
 DO  i = 1,3
   IF (fl(i) /= 0.0D0) CYCLE
   ipoint = ipoint + 1
   iarray(ipoint) = ip + i
 END DO
 GO TO 440
 
! THE SINGULARITY IS OF ORDER 3
 
 430 iorder = 3
 nwds   = 3
 ttlwds = 5
 iarray(3) = ip + 1
 iarray(4) = ip + 2
 iarray(5) = ip + 3
 
! WRITE IARRAY ON THE GPST FILE.
 
 440 IF (igpst == 1) GO TO 450
 igpst = 1
 CALL gopen (gpst,iz(ibuf2),1)
 450 nsing = nsing + 1
 CALL WRITE (gpst,iarray,ttlwds,0)
 460 GO TO igoto, (470,480)
 470 ASSIGN 480 TO igoto
 ip = ip + 3
 GO TO iback, (20,430)
 480 CONTINUE
 RETURN
END SUBROUTINE gpstg
