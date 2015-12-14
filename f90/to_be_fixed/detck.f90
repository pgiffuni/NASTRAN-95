SUBROUTINE detck (jarg,ifgpst,npvt)
     
!     COMMENTS FROM G.CHAN/UNISYS, 5/1991,
!     THIS ROUTINE WAS NAMED DETCKX BEFORE, WHICH HAD NOT BEEN TESTED.
!     THE ONE THAT USED TO BE DETCK APPEARS TO BE AN OLDER VERSION, AND
!     SHOULD BE REPLACED BY THIS ONE, IF THIS ONE WORKS
 
!     THIS ROUTINE GENERATES THE GRID POINT SINGULARITY TABLE BY
!     EXAMINING THE TRANSLATIONAL AND DIAGONAL 3 X 3 SUBMATRICES OF THE
!     KGG MATRIX.
!     IF JARG = 0, THE PIVOT POINT HAS ELEMENTS ATTACHED TO IT.
!     IF JARG =-1, THE PIVOT IS A SCALAR POINT AND NO ELEMENTS ARE
!                  CONNECTED TO IT.
!     IF JARG = 1, THE PIVOT POINT IS A GRID POINT AND NO ELEMENTS ARE
!                  CONNECTED TO IT.
 
 
 INTEGER, INTENT(IN)                      :: jarg
 INTEGER, INTENT(IN OUT)                  :: ifgpst
 INTEGER, INTENT(IN)                      :: npvt
 INTEGER :: tnwds,iarray(8),iz(1),back,NAME(2)
 DOUBLE PRECISION :: d,b,dz,fl,r,m,temp,fm,fr,det,const,dtol
 COMMON /ma1xx /  d(18),b(9),dz(1),fl(3),r(3),m(3)
 COMMON /system/  isys(69),tolel
 EQUIVALENCE      (iz(1),dz(1)),(iarray(1),iorder),(iarray(2),nwds)
 DATA    NAME  /  4HDETC,4HK   /,  neor  / 0 /
 
 dtol = tolel
 iarg = jarg
 IF (iarg < 0) THEN
   GO TO    10
 ELSE IF (iarg == 0) THEN
   GO TO    20
 ELSE
   GO TO    25
 END IF
 10 iorder = 1
 nwds   = 1
 iarray(3) = npvt
 CALL WRITE (ifgpst,iarray(1),3,neor)
 RETURN
 
!     AT THIS POINT, BOTH TRANSLATIONAL AND ROTATIONAL DIAGONAL 3X3 S
!     ARE STORED IN THE D ARRAY.  HENCE WE PROCESS THEM.
 
 20 CONTINUE
 25 CONTINUE
 ip = npvt - 1
 ASSIGN 450 TO igoto
 IF (iarg /= 1) GO TO 30
 ASSIGN 425 TO back
 GO TO 425
 30 ASSIGN 50 TO back
 DO  i = 1,9
   b(i) = d(i)
 END DO
 GO TO 90
 50 DO  i = 1,9
   b(i) = d(i+9)
 END DO
 
!     INSURE THE SYMMETRY OF THE B MATRIX
 
 IF (b(2) /= 0.0D0 .AND. b(4) /= 0.0D0) GO TO 65
 b(2) = 0.0D0
 b(4) = 0.0D0
 GO TO 70
 65 temp = (b(2) + b(4))/2.0D0
 b(2) = temp
 b(4) = temp
 70 IF (b(3) /= 0.0D0 .AND. b(7) /= 0.0D0) GO TO 75
 b(3) = 0.0D0
 b(7) = 0.0D0
 GO TO 80
 75 temp = (b(3) + b(7))/2.0D0
 b(3) = temp
 b(7) = temp
 80 IF (b(6) /= 0.0D0 .AND. b(8) /= 0.0D0) GO TO 85
 b(6) = 0.0D0
 b(8) = 0.0D0
 GO TO 90
 85 temp = (b(6) + b(8))/2.0D0
 
!     SCALE THE MATRIX BY DIVIDING EACH ELEMENT OF B BY THE LARGEST
!     ELEMENT. IF THE LARGEST ELEMENT IS NON-POSITIVE, THE SINGULARITY
!     IS OF ORDER 3.
 
 90 temp = b(1)
 DO  i = 2,9
   IF (b(i) > temp) temp = b(i)
 END DO
 IF (temp <= 0.0D0) GO TO 425
 DO  i = 1,9
   b(i) = b(i)/temp
 END DO
 
!     FIND THE SQUARES OF THE MAGNITUDES OF THE VECTORS OF THE ROWS OF
!     THE B MATRIX.
 
 iorder = 0
 j = 0
 DO  i = 1,9,3
   j = j + 1
   fl(j) = b(i)**2 + b(i+1)**2 + b(i+2)**2
   IF (fl(j) == 0.0D0) iorder = iorder + 1
 END DO
 IF (iorder == 2)  GO TO 410
 IF (iorder == 0)  GO TO 250
 
!     AT THIS POINT ONE AND ONLY ONE FL(I) IS ZERO.
 
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
 CALL mesage (-30,26,NAME)
 140 fm = b(5)*b(9) - b(6)*b(8)
 fr = DSQRT((b(5)**2 + b(6)**2)*(b(8)**2 + b(9)**2))
 GO TO 170
 150 fm = b(1)*b(9) - b(3)*b(7)
 fr = DSQRT((b(1)**2 + b(3)**2)*(b(7)**2 + b(9)**2))
 GO TO 170
 160 fm = b(1)*b(5) - b(2)*b(4)
 fr = DSQRT((b(1)**2 + b(2)**2)*(b(4)**2 + b(5)**2))
 170 IF (fm == 0.0D0) GO TO 175
 IF (fr <= 0.0D0 .OR. fm/fr >= dtol) GO TO 240
 
!     HERE WE HAVE THAT THE ORDER OF THE SINGULARITY IS 2.
 
 175 iorder = 2
 nwds   = 0
 tnwds  = 2
 SELECT CASE ( isave )
   CASE (    1)
     GO TO 180
   CASE (    2)
     GO TO 190
   CASE (    3)
     GO TO 200
 END SELECT
 180 k1   = 5
 k2   = 9
 inc1 = 1
 inc2 = 3
 inc3 = 2
 GO TO 210
 190 k1   = 1
 k2   = 9
 inc1 = 2
 inc2 = 3
 inc3 = 1
 GO TO 210
 200 k1   = 1
 k2   = 5
 inc1 = 3
 inc2 = 2
 inc3 = 1
 210 IF (b(k1) <= 0.0D0 .AND. b(k2) <= 0.0D0) GO TO 425
 IF (b(k1) <= 0.0D0) GO TO 220
 nwds      = 2
 tnwds     = 4
 iarray(3) = ip + inc1
 iarray(4) = ip + inc2
 ipoint    = 5
 GO TO 230
 220 ipoint = 3
 230 IF (b(k2) <= 0.0D0) GO TO 430
 nwds  = nwds  + 2
 tnwds = tnwds + 2
 iarray(ipoint  ) = ip + inc1
 iarray(ipoint+1) = ip + inc3
 GO TO 430
 
!     AT THIS POINT WE HAVE THAT ONE AND ONLY ONE FL IS ZERO BUT THAT
!     ORDER OF THE SINGULARITY IS 1.
 
 240 iorder    = 1
 nwds      = 1
 tnwds     = 3
 iarray(3) = ip + isave
 GO TO 430
 
!     AT STATEMENT NO. 250, WE HAVE THAT ALL THE FL(I) ARE .GT. 0.0D0,
!     SO THAT THE DETERMINANT, DET, OF B MUST BE COMPUTED.
 
 250 det = b(1)*(b(5)*b(9) - b(6)*b(8)) - b(2)*(b(4)*b(9) - b(6)*b(7))  &
     + b(3)*(b(4)*b(8) - b(5)*b(7))
 const = 0.05D0*dtol*fl(1)*fl(2)*fl(3)
 IF (det > const) GO TO 440
 
!     COMPUTE M(I) AND R(I)
 
 m(1) = b(5)*b(9) - b(6)*b(8)
 m(2) = b(1)*b(9) - b(3)*b(7)
 m(3) = b(1)*b(5) - b(2)*b(4)
 r(1) = DSQRT(b(5)**2 + b(6)**2) * DSQRT(b(8)**2 + b(9)**2)
 r(2) = DSQRT(b(1)**2 + b(3)**2) * DSQRT(b(7)**2 + b(9)**2)
 r(3) = DSQRT(b(1)**2 + b(2)**2) * DSQRT(b(4)**2 + b(5)**2)
 
!     FIND I1, J1, K1
!     SUCH THAT M(I1)/R(I1) .GE. M(J1)/R(J1) .GE. M(K1)/R(K1)
 
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
 
!     HERE THE SINGULARITY IS OF ORDER 2.
 
 nwds   = 0
 tnwds  = 2
 iorder = 2
 
!     FIND II, JJ, KK SUCH THAT B(II) .GE. B(JJ) .GE. B(KK)
 
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
 330 IF (b(ll) <= 0.0D0) GO TO 430
 nwds  = nwds  + 2
 tnwds = tnwds + 2
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
 370 iarray(ipoint  ) = ip + inc1
 iarray(ipoint+1) = ip + inc2
 ipoint = ipoint + 2
 kount  = kount  + 1
 IF (kount - 2 < 0) THEN
   GO TO   380
 ELSE IF (kount - 2 == 0) THEN
   GO TO   390
 ELSE
   GO TO   430
 END IF
 380 ll = jj
 GO TO 330
 390 ll = kk
 GO TO 330
 
!     AT THIS POINT THE SINGULARITY IS OF ORDER 1.
 
 400 iorder = 1
 nwds   = 1
 tnwds  = 3
 iarray(3) = ip + i1
 IF (m(j1) < r(j1)*dtol) GO TO 430
 nwds  = 2
 tnwds = 4
 iarray(4) = ip + j1
 IF (m(k1) < r(k1)*dtol) GO TO 430
 nwds  = 3
 tnwds = 5
 iarray(5) = ip + k1
 GO TO 430
 
!     AT THIS POINT 2 ROWS OF THE B MATRIX ARE IDENTICALLY ZERO.
 
 410 nwds   = 2
 tnwds  = 4
 ipoint = 2
 DO  i = 1,3
   IF (fl(i) /= 0.0D0) CYCLE
   ipoint = ipoint + 1
   iarray(ipoint) = ip + i
 END DO
 GO TO 430
 
!     THE SINGULARITY IS OF ORDER 3
 
 425 iorder = 3
 nwds   = 3
 tnwds  = 5
 iarray(3) = ip + 1
 iarray(4) = ip + 2
 iarray(5) = ip + 3
 
!     WRITE IARRAY ON THE GPST FILE.
 
 430 CALL WRITE (ifgpst,iarray(1),tnwds,neor)
 440 GO TO igoto, (450,460)
 450 ASSIGN 460 TO igoto
 ip = ip + 3
 GO TO back, (50,425)
 460 CONTINUE
 RETURN
END SUBROUTINE detck
