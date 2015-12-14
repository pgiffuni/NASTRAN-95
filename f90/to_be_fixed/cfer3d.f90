SUBROUTINE cfer3d (v1,v1l,v2,v2l,v3,v3l,v4,v4l,v5,v5l,zb,zc)
     
!     CFER3D IS A DOUBLE PRECISION ROUTINE (CALLED BY CFEER3) WHICH
!     PERFORMS THE TRIDIAGONAL REDUCTION FOR THE COMPLEX FEER METHOD
 
 
 DOUBLE PRECISION, INTENT(OUT)            :: v1(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v1l(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v2(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v2l(1)
 DOUBLE PRECISION, INTENT(IN)             :: v3(1)
 DOUBLE PRECISION, INTENT(IN)             :: v3l(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v4(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v4l(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v5(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v5l(1)
 REAL, INTENT(IN OUT)                     :: zb(1)
 REAL, INTENT(IN OUT)                     :: zc(1)
 LOGICAL :: sucess   ,no b     ,skip     ,again     , qpr      ,symmet
 INTEGER :: cdp      ,NAME(2)
 DOUBLE PRECISION :: zero     ,dsave(2)  , ss       ,lambda   ,d(4)     ,a(2)
 DIMENSION  s(8)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /  ufm      ,uwm
 COMMON  /feeraa/  ikmb(7,3),ilam(7)  ,iphi(7)  ,dudxx     ,  &
     iscr(11) ,dumaa(84),mcbvec(7)
 COMMON  /feerxc/  lambda(2),symmet   ,mreduc   ,nord      ,  &
     idiag    ,epsdum(2),northo   ,nord2     ,  &
     nord4    ,nordp1   ,nswp(2)  ,no b      ,  &
     it       ,ten2mt   ,tenmht   ,nstart    ,  &
     qpr      ,regdum(2),nzero    ,xcdum(3)  , numran
 COMMON  /system/  ksystm(65)
 COMMON  /unpakx/  iprc     ,ii       ,nn       ,incr
 COMMON  /packx /  itp1     ,itp2     ,iip      ,nnp       , incrp
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew    ,  &
     rew      ,norew    ,eofnrw   ,rsp       , rdp      ,csp      ,cdp      ,sqr
 EQUIVALENCE       (a(1),d(3))        ,(ksystm(2),nout)    , (d(1),s(1))
 DATA     zero  /  0.d0     /
 DATA     NAME  /  4HCFER   ,4H3D     /
 
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
 
!     V1,V2,V3,V4,V5  = AREAS OF OPEN CORE DESIGNATED BY SUBROUTINE
!                       CFEER3 AND USED INTERNALLY AS WORKING VECTORS,
!                       USUALLY RIGHT-HANDED
!     V1L,..,V5L      = SAME AS V1 THRU V5 BUT USUALLY LEFT-HANDED
!     RESTRICTION ...   LEFT-HANDED VECTOR MUST IMMEDIATELY FOLLOW
!                       CORRESPONDING RIGHT-HANDED VECTOR IN CORE
!     ZB,ZC           = REQUIRED GINO BUFFERS
 
!     DEFINITION OF INTERNAL PARAMETERS
 
!     A        = DIAGONAL ELEMENTS OF REDUCED TRIDIAGONAL MATRIX
!     D        = OFF-DIAG ELEMENTS OF REDUCED TRIDIAGONAL MATRIX
!     AGAIN    = LOGICAL INDICATOR FOR CYCLING THRU LOGIC AGAIN WHEN
!                NULL VECTOR TEST (D-BAR) FAILS
!     SKIP     = LOGICAL INDICATOR FOR AVOIDING REDUNDANT OPERATIONS
!     NORTHO   = TOTAL CURRENT NUMBER OF VECTOR PAIRS ON ORTHOGONAL
!                VECTOR FILE
!     NZERO    = NUMBER OF EIGENVECTOR PAIRS ON EIGENVECTOR FILE
!                (RESTART AND PRIOR NEIGHBORHOODS)
!     LANCOS   = LANCZOS ALGORITHM COUNTER
!     NSTART   = NUMBER OF INITIAL REORTHOGONALIZATION ATTEMPTS
 
 IF (qpr) WRITE (nout,8887)
 8887 FORMAT (1H1,50X,6HCFER3D, //)
 
!     SET PACK AND UNPACK CONSTANTS
 
 iprc = cdp
 incr = 1
 itp1 = iprc
 itp2 = itp1
 incrp= incr
 ii   = 1
 iip  = 1
 
!     NN AND NNP ARE SET LOCALLY
 
 CALL gopen (iscr(7),zb(1),wrtrew)
 CALL CLOSE (iscr(7),norew)
 IF (northo == 0) GO TO 20
 
!     LOAD AND RE-NORMALIZE ALL EXISTING VECTORS ON THE NASTRAN
!     EIGENVECTOR FILE (INCLUDES ANY RESTART VECTORS AND ALL VECTORS
!     OBTAINED IN PRIOR NEIGHBORHOODS). PACK THESE VECTORS ON
!     THE ORTHOGONAL VECTOR SCRATCH FILE.
 
 CALL OPEN (*170,iphi(1),zc(1),0)
 
!     LEFT-HAND VECTOR IS STORED IMMEDIATELY AFTER RIGHT-HAND VECTOR
 
 nnp   = nord2
 nord8 = 2*nord4
 DO  i = 1,northo
   IF (qpr) WRITE (nout,8802) i
   8802       FORMAT (1H ,13(10H----------),/,' ORTHOGONAL VECTOR',i3)
   CALL READ (*190,*5,iphi(1),v1(1),nord8+10,0,n3)
   GO TO 210
   5 IF (idiag == 0) GO TO 13
   DO   j = 1,nord4
     IF (v1(j) /= zero) GO TO 13
   END DO
   WRITE (nout,590) i
   13 CONTINUE
   IF (qpr) WRITE (nout,8803) (v1 (j),j=1,nord2)
   IF (qpr) WRITE (nout,8803) (v1l(j),j=1,nord2)
   8803        FORMAT (1H ,(1H ,4D25.16))
   CALL cfnor2 (v1(1),v1l(1),nord2,0,d(1))
   IF (idiag /= 0 .AND. nord2 <= 70) WRITE (nout,570) i,(v1(j),  &
       v1(j+1),v1l(j),v1l(j+1),j=1,nord2,2)
   CALL gopen (iscr(7),zb(1),wrt)
   CALL pack  (v1(1),iscr(7),mcbvec(1))
   CALL CLOSE (iscr(7),norew)
 END DO
 CALL CLOSE (iphi(1),norew)
 IF (idiag /= 0) WRITE (nout,580) northo,mcbvec
 
!     GENERATE INITIAL PSEUDO-RANDOM VECTORS
 
 20 n3 = 3*nord
 ij = 0
 ss = 1.d0
 nzero  = northo
 nstart = 0
 lancos = 0
 again  = .false.
 d(1) = zero
 d(2) = zero
 25 numran = numran + 1
 DO   i = 1,nord4
   ij = ij + 1
   ss = -ss
   IF (i > nord2) GO TO 28
   IF (i > nord ) GO TO 27
   jj = 2*i - 1
   GO TO 30
   27 jj = 2*(i-nord)
   GO TO 30
   28 IF (i > n3) GO TO 29
   jj = 2*i - 1 - nord2
   GO TO 30
   29 jj = 2*(i-n3) + nord2
   
!     THIS LOADS VALUES INTO V1 AND V1L
   
   30  v1(jj) = ss*(MOD(ij,3)+1)/(3.d0*  &
       (MOD(ij,13)+1)*(1+5.0D0*FLOAT(i)/nord))
 END DO
 IF (qpr) WRITE (nout,8844) (v1(i),i=1,nord4)
 8844 FORMAT (1H0,13(10H----------)/(1H ,4D25.16))
 IF (qpr) WRITE (nout,8845)
 8845 FORMAT (1H ,13(10H----------))
 
!     NORMALIZE RIGHT AND LEFT START VECTORS
 
 CALL cfnor2 (v1(1),v1l(1),nord2,0,d(1))
 
!     REORTHOGONALIZE START VECTORS W.R.T. RESTART AND
!     PRIOR-NEIGHBORHOOD VECTORS
 
 CALL cf2ort (sucess,10,ten2mt,nzero,lancos,  &
     v1(1),v1l(1),v5(1),v5l(1),v3(1),v3l(1),zb(1))
 IF (sucess) GO TO  40
 IF (again ) GO TO 160
 35 nstart = nstart + 1
 IF (nstart <= 2) GO TO 25
 WRITE (nout,600) uwm,lambda
 GO TO 450
 40 IF (again) GO TO 90
 
!     SWEEP START VECTORS CLEAN OF ZERO-ROOT EIGENVECTORS
 
 CALL cfe2ao (.false.,v1 (1),v2 (1),v3 (1),zb(1))
 CALL cfe2ao (.true .,v1l(1),v2l(1),v3l(1),zb(1))
 
!     NORMALIZE THE PURIFIED VECTOR AND OBTAIN D(1)
 
 CALL cfnor2 (v2(1),v2l(1),nord2,0,d(1))
 IF (nzero == 0 .OR. northo > nzero) GO TO 50
 
!     IF RESTART OR BEGINNING OF NEXT NEIGHBORHOOD, PERFORM
!     REORTHOGONALIZATION AND RENORMALIZATION
 
 CALL cf2ort (sucess,10,ten2mt,nzero,lancos,  &
     v2(1),v2l(1),v5(1),v5l(1),v3(1),v3l(1),zb(1))
 IF (.NOT.sucess) GO TO 35
 CALL cfnor2 (v2(1),v2l(1),nord2,0,d(1))
 
!     LOAD FIRST VECTORS TO ORTHOGONAL VECTOR FILE
 
 50 CALL gopen (iscr(7),zb(1),wrt)
 nnp = nord2
 CALL pack (v2(1),iscr(7),mcbvec(1))
 CALL CLOSE (iscr(7),norew)
 northo = northo + 1
 
!     COMMENCE LANCZOS ALGORITHM
 
!     INITIALIZE BY CREATING NULL VECTOR
 
 DO   i = 1,nord2
   v1 (i) = zero
   v1l(i) = zero
 END DO
 skip   =.false.
 
!     ENTER LANCZOS LOOP
 
 70 lancos = lancos + 1
 
!     GENERATE DIAGONAL ELEMENT OF REDUCED TRIDIAGONAL MATRIX
 
 IF (.NOT.skip) CALL cfe2ao (.false.,v2(1),v3(1),v5(1),zb(1))
 skip = .false.
 CALL cfnor2 (v3(1),v2l(1),nord2,1,a(1))
 
!     COMPUTE D-BAR
 
 CALL cfe2ao (.true.,v2l(1),v3l(1),v5(1),zb(1))
 DO  i = 1,nord2,2
   j = i + 1
   v4(i)  = v3(i)  - a(1)*v2(i)  + a(2)*v2(j) - d(1)*v1(i)  + d(2)*v1(j)
   v4(j)  = v3(j)  - a(1)*v2(j)  - a(2)*v2(i) - d(1)*v1(j)  - d(2)*v1(i)
   v4l(i) = v3l(i) - a(1)*v2l(i) + a(2)*v2l(j) - d(1)*v1l(i) + d(2)*v1l(j)
   v4l(j) = v3l(j) - a(1)*v2l(j) - a(2)*v2l(i) - d(1)*v1l(j) - d(2)*v1l(i)
 END DO
 CALL cfnor2 (v4(1),v4l(1),nord2,2,d(1))
 dsave(1) = d(1)
 dsave(2) = d(2)
 
!     TEST IF LANCZOS ALGORITHM FINISHED
 
 IF (lancos == mreduc) GO TO 150
 IF (.NOT.qpr) GO TO 85
 WRITE (nout,8845)
 WRITE (nout,8886) d
 8886 FORMAT (8H d-bar =,2D16.8,9X,3HA =,2D16.8)
 WRITE (nout,8844) (v4 (i),i=1,nord2)
 WRITE (nout,8844) (v4l(i),i=1,nord2)
 WRITE (nout,8845)
 85 CONTINUE
 
!     NULL VECTOR TEST
 
 IF (DSQRT(d(1)**2+d(2)**2) > DSQRT(a(1)**2+a(2)**2)*DBLE(tenmht)) GO TO 100
 IF (idiag /= 0) WRITE (nout,610) d
 again = .true.
 GO TO 25
 90 CALL cfe2ao (.false.,v1 (1),v4 (1),v3 (1),zb(1))
 CALL cfe2ao (.true .,v1l(1),v4l(1),v3l(1),zb(1))
 
!     PERFORM REORTHOGONALIZATION
 
 100 CALL cfnor2 (v4(1),v4l(1),nord2,0,d(1))
 CALL cf2ort (sucess,10,ten2mt,nzero,lancos,  &
     v4(1),v4l(1),v3(1),v3l(1),v5(1),v5l(1),zb(1))
 IF (.NOT.sucess) GO TO 160
 
!     NORMALIZE THE REORTHOGONALIZED VECTORS
 
 CALL cfnor2 (v4(1),v4l(1),nord2,0,d(1))
 
!     GENERATE OFF-DIAGONAL ELEMENT OF REDUCED TRIDIAGONAL MATRIX
 
 CALL cfe2ao (.false.,v4(1),v3(1),v5(1),zb(1))
 skip = .true.
 CALL cfnor2 (v3(1),v2l(1),nord2,1,d(1))
 IF (again) GO TO 105
 
!     NULL VECTOR TEST
 
 IF (DSQRT(d(1)**2+d(2)**2) <= DSQRT(a(1)**2+a(2)**2)*DBLE(tenmht)) GO TO 160
 GO TO 110
 105 again = .false.
 d(1) = zero
 d(2) = zero
 
!     TRANSFER TWO ELEMENTS TO REDUCED TRIDIAGONAL MATRIX FILE
 
 110 CALL WRITE (iscr(5),s(1),8,1)
 IF (idiag /= 0) WRITE (nout,560) lancos,d
 
!     LOAD CURRENT VECTORS TO ORTHOGONAL VECTOR FILE
 
 CALL gopen (iscr(7),zb(1),wrt)
 nnp = nord2
 CALL pack (v4(1),iscr(7),mcbvec(1))
 CALL CLOSE (iscr(7),norew)
 northo = northo + 1
 
!     TRANSFER (I+1)-VECTORS TO (I)-VECTORS AND CONTINUE LANCZOS LOOP
 
 DO   i = 1,nord2
   v1 (i) = v2 (i)
   v1l(i) = v2l(i)
   v2 (i) = v4 (i)
   v2l(i) = v4l(i)
 END DO
 GO TO 70
 
!     TRANSFER TWO ELEMENTS TO REDUCED TRIDIAGONAL MATRIX FILE
 
 150 IF (d(1) /= zero .OR. d(2) /= zero) GO TO 155
 d(1) = dsave(1)
 d(2) = dsave(2)
 155 CALL WRITE (iscr(5),s(1),8,1)
 IF (idiag /= 0) WRITE (nout,560) lancos,d
 GO TO 450
 160 mreduc = lancos
 WRITE (nout,500) uwm,mreduc,lambda
 IF (.NOT.again) GO TO 150
 d(1) = zero
 d(2) = zero
 GO TO 150
 170 i = -1
 180 CALL mesage (i,iphi(1),NAME)
 190 i = -2
 GO TO 180
 210 i = -8
 GO TO 180
 450 RETURN
 
 500 FORMAT (a25,' 3157',//5X,'FEER PROCESS MAY HAVE CALCULATED FEWER',  &
     ' ACCURATE MODES',i5,' THAN REQUESTED IN THE NEIGHBORHOOD', ' OF',2D14.6//)
 560 FORMAT (36H reduced tridiagonal matrix elements,5X,3HROW,i4, /10X,  &
     14HOFF-diagonal =,2D24.16, /14X,10HDIAGONAL =,2D24.16)
 570 FORMAT (18H0ORTHOGONAL vector,i4, /1H0,23X,5HRIGHT,56X,4HLEFT, //,  &
     2(1H ,2D25.16,10X,2D25.16))
 580 FORMAT (1H0,i10,32H orthogonal vector pairs on FILE,5X,12X,6I8,/)
 590 FORMAT (18H orthogonal vector,i4,8H is null)
 600 FORMAT (a25,' 3158',//5X,'NO ADDITIONAL MODES CAN BE FOUND BY ',  &
     'FEER IN THE NEIGHBORHOOD OF ',2D14.6,//)
 610 FORMAT (14H d-bar is null,10X,4D20.12)
END SUBROUTINE cfer3d
