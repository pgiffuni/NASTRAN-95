SUBROUTINE cf1ort (sucess,maxits,ten2mt,nzero,iortho,  &
        vr,vl,v1,v1l,v2,v2l,zb)
!*******
!     CF1ORT IS A SINGLE-PRECISION ROUTINE (CREATED FOR USE BY
!     THE COMPLEX FEER METHOD) WHICH PERFORMS THE
!     REORTHOGONALIZATION ALGORITHM
!*******
!     DEFINITION OF INPUT AND OUTPUT PARAMETERS
!*******
!     SUCESS   = LOGICAL INDICATOR FOR SUCCESSFUL REORTHOGONALIZATION
!                (OUTPUT)
!     MAXITS   = MAXIMUM NUMBER OF ALLOWED ITERATIONS (INPUT)
!     TEN2MT   = CONVERGENCE CRITERION
!     NZERO    = NUMBER OF ORTHOGONAL VECTOR PAIRS IN PRIOR
!                NEIGHBORHOODS INCLUDING RESTART
!     IORTHO   = NUMBER OF EXISTING ORTHOGONAL VECTOR PAIRS
!                IN CURRENT NEIGHBORHOOD
!     VR       = RIGHT-HANDED VECTOR TO BE REORTHOGONALIZED
!     VL       = LEFT -HANDED VECTOR TO BE REORTHOGONALIZED
!     V1,V1L,  = WORKING SPACE FOR FOUR VECTORS (V1L MUST
!     V2,V2L     FOLLOW V1 IN CORE)
!     ZB       = WORKING SPACE FOR ONE GINO BUFFER
!*******
 
 LOGICAL, INTENT(OUT)                     :: sucess
 INTEGER, INTENT(IN OUT)                  :: maxits
 REAL, INTENT(IN)                         :: ten2mt
 INTEGER, INTENT(IN)                      :: nzero
 INTEGER, INTENT(IN)                      :: iortho
 REAL, INTENT(IN OUT)                     :: vr(1)
 REAL, INTENT(IN OUT)                     :: vl(1)
 REAL, INTENT(IN)                         :: v1(1)
 REAL, INTENT(IN)                         :: v1l(1)
 REAL, INTENT(OUT)                        :: v2(1)
 REAL, INTENT(OUT)                        :: v2l(1)
 INTEGER, INTENT(IN OUT)                  :: zb(1)
 DIMENSION  , a(2)     ,otest(4)
 LOGICAL :: qpr      ,skip
 
 COMMON  /feeraa/  dumaa(42),iscr7
 COMMON  /feerxc/  dumxc(7) ,idiag    ,xcdum(3) ,nord2  &
     ,xcdum2(9),qpr      ,xcdum3(5),numort
 COMMON  /unpakx/  iprc     ,ii       ,nn       ,incr
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew  &
     ,rew      ,norew    ,eofnrw
 COMMON  /system/  ksys     ,nout
 
 mortho = nzero+iortho
 IF (mortho <= 0) GO TO 500
 IF (qpr) WRITE (nout,700)
 numort = numort + 1
 k = 0
 sucess = .false.
 nn = nord2
 critf = 100.*ten2mt**2
 DO   i = 1,nord2
   v2 (i) = vr(i)
   v2l(i) = vl(i)
 END DO
 CALL gopen (iscr7,zb(1),rdrew)
 8 DO   i = 1,4
   otest(i) = 0.
 END DO
 ll = 2
!*******
!     ENTER LOOP
!*******
 DO   i = 1,mortho
   IF (i == nzero+1) ll = 0
   IF (qpr) WRITE (nout,701) i
!     VALUES ARE UNPACKED INTO BOTH V1 AND V1L
   CALL unpack(*10,iscr7,v1(1))
   IF (.NOT.qpr) GO TO 20
   WRITE (nout,702) (v1 (j),j=1,nord2)
   WRITE (nout,702) (v1l(j),j=1,nord2)
   GO TO 20
   10 IF (idiag /= 0) WRITE (nout,710) i
   CYCLE
!*******
!     OBTAIN RIGHT-HAND INNER-PRODUCT TERM
!*******
   20 CALL cfnor1 (vr(1),v1l(1),nord2,1,a(1))
!*******
!     SUBTRACT OFF RIGHT-HAND INNER-PRODUCT TERM
!*******
   DO   j = 1,nord2,2
     l = j+1
     v2(j) = v2(j) - a(1)*v1(j) + a(2)*v1(l)
     v2(l) = v2(l) - a(1)*v1(l) - a(2)*v1(j)
   END DO
!*******
!     COMPUTE MAXIMUM RIGHT-HAND SQUARED-ERROR
!*******
   a(1) = a(1)**2+a(2)**2
   IF (otest(ll+1) < a(1)) otest(ll+1) = a(1)
!*******
!     OBTAIN LEFT-HAND INNER-PRODUCT TERM
!*******
   CALL cfnor1 (vl(1),v1(1),nord2,1,a(1))
!*******
!     SUBTRACT OFF LEFT-HAND INNER-PRODUCT TERM
!*******
   DO   j = 1,nord2,2
     l = j+1
     v2l(j) = v2l(j) - a(1)*v1l(j) + a(2)*v1l(l)
     v2l(l) = v2l(l) - a(1)*v1l(l) - a(2)*v1l(j)
   END DO
!*******
!     COMPUTE MAXIMUM LEFT-HAND SQUARED-ERROR
!*******
   a(1) = a(1)**2+a(2)**2
   IF (otest(ll+2) < a(1)) otest(ll+2) = a(1)
 END DO
 DO   i = 1,nord2
   vr(i) = v2 (i)
   vl(i) = v2l(i)
 END DO
 skip = .false.
 IF (.NOT.qpr) GO TO 91
 WRITE (nout,702) (vr(i),i=1,nord2)
 WRITE (nout,702) (vl(i),i=1,nord2)
!*******
!     TEST FOR CONVERGENCE
!*******
 91 IF (idiag /= 0) WRITE (nout,703) k,critf,otest
 IF (otest(1) <= critf .AND. otest(2) <= critf .AND.  &
     otest(3) <= critf .AND. otest(4) <= critf) GO TO 450
 IF (skip) GO TO 92
 IF (k /= 1.AND.k /= 3.AND.k /= 5) GO TO 92
 IF (idiag /= 0) WRITE (nout,720)
 critf = 100.*critf
 skip = .true.
 GO TO 91
 92 k = k + 1
 IF (k > maxits) GO TO 95
 CALL CLOSE (iscr7,eofnrw)
 CALL gopen (iscr7,zb(1),rdrew)
 GO TO 8
 95 CALL CLOSE (iscr7,norew)
 GO TO 600
 450 CALL CLOSE (iscr7,norew)
 500 sucess = .true.
 600 RETURN
 700 FORMAT(1H0,//26H begin reorthogonalization,//)
 701 FORMAT(1H ,13HUNPACK vector,i4)
 702 FORMAT(3H --,32(4H----),/(1H ,4E25.16))
 703 FORMAT(32H   reorthogonalization iteration,i3,  &
     9X,14HTARGET value =,e12.4,4X,8HERRORS =,4E12.4)
 710 FORMAT(18H orthogonal vector,i4,  &
     39H is null in reorthogonalization routine)
 720 FORMAT(52H   reorthogonalization tolerance temporarily relaxed)
END SUBROUTINE cf1ort
