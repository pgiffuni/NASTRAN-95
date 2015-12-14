SUBROUTINE ortck (x,mass,ibuf,num,ndim,gm,accum,eps)
     
!     ORTCK WILL GENERATE THE GENERALIZED MASS MATRIX FOR CLOSE ROOTS
!     AND MAKE THE EPSILON TEST TO DETERMINE IF THE VECTORS SHOULD BE
!     ORTHOGONALIZED
 
 
 REAL, INTENT(IN OUT)                     :: x(ndim,1)
 INTEGER, INTENT(IN)                      :: mass
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 INTEGER, INTENT(IN)                      :: num
 INTEGER, INTENT(IN)                      :: ndim
 REAL, INTENT(OUT)                        :: gm(num,1)
 DOUBLE PRECISION, INTENT(OUT)            :: accum(1)
 REAL, INTENT(IN OUT)                     :: eps
 
 DIMENSION  im(7)
!     COMMON   /DESCRP/  LENGTH    ,MAJOR(1)
 COMMON   /zntpkx/  z(4)      ,ii       ,ieol
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp
 
 iden  = 0
 im(1) = mass
 CALL rdtrl (im)
 IF (im(4) == 8) iden = 1
 IF (iden  == 1) GO TO 10
 CALL gopen (mass,ibuf,0)
 10 k = 1
 20 DO  i = 1,num
   DO  j = 1,num
     gm(i,j) = 0.
   END DO
 END DO
 DO  i = 1,ndim
   DO  j = 1,num
     accum(j) = 0.d0
   END DO
   IF (iden == 1) GO TO 80
   CALL intpk (*110,mass,0,rsp,0)
   50 IF (ieol == 1)GO TO 90
   CALL zntpki
   60 DO  j = 1,num
     accum(j) = accum(j) + z(1)*x(ii,j)
   END DO
   GO TO 50
   
!     IDENTITY
   
   80 ieol = 1
   ii   = i
   z(1) = 1.0
   GO TO 60
   90 DO  j = 1,num
     DO  m = 1,num
       gm(j,m) = gm(j,m) + accum(j)*x(i,m)
     END DO
   END DO
 END DO
 IF (iden == 1) GO TO 120
 CALL REWIND (mass)
 CALL skprec (mass,1)
 120 gm(1,1) = SQRT(gm(1,1))
 DO  i = 2,num
   gm(i,i) = SQRT(gm(i,i))
   ii = i - 1
   DO  j = 1,ii
     gm(i,j) = gm(i,j)/(gm(i,i)*gm(j,j))
   END DO
 END DO
 DO  i = 1,num
   DO  j = 1,ndim
     x(j,i) = x(j,i)/gm(i,i)
   END DO
 END DO
 j = 0
 150 DO  kk = 1,k
   IF (ABS(gm(k+1,kk)) < eps) CYCLE
   j = 1
   DO  i = 1,ndim
     x(i,k+1) = x(i,k+1) - gm(k+1,kk)*x(i,kk)
   END DO
 END DO
 k = k + 1
 IF (k >= num) GO TO 180
 IF (j ==   0) GO TO 150
 GO TO 20
 180 IF (iden == 1) GO TO 190
 CALL CLOSE (mass,rew)
 190 RETURN
END SUBROUTINE ortck
