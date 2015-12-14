SUBROUTINE sqrtm(a,ia,b,ib)
     
!     SCALED ARITHMETIC ROUTINES--SQUARE ROOT
 
 
 DOUBLE PRECISION, INTENT(OUT)            :: a
 INTEGER, INTENT(OUT)                     :: ia
 DOUBLE PRECISION, INTENT(IN)             :: b
 INTEGER, INTENT(IN)                      :: ib
 DIMENSION ipsw(1)
 
 DOUBLE PRECISION :: detsw(1)
 
 a = b
 ia = ib
 IF(MOD(ia,2) == 0) GO TO 10
 ia = ia-1
 a = a*10.0
 10 ia=ia/2
 a = DSQRT(DMAX1(a,0.d0))
 20 RETURN
 
!     DCALE OF DETERMINANT BY FACTORS OF 10
 
 ENTRY detm6(detsw,ipsw)
 IF(detsw(1) == 0.0D0) GO TO 20
 30 IF(DABS(detsw(1)) > 10.0D0) GO TO 50
 40 IF(DABS(detsw(1)) < 0.1D0) GO TO 60
 GO TO 20
 50 detsw(1) = detsw(1)*0.1D0
 ipsw(1) = ipsw(1)+1
 GO TO 30
 60 detsw(1) = detsw(1)*10.0D0
 ipsw(1) = ipsw(1)-1
 GO TO 40
END SUBROUTINE sqrtm
