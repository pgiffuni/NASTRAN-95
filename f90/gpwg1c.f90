SUBROUTINE gpwg1c  (b,e,eig,iflag)
     
!     DOUBLE PRECISION VERSION, BY G.CHAN/SPERRY    8/86
 
!     IFLAG=0 MEANS RUN OK
!     IFLAG=1 MEANS NO SOLUTION IN 20 ITERATIONS
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: b(3,3)
 DOUBLE PRECISION, INTENT(OUT)            :: e(3,3)
 DOUBLE PRECISION, INTENT(OUT)            :: eig(3)
 INTEGER, INTENT(OUT)                     :: iflag
 DOUBLE PRECISION :: ep(3,3), bp(3,3)
 DOUBLE PRECISION :: detb,epsil,bmax,r,s,c,t
 
 detb = 0.0D0
 DO  i = 1,3
   DO  j = 1,3
     detb = detb+b(i,j)*b(i,j)
   END DO
 END DO
 epsil = DSQRT(detb)*1.0D-5
 iflag =0
 ii = 1
 DO  i=1,3
   DO  j=1,3
     e(i,j) = 0.0D0
     IF (i == j) e(i,j) = 1.0D0
   END DO
 END DO
 IF (detb == 0.0D0) GO TO 100
 15 bmax = DMAX1(DABS(b(1,2)),DABS(b(1,3)),DABS(b(2,3)))
 IF (DABS(bmax)  < epsil) GO TO 100
 IF (bmax /= DABS(b(1,2))) GO TO 20
 i = 1
 j = 2
 k = 3
 GO TO 40
 20 IF (bmax /= DABS(b(1,3))) GO TO 30
 i = 1
 j = 3
 k = 2
 GO TO 40
 30 i = 2
 j = 3
 k = 1
 40 r = (b(j,j)-b(i,i))/b(i,j)
 IF (DABS(r) < 1.0D-6) GO TO 50
 IF (DABS(r) > 1.0D+6) GO TO 60
 t = DSQRT((r*r)/4.0D0+1.0D0)-0.5D0*r
 c = DSQRT(1.0D0+t*t)
 s = t/c
 c = 1.0D0/c
 GO TO 70
 50 s = DSQRT(.5D0)
 c = s
 GO TO 70
 60 s = 0.0D0
 c = 1.0D0
 70 bp(i,i) = b(i,i)*c*c+b(j,j)*s*s-2.0D0*b(i,j)*s*c
 bp(j,j) = b(i,i)*s*s+b(j,j)*c*c+2.0D0*b(i,j)*s*c
 bp(k,k) = b(k,k)
 bp(j,i) = 0.0D0
 bp(i,j) = 0.0D0
 bp(k,i) = b(i,k)*c-b(j,k)*s
 bp(i,k) = bp(k,i)
 bp(k,j) = b(j,k)*c+b(i,k)*s
 bp(j,k) = bp(k,j)
 ep(i,1) = e(i,1)*c-e(j,1)*s
 ep(j,1) = e(i,1)*s+e(j,1)*c
 ep(k,1) = e(k,1)
 ep(i,2) = e(i,2)*c-e(j,2)*s
 ep(j,2) = e(i,2)*s+e(j,2)*c
 ep(k,2) = e(k,2)
 ep(i,3) = e(i,3)*c-e(j,3)*s
 ep(j,3) = e(i,3)*s+e(j,3)*c
 ep(k,3) = e(k,3)
 DO   i=1,3
   DO   j=1,3
     b(i,j) = bp(i,j)
     e(i,j) = ep(i,j)
   END DO
 END DO
 IF (ii >= 21) GO TO 90
 ii = ii+1
 GO TO 15
 90 iflag=1
 GO TO 120
 100 DO   i=1,3
   eig(i) = b(i,i)
 END DO
 120 RETURN
END SUBROUTINE gpwg1c
