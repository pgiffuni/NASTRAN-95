SUBROUTINE q4nrms (bgpdt,gpnorm,iorder,iflag)
!     &    ENTRY Q4NRMD (BGPDT,GPNORM,IORDER,IFLAG)
 
!*****
!     COMPUTES UNIT NORMAL VECTORS FOR QUAD4 GRID POINTS.
!*****
 
 REAL, INTENT(IN)                         :: bgpdt(4,4)
 REAL, INTENT(OUT)                        :: gpnorm(4,4)
 INTEGER, INTENT(IN)                      :: iorder(4)
 INTEGER, INTENT(OUT)                     :: iflag
 
 REAL :: shp(4), sshp(4,2),  v(3,3), tshp(4), tsshp(4,2),  &
     axi(4),aeta(4),eta,vmag,xi
 DOUBLE PRECISION :: dshp(4),dsshp(4,2),dv(3,3),tdshp(4),tdsshp(4,2),  &
     adi(4),aetd(4),etd,dmag,di
 DATA   axi     / -1.0  ,  1.0  , 1.0  , -1.0   /
 DATA   aeta    / -1.0  , -1.0  , 1.0  ,  1.0   /
 DATA   adi     / -1.0D0,  1.0D0, 1.0D0, -1.0D0 /
 DATA   aetd    / -1.0D0, -1.0D0, 1.0D0,  1.0D0 /
 
!     SINGLE PRECISION SECTION -
!*****
!     COMPUTE SHAPE FUNCTION DERIVATIVES
!*****
 iflag = 0
 DO  ii=1,4
   io  = iorder(ii)
   xi  = axi(io)
   eta = aeta(io)
   CALL q4shps (xi,eta,shp,sshp)
!*****
!     SORT THE SHAPE FUNCTIONS
!*****
   DO  i=1,4
     tshp(i) = shp(i)
     DO  j=1,2
       tsshp(i,j) = sshp(i,j)
     END DO
   END DO
   
   DO  ik=1,4
     i = iorder(ik)
     shp(ik) = tshp(i)
     DO  j=1,2
       sshp(ik,j) = tsshp(i,j)
     END DO
   END DO
!*****
!     COMPUTE VECTOR
!*****
   DO  i=1,2
     DO  j=1,3
       v(i,j) = 0.0
       j1 = j + 1
       DO  k=1,4
         v(i,j) = v(i,j) + sshp(k,i)*bgpdt(j1,k)
       END DO
     END DO
   END DO
   
   v(3,1) = v(1,2)*v(2,3) - v(2,2)*v(1,3)
   v(3,2) = v(1,3)*v(2,1) - v(2,3)*v(1,1)
   v(3,3) = v(1,1)*v(2,2) - v(2,1)*v(1,2)
   vmag   = v(3,1)**2+v(3,2)**2+v(3,3)**2
   
   IF (vmag > 1.0E-11) GO TO 40
   iflag = 1
   GO TO 200
   
   40 vmag =  SQRT(vmag)
   gpnorm(2,ii) = v(3,1)/vmag
   gpnorm(3,ii) = v(3,2)/vmag
   gpnorm(4,ii) = v(3,3)/vmag
 END DO
 GO TO 200
 
 ENTRY q4nrmd (bgpdt,gpnorm,iorder,iflag)
!     =======================================
 
!     DOUBLE PRECISION SECTION -
 
!*****
!     COMPUTE SHAPE FUNCTION DERIVATIVES
!*****
 iflag = 0
 DO  ii=1,4
   io = iorder(ii)
   di = adi(io)
   etd = aetd(io)
   CALL q4shpd (di,etd,dshp,dsshp)
   
!     SORT THE SHAPE FUNCTIONS
   
   DO  i=1,4
     tdshp(i) = dshp(i)
     DO  j=1,2
       tdsshp(i,j) = dsshp(i,j)
     END DO
   END DO
   
   DO  ik=1,4
     i = iorder(ik)
     dshp(ik) = tdshp(i)
     DO  j=1,2
       dsshp(ik,j) = tdsshp(i,j)
     END DO
   END DO
!*****
!     COMPUTE VECTOR
!*****
   DO  i=1,2
     DO  j=1,3
       dv(i,j) = 0.0D0
       j1 = j + 1
       DO  k=1,4
         dv(i,j) = dv(i,j) + dsshp(k,i)*bgpdt(j1,k)
       END DO
     END DO
   END DO
   
   dv(3,1) = dv(1,2)*dv(2,3) - dv(2,2)*dv(1,3)
   dv(3,2) = dv(1,3)*dv(2,1) - dv(2,3)*dv(1,1)
   dv(3,3) = dv(1,1)*dv(2,2) - dv(2,1)*dv(1,2)
   dmag    = dv(3,1)**2+dv(3,2)**2+dv(3,3)**2
   
   IF (dmag > 1.0D-11) GO TO 140
   iflag = 1
   GO TO 200
   
   140 dmag = DSQRT(dmag)
   gpnorm(2,ii) = dv(3,1)/dmag
   gpnorm(3,ii) = dv(3,2)/dmag
   gpnorm(4,ii) = dv(3,3)/dmag
 END DO
 
 200 RETURN
END SUBROUTINE q4nrms
