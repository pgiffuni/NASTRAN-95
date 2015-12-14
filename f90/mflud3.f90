SUBROUTINE mflud3
!*****
!     THIS ROUTINE GENERATES THE PSUEDO   MASS    MATRIX TERMS
!     FOR THE TRIANGULAR FLUID ELEMENT
!*****
!     THE ECPT DATA IS THE FOLLOWING
 
!         FIELD         SYMBOL
!           1             ID
!           2             SIL1
!           3             SIL2
!           4             SIL3
!           5             RHO
!           6             BULK
!           7             N
!           8             CSF
!           9             R1
!           10            Z1
!           11            -
!           12            CSF
!           13            R2
!           14            Z2
!           15            -
!           16            CSF
!           17            R3
!           18            Z3
!           19            -
!           20            -
!****
 DOUBLE PRECISION :: r         ,z        ,are2 ,piab     ,emass
!*****
 INTEGER :: necpt(100)
 COMMON/sma2cl/ dum(2),npvt
 COMMON/sma2io/ dum1(10),ifmgg
 COMMON  /sma2et/ ecpt(100)
 COMMON/sma2dp/          r(3)      ,z(3)     ,are2 ,piab      ,emass    ,jp  &
     ,ir        ,jpvt     ,igrid
 EQUIVALENCE        (ecpt(1),necpt(1))
!*****
!*****
 
 IF(ecpt(6) == 0.0) RETURN
!*****
!     STORE THE POINT LOCATIONS AND FIND THE PIVOT POINT
!*****
 jp =0
 DO  i=1,3
   ir =  9 + 4*(i-1)
   r(i)  = ecpt(ir)
   IF(ecpt(ir) <= 0.0) GO TO 1000
   z(i)  = ecpt(ir+1)
   IF( npvt /= necpt(i+1)) CYCLE
   jp = i
 END DO
 IF( jp == 0) GO TO 1000
 are2=DABS((r(2) -r(1))*(z(3)-z(1)) - (r(3)-r(1))*(z(2)-z(1)) )
 piab = 2.617994D-2 *are2 / DBLE( ecpt(6) )
 IF (necpt(7) == 0) piab = piab*2.0D0
 jpvt = npvt
 DO   i = 1,3
   igrid = necpt(i+1)
   emass = piab*( r(1)+r(2)+r(3) +r(jp) +r(i))
   IF (i == jp) emass = emass*2.0D0
   CALL sma2b ( emass, igrid,jpvt,ifmgg,0.0D0)
 END DO
 1000 RETURN
END SUBROUTINE mflud3
