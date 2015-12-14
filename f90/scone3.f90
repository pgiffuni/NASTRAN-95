SUBROUTINE scone3( again )
     
 
 LOGICAL, INTENT(IN OUT)                  :: again
 REAL :: iii
 
 INTEGER :: iforce(25), istres(100), elemid, iblock(9,14)
 
 
 
 COMMON /sdr2x7/ dum(5),iii,zoff(2),dum2(92),stress(100),force(25)
 
 COMMON /sdr2x8/ vec(8), sum(8), sig(3), sig1, sig2, sig12, temp,  &
     delta, theta, npoint, zoveri, ipt, BLOCK(9,14), elhar, elemid,  &
     harm, n, sinphi, conphi, nphi, nangle
 
 EQUIVALENCE( istres(1), stress(1) )
 EQUIVALENCE( iforce(1), force (1) )
 EQUIVALENCE( iblock(1,1),BLOCK(1,1) )
 
 IF( again ) GO TO 10
 again = .true.
 nangle = 0
 10 nangle = nangle + 1
!*****
!     OUTPUT FORCES FOR THIS ANGLE
!*****
 iforce(1) = elemid
 force(2) = BLOCK(1,nangle)
 force(3) = BLOCK(5,nangle)
 force(4) = BLOCK(6,nangle)
 force(5) = BLOCK(7,nangle)
 force(6) = BLOCK(8,nangle)
 force(7) = BLOCK(9,nangle)
!*****
! COMPUTE AND OUTPUT STRESSES AND PRINCIPAL STRESSES
!*****
 istres(1) = elemid
 stress(2) = BLOCK(1,nangle)
 DO  i = 1,2
   zoveri=0.0
   IF (iii /= 0.0) zoveri=zoff(i)/iii
   DO  j = 1,3
     sig(j) = BLOCK(j+1,nangle) + BLOCK(j+4,nangle) * zoveri
   END DO
   temp = sig(1) - sig(2)
   sig12 = SQRT( (temp*0.50E0)**2 + sig(3)**2 )
   delta = ( sig(1) + sig(2) ) * 0.50E0
   sig1 = delta + sig12
   sig2 = delta - sig12
   delta = 2.0E0 * sig(3)
   IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15 )GO TO 50
   theta = ATAN2( delta, temp ) * 28.6478898E0
   GO TO 60
   50 theta = 0.0E0
   60 ipt = 8*i-6
   stress(ipt+1) = zoff(i)
   stress(ipt+2) = sig(1)
   stress(ipt+3) = sig(2)
   stress(ipt+4) = sig(3)
   stress(ipt+5) = theta
   stress(ipt+6) = sig1
   stress(ipt+7) = sig2
   stress(ipt+8) = sig12
 END DO
!*****
! SET AGAIN .FALSE. IF SDR2E IS NOT TO CALL THIS ROUTINE AGAIN FOR THIS
! ELEMENT.. E.G. ALL THE ANGLES DESIRED HAVE BEEN PROCESSED...
!*****
 IF( nangle == 14 ) GO TO 100
 IF( iblock(1,nangle+1) == 1 ) GO TO 100
 RETURN
 100 again = .false.
 RETURN
END SUBROUTINE scone3
