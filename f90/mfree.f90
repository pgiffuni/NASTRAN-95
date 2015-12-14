SUBROUTINE mfree
!     THIS ROUTINE GENERATES MASS TERMS FOR THE INTERNALLY CREATED
!     ELEMENT WHICH DESCRIBES FREE SURFACE EFFECTS
!*****
!     THE ECPT DATA IS
!         NO.       DESCRIPTION
!         1         EL ID
!         2         SIL 1
!         3         SIL 2
!         4         GAMMA
!         5         N
!         6         0
!         7         R1
!         8         Z1
!         9         -
!         10        0
!         11        R2
!         12        Z2
!         13        -
 DOUBLE PRECISION :: rp             ,rn ,dr             ,ct  &
     ,em
 
 INTEGER :: necpt(100)
 
 COMMON /sma2dp/    rp             ,rn ,dr             ,ct  &
     ,em
 COMMON /sma2io/    io(36)
 
 COMMON /sma2cl/    dum(2)        ,npvt
 
 COMMON /sma2et/    ecpt(100)
 EQUIVALENCE   (necpt(1),ecpt(1))
 ifile = io(11)
 IF (ecpt(4) == 0.0) GO TO 1100
 IF(necpt(2) == necpt(3)) GO TO 500
 dr = ecpt(11) - ecpt(7)
 IF(npvt == necpt(2)) GO TO 20
 IF(npvt /= necpt(3)) GO TO 1000
 
 rp = ecpt(11)
 rn = ecpt(7)
 ip = necpt(3)
 in = necpt(2)
 
 GO TO 50
 
 20 rp = ecpt(7)
 rn = ecpt(11)
 ip =necpt(2)
 in =necpt(3)
 50 ct = (0.2617994D0/ecpt(4)) * dr
 IF( necpt(5) == 0) ct = 2.0D0 * ct
 em = ct * (3.0D0*rp +rn)
 CALL sma2b (em,ip,ip,ifile,0.0D0)
 em = ct * ( rp +rn)
 CALL sma2b (em,in,ip,ifile,0.0D0)
 GO TO 1100
 
!      CASE OF CENTER ELEMENT CONNECTED TO ONE POINT
 
 500 IF(necpt(2) /= npvt) GO TO 1000
 ct = 1.5707963D0 / DBLE( ecpt(4) )
 rp = ecpt (7)
 IF( necpt(5) <= 0 ) GO TO 510
 rn =necpt(5)
 ct  = ct/ (2.0D0*rn +2.0D0)
 510 em = ct*rp**2
 ip = npvt
 CALL sma2b( em,ip,ip,ifile,0.0D0)
 1000 RETURN
 1100 RETURN
END SUBROUTINE mfree
