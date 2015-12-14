SUBROUTINE solve1(a1,r1,rp,xi,lam2,lam3,lam4,cont)
     
!     ROUTINE TO SOLVE FOR LAMBDAS AS FNCTS. OF XI
 
 
 
 
 REAL, INTENT(IN)                         :: a1
 REAL, INTENT(IN)                         :: r1
 REAL, INTENT(IN)                         :: rp
 REAL, INTENT(IN)                         :: xi
 REAL, INTENT(OUT)                        :: lam2
 REAL, INTENT(OUT)                        :: lam3
 REAL, INTENT(OUT)                        :: lam4
 REAL, INTENT(IN OUT)                     :: cont
 
 
 IF (rp == 0.0) GO TO 20
 
 sum = a1 + xi / rp
 sinsum = SIN(sum)
 bb = r1 - rp * (SIN(a1) - sinsum)
 rt = 0.0E0
 IF( sinsum /= 0.0E0 ) rt = bb / sinsum
 psi1 = COS(sum)
 psi2 = -sinsum / rp
 
!     CHECK FOR SHELL CAP CASE
 IF ( a1 /= 0.0 )  GO TO 40
 lam2  = 0.0E0
 IF( bb /= 0.0E0 ) lam2  = psi1 / bb
 lam3  =  1.0 / rp
 lam4  = -1.0 / rp**2
 GO TO 50
 
!     ALF1 = ALF2
 
 20 sina = SIN(a1)
 cosa = COS(a1)
 bb = r1 + xi * cosa
 rt = 0.0E0
 IF( sina /= 0.0E0 ) rt = bb / sina
 psi1 = cosa
 psi2 = 0.0
 
 40 lam2 = 0.0E0
 IF( bb /= 0.0E0 ) lam2 = psi1 / bb
 lam3 = 0.0E0
 IF( rt /= 0.0E0 ) lam3 = 1.0E0 / rt
 lam4 = 0.0E0
 IF( bb /= 0.0E0 ) lam4 = psi2 / bb
 
 50 CONTINUE
 RETURN
END SUBROUTINE solve1
