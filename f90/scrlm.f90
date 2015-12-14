SUBROUTINE scrlm (scurl, xxi, e, h, cont, rp, alf1, r1, lam1,hf)
     
! THIS SUBROUTINE COMPUTES THE STRESS MATRIX IN FIELD COORDINATES
! FOR THE TOROIDAL RING ELEMENT
 
 
! NOTE THE DOUBLE SUBSCRIPTING USED IN THE SCRLM SUBROUTINE IS
! COMPATIBLE WITH THE CALLING PROGRAM. THE SEL ARRAY WILL RETURN WITH
! THE STRESS MATRIX TRANSPOSED (10X15, STORED ROWWISE) BUT IN THE SCRLM
! SUBROUTINE THE STRESS MATRIX IS COMPUTED AS A DOUBLY SUBSCRIPTED
! 15X10 ARRAY (STORED COLUMNWISE).
 
 
 REAL, INTENT(OUT)                        :: scurl(15,10)
 REAL, INTENT(IN)                         :: xxi(3)
 REAL, INTENT(IN)                         :: e(2,2)
 REAL, INTENT(IN)                         :: h
 REAL, INTENT(IN OUT)                     :: cont
 REAL, INTENT(IN OUT)                     :: rp
 REAL, INTENT(IN)                         :: alf1
 REAL, INTENT(IN OUT)                     :: r1
 REAL, INTENT(IN)                         :: lam1
 REAL, INTENT(IN)                         :: hf
 
 REAL :: lam2 , lam3 ,lam4
!     ------------------------------------------------------------------
 
 sc1 = h
 sc2 =hf**3 / 12.0
 jj = 1
 kk = 3
 ll = 5
 
 DO  i = 1,3
   xx1 = xxi(i)
   xx2 = xx1 * xx1
   xx3 = xx2 * xx1
   xx4 = xx3 * xx1
   xx5 = xx4 * xx1
   CALL solve1(alf1,r1,rp,xx1,lam2,lam3,lam4,cont)
   DO  j = 1,2
     scurl(jj, 1) = lam2 * e(j,2)
     scurl(jj, 2) = scurl(jj,1) * xx1 + e(1,j)
     scurl(jj, 3) = scurl(jj,1) * xx2 + e(1,j) * 2.0 * xx1
     scurl(jj, 4) = scurl(jj,1) * xx3 + e(1,j) * 3.0 * xx2
     scurl(jj, 5) = lam1 * e(1,j) + lam3 * e(j,2)
     scurl(jj, 6) = scurl(jj,5) * xx1
     scurl(jj, 7) = scurl(jj,5) * xx2
     scurl(jj, 8) = scurl(jj,5) * xx3
     scurl(jj, 9) = scurl(jj,5) * xx4
     scurl(jj,10) = scurl(jj,5) * xx5
     jj = jj + 1
   END DO
   jj = jj + 3
   DO  k = 1,2
     scurl (kk,1) = 0.0
     scurl (kk,2) = 0.0
     scurl (kk,3) = 0.0
     scurl (kk,4) = 0.0
     scurl(kk, 5) = 0.0
     scurl(kk, 6) = -lam2 * e(k,2)
     scurl(kk, 7) = scurl(kk,6) * 2.0 * xx1 - e(1,k) *  2.0
     scurl(kk, 8) = scurl(kk,6) * 3.0 * xx2 - e(1,k) *  6.0 * xx1
     scurl(kk, 9) = scurl(kk,6) * 4.0 * xx3 - e(1,k) * 12.0 * xx2
     scurl(kk,10) = scurl(kk,6) * 5.0 * xx4 - e(1,k) * 20.0 * xx3
     kk = kk + 1
   END DO
   kk = kk + 3
   el = e(1,1) * lam2
   ell = el * lam1
   eel = e(1,1) * lam1
   scurl (ll,1) = 0.0
   scurl (ll,2) = 0.0
   scurl (ll,3) = 0.0
   scurl (ll,4) = 0.0
   scurl (ll,5) = 0.0
   scurl(ll, 6) = lam2**2 * e(2,2) - lam4 * e(1,2)
   scurl(ll, 7) = scurl(ll,6) * 2.0 * xx1 - 2.0 * el
   scurl(ll, 8) = scurl(ll,6) * 3.0 * xx2 - 6.0 * (el * xx1 + e(1,1))
   scurl(ll, 9) = scurl(ll,6) * 4.0 * xx3 - 12.0 * el * xx2 - 24.0 *  &
       e(1,1) * xx1
   scurl(ll,10) = scurl(ll,6) * 5.0 * xx4 - 20.0 * el * xx3 - 60.0 *  &
       e(1,1) * xx2
   ll = ll + 5
 END DO
 
!     ADJUSTMENT FOR SHELL CAP CASE
 IF ( alf1 /= 0.0 )  GO TO 198
 scurl (1,2) =  e(1,2) + e(1,1)
 scurl (2,2) =  e(2,2) + e(1,2)
 scurl (3,7) = -2. * (e(1,2) + e(1,1) )
 scurl (4,7) =  2. * (e(2,2) + e(1,2) )
 scurl (5,8) =  3. * (e(2,2) - 4.*e(1,1) )
 198 DO  j = 1,15,5
   DO  i = 1,10
     scurl(j  ,i) = scurl(j  ,i) * sc1
     scurl(j+1,i) = scurl(j+1,i) * sc1
     scurl(j+2,i) = scurl(j+2,i) * sc2
     scurl(j+3,i) = scurl(j+3,i) * sc2
     scurl(j+4,i) = scurl(j+4,i) * sc2
   END DO
 END DO
 RETURN
END SUBROUTINE scrlm
