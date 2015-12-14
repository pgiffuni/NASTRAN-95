SUBROUTINE fcurl (fmeo, fme1, ffeo, ffe1, yi, s, lam1)
!     ------------------------------------------------------------------
 
 REAL, INTENT(OUT)                        :: fmeo(10,2)
 REAL, INTENT(OUT)                        :: fme1(10,2)
 REAL, INTENT(OUT)                        :: ffeo(10,2)
 REAL, INTENT(OUT)                        :: ffe1(10,2)
 REAL, INTENT(IN)                         :: yi(6, 7)
 REAL, INTENT(IN)                         :: s
 REAL, INTENT(IN)                         :: lam1
 
 
 
 
 fmeo( 1,1) = 0.0
 fmeo( 2,1) = yi(1,1)
 fmeo(3,1)  =  yi(1,2)  *  2.0
 fmeo(4,1)  =  yi(1,3)  *  3.0
 fmeo( 5,1) = yi(1,1) * lam1
 fmeo( 6,1) = yi(1,2) * lam1
 fmeo( 7,1) = yi(1,3) * lam1
 fmeo( 8,1) = yi(1,4) * lam1
 fmeo( 9,1) = yi(1,5) * lam1
 fmeo(10,1) = yi(1,6) * lam1
 
 fmeo( 1,2) = yi(4,1)
 fmeo( 2,2) = yi(4,2)
 fmeo( 3,2) = yi(4,3)
 fmeo( 4,2) = yi(4,4)
 fmeo( 5,2) = yi(2,1)
 fmeo( 6,2) = yi(2,2)
 fmeo( 7,2) = yi(2,3)
 fmeo( 8,2) = yi(2,4)
 fmeo( 9,2) = yi(2,5)
 fmeo(10,2) = yi(2,6)
 
 s1 = 1.0 / s
 fme1( 1,1) = 0.0
 fme1( 2,1) = s1 * yi(1,2)
 fme1(3,1)  =  s1  *  yi(1,3) * 2.0
 fme1(4,1)  =  s1  *  yi(1,4) * 3.0
 fme1( 5,1) = s1 * yi(1,2) * lam1
 fme1( 6,1) = s1 * yi(1,3) * lam1
 fme1( 7,1) = s1 * yi(1,4) * lam1
 fme1( 8,1) = s1 * yi(1,5) * lam1
 fme1( 9,1) = s1 * yi(1,6) * lam1
 fme1(10,1) = s1 * yi(1,7) * lam1
 
 fme1( 1,2) = s1 * yi(4,2)
 fme1( 2,2) = s1 * yi(4,3)
 fme1( 3,2) = s1 * yi(4,4)
 fme1( 4,2) = s1 * yi(4,5)
 fme1( 5,2) = s1 * yi(2,2)
 fme1( 6,2) = s1 * yi(2,3)
 fme1( 7,2) = s1 * yi(2,4)
 fme1( 8,2) = s1 * yi(2,5)
 fme1( 9,2) = s1 * yi(2,6)
 fme1(10,2) = s1 * yi(2,7)
 
 ffeo( 1,1) = 0.0
 ffeo (2,1) = 0.0
 ffeo (3,1) = 0.0
 ffeo (4,1) = 0.0
 ffeo( 5,1) = 0.0
 ffeo( 6,1) = 0.0
 ffeo( 7,1) = - 2.0 * yi(1,1)
 ffeo( 8,1) = - 6.0 * yi(1,2)
 ffeo( 9,1) = -12.0 * yi(1,3)
 ffeo(10,1) = -20.0 * yi(1,4)
 
 ffeo (1,2) = 0.0
 ffeo (2,2) = 0.0
 ffeo (3,2) = 0.0
 ffeo (4,2) = 0.0
 ffeo( 5,2) = 0.0
 ffeo( 6,2) = -yi(4,1)
 ffeo( 7,2) = -2.0 * yi(4,2)
 ffeo( 8,2) = -3.0 * yi(4,3)
 ffeo( 9,2) = -4.0 * yi(4,4)
 ffeo(10,2) = -5.0 * yi(4,5)
 
 ffe1( 1,1) = 0.0
 ffe1 (2,1) = 0.0
 ffe1 (3,1) = 0.0
 ffe1 (4,1) = 0.0
 ffe1( 5,1) = 0.0
 ffe1( 6,1) = 0.0
 ffe1( 7,1) = -s1 *  2.0 * yi(1,2)
 ffe1( 8,1) = -s1 *  6.0 * yi(1,3)
 ffe1( 9,1) = -s1 * 12.0 * yi(1,4)
 ffe1(10,1) = -s1 * 20.0 * yi(1,5)
 
 ffe1 (1,2) = 0.0
 ffe1 (2,2) = 0.0
 ffe1 (3,2) = 0.0
 ffe1 (4,2) = 0.0
 ffe1( 5,2) = 0.0
 ffe1( 6,2) = -s1 * yi(4,2)
 ffe1( 7,2) = -s1 * 2.0 * yi(4,3)
 ffe1( 8,2) = -s1 * 3.0 * yi(4,4)
 ffe1( 9,2) = -s1 * 4.0 * yi(4,5)
 ffe1(10,2) = -s1 * 5.0 * yi(4,6)
 RETURN
END SUBROUTINE fcurl
