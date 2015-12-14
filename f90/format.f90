SUBROUTINE FORMAT (a,n1x,n2x,n3x,l1x,l2x)
     
!  $MIXED_FORMATS
 
 
 REAL, INTENT(IN OUT)                     :: a(6)
 INTEGER, INTENT(IN)                      :: n1x
 INTEGER, INTENT(IN)                      :: n2x
 INTEGER, INTENT(IN)                      :: n3x
 INTEGER, INTENT(IN)                      :: l1x
 INTEGER, INTENT(IN)                      :: l2x
 REAL :: ff(6),f(6,6,2),f11(6),f21(6),f31(6),f41(6),f51(6),  &
     f61(6),f12(6),f22(6),f32(6),f42(6),f52(6),f62(6)
 
 COMMON /system/ skip,mo
 EQUIVALENCE     (f11(1),f(1,1,1)), (f21(1),f(1,2,1)),  &
     (f31(1),f(1,3,1)), (f41(1),f(1,4,1)),  &
     (f51(1),f(1,5,1)), (f61(1),f(1,6,1)),  &
     (f12(1),f(1,1,2)), (f22(1),f(1,2,2)),  &
     (f32(1),f(1,3,2)), (f42(1),f(1,4,2)), (f52(1),f(1,5,2)), (f62(1),f(1,6,2))
 DATA    f11   / 4H(i5, ,4H49X, ,4H1P,1 ,4HE19. ,4H6,i0 ,4H58)  /,  &
     f21   / 4H(i5, ,4H40X, ,4H1P,2 ,4HE19. ,4H6,i0 ,4H48)  /,  &
     f31   / 4H(i5, ,4H30X, ,4H1P,0 ,4HE19. ,4H6,i0 ,4H39)  /,  &
     f41   / 4H(i5, ,4H21X, ,4H1P,4 ,4HE19. ,4H6,i0 ,4H29)  /,  &
     f51   / 4H(i5, ,4H11X, ,4H1P,5 ,4HE19. ,4H6,i0 ,4H20)  /,  &
     f61   / 4H(i5, ,4H02X, ,4H1P,6 ,4HE19. ,4H6,i0 ,4H10)  /
 DATA    f12   / 4H(i5, ,4H02X, ,4H1P,1 ,4HE19. ,4H6,i1 ,4H05)  /,  &
     f22   / 4H(i5, ,4H02X, ,4H1P,2 ,4HE19. ,4H6,i0 ,4H86)  /,  &
     f32   / 4H(i5, ,4H02X, ,4H1P,3 ,4HE19. ,4H6,i0 ,4H67)  /,  &
     f42   / 4H(i5, ,4H02X, ,4H1P,4 ,4HE19. ,4H6,i0 ,4H48)  /,  &
     f52   / 4H(i5, ,4H02X, ,4H1P,5 ,4HE19. ,4H6,i0 ,4H29)  /,  &
     f62   / 4H(i5, ,4H02X, ,4H1P,6 ,4HE19. ,4H6,i0 ,4H10)  /
 
 n1 = n1x
 n2 = n2x
 n3 = n3x
 l1 = l1x
 l2 = l2x
 n  = (n2-n1+n3)/n3
 IF (n <= 0) GO TO 20
 IF (n > 6) n = 6
 l  = 2
 IF (l1 <= 0 .OR. l2 <= 0) l = 1
 DO  i = 1,6
   ff(i) = f(i,n,l)
 END DO
 l1 = IABS(l1)
 l2 = IABS(l2)
 WRITE (mo,ff,ERR=20) l1,(a(i),i=n1,n2,n3),l2
 
 20  RETURN
END SUBROUTINE FORMAT
