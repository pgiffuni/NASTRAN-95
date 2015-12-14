COMPLEX FUNCTION sumphi(ixr,iyr,nd1,ndn,capphi,dss,n,m,asym)
     
!     FUNCTION TO COMPUTE SUM OF CAPPHI-DELTA SOURCE STENGTH PRODUCT
 
 
 INTEGER, INTENT(IN)                      :: ixr
 INTEGER, INTENT(IN)                      :: iyr
 INTEGER, INTENT(IN OUT)                  :: nd1(1)
 INTEGER, INTENT(IN)                      :: ndn(1)
 COMPLEX, INTENT(IN)                      :: capphi(1)
 COMPLEX, INTENT(IN)                      :: dss(n,m)
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN OUT)                  :: m
 LOGICAL, INTENT(IN)                      :: asym
 
 
 
 
 sumphi  =  ( 0.0 , 0.0 )
 IF ( ixr == 0 )   RETURN
 DO      i = 1 , ixr
   ixs  =  i - 1
   ip  =  ixr - ixs
   ltot  =  2 * ip + 1
   iphi  =  ( ip * ( ip + 1 ) ) / 2
   iys  =  iyr - ixr + ixs
   DO      l = 1 , ltot
     IF ( asym .AND. iys == 0 )   GO TO  200
     j  =  IABS ( iys ) + 1
     IF ( .NOT. ( i >= (nd1(j)) .AND. i <= ndn(j) ) )   GO TO 200
     s  =  1.0
     IF ( asym .AND. iys < 0 )   s  =  -s
     ijphi  =  iphi + 1 + IABS ( iyr - iys )
     sumphi  =  sumphi + s * capphi(ijphi) * dss(i,j)
     200  iys  =  iys + 1
   END DO
 END DO
 RETURN
END FUNCTION sumphi
