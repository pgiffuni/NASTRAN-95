SUBROUTINE ifp1f (*,iword,ii)
     
!     FINDS FIRST 4 NON-BLANK CHARACTERS
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(OUT)                     :: iword
 INTEGER, INTENT(OUT)                     :: ii
 DIMENSION       core(1),corey(401)
 COMMON /zzzzzz/ corex(1)
 COMMON /ifp1a / skip1(4),ncpw4,skip2(4),izzzbb,skip3(3),iben
 EQUIVALENCE     (corex(1),corey(1)), (core(1),corey(401))
 
 iword = izzzbb
 l  = 1
 ii = 0
 DO  i = 1,18
   DO  j = 1,ncpw4
     k = khrfn1(izzzbb,1,core(i),j)
     IF (k == iben) CYCLE
     IF (ii == 0) ii = i
     iword = khrfn1(iword,l,k,1)
     l = l + 1
     IF (l > ncpw4) GO TO 20
   END DO
 END DO
 RETURN 1
 20   RETURN
END SUBROUTINE ifp1f
