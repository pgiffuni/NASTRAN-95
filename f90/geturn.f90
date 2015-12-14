SUBROUTINE geturn ( namfil )
     INCLUDE 'DSIOF.COM'
 
 INTEGER, INTENT(IN OUT)                  :: namfil
 COMMON / xfist  / ifstmx, ifstca, ifist( 100 )
 COMMON / xfiat  / ifatuf, ifatmx, ifatca, ifiat( 640 )
 COMMON / xxfiat / ixfiat( 19 )
 INTEGER*2         iunit
 COMMON / dsunit / iunit( 220 )
 DATA  mask / '00007FFF'x /
 
 IF ( namfil == lasnam .AND. ifilex /= 0 ) GO TO 20
 ifilex = 0
 lim = 2 * ifstca - 1
 DO  ifst = 1, lim, 2
   IF ( namfil /= ifist( ifst ) ) CYCLE
   IF ( namfil >= 101 .AND. namfil <= 320 ) GO TO 10
   IF ( ifist( ifst + 1 ) > 0 ) GO TO 5
   ifilex = ixfiat( IABS( ifist( ifst+1 ) ) + 1 )
   IF (ifilex <= maxpri) GO TO 20
   2    ifilex = 0
   GO TO 200
   5    ifilex = IAND( ifiat( ifist( ifst+1 ) - 2 ), mask )
   GO TO 20
   10   ifilex = IAND( ifiat( ifist( ifst+1 ) - 2 ), mask )
   IF (ifilex > maxpri) GO TO 2
   iunit( namfil-100 ) = ifilex
   GO TO 20
 END DO
 GO TO 200
 20   iprvop = fcb( 1, ifilex )
 IF ( iprvop == 2 ) iprvop = 0
 nlr    = fcb( 3, ifilex )
 nblock = fcb( 4, ifilex )
 lasnam = namfil
 200   CONTINUE
 RETURN
END SUBROUTINE geturn
