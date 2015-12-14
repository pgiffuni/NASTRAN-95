SUBROUTINE smcout ( zi, zr, zd, zrs, zrd )
     
! SMCOUT DOES THE FINAL DIVISION OF THE TERMS OF THE PIVOTAL COLUMN
! AND WRITES THE COLUMN DATA TO THE LOWER TRIANGULAR MATRIX.
 
! THE FOLLOWING CALCULATIONS ARE DONE IN SUBROUTINE SMC2-RS,RD,CS,CD
 
!      do 100 k = 1,n
!         do 10  i = k,n
!         temp = 0.
!         do 5  l = 1,k-1
!            temp = temp + a(i,l)*a(k,l) / a(l,l)
!    5       continue
!         a(i,k) = a(i,k) - temp
!   10    continue
 
!            THE FOLLOWING LAST COMPUTATION TAKES PLACE IN THIS SUBROUTINE.
!            THE RESULTS OF THE DIVISION ARE WRITTEN TO THE OUTPUT FILE BUT
!            THE RESULTS OF THE ABOVE (WITHOUT THE DIVISION BELOW) IS
!            MAINTAINED IN MEMORY FOR REMAINING COLUMN COMPUTATIONS.
 
!         do 11 j = k+1,n
!           a(k,j) = a(j,k) / a( k,k )
!   11      continue
!  100 continue
 
!  THE FINAL COMPUTATIONS ARE WRITTEN TO THE LLL MATRIX USING PUTSTR/ENDPUT.
 
 
 INTEGER, INTENT(IN)                      :: zi(10)
 REAL, INTENT(IN OUT)                     :: zr(10)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(10)
 REAL, INTENT(OUT)                        :: zrs(10)
 DOUBLE PRECISION, INTENT(OUT)            :: zrd(10)
 
 REAL :: minds
 DOUBLE PRECISION :: dakk2, dakkr, dakki, dakk, xnd(10)
 DOUBLE PRECISION :: dr
 INCLUDE  'SMCOMX.COM'
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON / xmssg  / ufm  , uwm  , uim  , sfm
 COMMON / zzzzzz / xns(10)
 EQUIVALENCE      ( xns, xnd )
 EQUIVALENCE      ( minds, mindd ), (dsr, ddr ), (dsc, ddc )
 DATA               rzero / 1.0E-10  /
!      PRINT *,' SMCOUT-ENTER,KCOL=',KCOL
 kdir  = ( kcol-1 ) * 4 + 1
 kmidx = zi( kdir )
 kridx = kmidx + 4
 km2   = zi( kmidx+1 )
 kridn = kridx + km2
 nrows = 0
 DO  i = 1, km2, 2
   nrows  = nrows + zi( kridx + i )
 END DO
 kvidx = kridx + km2
 SELECT CASE ( ktype )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO  2000
   CASE (    3)
     GO TO  3000
   CASE (    4)
     GO TO  4000
 END SELECT
 
! DO DIVISION IN REAL SINGLE PRECISION AND COMPUTE THE  DETERMINANT DDS.
! CHECK FOR THE SMALLEST VALUE OF ANY DIAGONAL ELEMENT ("MINDS")
 
 1000  CONTINUE
 akk          =  zr( kvidx )
!WKBI  7/95 SPR95005
 IF ( akk == 0.0 ) GO TO 7002
 1010  IF ( ABS( dsr ) < 10. ) GO TO 1020
 dsr = dsr / 10.
 power = power + 1
 GO TO 1010
 1020  IF ( ABS( dsr ) > 0.1 ) GO TO 1030
 dsr = dsr * 10.
 power = power - 1
 GO TO 1020
 1030  dsr    = dsr * akk
 minds  = AMIN1 ( ABS(akk), minds )
 IF ( chlsky == 0 ) GO TO 1040
 IF ( akk <= 0 ) GO TO 7001
 akk = SQRT( akk )
 1040  IF ( akk < 0. ) sturm = sturm + 1
!WKBD  7/95 SPR95005
!      IF ( AKK .EQ. 0. ) GO TO 7002
 1050  zrs( 1 ) = akk
!WKBR 7/95 SPR95005
!      AKK      = 1. / AKK
 akk      = -1. / akk
 DO  i = 2, nrows
!WKBR 7/95 SPR95005
!      ZRS( I ) = -1.0 * ZR( KVIDX + I - 1 ) * AKK
   zrs( i ) =  zr( kvidx + i - 1 ) * akk
 END DO
 GO TO 5000
 
! DO DIVISION IN REAL DOUBLE PRECISION AND COMPUTE THE  DETERMINANT DDR.
! CHECK FOR THE SMALLEST VALUE OF ANY DIAGONAL ELEMENT ("MINDD")
 
 2000  CONTINUE
 kvidx = ( kvidx/2 ) + 1
 dakk  =  zd( kvidx )
!WKBI  7/95 SPR95005
 IF ( dakk == 0.d0 ) GO TO 7002
 2010  IF ( DABS( ddr ) < 10.d0 ) GO TO 2020
 ddr = ddr / 10.d0
 power = power + 1
 GO TO 2010
 2020  IF ( DABS( ddr ) > 0.1 ) GO TO 2030
 ddr   = ddr * 10.d0
 power = power - 1
 GO TO 2020
 2030  ddr   = ddr * dakk
 mindd = DMIN1 ( DABS(dakk), mindd )
 IF ( chlsky == 0 ) GO TO 2040
 IF ( dakk <= 0 ) GO TO 7001
 dakk = DSQRT( dakk )
 2040  IF ( dakk < 0.d0 ) sturm = sturm + 1
!WKBD 7/95 SPR95005
!      IF ( DAKK .EQ. 0.D0 ) GO TO 7002
 2050  zrd( 1 ) = dakk
!WKBR 7/95 SPR95005
!      DAKK        = 1.D0 / DAKK
 dakk        = -1.d0 / dakk
 DO  i = 2, nrows
!WKBR 7/95 SPR95005
!      ZRD( I ) = -1.0D0 * ZD( KVIDX + I - 1 ) * DAKK
   zrd( i ) = zd( kvidx + i - 1 ) * dakk
 END DO
 GO TO 5000
 
! DO DIVISION IN COMPLEX SINGLE PRECISION AND COMPUTE THE DETERMINANT
! DSR AND DSC.
! CHECK FOR THE SMALLEST VALUE OF ANY DIAGONAL ELEMENT ("MINDS")
 
 3000  CONTINUE
!   (A+Bi) / (C+Di) = (AC + DB + ( CB-AD )i ) / (C**2 + D**2)
 akkr  = zr( kvidx   )
 akki  = zr( kvidx+1 )
 akk2  = akkr*akkr + akki*akki
!WKBI  7/95 SPR95005
 IF ( akk2 == 0. ) GO TO 7002
 3010  IF ( ABS( dsr**2 + dsc**2 ) < 10. ) GO TO 3020
 dsr   = dsr / 10.
 dsc   = dsc / 10.
 power = power + 1
 GO TO 3010
 3020  IF ( ABS( dsr**2 + dsc**2 ) > 0.1 ) GO TO 3030
 dsr   = dsr * 10.
 dsc   = dsc * 10.
 power = power - 1
 GO TO 3020
 3030  rs    = dsr*akkr - dsc*akki
 dsc   = dsr*akki + dsc*akkr
 dsr   = rs
 minds = AMIN1 ( ABS(akk2), minds )
 IF ( chlsky == 0 ) GO TO 3040
 IF ( akk2 <= 0 ) GO TO 7001
 akk2 = SQRT( akk2 )
 3040  IF ( akkr < 0. ) sturm = sturm + 1
!WKBD  7/95 SPR95005
!      IF ( AKK2 .EQ. 0. ) GO TO 7002
 3050  zrs( 1 ) = akkr
 zrs( 2 ) = akki
 nrowm   = nrows * 2 - 1
!WKBR 7/95 SPR95005
!      AKK2    = 1. / AKK2
 akk2    = -1. / akk2
 kvidx   = kvidx + 1
 DO  i = 2, nrowm, 2
!WKBDB 7/95 SPR95005
!      ZRS( I+1 ) =  -1.0 * ( ZR( KVIDX+I-1   ) * AKKR  +
!     &                         ZR( KVIDX+I   ) * AKKI  ) * AKK2
!      ZRS( I+2 ) =  -1.0 * ( ZR( KVIDX+I     ) * AKKR  -
!     &                         ZR( KVIDX+I-1 ) * AKKI  ) * AKK2
!WKBDE 7/95 SPR95005
!WKBIB 7/95 SPR95005
   zrs( i+1 ) =  ( zr( kvidx+i-1 ) * akkr  + zr( kvidx+i   ) * akki  ) * akk2
   zrs( i+2 ) =  ( zr( kvidx+i   ) * akkr  - zr( kvidx+i-1 ) * akki  ) * akk2
!WKBIE 7/95 SPR95005
 END DO
 GO TO 5000
 
! DO DIVISION IN COMPLEX DOUBLE PRECISION AND COMPUTE THE DETERMINANT
! DDR AND DDC.
! CHECK FOR THE SMALLEST VALUE OF ANY DIAGONAL ELEMENT ("MINDD")
 
 4000  CONTINUE
 kvidx = ( kvidx/2 ) + 1
 dakkr = zd( kvidx   )
 dakki = zd( kvidx+1 )
 dakk2 = dakkr*dakkr + dakki*dakki
!WKBI  7/95 SPR95005
 IF ( dakk2 == 0. ) GO TO 7002
 4010  IF ( DABS( ddr**2 + ddc**2 ) < 10.d0 ) GO TO 4020
 ddr   = ddr / 10.
 ddc   = ddc / 10.
 power = power + 1
 GO TO 4010
 4020  IF ( DABS( ddr**2 + ddc**2 ) > 0.1D0 ) GO TO 4030
 ddr   = ddr * 10.
 ddc   = ddc * 10.
 power = power - 1
 GO TO 4020
 4030  dr    = ddr*dakkr - ddc*dakki
 ddc   = ddr*dakki + ddc*dakkr
 ddr   = dr
 mindd = DMIN1 ( DABS(dakk2), mindd )
 IF ( chlsky == 0 ) GO TO 4040
 IF ( dakk2 <= 0 ) GO TO 7001
 dakk2 = DSQRT( dakk2 )
 4040  IF ( dakkr < 0. ) sturm = sturm + 1
!WKBD  7/95 SPR95005
!      IF ( DAKK2 .EQ. 0. ) GO TO 7002
 4050  zrd( 1 ) = dakkr
 zrd( 2 ) = dakki
 nrowm1  = nrows * 2 - 1
!WKBR 7/95 SPR95005
!      DAKK2   = 1.D0 / (DAKK2 )
 dakk2   = -1.d0 / (dakk2 )
 kvidx   = kvidx + 1
 DO  i = 2, nrowm1, 2
!WKBDB 7/95 SPR95005
!      ZRD( I+1 ) = -1.0D0 * ( ZD( KVIDX+I-1 ) * DAKKR  +
!     &                          ZD( KVIDX+I ) * DAKKI  ) * DAKK2
!      ZRD( I+2 ) = -1.0D0 * ( ZD( KVIDX+I )   * DAKKR  -
!     &                        ZD( KVIDX+I-1 ) * DAKKI  ) * DAKK2
!WKBDE 7/95 SPR95005
!WKBIB 7/95 SPR95005
   zrd( i+1 ) = ( zd( kvidx+i-1 ) * dakkr  +  &
       zd( kvidx+i   ) * dakki  ) * dakk2
   zrd( i+2 ) = ( zd( kvidx+i   ) * dakkr  -  &
       zd( kvidx+i-1 ) * dakki  ) * dakk2
!WKBIE 7/95 SPR95005
 END DO
 GO TO 5000
 
! NOW WRITE THE COLUMN OUT TO THE OUTPUT MATRIX
 
 5000  CONTINUE
 itwrds      = 0
 moblk( 8 )  = -1
 moblk( 12 ) = kcol
 krow  = zi( kridx )
 krows = zi( kridx+ 1)
 IF ( ktype <= 2 ) nwds = 1
 IF ( ktype > 2 ) nwds = 2
 iol = 1
 5050  CALL putstr ( moblk )
 moblk( 4 ) = krow
 moblk( 7 ) = MIN0 ( krows, moblk(6) )
 jstr = moblk( 5 )
 nstr = jstr + (moblk( 7 ) - 1 ) * nwds
 IF ( ktype >= 3 ) nstr = nstr + 1
 IF ( kprec == 2 ) GO TO 5200
 
! MOVE REAL SINGLE AND SINGLE COMPLEX VALUES INTO BUFFER
 
 5100  DO  jj = jstr, nstr
   xns( jj ) = zrs( iol  )
   iol = iol + 1
 END DO
 itwrds = nstr - jstr + 1
 GO TO 5500
 
! MOVE REAL DOUBLE AND DOUBLE COMPLEX VALUES INTO BUFFER
 
 5200  DO  jj = jstr, nstr
   xnd( jj ) = zrd( iol )
!      PRINT *,' SMCOUT,ROW,NUM,TERM=',MOBLK(4),MOBLK(7),XND(JJ)
   iol = iol + 1
 END DO
 itwrds = ( nstr-jstr+1 ) * 2
 5500  CONTINUE
 
! CHECK TO SEE IF ALL CONSECUTIVE ROWS CAN BE STORED IN THE BUFFER
! I.E., ARE THERE ENOUGH WORDS IN THE AVAILABLE STRING
 
 IF ( moblk( 7 ) == krows ) GO TO 5600
 istore = moblk( 7 )
 krows  = krows - istore
 krow   = krow  + istore
 CALL endput ( moblk )
 GO TO 5050
 
! ALL OF THE CURRENT CONSECUTIVE ROWS WERE STORED IN THE BUFFER.
! GO AND GET THE NEXT SET OF CONSECUTIVE ROWS, IF ANY EXIST.
 
 5600  kridx = kridx + 2
 IF ( kridx >= kridn ) GO TO 7000
 CALL endput ( moblk )
 krow  = zi( kridx )
 krows = zi( kridx+1 )
 GO TO 5050
 
! ALL ROWS OF THIS COLUMN HAVE BEEN STORED, CLOSE OUT THE COLUMN
 
 7000  moblk( 8 ) = 1
 CALL endput ( moblk )
 GO TO 7777
 7001  WRITE ( nout, 9001 ) ufm, kcol
 9001  FORMAT(a23,' 3181, ATTEMPT TO PERFORM CHOLESKY DECOMPOSITION'  &
     ,' ON A NEGATIVE DEFINITE MATRIX IN SUBROUTINE SMCOMP.'  &
     ,/,' NEGATIVE DIAGONAL TERM FOUND ON COLUMN ',i6)
 ierror = 4
 CALL mesage ( -61, 0, 0 )
 7002  WRITE ( nout, 9002 ) uwm, kcol, rzero
 9002  FORMAT(a25,' 2396, SMCOMP COMPUTED A ZERO ON THE DIAGONAL '  &
     ,'DURING DECOMPOSITION OF ROW NUMBER ',i6,'.',/  &
     ,' USE OF DIAG 22 OUTPUT SHOULD PERMIT YOU TO CORRELATE THE'  &
     ,' ROW WITH A MODEL D.O.F.',/,' A VALUE OF ',e13.6  &
     ,' WILL BE USED IN PLACE OF THE ZERO, HOWEVER',/  &
     ,' THE ACCURACY OF THE DECOMPOSITION MAY BE IN DOUBT.')
 akk   = rzero
 dakk  = rzero
 akkr  = rzero
 akki  = rzero
 dakkr = rzero
 dakki = rzero
 akk2  = akkr*akkr   + akki*akki
 dakk2 = dakkr*dakkr + dakki*dakki
!WKBIB 7/95 SPR95005
 SELECT CASE ( ktype )
   CASE (    1)
     GO TO  7010
   CASE (    2)
     GO TO  7020
   CASE (    3)
     GO TO  7030
   CASE (    4)
     GO TO  7040
 END SELECT
 7010  CONTINUE
 zr( kvidx ) = akk
 GO TO 1010
 7020  CONTINUE
 zd( kvidx ) = dakk
 GO TO 2010
 7030  CONTINUE
 zr( kvidx   ) = akkr
 zr( kvidx+1 ) = akki
 GO TO 3010
 7040  CONTINUE
 zd( kvidx   ) = dakkr
 zd( kvidx+1 ) = dakki
 GO TO 4010
!WKBIE 7/95 SPR95005
!WKBD  7/95 SPR95005
!      GO TO ( 1050, 2050, 3050, 4050 ), KTYPE
 7777  lll( 6 ) = MAX0( lll(6), itwrds )
 lll( 7 ) = lll( 7 ) + itwrds
 RETURN
END SUBROUTINE smcout
