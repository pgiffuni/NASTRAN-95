SUBROUTINE dssize ( namfil, ncols, nterms, nstrgs, nwdtrm )
     
!   DSSIZE DETERMINES THE SIZE OF A GIVEN MATRIX FILE
!      NCOLS  = NUMBER OF COLUMNS
!      NTERMS = TOTAL NUMBER OF NON-ZERO TERMS IN MATRIX
!      NSTRGS = TOTAL NUMBER OF STRINGS OF CONSECUTIVE TERMS IN MATRIX
!      NWDTRM = NUMBER OF WORDS PER TERM
 
 INCLUDE 'DSIOF.COM'
 
 INTEGER, INTENT(IN)                      :: namfil
 INTEGER, INTENT(OUT)                     :: ncols
 INTEGER, INTENT(OUT)                     :: nterms
 INTEGER, INTENT(OUT)                     :: nstrgs
 INTEGER, INTENT(OUT)                     :: nwdtrm
 COMMON / zzzzzz / mem( 4 )
 INTEGER :: mcb(7)
 
 CALL geturn( namfil )
 IF ( ifilex == 0 ) GO TO 701
 mcb( 1 ) = namfil
 CALL rdtrl ( mcb )
 ncols  = mcb( 2 )
 nstrgs = fcb( 16, ifilex )
 nterms = fcb( 17, ifilex )
 nwdtrm = 2
 IF ( mcb( 5 ) == 1 ) nwdtrm = 1
 IF ( mcb( 5 ) == 4 ) nwdtrm = 4
 GO TO 777
 701   nterms = 0
 nstrgs = 0
 ncols  = 0
 nwdtrm = 0
 777   CONTINUE
 RETURN
END SUBROUTINE dssize
