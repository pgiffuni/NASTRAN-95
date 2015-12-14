SUBROUTINE smcph2 ( zi, zr, zd )
     
! SMCPH2 PERFORMS THE ACTUAL DECOMPOSITION OF THE MATRIX THAT WAS
! SETUP IN MEMORY AND/OR THE SPILL BY SMCPH1.
! SEE SMCPH1 FOR THE DEFINITION OF SOME OF THE VARIABLES IN /SMCOMX/
 
 
 INTEGER, INTENT(IN OUT)                  :: zi(10)
 REAL, INTENT(IN OUT)                     :: zr(10)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(10)
 
 DOUBLE PRECISION :: xnd(10)
 INTEGER :: itemp(4)
 INTEGER :: prc     ,words   ,rlcmpx  ,NAME(2)
 CHARACTER (LEN=4) :: cname(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 INCLUDE           'SMCOMX.COM'
 COMMON  /logout/  lout
 COMMON  /xmssg /  ufm     ,uwm     ,uim     ,sfm
 COMMON  /names /  rdnrw   ,rdrew   ,wrt     ,wrtrew  ,rew  &
     ,                 norew   ,eofnrw  ,rsp     ,rdp     ,csp  &
     ,                 cdp     ,sqr     ,rect    ,diag    ,lowtri  &
     ,                 uprtri  ,sym
 COMMON  /zzzzzz/  xns(10)
 COMMON  /TYPE  /  prc(2)  , words(4), rlcmpx(4)
 EQUIVALENCE       ( xns, xnd )
 EQUIVALENCE       ( mblk(6), mterms ), (mblk(5), mstr )
 EQUIVALENCE       ( mblk(4), mrow   ), (mblk(2), mtype)
 EQUIVALENCE       ( NAME   , cname  )
 
!   open core is allocated as follows for the decomposition
 
!       -------------------------------
!       zi(1)
!       Directory  (4,n) , n=number of columns of matrix
!                  (1,i) = index to active rows and terms within memory
!                  (2,i) = first column data needed for this pivot
!                  (3,i) = last pivot column to use this data
!                          (also, the last active row of this column)
!                  (4,i) = savpos position pointer for data spilled to a
!                          scratch file
!       -------------------------------
!       zi(nar)
!       Area for storage of row numbers used for previous column of
!       decomposition (length=MAXNAR+2)
!       -------------------------------
!       zi(ispill)
!       Area to read data from spill file (length =MXRECL+4)
!       This area is not needed if no columns written to spill file
!       -------------------------------
!       zi(ILSROW)
!       Area for storage of last non-zero row term for a given column
!       (length=MAXNCOL)
!       -------------------------------
!       zi(ioloop)
!       Values for outer loop terms in all row computations in the
!       current pivotal column.
 
!         temp = temp + a(i,j) * a(k,j) / a(j,j)
!                                ===============
!            i = row of column
!            k = pivotal column being processed
!            j = 1, k-1
!          a(i,k) = a(i,k) - temp
!          (Note, length is 2*MAXNCOL)
!          MAXNCOL = maximum number of columns referenced by any
!                    pivotal column
!       -------------------------------
!       zi(iiloop)
!       Values for inner loop terms in each row computation
!         temp = temp + a(i,j) * a(k,j) / a(j,j)
!                       ======
!            i = row of column
!            k = pivotal column being processed
!            j = 1, k-1
!          a(i,k) = a(i,k) - temp
!          (Note, length is  MAXNCOL*MAXNAC)
!          MAXNAC = MAXIMUM NUMBER OF ACTIVE ROWS FOR ANY GIVEN COLUMN
!       -------------------------------
!       zi(iwork)
!       Temporary storage for storing "temp" values for each row (see
!       "temp" in above equation for zi(iiloop) )
!       -------------------------------
!       zi(idbase)
!       Memory for rows and terms of columns as pointed to by directory
!       in the first part of open core.  This data is loaded from the
!       bottom up to allow for better management of open core.
!       The format for the storage of this data is as follows:
!          (index from directory above points to the first word of the
!           data that follows)
!               1.  Column number
!               2.  Length of active row section (m*2), m=number of
!                   repeats of contents of words 5 and 6 below.
!               3.  Total number of words in this block of allocation
!               4.  Length of values section
!               5.  row number
!               6.  number of consecutive values beginning at this row
!                   (words 5 and 6 repeat m times)
!           5+2*m.  value for first row
!     5+2*m+iprec.  next row value (iprec=1,2,4)
!   5+2*m+iprec*l.  last row value for column (l=total values)
!       -------------------------------
!       zi(ibuf2)
!       Buffer for spill file if all column values can not be kept in memory
!       -------------------------------
!       zi(ibuf1)
!       Buffer for input matrix file to be decomposed
!       -------------------------------
 
!      CALL AUDIT ('SMCPH2  ',1 )
 nar    = ncol*4 + 1
 ilsrow = nar + maxnar + 2
 mspill = 0
 ispill = 0
 IF ( nspill == 0 ) GO TO 5
 ispill = nar + maxnar + 1
 IF ( MOD( ispill,2 ) == 0 ) ispill = ispill + 1
 ilsrow = ispill + mxrecl + 4
 5     CONTINUE
 ioloop = ilsrow + maxncol + 2
 IF ( MOD( ioloop,2 ) == 0 ) ioloop = ioloop + 1
 iiloop = ioloop + 2*maxncol*ivwrds
 iwork  = iiloop + maxncol*maxnac*ivwrds
 itotal = iwork  + maxnac*ivwrds
 inddir = lascol * 4 - 3
 idbase = zi( inddir )
 memfre = 0
 memlas = 0
 memlck = 0
 memcol1= 1
 xncol  = ncol
 xspill = nspill
 xfact  = xncol / ( xncol-nspill )
 
!  MORE = ESTIMATED NUMBER OF WORDS NEEDED FOR STORING ALL OF MATRIX
!  IADJ = WORDS OF COLUMN DATA THAT WILL NEED TO BE WRITTEN TO THE SPILL
!         FILE TO ALLOW FOR "ITOTAL" WORDS FOR THE PHASE II ARRAYS.
 
 more   = xfact * ( lcore - idbase )
 iadj   = 0
 IF ( itotal > idbase ) iadj = itotal - idbase
 maxmem = itotal + more + iadj
 CALL sswtch ( 45, l45 )
 IF ( nspill == 0 .AND. l45 == 0 ) GO TO 10
 WRITE ( lout, 8002 ) maxmem, lcore
 8002  FORMAT(  &
     7X,' ESTIMATED OPEN CORE NEEDED TO ELIMINATE USE OF SPILL=',i8  &
     ,/, 7X,' OPEN CORE AVAILABLE FOR THIS DECOMPOSITION          =',i8 )
 
! TEST TO BE SURE THAT AT LEAST HALF OF THE MEMORY IS AVAILABLE.
! IF NOT, USE OLD METHOD INSTEAD OF THIS ONE.
 
 xmaxmem = maxmem
 xcore   = lcore
 percnt  = xcore / xmaxmem
 IF ( percnt < .5 ) GO TO 7008
 
! CHECK TO SEE IF ENOUGH OPEN CORE FOR INNER AND OUTER LOOP VALUES
 
 10    IF ( itotal < idbase ) GO TO 500
 
! NEED MORE OPEN CORE FOR LOOP AREAS.  WRITE COLUMN DATA TO SPILL FILE.
! IF COLUMNS WERE WRITTEN TO SPILL FILE FROM SMCPH1, THEN FILE WILL
! STILL BE OPEN.  IF NOT, MUST ALLOW FOR SPILL AREA IN OPEN CORE AND
! RE-ADJUST THE OPEN CORE POINTERS.
 
 nextra = 0
 IF ( opnscr ) GO TO 20
 opnscr = .true.
 CALL OPEN ( *7003, iscr1, zi( ibuf2 ), wrtrew )
 ispill = nar + maxnar + 1
 IF ( MOD( ispill,2 ) == 0 ) ispill = ispill + 1
 ilsrow = ispill + mxrecl + 4
 ioloop = ilsrow + maxncol + 2
 IF ( MOD( ioloop,2 ) == 0 ) ioloop = ioloop + 1
 iiloop = ioloop + 2*maxncol*ivwrds
 iwork  = iiloop + maxncol*maxnac*ivwrds
 itotal = iwork  + maxnac*ivwrds
 20    CONTINUE
 
! WRITE THE LAST COLUMN OF DATA CURRENTLY IN MEMORY TO THE SPILL FILE
 
 INDEX  = zi( inddir )
 irval  = INDEX + 4
 nrvals = zi( INDEX+1 )
 nterms = zi( INDEX+3 )
 ivval  = irval + nrvals
 itemp( 1 ) = zi( INDEX )
 itemp( 2 ) = nrvals
 itemp( 3 ) = 0
 itemp( 4 ) = nterms
!      PRINT *,' SMCPH2 CALLING WRITE FOR ITEMP,NRVALS,NTERMS,IVWRDS'
!      PRINT *,                           ITEMP,NRVALS,NTERMS,IVWRDS
 CALL WRITE ( iscr1, itemp, 4, 0 )
 CALL savpos( iscr1, kpos )
 CALL WRITE ( iscr1, zr( irval ), nrvals, 0 )
 CALL WRITE ( iscr1, zr( ivval ), nterms*ivwrds, 1 )
 zi( inddir   ) = 0
 zi( inddir+3 ) = kpos
 50    inddir = inddir - 4
 IF ( inddir <= 0 ) GO TO 7008
 IF ( zi ( inddir ) == 0 ) GO TO 50
 
! RESET IDBASE TO INDICATE THE LAST COLUMN OF DATA IN MEMORY
 
 idbase = zi( inddir )
 mspill = mspill + 1
 GO TO 10
 
! OPEN THE OUTPUT FILE
 
 500   CONTINUE
 left   = idbase - itotal
 
! DETERMINE HOW MANY MORE COLUMNS OF THE INNER LOOP AREA AND
! EXTRA TERMS OF THE OUTER LOOP AREA ARE AVAILABLE
!   NEXTRA = NUMBER OF EXTRA COLUMNS AVAILABLE IN THE INNER LOOP AREA
!          = NUMBER OF EXTRA COLUMNS AVAILABLE IN THE OUTER LOOP AREA
!            (INNER LOOP AREA SIZE = MAXNAC * ( MAXNCOL + NEXTRA ) )
!            (OUTER LOOP AREA SIZE = 2      * ( MAXNCOL + NEXTRA ) )
!          = NUMBER OF EXTRA ROWS IN THE "ILSROW" ARRAY (MAXNCOL+NEXTRA)
!  (Note: for each column added, we need the following:
!           for array ILSROW:                1
!           to insure double word boundary:  1
!           for outer loop:                  2*IVWRDS
!           for inner loop:                  MAXNAC*IVWRDS
!           ( must allow for temp array size:    MAXNAC*IVWRDS
 need   = 2 + 2*ivwrds + maxnac*ivwrds
 nextra = ( left - 2 - (maxnac*ivwrds) ) / need
!      PRINT *,' LEFT,NEED,NEXTRA=',LEFT,NEED,NEXTRA
 IF ( nextra == 0 ) GO TO 505
 ioloop = ilsrow +            ( maxncol+nextra ) + 2
 IF ( MOD( ioloop,2 ) == 0 ) ioloop = ioloop + 1
 iiloop = ioloop + ( 2      * ( maxncol+nextra ) ) * ivwrds
 iwork  = iiloop + ( maxnac * ( maxncol+nextra ) ) * ivwrds
 itotal = iwork  +            ( maxnac           ) * ivwrds
 505   IF ( kprec == 2 ) ioloop = ioloop / 2 + 1
 IF ( kprec == 2 ) iiloop = iiloop / 2 + 1
 IF ( kprec == 2 ) iwork  = iwork  / 2 + 1
 nvterm = 1
 IF ( ktype >= 3 ) nvterm = 2
 IF ( mspill /= 0 ) WRITE ( lout, 8001 ) mspill
 8001  FORMAT(8X,'ADDITIONAL COLUMNS WRITTEN TO SPILL '  &
     ,'FOR PHASE II PROCESSING =',i6)
 IF ( .NOT. opnscr ) GO TO 510
 CALL CLOSE ( iscr1, 1 )
 CALL OPEN  ( *7002, iscr1, zi( ibuf2 ), rdrew )
 510   CONTINUE
 CALL OPEN ( *7001, lll, zi( ibuf1 ), wrtrew )
 CALL fname ( lll, NAME )
 CALL WRITE ( lll, NAME, 2, 1 )
 
! DO THE DECOMPOSITION NOW
 
!      CALL AUDIT ( 'SMC2RD  ', 1 )
!      PRINT *,' IILOOP,IOLOOP,NAR,ILSROW,NEXTRA,IDBASE,IWORK,ISPILL'
!      PRINT *,  IILOOP,IOLOOP,NAR,ILSROW,NEXTRA,IDBASE,IWORK,ISPILL
 SELECT CASE ( ktype )
   CASE (    1)
     GO TO  1000
   CASE (    2)
     GO TO  2000
   CASE (    3)
     GO TO  3000
   CASE (    4)
     GO TO  4000
 END SELECT
 1000  CONTINUE
 CALL smc2rs ( zi, zr, zr( iiloop ), zr( ioloop ), zi( nar )  &
     ,    zi( ilsrow ), zr( iwork ), maxnac, maxncol+nextra, maxnar )
 GO TO 5000
 2000  CONTINUE
 CALL smc2rd ( zi, zd, zd( iiloop ), zd( ioloop ), zi( nar )  &
     ,    zi( ilsrow ), zd( iwork ), maxnac, maxncol+nextra, maxnar )
 GO TO 5000
 3000  CONTINUE
!      PRINT *,' CALLING SMC2CS'
 CALL smc2cs ( zi, zr, zd( iiloop ), zd( ioloop ), zi( nar )  &
     ,    zi( ilsrow ), zd( iwork ), maxnac, maxncol+nextra, maxnar )
 GO TO 5000
 4000  CONTINUE
!      PRINT *,' CALLING SMC2CD'
 CALL smc2cd ( zi, zd, zd( iiloop ), zd( ioloop ), zi( nar )  &
     ,    zi( ilsrow ), zd( iwork ), maxnac, maxncol+nextra, maxnar )
 GO TO 5000
 5000  CONTINUE
!      CALL AUDIT ( 'SMC2RD  ', 2 )
 CALL CLOSE ( lll  , 1 )
 CALL CLOSE ( iscr1, 1 )
 GO TO 7777
 7001  CONTINUE
 CALL fname ( lll, NAME )
 ierror = 2
 WRITE ( nout, 9001 ) ufm, lll(1), cname
 9001  FORMAT(1X,a23,/,' SMCPH2 UNABLE TO OPEN FILE ',i4,' ;FILE NAME ='  &
     ,  2A4 )
 GO TO 7100
 7002  CALL fname ( iscr1, NAME )
 WRITE ( nout, 9001 ) ufm, iscr1, cname
 ierror = 3
 GO TO 7100
 7003  CONTINUE
 ierror = 2
 CALL fname ( iscr1, NAME )
 WRITE ( nout, 9001 ) ufm, iscr1, cname
 GO TO 7100
 7008  CONTINUE
 CALL fname ( lll, NAME )
 minum = (.5 * xmaxmem ) - lcore
 WRITE ( lout, 9008 ) ncol, minum
 9008  FORMAT(8X,'INSUFFICIENT OPEN CORE FOR DECOMPOSITION WITH NEW'  &
     ,' METHOD' ,/,    8X,'TOTAL NUMBER OF COLUMNS IN MATRIX =',i8  &
     ,/,    8X,'SUGGESTED ADDITIONAL OPEN CORE IS =',i8)
 CALL CLOSE ( iscr1, 1 )
!      CALL MESAGE ( -8, 0, 0 )
 ierror = 1
 GO TO 7777
 7100  CALL mesage ( -61, 0, 0 )
 7777  CONTINUE
!      CALL AUDIT ( 'SMCPH2  ',2)
!      CALL AUDIT ( 'END     ',1)
!      IF ( NCOL .NE. 0 ) STOP
 RETURN
END SUBROUTINE smcph2
