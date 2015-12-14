SUBROUTINE smcph1 ( zi, zr, zd )
     
 INTEGER, INTENT(OUT)                     :: zi(4)
 REAL, INTENT(OUT)                        :: zr(4)
 DOUBLE PRECISION, INTENT(OUT)            :: zd(4)
 REAL :: minds
 LOGICAL :: frstval
 INTEGER :: itemp(4)
 INTEGER :: prc     ,words   ,rlcmpx  ,NAME(2), rew
 DOUBLE PRECISION :: xnd(10)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 CHARACTER (LEN=4) :: cname(2)
 CHARACTER (LEN=14) :: ctype(4)
 
!  KTYPE   = TYPE (1-RS,2-RD,3-CS,4-CD) OF LOWER TRIANGULAR MATRIX
!  KPREC   = PRECISION (1-SINGL, 2-DOUBL) OF LOWER TRIANGULAR MATRIX
!  MAXROW  = HIGHEST ROW NUMBER REFERENCED THUS FAR IN PROCESSING
!            A GIVEN COLUMN - NEEDED TO DETERMINE CREATED TERMS DURING
!            DECOMPOSITION
!  MAXINLOP= MAXIMUM TERMS FOR ANY GIVEN INNER LOOP
!  MAXNCOL = MAXIMUM NUMBER OF COLUMNS REFERENCED BY ANY GIVEN COLUMN
!  LASCOL  = LAST COLUMN NUMBER OF MATRIX TO BE DECOMPOSED
!  NEXCOL  = FIRST NON-ZERO TERM IN CURRENT PIVOT COLUMN BELOW DIAGONAL
!            USED TO DETERMINE THE NEXT PIVOT COLUMN WHERE THE ROW
!            WILL BE NEEDED.
!  ICURCOL = CURRENT COLUMN BEING PROCESSED
!  MXRECL  = MAXIMUM SIZE IN WORDS OF ANY ONE RECORD WRITTEN TO THE
!            SPILL FILE
!  NSPILL  = NUMBER OF COLUMNS WRITTEN TO THE SPILL FILE
 
 INCLUDE  'SMCOMX.COM'
 COMMON  /xmssg /  ufm     ,uwm     ,uim     ,sfm
 COMMON  /ntime /  nitems  ,tmio    ,tmbpak  ,tmipak  ,tmpak  &
     ,                 tmupak  ,tmgstr  ,tmpstr  ,tmt(4)  ,tml(4)
 COMMON  /names /  rdnrw   ,rdrew   ,wrt     ,wrtrew  ,rew  &
     ,                 norew   ,eofnrw  ,rsp     ,rdp     ,csp  &
     ,                 cdp     ,sqr     ,rect    ,diag    ,lowtri  &
     ,                 uprtri  ,sym
 COMMON  /zzzzzz/  xns(10)
 COMMON  /TYPE  /  prc(2)  , words(4), rlcmpx(4)
 COMMON  /logout/  lout
 EQUIVALENCE       ( ddr     , dsr    ), (ddc    , dsc  )
 EQUIVALENCE       ( mindd   , minds  ), (xns    , xnd  )
 EQUIVALENCE       ( mblk(6) , mterms ), (mblk(5), mstr )
 EQUIVALENCE       ( mblk(4) , mrow   ), (mblk(2), mtype)
 EQUIVALENCE       ( cname   , NAME   )
 DATA              ctype / 'REAL SINGLE   ', 'REAL DOUBLE   '  &
     ,                         'COMPLEX SINGLE', 'COMPLEX DOUBLE' /
 
!   open core is allocated as follows for phase1 of the decomposition
 
!       -------------------------------
!       zi(1) -  Beginning of directory for in-memory column data
!       Directory (4,n) , n=number of columns of matrix
!                  (1,i) = index to active rows and terms within memory
!                  (2,i) = first column data needed for this pivot
!                  (3,i) = last pivot column to use this data
!                  (4,i) = savpos position pointer for data spilled to a
!                          scratch file
!       -------------------------------
!       zi(iacrow) - Beginning of active row vector.
!       Vector for determining active rows for each column, n words
!       Each row value will define the next column where the row value
!       is next needed for calculation of the lll matrix.
!       -------------------------------
!       zi(IRVAL) - Stagging area for storing data
!       Defines the values in the next section of open core, 2*n
!                  (1,i) = row number
!                  (2,i) = number of consecutive terms beginning at row
!       This section and the next section are staging areas for storing
!       of rows and row values of columns to be pointed to by the directory
!       in the first part of open core.
!       -------------------------------
!       zi(IVVAL)
!       Row values of column as defined by previous section, n*iprec words
!       -------------------------------
!       zi(idbase)
!       Memory for rows and terms of columns as pointed to by directory
!       in the first part of open core.  This data is loaded from the
!       bottom up to allow for better management of open core in
!       subroutine smcph2 which is called after this subroutine.
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
!           5+2*m.  row value for first row
!     5+2*m+iprec.  next row value (iprec=1,2,4)
!   5+2*m+iprec*l.  last row value for column (l=total values)
!       -------------------------------
!       zi(ibuf2)
!       Buffer for spill file if all column values can not be kept in memory
!       -------------------------------
!       zi(ibuf1)
!       Buffer for input matrix file to be decomposed
!       -------------------------------
 
!      CALL AUDIT ( 'BEGIN   ', 1 )
!      CALL AUDIT ( 'SMCPH1  ', 1 )
 CALL fname ( mcb, NAME )
 ncol   = mcb( 2 )
 memcoln= 0
 mxrecl = 0
 maxnac = 0
 maxnar = 0
 iprec  = prc  ( mcb( 5 ) )
 ktype  = mcb( 5 )
 IF ( isprec == 2 .AND. ktype == 1 ) ktype = 2
 IF ( isprec == 2 .AND. ktype == 3 ) ktype = 4
 IF ( isprec == 1 .AND. ktype == 2 ) ktype = 1
 IF ( isprec == 1 .AND. ktype == 4 ) ktype = 3
 IF ( ktype  == 1 .OR.  ktype == 3 ) kprec = 1
 IF ( ktype  == 2 .OR.  ktype == 4 ) kprec = 2
 ivwrds = words( ktype )
 iacrow = 4*ncol + 1
 irval  = iacrow + ncol
 ivval  = irval  + 2*ncol
 
! ENSURE THAT IVVAL IS ON A DOUBLE WORD BOUNDARY
 
 IF ( MOD( ivval,2 )  == 0 ) ivval = ivval + 1
 idbase = ivval  + ivwrds*ncol
 
! ENSURE THAT IDBASE IS ON A DOUBLE WORD BOUNDARY
 
 IF ( MOD( idbase,2 ) == 0 ) idbase = idbase + 1
 IF ( lcore < (idbase + 2*isysbf) ) GO TO 7001
 ibuf1  = lcore - isysbf
 ibuf2  = ibuf1 - isysbf
 idbmax = ibuf2 - 1
 
! ENSURE THAT IDBMAX IS ON A DOUBLE WORD BOUNDARY
 
 IF ( MOD( idbmax,2 ) == 0 ) idbmax = idbmax - 1
 idbind = idbmax
 CALL OPEN   ( *7002, mcb, zi(ibuf1), rdrew )
 CALL skprec ( mcb, 1 )
 mblk(1) = mcb( 1 )
 lll(2)  = mcb( 2 )
 lll(3)  = mcb( 2 )
 lll(4)  = 4
 lll(5)  = ktype
 lll(6)  = 0
 lll(7)  = lshift( 1, nbpw-2 - (nbpw-32) )
 icurcol = 1
 opnscr  = .false.
 nspill  = 0
 maxrow  = 0
 maxinlop= 0
 maxncol = 0
 lascol  = 0
 power   = 0
 IF ( kprec /= 2 ) GO TO 5
 ddr     = 1.0D0
 ddc     = 0.0D0
 mindd   = 1.0D+25
 GO TO 8
 5     dsr     = 1.0
 dcr     = 0.0
 minds   = 1.0E+25
 8     CONTINUE
 moblk( 1 ) = lll( 1 )
 moblk( 2 ) = ktype
 moblk( 3 ) = 1
 
! ZERO OUT THE ACTIVE COLUMN VECTOR
 
 DO  i = 1, ncol
   zi( iacrow + i - 1 ) = 0
 END DO
 LEN = ncol*4
 
! ZERO OUT THE DIRECTORY
 
 DO  i = 1, LEN
   zi( i )   = 0
 END DO
 50    CONTINUE
 nterms  = 0
 frstval = .true.
 indexr  = irval
 indexv  = ivval
 indexvd = ( indexv / 2 ) + 1
 mblk(8) = -1
 nexcol  = 0
 inddir = (icurcol-1)*4 + 1
 100   CALL getstr ( *1000, mblk )
 IF ( icurcol <= ( mrow+mterms-1) ) GO TO 120
 
! ALL ROW TERMS ARE BEFORE CURRENT PIVOT COLUMN; SKIP THESE TERMS
! AND GET NEXT STRING.
! CHECK TO SEE IF THIS IS THE FIRST TERM OF THE PIVOT COLUMN
 
 IF ( zi( inddir + 1 ) == 0 ) zi( inddir + 1 ) = mrow
 CALL endget ( mblk )
 GO TO 100
 
! SAVE SOME OR ALL OF THE TERMS IN THIS STRING
! IF THIS IS NOT THE FIRST STRING TO PROCESS, THEN SAVE ALL VALUES
 
 120   iskip = 0
 IF ( .NOT. frstval ) GO TO 140
 
! CHECK IF THIS IS THE FIRST TERM OF THE PIVOT COLUMN
 
 IF ( zi( inddir + 1 ) == 0 ) zi( inddir + 1 ) = mrow
 
! OTHERWISE, CHECK IF ALL TERMS OR ONLY SOME ARE TO BE SAVED
 
 frstval = .false.
 IF ( icurcol == mrow ) GO TO 130
 
! CHECK FOR ZERO ON THE DIAGONAL
 
 IF ( icurcol < mrow ) GO TO 7004
 
! SKIP ALL TERMS BEFORE THE CURRENT PIVOT COLUMN
 
 iskip = icurcol - mrow
 zi( indexr   ) = icurcol
 nterms         = mterms - iskip
 zi( indexr+1 ) = nterms
 IF ( ( mterms-iskip ) > 1 ) nexcol = icurcol + 1
 GO TO 200
 130   CONTINUE
 IF ( mterms > 1 ) nexcol = mrow + 1
 zi(indexr   ) = mrow
 zi(indexr+1 ) = mterms
 nterms        = mterms
 GO TO 200
 140   CONTINUE
 
! CHECK TO SEE IF CURRENT STRING IS AN EXTENSION OF PREVIOUS STRING
 
 IF ( (zi( indexr )+zi( indexr+1 ) ) == mrow ) GO TO 170
 
! NO, MUST CREATE NEW POINTER FOR VALUES
! BUT FIRST, CHECK FOR PROVIDING FOR COMPUTED TERMS OF
! PREVIOUS PIVOT COLUMNS
 
 irow1 = zi( indexr ) + zi( indexr+1 )
 irown = mrow - 1
 irflag = 1
 GO TO 6000
 
! NOW CHECK IF THE ADDED TERMS ARE PART OF SAME STRING AS THAT JUST
! GOTTEN FROM GETSTR CALL
 
 150   IF ( ( zi(indexr) + zi(indexr+1) ) == mrow ) GO TO 170
 
! NEW STRING TO BE DEFINED FOR THE CURRENT TERMS FROM GETSTR
 
 160   indexr = indexr + 2
 zi(indexr   ) = mrow
 zi(indexr+1 ) = mterms
 nterms        = nterms + mterms
 IF ( nexcol == 0 ) nexcol = mrow
 GO TO 200
 170   CONTINUE
 
! TERMS ARE AN EXTENSION OF EXISTING DATA, CHANGE THE NUMBER OF TERMS
 
 zi( indexr+1 ) = zi( indexr+1 ) + mterms
 nterms         = nterms + mterms
 IF ( nexcol == 0 ) nexcol = mrow
 200   CALL smcrtr ( zr, zd )
 
! SET ACTIVE COLUMN ROW NUMBERS FOR POSSIBLE EXPANDED TERMS
 
 irow1 = mrow + iskip
 irown = irow1 + mterms - 1 - iskip
 DO  k = irow1, irown
   zi( iacrow + k - 1 ) = nexcol
 END DO
 
! GO AND GET ADDITIONAL STRINGS IF ANY
 
 CALL endget ( mblk )
 GO TO 100
 1000  CONTINUE
 
! END OF READING CURRENT COLUMN, CHECK IF DIAGONAL TERM FOUND
 
!      PRINT *,' SMCPH1,ICURCOL,NEXCOL,MAXROW=',ICURCOL,NEXCOL,MAXROW
 IF ( frstval ) GO TO 7004
 
! SEE IF ANY COMPUTED TERMS FROM PREVIOUS PIVOT COLUMNS ARE TO BE
! ADDED ONTO THE END OF THE CURRENT ACTIVE ROWS FOR THIS COLUMN
 
 lrow  = zi( indexr ) + zi( indexr+1 ) - 1
 IF ( lrow > maxrow ) maxrow = lrow
 irow1 = lrow + 1
 irown = maxrow
 irflag = 2
!      PRINT *,' B1050,ICURCOL,IROWN,IROW1=',ICURCOL,IROWN,IROW1
 IF ( irown >= irow1 ) GO TO 6000
 
! SET UP DIRECTORY AND SAVE DATA EITHER
! IN MEMORY OR ON SPILL FILE
 
 1050 CONTINUE
 
! RECOMPUTE LROW IN CASE NEW TERMS WERE ADDED FROM PREVIOUS PIVOT COLUMNS
 
 lrow   = zi( indexr ) + zi( indexr+1 ) - 1
 
! INDEXR POINTS TO CURRENT DIRECTORY ENTRY BUT INDEXV POINTS TO NEXT
! AVAILABLE POSITION FOR STORING TERMS
 
 nrvals = indexr - irval + 2
 nvvals = indexv - ivval
 nwords = nrvals + nvvals + 4
 
! SAVE DATA IN MEMORY AND SET DIRECTORY ACCORDINGLY
 
 itest = idbind - nwords + 1
 
! MAKE SURE ITEST IS ON DOUBLE WORD BOUNDARY
 
 IF ( MOD( itest,2 ) == 0 ) itest = itest - 1
 
! CHECK TO SEE IF THERE IS SUFFICIENT MEMORY
 
 maxnar = MAX0( nrvals, maxnar )
 IF ( itest < idbase ) GO TO 1800
 idbind = itest
 zi( inddir     ) = idbind
 zi( inddir + 3 ) = 0
 zi( idbind     ) = icurcol
 zi( idbind + 1 ) = nrvals
 zi( idbind + 2 ) = nwords
 zi( idbind + 3 ) = nterms
 idbind = idbind + 3
 DO  k = 1, nrvals
   zi( idbind + k ) = zi( irval + k - 1 )
 END DO
 idbind = idbind + nrvals
 IF ( kprec == 2 ) GO TO 1300
 DO  k = 1, nvvals
   zi( idbind + k ) = zi( ivval + k - 1 )
 END DO
 GO TO 1400
 1300  indxv = idbind / 2
 nv    = nvvals / 2
 ivd   = ivval  / 2
 DO  k = 1, nv
   zd( indxv+k ) = zd( ivd+k )
 END DO
 1400  CONTINUE
 idbind = idbind + nvvals - nwords
 lascol = icurcol
 memcoln= icurcol
 itest  = nrvals + nvvals + 4
 IF ( itest > mxrecl ) mxrecl = itest
 GO TO 2000
 1800  CONTINUE
 IF ( opnscr ) GO TO 1810
 opnscr = .true.
 CALL OPEN ( *7003, iscr1, zi(ibuf2), wrtrew )
 
! NO MORE MEMORY, SAVE COLUMN DATA TO SPILL FILE, KEEP RECORD POSITION
 
 1810  itemp( 1 ) = icurcol
 itemp( 2 ) = nrvals
 itemp( 3 ) = 0
 itemp( 4 ) = nterms
 CALL WRITE ( iscr1, itemp, 4, 0 )
 CALL savpos( iscr1, kpos )
 CALL WRITE ( iscr1, zi( irval ), indexr-irval+2, 0 )
 CALL WRITE ( iscr1, zi( ivval ), indexv-ivval+2, 1 )
 zi( inddir   ) = 0
 zi( inddir+3 ) = kpos
 nspill = nspill + 1
 itest  = nrvals + nvvals + 4
 IF ( itest > mxrecl ) mxrecl = itest
 2000  CONTINUE
 lrow   = zi( indexr ) + zi( indexr+1 ) - 1
 
! SAVE LAST PIVOT COLUMN FOR WHICH DATA IN THIS COLUMN IS USED
 
 IF ( nterms > maxnac ) maxnac = nterms
 zi( inddir+2 ) = lrow
 ifirstc = zi( inddir+1 )
 maxtes  = ( icurcol - ifirstc + 1 )
 IF ( maxtes > maxncol  ) maxncol = maxtes
 maxtes  = nterms * ( icurcol - ifirstc )
 IF ( maxtes > maxinlop ) maxinlop = maxtes
 
! CHECK TO DETERMINE IF ALL COLUMNS HAVE BEEN PROCESSED
 
 IF ( icurcol >= ncol ) GO TO 7777
 
! CHECK IF ONLY ONE TERM IN THIS COLUMN
 
 IF ( nexcol /= 0 ) GO TO 2005
 
! MUST FIND FIRST NON-ZERO TERM FOLLOWING THE CURRENT PIVOT
 
 DO  k = icurcol+1, ncol
   IF ( zi( iacrow+k ) /= icurcol ) CYCLE
   nexcol = k
   GO TO 2005
 END DO
 WRITE ( nout, 9901 ) icurcol
 GO TO 2030
 9901  FORMAT(' SYMMETRIC DECOMPOSITION FOUND NO TERMS BEING '  &
     ,' CONNECTED TO DIAGONAL ON COLUMN ',i6)
 2005  CONTINUE
!      PRINT *,' AFTER 2005,ICURCOL,NEXCOL=',ICURCOL,NEXCOL
 
! UPDATE ACTIVE ROWS IN COLUMN VECTOR FOR ALL TERMS OF THIS COLUMN
 
 LEN = irval + nrvals - 1
 DO  k = irval, LEN, 2
   irow = zi( k )
   nrow = irow + zi( k+1 ) - 1
   DO  l = irow, nrow
     zi( iacrow + l - 1 ) = nexcol
   END DO
 END DO
 DO  l = icurcol+1, maxrow
   IF ( zi( iacrow+l-1 ) == icurcol ) zi( iacrow+l-1) = nexcol
 END DO
 2030  CONTINUE
 
! END OF CURRENT COLUMN, PREPARE FOR NEXT COLUMN
 
!      write ( nout, 901 ) icurcol
!901   format(20x,' Active rows after processing column ',i10)
!      do 2040 l = 1, ncol
!      write ( nout, 902 ) l, zi(iacrow+l-1)
!902   format(' Row, next reference =',2i7)
!2040  continue
!      write ( nout, 903 )
!903   format(20x, ' Directory',/,
!     &' Column  Memory Index   First Used    Last Used    Savpos')
!      do 2050 l = 1, ncol
!      ind = ( l-1 ) * 4 + 1
!      write ( nout, 904 ) l, zi(ind), zi(ind+1), zi(ind+2), zi(ind+3)
!904   format( i7, i14, i13, i13, i9)
!2050  continue
 icurcol = icurcol + 1
 inddir  = (icurcol-1)*4 + 1
 GO TO 50
 
! THE FOLLOWING IS AN INTERNAL ROUTINE TO ADD COMPUTED TERMS RESULTING
! FROM THE PROCESSING OF PREVIOUS PIVOT COLUMNS INTO THE CURRENT ACTIVE
! ROWS FOR THE CURRENT COLUMN
 
 6000  CONTINUE
 DO  k = irow1, irown
   IF ( zi(iacrow + k - 1 ) < icurcol ) CYCLE
   IF ( nexcol == 0 ) nexcol = k
   
! NEED TO ADD THIS TERM TO THE ACTIVE ROWS
! CHECK TO SEE IF THIS TERM IS AN EXTENSION OF CURRENT TERMS
   
   IF ( (zi( indexr ) + zi( indexr+1 ) ) == k ) GO TO 6010
   
! NO, NEED TO CREATE ANOTHER POINTER
   
   indexr = indexr + 2
   zi( indexr )   = k
   zi( indexr+1 ) = 1
   nterms         = nterms +1
   GO TO 6020
   6010  CONTINUE
   
! JUST ADD TO THE NUMBER OF CONSECUTIVE VALUES FOR CURRENT ROW
   
   zi( indexr+1 ) = zi( indexr+1 ) + 1
   nterms         = nterms + 1
   6020  CONTINUE
   
! NOW, ZERO OUT ROW VALUE
   
   SELECT CASE ( ktype )
     CASE (    1)
       GO TO  6030
     CASE (    2)
       GO TO  6040
     CASE (    3)
       GO TO  6050
     CASE (    4)
       GO TO  6060
   END SELECT
   
! TYPE IS REAL SP
   
   6030  zr( indexv ) = 0.
   indexv = indexv + 1
   CYCLE
   
! TYPE IS REAL DP
   
   6040  zd( indexvd ) = 0.d0
   indexvd = indexvd + 1
   indexv  = indexv  + 2
   CYCLE
   
! TYPE IS COMPLEX SP
   
   6050  zr( indexv   ) = 0.
   zr( indexv+1 ) = 0.
   indexv = indexv + 2
   CYCLE
   
! TYPE IS COMPLEX DP
   
   6060  zd( indexvd   ) = 0.d0
   zd( indexvd+1 ) = 0.d0
   indexvd = indexvd + 2
   indexv  = indexv  + 4
   CYCLE
 END DO
 SELECT CASE ( irflag )
   CASE (    1)
     GO TO  150
   CASE (    2)
     GO TO  1050
 END SELECT
 
! INSUFFICIENT MEMORY
 
 7001  minmum = ncol*7 + 2*ncol*ivwrds + 2*isysbf
 WRITE ( nout, 9001 ) ufm, mcb(1), cname, ncol, ktype  &
     ,                    lcore, minmum
 9001  FORMAT(1X,a23,/,' INSUFFICIENT MEMORY TO DECOMPOSE MATRIX IN '  &
     ,i4,' FILE NAME=',2A4  &
     ,/,' NUMBER OF COLUMNS=',i7,' TYPE=',i2,' MEMORY AVAILABLE =',i10  &
     ,/,' MINIMUM REQUIRED IS =',i10)
!      CALL MESAGE ( -8, 0, 0 )
 ierror = 1
 GO TO 7777
 7002  CALL fname ( mcb, NAME )
 WRITE ( nout, 9002 ) ufm, mcb(1), cname
 9002  FORMAT(1X, a23, /,' SMCPH1 UNABLE TO OPEN FILE ',i4,' NAME= ',2A4)
 ierror = 2
 CALL mesage ( -61, 0, 0 )
 7003  CALL fname ( iscr1, NAME )
 WRITE ( nout, 9003 ) ufm, iscr1, cname
 9003  FORMAT(1X, a23, /,' SMCPH1 UNABLE TO OPEN FILE ',i4,' NAME= ',2A4)
 ierror = 2
 CALL mesage ( -61, 0, 0 )
 
! ZERO ON DIAGONAL, TERMINATE DECOMPOSITION BUT FIRST SCAN REST OF
! MATRIX TO DETERMINE OTHER COLUMNS WITH ZERO DIAGONALS.
 
 7004  CONTINUE
 ierror = 7
 izeros = 1
 indexz = 0
 IF ( frstval ) GO TO 7020
 CALL endget ( mblk )
 7010  CALL skprec ( mblk, 1 )
 7020  indexz = indexz + 1
 zi ( indexz ) = icurcol
 7025  icurcol = icurcol + 1
 IF ( icurcol > ncol ) GO TO 7050
 mblk( 8 ) = -1
 7030  CALL getstr ( *7020, mblk )
 CALL endget ( mblk )
 IF ( icurcol >= mrow .AND. icurcol <= mrow+mterms-1)GO TO 7040
 IF ( mrow    > icurcol ) GO TO 7010
 GO TO 7030
 7040  CALL skprec ( mblk, 1 )
 GO TO 7025
 7050  CALL CLOSE ( mcb  , rew )
 CALL CLOSE ( iscr1, rew )
 WRITE ( nout, 9050 ) ufm, cname, (zi(k),k=1,indexz)
 9050  FORMAT(a23,' 3097, SYMMETRIC DECOMPOSITION OF DATA BLOCK ',2A4  &
     ,      ' ABORTED BECAUSE THE FOLLOWING COLUMNS ARE SINGULAR -'  &
     ,/,(5X,20I6,/))
 RETURN
 7777  CONTINUE
!      CALL SMCHLP
!      CALL SMCDMP ( ZI, ZR, ZD )
 7778  CONTINUE
 CALL CLOSE ( mcb, rew )
 itwrds  = idbmax - idbind
 itcols  = ncol   - nspill
 CALL sswtch ( 45, l45 )
 IF ( l45 == 0 ) GO TO 7779
 WRITE ( lout, 9004 ) itcols , nspill , maxnac, maxncol  &
     ,                    maxinlop, itwrds
 9004  FORMAT(/ ,14X,' STATISTICS FOR SYMMETRIC DECOMPOSITION OF FILE ',/  &
     ,/, 7X,' COLUMNS CONTAINED IN MEMORY                         =',i8  &
     ,/, 7X,' COLUMNS WRITTEN TO SPILL FILE                       =',i8  &
     ,/, 7X,' MAX. NO. OF ACTIVE ROWS FOR ANY ACTIVE COLUMN       =',i8  &
     ,/, 7X,' MAX. NUMBER OF COLUMNS REFERENCED BY A PIVOT COLUMN =',i8  &
     ,/, 7X,' MAX. TERMS FOR ANY GIVEN INNER LOOP                 =',i8  &
     ,/, 7X,' TOTAL WORDS IN OPEN CORE USED FOR COLUMN DATA       =',i8 )
 WRITE ( lout, 9005 ) 'INPUT ', cname, ctype( mcb( 5 ) )
 CALL fname ( lll, NAME )
 WRITE ( lout, 9005 ) 'OUTPUT', cname, ctype( ktype )
 9005  FORMAT( 8X, a6,' FILE: ',2A4     ,'      DATA TYPE= ',a14 )
!      CALL AUDIT( 'SMCPH1  ', 2 )
 7779  CONTINUE
 RETURN
END SUBROUTINE smcph1
