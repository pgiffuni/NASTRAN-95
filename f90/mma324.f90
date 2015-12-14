SUBROUTINE mma324 ( zi, zd )
     
!     MMA324 PERFORMS THE MATRIX OPERATION IN COMPLEX DOUBLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA324 USES METHOD 32 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  CALL MMARM1 TO PACK AS MANY COLUMNS OF "A" INTO MEMORY
!           AS POSSIBLE LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  ADD EACH ROW TERM OF "C" TO "D" MATRIX COLUMN
!       4.  CALL MMARC1,2,3,4 TO READ COLUMNS OF "B" AND "C".
 
 
 INTEGER, INTENT(IN)                      :: zi(2)
 DOUBLE PRECISION, INTENT(IN)             :: zd(2)
 INTEGER :: t
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: ofile      ,filea    ,fileb ,filec , filed
 
 DOUBLE COMPLEX    cdtemp    ,cd
 INCLUDE           'MMACOM.COM'
 COMMON / names  / rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 COMMON / TYPE   / iprc(2)    ,nwords(4),irc(4)
 COMMON / mpyadx / filea(7)   ,fileb(7) ,filec(7)  &
     ,                 filed(7)   ,nz       ,t     ,signab,signc ,prec1  &
     ,                 scrtch     ,time
 COMMON / system / ksystm(152)
 COMMON / zblpkx / d(4)       ,kdrow
 COMMON / unpakx / typeu      ,iurow1   ,iurown, incru
 COMMON / packx  / typei      ,typep    ,iprow1, iprown , incrp
 EQUIVALENCE       ( d(1)     ,cd    )
 EQUIVALENCE       (ksystm( 1),sysbuf)  , (ksystm( 2),iwr   )
 EQUIVALENCE       (filea(2)  ,nac   )  , (filea(3)  ,nar   )  &
     ,                 (filea(4)  ,naform)  , (filea(5)  ,natype)  &
     ,                 (filea(6)  ,nanzwd)  , (filea(7)  ,nadens)
 EQUIVALENCE       (fileb(2)  ,nbc   )  , (fileb(3)  ,nbr   )  &
     ,                 (fileb(4)  ,nbform)  , (fileb(5)  ,nbtype)  &
     ,                 (fileb(6)  ,nbnzwd)  , (fileb(7)  ,nbdens)
 EQUIVALENCE       (filec(2)  ,ncc   )  , (filec(3)  ,ncr   )  &
     ,                 (filec(4)  ,ncform)  , (filec(5)  ,nctype)  &
     ,                 (filec(6)  ,ncnzwd)  , (filec(7)  ,ncdens)
 EQUIVALENCE       (filed(2)  ,ndc   )  , (filed(3)  ,ndr   )  &
     ,                 (filed(4)  ,ndform)  , (filed(5)  ,ndtype)  &
     ,                 (filed(6)  ,ndnzwd)  , (filed(7)  ,nddens)
 
!   OPEN CORE ALLOCATION AS FOLLOWS:
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "B" MATRIX IN COMPACT FORM
!     Z( IDX      ) = ARRAY FOR ONE COLUMN OF "D" MATRIX
!     Z( IAX      ) = ARRAY FOR MULTIPLE COLUMNS OF "A" MATRIX
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE
 
 filed( 2 ) = 0
 filed( 6 ) = 0
 filed( 7 ) = 0
 idrow = ibrow
 DO  ii = 1, nbc
   CALL bldpk ( ndtype, ndtype, ofile, 0, 0 )
!      PRINT *,' PROCESSING B MATRIX COLUMN, II=',II
   
! READ A COLUMN FROM THE "B" MATRIX
   
   SIGN   = 1
   irfile = fileb( 1 )
   CALL mmarc4 ( zi, zd )
   lasindb = lasind
   
! NOW READ "C", OR SCRATCH FILE WITH INTERMEDIATE RESULTS.
! IF NO "C" FILE AND THIS IS THE FIRST PASS, INITIALIZE "D" COLUMN TO ZERO.
   
   IF ( ifile == 0 ) GO TO 950
   IF ( ipass == 1 ) SIGN = signc
   irfile = ifile
   
! READ A COLUMN FROM THE "C" MATRIX
   
   CALL mmarc4 ( zi( idx ), zd( idx2+1 ) )
   lasindc = lasind + idx - 1
   950   CONTINUE
   
! CHECK IF COLUMN OF "B" IS NULL
   
   IF ( zi( 1 ) /= 0 ) GO TO 1000
   IF ( ifile /= 0 ) GO TO 960
   951   cd = ( 0.0D0, 0.0D0)
   kdrow = 1
   CALL zblpki
   GO TO 55000
   960   IF ( zi( idx ) == 0 ) GO TO 951
   indxc = idx
   961   IF ( indxc >= lasindc ) GO TO 55000
   irowc1 = zi( indxc )
   icrows = zi( indxc+1 )
   irowcn = irowc1 + icrows - 1
   indxcv = ( indxc+3 ) / 2
   DO  i = irowc1, irowcn
     cd = DCMPLX( zd( indxcv ), zd( indxcv+1 ) )
     kdrow  = i
     CALL zblpki
     indxcv = indxcv + 2
   END DO
   indxc  = indxc + 2 + icrows*nwdd
   GO TO 961
   1000  CONTINUE
   irowb1 = zi( 1 )
   irows  = zi( 2 )
   irowbn = irowb1 + irows - 1
   indxb  = 1
   indxa  = iax
   indxc  = idx
   IF ( ifile == 0 .OR. indxc >= lasindc ) GO TO 9000
   irowc1 = zi( indxc )
   icrows = zi( indxc+1 )
   irowcn = irowc1 + icrows - 1
   1010  CONTINUE
   
! CHECK TO ADD TERMS FROM "C" OR INTERIM SCRATCH FILE BEFORE CURRENT ROW
   
   IF ( idrow == 0 .OR. irowc1 > idrow ) GO TO 9000
   4000  CONTINUE
   irown = idrow
   IF ( irowcn < idrow ) irown = irowcn
   5000  CONTINUE
   indxcv = ( indxc+3 ) / 2
   nrows = irown - irowc1 + 1
   DO  i = 1, nrows
     kdrow = irowc1 + i - 1
     cd = DCMPLX( zd( indxcv ), zd( indxcv+1 ) )
     indxcv = indxcv + 2
     CALL zblpki
   END DO
   IF ( irowcn >= idrow ) GO TO 9000
   indxc = indxc + 2 + icrows*nwdd
   IF ( indxc >= lasindc ) GO TO 9000
   irowc1 = zi ( indxc )
   icrows = zi ( indxc+1 )
   irowcn = irowc1 + icrows - 1
   GO TO 4000
   9000  CONTINUE
   
! CHECK FOR NULL COLUMN FROM "B" MATRIX
   
   IF ( irowb1 == 0 ) GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
! COMPLEX DOUBLE PRECISION
   10000 CONTINUE
   DO  i = 1, ncolpp
     cd  = ( 0.0D0, 0.0D0 )
     kdrow  = idrow + i
     icola  = ibrow + i
     IF ( icola /= IABS( zi( indxa ) ) ) GO TO 70001
     indxal = zi( indxa+1 ) + iax - 1
     indxa  = indxa + 2
     indxb  = 1
     11000 IF ( indxb >= lasindb ) GO TO 14500
     irowb1 = zi( indxb )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     indxbs = indxb
     indxb  = indxb + 2 + irows*nwdd
     12000 CONTINUE
     IF ( indxa >= indxal ) GO TO 14500
     irowa1 = zi( indxa )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     IF ( irowbn < irowa1 ) GO TO 11000
     IF ( irowan < irowb1 ) GO TO 14200
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     indxav = ( ( indxa +3 ) / 2 ) + 2*( irow1 - irowa1 ) - 1
     indxbv = ( ( indxbs+3 ) / 2 ) + 2*( irow1 - irowb1 ) - 1
     cdtemp = ( 0.0D0, 0.0D0)
     kcnt   = ( irown - irow1 ) * 2 + 1
     DO  k = 1, kcnt, 2
       cdtemp = cdtemp + DCMPLX( zd( indxav+k), zd( indxav+k+1 ) ) *  &
           DCMPLX( zd( indxbv+k), zd( indxbv+k+1 ) )
     END DO
     cd = cd + cdtemp
     IF ( irowan > irowbn ) GO TO 11000
     14200 CONTINUE
     indxa  = indxa + 2 + ntms*nwdd
     GO TO 12000
     14500 indxa  = indxal
     14510 IF ( indxc >= lasindc .OR. ifile == 0 ) GO TO 14600
     IF ( kdrow < irowc1 ) GO TO 14600
     IF ( kdrow > irowcn ) GO TO 14550
     indxcv  = ( indxc+3 ) / 2 + 2*( kdrow - irowc1 )
     cd = cd + DCMPLX( zd( indxcv ), zd( indxcv+1 ) )
     GO TO 14600
     14550 indxc   = indxc + 2 + icrows*nwdd
     IF ( indxc >= lasindc ) GO TO 14600
     irowc1  = zi( indxc )
     icrows  = zi( indxc+1 )
     irowcn  = irowc1 + icrows - 1
     GO TO 14510
     14600 CONTINUE
     CALL zblpki
   END DO
   50000 CONTINUE
   IF ( kdrow == ndr .OR. ifile == 0 .OR. indxc >= lasindc ) GO TO 55000
   
! ADD REMAINING TERMS FROM EITHER THE "C" MATRIX OR INTERIM SCRATCH MATRIX
   
   irow1 = kdrow + 1
   50100 CONTINUE
   indxcv = ( indxc+3 ) / 2
   IF ( irow1 < irowc1 ) GO TO 51000
   IF ( irow1 <= irowcn ) GO TO 50900
   indxc   = indxc + 2 + icrows*nwdd
   IF ( indxc >= lasindc ) GO TO 55000
   irowc1  = zi( indxc )
   icrows  = zi( indxc+1 )
   irowcn  = irowc1 + icrows - 1
   GO TO 50100
   50900 CONTINUE
   indxcv = ( indxc+3 ) / 2 + 2*( irow1 - irowc1 )
   irowc1 = irow1
   51000 CONTINUE
   nrows = irowcn - irowc1 + 1
   DO  k = 1, nrows
     kdrow = irowc1 + k - 1
     cd = DCMPLX( zd( indxcv ), zd( indxcv+1 ) )
     indxcv = indxcv + 2
     CALL zblpki
   END DO
   indxc   = indxc + 2 + icrows*nwdd
   IF ( indxc >= lasindc ) GO TO 55000
   irowc1  = zi( indxc )
   icrows  = zi( indxc+1 )
   irowcn  = irowc1 + icrows - 1
   indxcv  = ( indxc+3 ) / 2
   GO TO 51000
   55000 CONTINUE
   CALL bldpkn ( ofile, 0, filed )
! END OF PROCESSING THIS COLUMN FOR THIS PASS
 END DO
 GO TO 70000
 70001 CONTINUE
 WRITE ( iwr, 90001 ) icola, zi( indxa ), iax, indxa
 90001 FORMAT(' UNEXPECTED COLUMN FOUND IN PROCESSING MATRIX A'  &
     ,/,' COLUMN EXPECTED:',i6 ,/,' COLUND FOUND   :',i6  &
     ,/,' IAX =',i7,' INDXA=',i7 )
 CALL mesage ( -61, 0, 0 )
 70000 RETURN
END SUBROUTINE mma324
