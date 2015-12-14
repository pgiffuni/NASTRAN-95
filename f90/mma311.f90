SUBROUTINE mma311 ( zi, zr )
     
!     MMA311 PERFORMS THE MATRIX OPERATION IN REAL SINGLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA311 USES METHOD 31 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  CALL MMARM1 TO PACK AS MANY COLUMNS OF "A" INTO MEMORY
!           AS POSSIBLE LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  UNPACK COLUMNS OF "C" MATRIX BUT USE GETSTR (MMARC1,2,3,4)
!           TO READ COLUMNS OF "B".
 
 
 INTEGER, INTENT(IN)                      :: zi(2)
 REAL, INTENT(OUT)                        :: zr(2)
 INTEGER :: t
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: ofile      ,filea    ,fileb ,filec , filed
 
 INCLUDE           'MMACOM.COM'
 COMMON / names  / rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 COMMON / TYPE   / iprc(2)    ,nwords(4),irc(4)
 COMMON / mpyadx / filea(7)   ,fileb(7) ,filec(7)  &
     ,                 filed(7)   ,nz       ,t     ,signab,signc ,prec1  &
     ,                 scrtch     ,time
 COMMON / system / ksystm(152)
 COMMON / unpakx / typeu      ,iurow1   ,iurown, incru
 COMMON / packx  / typei      ,typep    ,iprow1, iprown , incrp
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
 
 irfile = fileb( 1 )
 DO  ii = 1, nbc
!      PRINT *,' PROCESSING COLUMN=',II
   
! READ A COLUMN FROM THE "B" MATRIX
   
   CALL mmarc1 ( zi, zr )
   
! NOW READ "C", OR SCRATCH FILE WITH INTERMEDIATE RESULTS.
! IF NO "C" FILE AND THIS IS THE FIRST PASS, INITIALIZE "D" COLUMN TO ZERO.
   
   IF ( ifile == 0 ) GO TO 950
   iurow1 = 1
   iurown = ndr
   typeu  = ndtype
   IF ( ipass == 1 ) typeu = ndtype * signc
   CALL unpack (*950, ifile, zr( idx ) )
   GO TO 980
   950   CONTINUE
   DO  j = 1, ndr
     zr( idx+j-1 ) = 0
   END DO
   980   CONTINUE
   
! CHECK IF COLUMN OF "B" IS NULL
   
   irowb1 = zi( 1 )
   irows  = zi( 2 )
   irowbn = irowb1 + irows - 1
   indxb  = 1
   indxa  = iax
   
! CHECK FOR NULL COLUMN FROM "B" MATRIX
   
   IF ( irowb1 == 0 ) GO TO 50000
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A * B  +  C )
   
! SINGLE PRECISION
   1000  CONTINUE
   DO  i = 1, ncolpp
     indxal = zi( indxa + 1 ) + iax - 1
     icola  = ibrow+i
     IF ( icola /= IABS( zi( indxa ) ) ) GO TO 70001
     indxa  = indxa + 2
     1100  CONTINUE
     IF ( icola < irowb1 ) GO TO 1450
     IF ( icola <= irowbn ) GO TO 1200
     indxb  = indxb + 2 + irows
     IF ( indxb > lasind ) GO TO 50000
     irowb1 = zi( indxb )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     GO TO 1100
     1200  CONTINUE
     indxbv = icola - irowb1 + indxb + 2
     IF ( zr( indxbv ) == 0. ) GO TO 1450
     1300  CONTINUE
     IF ( indxa >= indxal ) GO TO 1450
     irowa1 = zi( indxa )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     indxav = indxa + 2 - irowa1
     DO  k = irowa1, irowan
       zr( idx+k-1 ) = zr( idx+k-1 ) +  zr( indxav+k ) * zr( indxbv )
     END DO
     indxa  = indxa + 2 + ntms
     GO TO 1300
     1450  indxa  = indxal
   END DO
   GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
   idrow = ibrow
   idxx  = idx + idrow - 1
! SINGLE PRECISION
   10000 CONTINUE
   DO  i = 1, ncolpp
     icola  = ibrow + i
     IF ( icola /= IABS( zi( indxa ) ) ) GO TO 70001
     indxal = zi( indxa+1 ) + iax - 1
     indxa  = indxa + 2
     indxb  = 1
     11000 IF ( indxb >= lasind ) GO TO 14500
     irowb1 = zi( indxb )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     indxbv = indxb + 2 - irowb1
     indxb  = indxb + 2 + irows
     12000 CONTINUE
     IF ( indxa >= indxal ) GO TO 14500
     irowa1 = zi( indxa )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     IF ( irowbn < irowa1 ) GO TO 11000
     IF ( irowan < irowb1 ) GO TO 14200
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     indxav = indxa + 2 - irowa1
     DO  k = irow1, irown
       zr( idxx+i ) = zr( idxx+i ) +  zr( indxav+k ) * zr( indxbv+k )
     END DO
     IF ( irowan > irowbn ) GO TO 11000
     14200 CONTINUE
     indxa  = indxa + 2 + ntms
     GO TO 12000
     14500 indxa  = indxal
   END DO
   GO TO 50000
! END OF PROCESSING THIS COLUMN FOR THIS PASS
   50000 CONTINUE
!  NOW SAVE COLUMN
   CALL pack ( zr( idx ), ofile, filed )
 END DO
 GO TO 70000
 70001 CONTINUE
 WRITE ( iwr, 90001 ) icola, zi( indxa ), iax, indxa
 90001 FORMAT(' UNEXPECTED COLUMN FOUND IN PROCESSING MATRIX A'  &
     ,/,' COLUMN EXPECTED:',i6 ,/,' COLUND FOUND   :',i6  &
     ,/,' IAX =',i7,' INDXA=',i7 )
 CALL mesage ( -61, 0, 0 )
 70000 RETURN
END SUBROUTINE mma311
