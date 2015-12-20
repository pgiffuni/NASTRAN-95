SUBROUTINE mma314 ( zi, zd, zdc )
     
!     MMA314 PERFORMS THE MATRIX OPERATION IN COMPLEX DOUBLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA314 USES METHOD 31 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  CALL MMARM1 TO PACK AS MANY COLUMNS OF "A" INTO MEMORY
!           AS POSSIBLE LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  UNPACK COLUMNS OF "C" MATRIX BUT USE GETSTR (MMARC1,2,3,4)
!           TO READ COLUMNS OF "B".
 
 
 INTEGER, INTENT(IN)                      :: zi(2)
 DOUBLE PRECISION, INTENT(IN)             :: zd(2)
 INTEGER :: t
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: ofile      ,filea    ,fileb ,filec , filed
 
 DOUBLE COMPLEX    zdc(2)
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
   
   CALL mmarc4 ( zi, zd )
   
! NOW READ "C", OR SCRATCH FILE WITH INTERMEDIATE RESULTS.
! IF NO "C" FILE AND THIS IS THE FIRST PASS, INITIALIZE "D" COLUMN TO ZERO.
   
   IF ( ifile == 0 ) GO TO 950
   iurow1 = 1
   iurown = ndr
   typeu  = ndtype
   IF ( ipass == 1 ) typeu = ndtype * signc
   CALL unpack (*950, ifile, zdc( idx4+1 ) )
   GO TO 980
   950   CONTINUE
   DO  j = 1, ndr
     zdc( idx4+j ) = 0
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
   
! DOUBLE PRECISION
   1000  CONTINUE
   DO  i = 1, ncolpp
     indxal = zi( indxa + 1 ) + iax - 1
     icola  = ibrow+i
     IF ( icola /= IABS( zi( indxa ) ) ) GO TO 70001
     indxa  = indxa + 2
     1100  CONTINUE
     IF ( icola < irowb1 ) GO TO 1450
     IF ( icola <= irowbn ) GO TO 1200
     indxb  = indxb + 2 + irows*nwdd
     IF ( indxb > lasind ) GO TO 50000
     irowb1 = zi( indxb )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     GO TO 1100
     1200  CONTINUE
     indxbv = 2 * ( icola - irowb1 ) + ( indxb + 3 ) / 2
     IF ( zd( indxbv ) == 0. ) GO TO 1450
     1300  CONTINUE
     IF ( indxa >= indxal ) GO TO 1450
     irowa1 = zi( indxa )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     indxav = ( (indxa+3 ) / 2 )
     DO  k = irowa1, irowan
       zdc( idx4+k ) = zdc( idx4+k ) + DCMPLX( zd(indxav), zd(indxav+1) ) *  &
           DCMPLX( zd(indxbv), zd(indxbv+1) )
       indxav = indxav + 2
     END DO
     indxa  = indxa + 2 + ntms*nwdd
     GO TO 1300
     1450  indxa  = indxal
   END DO
   GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
   idrow = ibrow
   idxx  = idx4 + idrow
! DOUBLE PRECISION
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
     kcnt   = ( irown - irow1 ) * 2 + 1
     DO  k = 1, kcnt, 2
       zdc( idxx+i ) = zdc( idxx+i ) +  &
           DCMPLX(  zd( indxav+k ), zd( indxav+k+1 ) ) *  &
           DCMPLX(  zd( indxbv+k ), zd( indxbv+k+1 ) )
     END DO
     IF ( irowan > irowbn ) GO TO 11000
     14200 CONTINUE
     indxa  = indxa + 2 + ntms*nwdd
     GO TO 12000
     14500 indxa  = indxal
   END DO
   GO TO 50000
! END OF PROCESSING THIS COLUMN FOR THIS PASS
   50000 CONTINUE
!  NOW SAVE COLUMN
   CALL pack ( zdc( idx4+1 ), ofile, filed )
 END DO
 GO TO 70000
 70001 CONTINUE
 WRITE ( iwr, 90001 ) icola, zi( indxa ), iax, indxa
 90001 FORMAT(' UNEXPECTED COLUMN FOUND IN PROCESSING MATRIX A'  &
     ,/,' COLUMN EXPECTED:',i6 ,/,' COLUND FOUND   :',i6  &
     ,/,' IAX =',i7,' INDXA=',i7 )
 CALL mesage ( -61, 0, 0 )
 70000 RETURN
END SUBROUTINE mma314
