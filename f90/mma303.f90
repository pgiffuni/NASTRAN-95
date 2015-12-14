SUBROUTINE mma303 ( zi, zc )
     
!     MMA303 PERFORMS THE MATRIX OPERATION USING METHOD 30 AND
!       COMPLEX SINGLE PRECISION
 
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA303 USES METHOD 30 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  CALL 'MMARM1' TO PACK AS MANY COLUMNS OF "A" INTO MEMORY
!           AS POSSIBLE LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  USE UNPACK TO READ MATRICES "B" AND "C".
 
 
 INTEGER, INTENT(IN)                      :: zi(2)
 COMPLEX, INTENT(OUT)                     :: zc(2)
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
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "B" MATRIX
!     Z( IDX      ) = ARRAY FOR ONE COLUMN OF "D" MATRIX
!     Z( IAX      ) = ARRAY FOR MULTIPLE COLUMNS OF "A" MATRIX
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE
 
! PROCESS ALL OF THE COLUMNS OF "B";  ADD "C" DATA ON FIRST PASS
 DO  ii = 1, nbc
   iurow1 = -1
   typeu  = ndtype
   CALL unpack ( *930, fileb, zc( 1 ) )
   irowb1 = iurow1
   irowbn = iurown
   GO TO 940
   930   irowb1 = 0
   irowbn = 0
   940   CONTINUE
   IF ( ifile == 0 ) GO TO 950
   iurow1 = 1
   iurown = ndr
   typeu  = ndtype
   IF ( ipass == 1 ) typeu = ndtype * signc
   CALL unpack (*950, ifile, zc( idx2+1 ) )
   GO TO 980
   950   CONTINUE
   DO  j = 1, ndr
     zc( idx2+j ) = (0.0,0.0)
   END DO
   980   CONTINUE
   
! CHECK IF COLUMN OF "B" IS NULL
   
   indxa  = iax
   IF ( irowb1 == 0 ) GO TO 50000
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A * B  +  C )
   
! COMPLEX SINGLE PRECISION
   1000  CONTINUE
   DO  i = 1, ncolpp
     indxal = zi( indxa+1 ) + iax - 1
     icola  = ibrow+i
     IF ( icola < irowb1 .OR. icola > irowbn ) GO TO 1450
     ibrowi = icola - irowb1 + 1
     IF ( zc( ibrowi ) == 0. ) GO TO 1450
     IF ( icola /= IABS( zi( indxa ) ) ) GO TO 70001
     indxa  = indxa+2
     1100  CONTINUE
     IF ( indxa >= indxal ) GO TO 1450
     irowa1 = zi( indxa )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     indxav = ( ( indxa+3 ) / 2 ) - irowa1
     DO  k = irowa1, irowan
       
!         D = C + A*B
       
       zc( idx2+k ) = zc( idx2+k ) +  zc( indxav+k ) * zc( ibrowi )
     END DO
     indxa  = indxa + 2 + ntms*2
     GO TO 1100
     1450  indxa  = indxal
   END DO
   GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
! COMPLEX SINGLE PRECISION
   10000 CONTINUE
   indxb  = 1 - irowb1
   idxx   = idx2 + ibrow
   DO  i = 1, ncolpp
     icola  = ibrow + i
     IF ( icola /= IABS( zi( indxa ) ) ) GO TO 70001
     indxal = zi( indxa+1 ) + iax - 1
     indxa  = indxa+2
     11000 CONTINUE
     IF ( indxa >= indxal ) GO TO 14500
     irowa1 = zi( indxa )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) GO TO 14100
     indxav = ( ( indxa + 3 ) / 2 ) - irowa1
     
!         D = C + A*B
     
     DO  k = irow1, irown
       zc( idxx+i ) = zc( idxx+i ) +  zc( indxav+k ) * zc( indxb+k )
     END DO
     14100 CONTINUE
     indxa  = indxa + 2 + ntms*2
     GO TO 11000
     14500 indxa  = indxal
   END DO
   GO TO 50000
! END OF PROCESSING THIS COLUMN FOR THIS PASS
   50000 CONTINUE
!  NOW SAVE COLUMN
   CALL pack ( zc( idx2+1 ), ofile, filed )
 END DO
 GO TO 70000
 70001 CONTINUE
 WRITE ( iwr, 90001 ) icola, zi( indxa ), iax, indxa
 90001 FORMAT(' UNEXPECTED COLUMN FOUND IN PROCESSING MATRIX A'  &
     ,/,' COLUMN EXPECTED:',i6 ,/,' COLUMN FOUND   :',i6  &
     ,/,' IAX =',i7,'  INDXA=',i7 )
 CALL mesage ( -61, 0, 0 )
 70000 RETURN
END SUBROUTINE mma303

