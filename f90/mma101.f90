SUBROUTINE mma101 ( zi, zr )
     
!     MMA101 PERFORMS THE MATRIX OPERATION USING METHOD 10 AND
!       REAL SINGLE PRECISION
 
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA10 USES METHOD 10 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  UNPACK AS MANY COLUMNS OF "A" INTO MEMORY AS POSSIBLE
!           LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  USE UNPACK TO READ MATRICES "B" AND "C".
 
!     MEMORY FOR EACH COLUMN OF "A" IS AS FOLLOWS:
!         Z(1)   = FIRST NON-ZERO ROW NUMBER FOR COLUMN
!         Z(2)   = LAST NON-ZERO ROW NUMBER FOR COLUMN
!         Z(3-N) = VALUES OF NON-ZERO ROWS
 
 
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
 EQUIVALENCE       (ksystm( 1),sysbuf)  , (ksystm( 2),nout  )
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
 
! PROCESS ALL OF THE COLUMNS OF "B", ADD "C" DATA ON FIRST PASS
 DO  ii = 1, nbc
   iurow1 = -1
   typeu  = ndtype * signab
   CALL unpack ( *930, fileb, zr( 1 ) )
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
   CALL unpack (*950, ifile, zr( idx ) )
   GO TO 980
   950   CONTINUE
   DO  j = 1, ndr
     zr( idx+j-1 ) = 0
   END DO
   980   CONTINUE
   nwddnar = nwdd*nar
   
! CHECK IF COLUMN OF "B" IS NULL
   
   IF ( irowb1 == 0 ) GO TO 50000
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A * B  +  C )
   
! SINGLE PRECISION
   1000  CONTINUE
   DO  i = 1, ncolpp
     ibrowi = ibrow+i
     IF ( ibrowi < irowb1 .OR. ibrowi > irowbn ) CYCLE
     ibrowi = ibrowi - irowb1 + 1
     IF ( zr( ibrowi ) == 0. ) CYCLE
     indxa  = iax + 2*i + ( i-1 )*nar
     irowa1 = zi( indxa-2 )
     IF ( irowa1 == 0 ) CYCLE
     irowan = zi( indxa-1 )
     indxa  = indxa - irowa1
     DO  k = irowa1, irowan
       zr( idx+k-1 ) = zr( idx+k-1 ) +  zr( indxa+k ) * zr( ibrowi )
     END DO
   END DO
   GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
   idrow = ibrow
! SINGLE PRECISION
   10000 CONTINUE
   DO  i = 1, ncolpp
     indxa  = iax + 2*i + ( i-1 )*nar
     irowa1 = zi( indxa-2 )
     IF ( irowa1 == 0 ) CYCLE
     irowan = zi( indxa-1 )
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) CYCLE
     indxa  = indxa - irowa1
     indxb  = 1 - irowb1
     idxx   = idx + idrow - 1
     DO  k = irow1, irown
       zr( idxx+i ) = zr( idxx+i ) +  zr( indxa+k ) * zr( indxb+k )
     END DO
   END DO
   GO TO 50000
! END OF PROCESSING THIS COLUMN FOR THIS PASS
   50000 CONTINUE
!  NOW SAVE COLUMN
   CALL pack ( zr( idx ), ofile, filed )
 END DO
 RETURN
END SUBROUTINE mma101

