SUBROUTINE mma111 ( zi, zr )
     
!     MMA111 PERFORMS THE MATRIX OPERATION IN REAL SINGLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA111 USES METHOD 11 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  UNPACK AS MANY COLUMNS OF "A" INTO MEMORY AS POSSIBLE
!           LEAVING SPACE FOR ONE COLUMN OF "B" AND "D".
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  UNPACK COLUMNS OF "C" MATRIX BUT USE GETSTR (MMARC1,2,3,4)
!           TO READ COLUMNS OF "B".
 
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
 
 DO  ii = 1, nbc
   
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
   nwddnar = nwdd*nar
   
! CHECK IF COLUMN OF "B" IS NULL
   
   irowb1 = zi( 1 )
   irows  = zi( 2 )
   irowbn = irowb1 + irows - 1
   indx   = 1
   
! CHECK FOR NULL COLUMN ON "B"
   
   IF ( irowb1 == 0 ) GO TO 50000
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A * B  +  C )
   
! SINGLE PRECISION
   1000  CONTINUE
   DO  i = 1, ncolpp
     ibrowi = ibrow+i
     indxa  = iax + 2*i + ( i-1 )*nar
     irowa1 = zi( indxa-2 )
     IF ( irowa1 == 0 ) CYCLE
     irowan = zi( indxa-1 )
     indxa  = indxa - irowa1
     1100  CONTINUE
     IF ( ibrowi < irowb1 ) CYCLE
     IF ( ibrowi <= irowbn ) GO TO 1200
     indx   = indx + 2 + irows
     IF ( indx > lasind ) GO TO 50000
     irowb1 = zi( indx )
     irows  = zi( indx+1 )
     irowbn = irowb1 + irows - 1
     GO TO 1100
     1200  CONTINUE
     indxv  = ibrowi - irowb1 + indx + 2
     IF ( zr( indxv ) == 0. ) CYCLE
     DO  k = irowa1, irowan
       zr( idx+k-1 ) = zr( idx+k-1 ) +  zr( indxa+k ) * zr( indxv )
     END DO
   END DO
   GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
   idrow = ibrow
! SINGLE PRECISION
   10000 CONTINUE
   DO  i = 1, ncolpp
     indx   = 1
     indxa  = iax + 2*i + ( i-1 )*nar
     irowa1 = zi( indxa-2 )
     IF ( irowa1 == 0 ) CYCLE
     irowan = zi( indxa-1 )
     indxav = indxa - irowa1
     11000 IF ( indx >= lasind ) CYCLE
     irowb1 = zi( indx )
     irows  = zi( indx+1 )
     irowbn = irowb1 + irows - 1
     indxv  = indx+2
     indx   = indx + 2 + irows
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) GO TO 11000
     indxb  = indxv - irowb1
     idxx   = idx + idrow - 1
     DO  k = irow1, irown
       zr( idxx+i ) = zr( idxx+i ) +  zr( indxav+k ) * zr( indxb+k )
     END DO
     GO TO 11000
   END DO
   GO TO 50000
! END OF PROCESSING THIS COLUMN FOR THIS PASS
   50000 CONTINUE
!  NOW SAVE COLUMN
   CALL pack ( zr( idx ), ofile, filed )
 END DO
 RETURN
END SUBROUTINE mma111
