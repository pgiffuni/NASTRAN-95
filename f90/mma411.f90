SUBROUTINE mma411 ( zi, zr )
     
!     MMA411 PERFORMS THE MATRIX OPERATION USING METHOD 41
!       IN REAL SINGLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA411 USES METHOD 41 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  READ AS MANY COLUMNS OF "B" INTO MEMORY AS POSSIBLE
!           INTO MEMORY IN COMPACT FORM LEAVING SPACE FOR A FULL
!           COLUMN OF "D" FOR EVERY COLUMN "B" READ.  SEE SUBROUTINES
!           MMARM1,2,3,4 FOR FORMAT OF COMPACT FORM.
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  CALL UNPACK TO READ MATRICES "C".
!       5.  CALL MMARC1,2,3,4 TO READ COLUMNS OF MATRIX "A".
 
 
 INTEGER, INTENT(IN)                      :: zi(2)
 REAL, INTENT(IN OUT)                     :: zr(2)
 INTEGER :: t
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: filea      ,fileb ,filec , filed
 
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
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "A" MATRIX
!     Z( IBX      ) = ARRAY FOR MULTIPLE COLUMNS OF "B" MATRIX
!                     (STORED IN COMPACT FORM)
!     Z( IDX      ) = ARRAY FOR MULTIPLE COLUMNS OF "D" MATRIX
!                     (FULL COLUMN SPACE ALLOCATION)
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE
 
 
! PROCESS ALL OF THE COLUMNS OF "A"
 
 irfile = filea( 1 )
 SIGN   = 1
 DO  ii = 1, nac
!      print *,' processing column of a, ii=',ii
   
! READ A COLUMN FROM THE "A" MATRIX
   
   CALL mmarc1 ( zi, zr )
   
! CHECK FOR NULL COLUMN FROM THE "A" MATRIX
   
   IF ( zi ( 1 ) == 0 ) CYCLE
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A*B + C )
   
! SINGLE PRECISION
   1000  CONTINUE
   indxb  = ibx
   DO  i = 1, ncolpp
     icolb  = ibrow + i
     IF ( icolb /= IABS( zi( indxb ) ) ) GO TO 70001
     indxbl = zi( indxb+1 ) + ibx - 1
     indxb  = indxb + 2
     indxd  = ( idx + ( i-1 )*ndr ) - 1
     1100  CONTINUE
     IF ( indxb >= indxbl ) GO TO 1450
     irowb1 = zi( indxb   )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     IF ( ii > irowbn ) GO TO 1410
     IF ( ii < irowb1 ) GO TO 1450
     indxbv = indxb + 2 + ii - irowb1
     IF ( zr( indxbv ) == 0.0 ) GO TO 1450
     indxa  = 1
     1200  CONTINUE
     irowa1 = zi( indxa   )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     indxav = indxa + 2 - irowa1
     DO  k = irowa1, irowan
       zr( indxd+k ) = zr( indxd+k ) + zr( indxav+k ) * zr( indxbv )
     END DO
     indxa  = indxa + 2 + ntms
     IF ( indxa < lasind ) GO TO 1200
     GO TO 1450
     1410  CONTINUE
     indxb  = indxb + 2 + irows
     GO TO 1100
     1450  indxb  = indxbl
   END DO
   CYCLE
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
! SINGLE PRECISION
   10000 CONTINUE
   indxb  = ibx
   DO  i = 1, ncolpp
     indxa  = 1
     icolb  = ibrow + i
     IF ( icolb /= IABS( zi( indxb ) ) ) GO TO 70001
     indxbl = zi( indxb+1 ) + ibx - 1
     indxb  = indxb + 2
     indxd  = ( idx + ( i-1 )*ndr ) + ii - 1
     11000 CONTINUE
     IF ( indxb >= indxbl ) GO TO 14500
     irowb1 = zi( indxb )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     indxbv = indxb + 2 - irowb1
     indxb  = indxb + 2 + irows
     12000 CONTINUE
     irowa1 = zi( indxa   )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     IF ( irowbn < irowa1 ) GO TO 11000
     IF ( irowan < irowb1 ) GO TO 14200
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) GO TO 11000
     indxav = indxa + 2 - irowa1
     DO  k = irow1, irown
       zr( indxd ) = zr( indxd ) + zr( indxav+k ) * zr( indxbv+k )
     END DO
     IF ( irowan > irowbn ) GO TO 11000
     14200 CONTINUE
     indxa  = indxa + 2 + ntms
     IF ( indxa >= lasind ) GO TO 14500
     GO TO 12000
     14500 CONTINUE
     indxb  = indxbl
   END DO
! END OF PROCESSING THIS COLUMN OF "A" FOR THIS PASS
 END DO
!  NOW SAVE COLUMNS COMPLETED
 DO  k = 1, ncolpp
   indx = idx + ( k-1 ) * ndr
   CALL pack ( zr( indx ), filed, filed )
 END DO
 GO TO 70000
 70001 WRITE ( iwr, 90001 ) icolb, zi( indxb ), ibx, indxb
 90001 FORMAT(' UNEXPECTED COLUMN FOUND IN PROCESSING MATRIX B'  &
     ,/,' COLUMN EXPECTED:',i6 ,/,' COLUMN FOUND   :',i6  &
     ,/,' IBX =',i7,'  INDXB =',i7 )
 CALL mesage ( -61, 0, 0 )
 70000 RETURN
END SUBROUTINE mma411

