SUBROUTINE mma414 ( zi, zd, zdc )
     
!     MMA414 PERFORMS THE MATRIX OPERATION USING METHOD 41
!       IN COMPLEX DOUBLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA414 USES METHOD 41 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  READ AS MANY COLUMNS OF "B" INTO MEMORY AS POSSIBLE
!           INTO MEMORY IN COMPACT FORM LEAVING SPACE FOR A FULL
!           COLUMN OF "D" FOR EVERY COLUMN "B" READ.  SEE SUBROUTINES
!           MMARM1,2,3,4 FOR FORMAT OF COMPACT FORM.
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  CALL UNPACK TO READ MATRICES "C".
!       5.  CALL MMARC1,2,3,4 TO READ COLUMNS OF MATRIX "A".
 
 
 INTEGER, INTENT(IN)                      :: zi(2)
 DOUBLE PRECISION, INTENT(IN)             :: zd(2)
 REAL, INTENT(OUT)                        :: zdc
 INTEGER :: t
 INTEGER :: typei      ,typep    ,typeu ,signab, signc
 INTEGER :: rd         ,rdrew    ,wrt   ,wrtrew, clsrew,cls
 INTEGER :: filea      ,fileb ,filec , filed
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
   
! READ A COLUMN FROM THE "A" MATRIX
   
   CALL mmarc4 ( zi, zd )
   
! CHECK IF COLUMN FROM "A" MATRIX IS NULL
   
   IF ( zi ( 1 ) == 0 ) CYCLE
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A*B + C )
   
! COMPLEX DOUBLE PRECISION
   1000  CONTINUE
   indxb  = ibx
   DO  i = 1, ncolpp
     icolb  = ibrow + i
     IF ( icolb /= IABS( zi( indxb ) ) ) GO TO 70001
     indxbl = zi( indxb+1 ) + ibx - 1
     indxb  = indxb + 2
     indxd  = ( idx4 + ( i-1 )*ndr )
     1100  CONTINUE
     IF ( indxb >= indxbl ) GO TO 1450
     irowb1 = zi( indxb   )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     IF ( ii > irowbn ) GO TO 1410
     IF ( ii < irowb1 ) GO TO 1450
     indxbv = ( ( indxb + 3 ) / 2 ) + 2*( ii - irowb1 )
     IF ( zd( indxbv   ) == 0.0D0 .AND. zd( indxbv+1 ) == 0.0D0 ) GO TO 1450
     indxa  = 1
     1200  CONTINUE
     irowa1 = zi( indxa   )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     indxav = ( ( indxa + 3 ) / 2 )
     DO  k = irowa1, irowan
       zdc( indxd+k ) = zdc( indxd+k ) +  &
           DCMPLX( zd( indxav ), zd( indxav+1 ) ) *  &
           DCMPLX( zd( indxbv ), zd( indxbv+1 ) )
       indxav = indxav + 2
     END DO
     indxa  = indxa + 2 + ntms*nwdd
     IF ( indxa < lasind ) GO TO 1200
     GO TO 1450
     1410  CONTINUE
     indxb  = indxb + 2 + irows*nwdd
     GO TO 1100
     1450  indxb  = indxbl
   END DO
   CYCLE
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
! COMPLEX DOUBLE PRECISION
   10000 CONTINUE
   indxb  = ibx
   DO  i = 1, ncolpp
     indxa  = 1
     icolb  = ibrow + i
     IF ( icolb /= IABS( zi( indxb ) ) ) GO TO 70001
     indxbl = zi( indxb+1 ) + ibx - 1
     indxb  = indxb + 2
     indxd  = ( idx4 + ( i-1 )*ndr ) + ii
     11000 CONTINUE
     IF ( indxb >= indxbl ) GO TO 14500
     irowb1 = zi( indxb )
     irows  = zi( indxb+1 )
     irowbn = irowb1 + irows - 1
     indxbs = indxb
     indxb  = indxb + 2 + irows*nwdd
     12000 CONTINUE
     irowa1 = zi( indxa   )
     ntms   = zi( indxa+1 )
     irowan = irowa1 + ntms - 1
     IF ( irowbn < irowa1 ) GO TO 11000
     IF ( irowan < irowb1 ) GO TO 14200
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) GO TO 11000
     indxbv = ( ( indxbs+ 3 ) / 2 ) + 2*( irow1 - irowb1 ) - 1
     indxav = ( ( indxa + 3 ) / 2 ) + 2*( irow1 - irowa1 ) - 1
     icnt   = 2*( irown - irow1 ) + 1
     DO  k = 1, icnt, 2
       zdc( indxd ) = zdc( indxd ) +  &
           DCMPLX( zd( indxav+k ), zd( indxav+k+1 ) ) *  &
           DCMPLX( zd( indxbv+k ), zd( indxbv+k+1 ) )
     END DO
     IF ( irowan > irowbn ) GO TO 11000
     14200 CONTINUE
     indxa  = indxa + 2 + ntms*nwdd
     IF ( indxa >= lasind ) GO TO 14500
     GO TO 12000
     14500 CONTINUE
     indxb  = indxbl
   END DO
! END OF PROCESSING THIS COLUMN OF "A" FOR THIS PASS
 END DO
!  NOW SAVE COLUMNS COMPLETED
 DO  k = 1, ncolpp
   indx = idx4 + ( k-1 ) * ndr + 1
   CALL pack ( zdc( indx ), filed, filed )
 END DO
 GO TO 70000
 70001 WRITE ( iwr, 90001 ) icolb, zi( indxb ), ibx, indxb
 90001 FORMAT(' UNEXPECTED COLUMN FOUND IN PROCESSING MATRIX B'  &
     ,/,' COLUMN EXPECTED:',i6 ,/,' COLUMN FOUND   :',i6  &
     ,/,' IBX =',i7,'  INDXB =',i7 )
 CALL mesage ( -61, 0, 0 )
 70000 RETURN
END SUBROUTINE mma414

