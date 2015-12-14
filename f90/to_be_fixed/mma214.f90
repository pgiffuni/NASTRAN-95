SUBROUTINE mma214 ( zi, zd, zdc )
     
!     MMA214 PERFORMS THE MATRIX OPERATION USING METHOD 21
!       IN COMPLEX DOUBLE PRECISION
 
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA214 USES METHOD 21 WHICH IS AS FOLLOWS:
!       1.  THIS IS FOR "A" NON-TRANSPOSED AND TRANSPOSED
!       2.  UNPACK AS MANY COLUMNS OF "B" INTO MEMORY AS POSSIBLE
!           LEAVING SPACE FOR A COLUMN OF "D" FOR EVERY COLUMN "B" READ.
!       3.  INITIALIZE EACH COLUMN OF "D" WITH THE DATA FROM "C".
!       4.  USE UNPACK TO READ MATRICES "B" AND "C".
 
!     MEMORY FOR EACH COLUMN OF "B" IS AS FOLLOWS:
!         Z(1)   = FIRST NON-ZERO ROW NUMBER FOR COLUMN
!         Z(2)   = LAST NON-ZERO ROW NUMBER FOR COLUMN
!         Z(3-N) = VALUES OF NON-ZERO ROWS
 
 
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
!     Z( 1        ) = ARRAY FOR ONE COLUMN OF "A" MATRIX
!     Z( IDX      ) = ARRAY FOR MULTIPLE COLUMNS OF "D" MATRIX
!     Z( IBX      ) = ARRAY FOR MULTIPLE COLUMNS OF "B" MATRIX
!        THROUGH
!     Z( LASMEM   )
!     Z( IBUF4    ) = BUFFER FOR "D" FILE
!     Z( IBUF3    ) = BUFFER FOR "C" FILE
!     Z( IBUF2    ) = BUFFER FOR "B" FILE
!     Z( IBUF1    ) = BUFFER FOR "A" FILE
!     Z( NZ       ) = END OF OPEN CORE THAT IS AVAILABLE
 
 
! PROCESS ALL OF THE COLUMNS OF "A"
 
 DO  ii = 1, nac
   
! READ A COLUMN FROM THE "A" MATRIX
   
   CALL mmarc4 ( zi, zd )
   
! CHECK IF COLUMN FROM "A" IS NULL
   
   IF ( zi( 1 ) == 0 ) CYCLE
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A*B + C )
   
! COMLEX DOUBLE PRECISION
   4000  CONTINUE
   DO  i = 1, ncolpp
     indx   = 1
     irowa1 = zi ( indx )
     indxb  = ibx + 2*i + ( i-1 )*nwddnbr
     irowb1 = zi( indxb-2 )
     irowbn = zi( indxb-1 )
     IF ( ii < irowb1 .OR. ii > irowbn ) CYCLE
     indxb  = ( ( indxb+1 ) / 2 ) + 2*( ii - irowb1 )
     IF ( zd( indxb ) == 0.0D0 .AND. zd( indxb+1 ) == 0.0D0 ) CYCLE
     indxd  = idx4 + ( i-1 )*ndr
     4100  CONTINUE
     irows  = zi( indx+1 )
     irowan = irowa1 + irows - 1
     indxa  = (indx+3)/2
     DO  k = irowa1, irowan
       zdc( indxd+k ) = zdc( indxd+k ) +  &
           DCMPLX( zd( indxa ), zd( indxa+1 ) ) *  &
           DCMPLX( zd( indxb ), zd( indxb+1 ) )
       indxa  = indxa + 2
     END DO
     indx   = indx + 2 + irows*nwdd
     IF ( indx >= lasind ) CYCLE
     irowa1 = zi( indx )
     GO TO 4100
   END DO
   CYCLE
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
! COMLEX DOUBLE PRECISION
   40000 CONTINUE
   DO  i = 1, ncolpp
     indxb  = ibx + 2*i + ( i-1 )*nwddnbr
     irowb1 = zi( indxb-2 )
     IF ( irowb1 == 0 ) CYCLE
     irowbn = zi( indxb-1 )
     indx   = 1
     irowa1 = zi( indx )
     indxd  = idx4 + ( i-1 )*ndr + ii
     41000 CONTINUE
     irows  = zi( indx+1 )
     irowan = irowa1 + irows - 1
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) GO TO 44100
     indxbv = ( ( ( indxb+1 ) / 2 ) + 2*( irow1 - irowb1 ) ) - 1
     indxa  = ( ( ( indx+3  ) / 2 ) + 2*( irow1 - irowa1 ) ) - 1
     kcnt   = ( irown - irow1 ) * 2 + 1
     DO  k = 1, kcnt, 2
       zdc( indxd ) = zdc( indxd ) +  &
           DCMPLX( zd( indxa +k ), zd( indxa +k+1 ) ) *  &
           DCMPLX( zd( indxbv+k ), zd( indxbv+k+1 ) )
     END DO
     44100 CONTINUE
     indx   = indx + 2 + irows*nwdd
     IF ( indx >= lasind ) CYCLE
     irowa1 = zi( indx )
     GO TO 41000
   END DO
! END OF PROCESSING THIS COLUMN OF "A" FOR THIS PASS
 END DO
!  NOW SAVE COLUMNS COMPLETED
 DO  k = 1, ncolpp
   indx = idx4 + ( k-1 ) * ndr
   CALL pack ( zdc( indx+1 ), filed, filed )
 END DO
 RETURN
END SUBROUTINE mma214

