SUBROUTINE mma114 ( zi, zd, zdc )
     
!     MMA114 PERFORMS THE MATRIX OPERATION IN COMPLEX DOUBLE PRECISION
!       (+/-)A(T & NT) * B (+/-)C = D
 
!     MMA114 USES METHOD 10 WHICH IS AS FOLLOWS:
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
 DOUBLE PRECISION, INTENT(IN)             :: zd(2)
 REAL, INTENT(OUT)                        :: zdc
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
   
   CALL mmarc4 ( zi, zd )
   
! NOW READ "C", OR SCRATCH FILE WITH INTERMEDIATE RESULTS.
! IF NO "C" FILE AND THIS IS THE FIRST PASS, INITIALIZE "D" COLUMN TO ZERO
   
   IF ( ifile == 0 ) GO TO 950
   iurow1 = 1
   iurown = ndr
   typeu  = ndtype
   IF ( ipass == 1 ) typeu = ndtype * signc
   CALL unpack (*950, ifile, zdc( idx4+1 ) )
   GO TO 980
   950   CONTINUE
   DO  j = 1, ndr
     zdc( idx4+j ) = (0.0,0.0)
   END DO
   980   CONTINUE
   nwddnar = nwdd*nar
   
! CHECK IF COLUMN OF "B" IS NULL
   
   irowb1  = zi( 1 )
   irows   = zi( 2 )
   irowbn  = irowb1 + irows - 1
   indx    = 1
   
! CHECK FOR NULL COLUMN FROM THE "B" MATRIX
   
   IF ( irowb1 == 0 ) GO TO 50000
   IF ( t /= 0 ) GO TO 5000
   
! "A" NON-TRANSPOSE CASE    ( A * B  +  C )
   
! COMLEX DOUBLE PRECISION
   4000  CONTINUE
   DO  i = 1, ncolpp
     ibrowi = ibrow+i
     indxa  = iax + 2*i + ( i-1 )*nwddnar
     irowa1 = zi( indxa-2 )
     IF ( irowa1 == 0 ) CYCLE
     irowan = zi( indxa-1 )
     indxa1 = ( ( indxa+1 ) / 2 ) - 2
     4100  CONTINUE
     IF ( ibrowi < irowb1 ) CYCLE
     IF ( ibrowi <= irowbn ) GO TO 4200
     indx   = indx + 2 + irows*nwdd
     IF ( indx >= lasind ) GO TO 50000
     irowb1 = zi( indx )
     irows  = zi( indx+1 )
     irowbn = irowb1 + irows - 1
     GO TO 4100
     4200  CONTINUE
     indxv  =  2*( ibrowi - irowb1 ) + ( ( indx+3 ) / 2 )
     IF (   zd( indxv  ) == 0.d0 .AND. zd( indxv+1) == 0.d0 ) CYCLE
     indxav = indxa1
     DO  k = irowa1, irowan
       indxav = indxav + 2
       zdc( idx4+k ) = zdc( idx4+k ) +  &
           DCMPLX( zd(indxav ), zd(indxav +1) ) *  &
           DCMPLX( zd(indxv  ), zd(indxv+1  ) )
     END DO
   END DO
   GO TO 50000
   
!  TRANSPOSE CASE ( A(T) * B + C )
   
   5000  CONTINUE
   idrow = ibrow
! COMLEX DOUBLE PRECISION
   40000 CONTINUE
   DO  i = 1, ncolpp
     indx   = 1
     indxa1 = iax + 2*i + ( i-1 )*nwddnar
     irowa1 = zi( indxa1-2 )
     IF ( irowa1 == 0 ) CYCLE
     irowan = zi( indxa1-1 )
     41000 IF ( indx >= lasind ) CYCLE
     irowb1 = zi( indx )
     irows  = zi( indx+1 )
     irowbn = irowb1 + irows - 1
     indxv  = ( indx+3 ) / 2
     indx   = indx + 2 + irows*nwdd
     irow1  = MAX0( irowa1, irowb1 )
     irown  = MIN0( irowan, irowbn )
     IF ( irown < irow1 ) GO TO 41000
     indxa  = ( ( indxa1+1 ) / 2 ) + 2*( irow1 - irowa1 ) - 1
     idx4x  = idx4 + idrow
     indxb  = indxv + 2*( irow1 - irowb1 ) - 1
     kcnt   = ( irown-irow1 ) * 2  + 1
     DO  k = 1, kcnt, 2
       zdc( idx4x+i ) = zdc( idx4x+i ) +  &
           DCMPLX( zd(indxa+k ), zd(indxa+k+1) ) *  &
           DCMPLX( zd(indxb+k ), zd(indxb+k+1) )
     END DO
     GO TO 41000
   END DO
   GO TO 50000
! END OF PROCESSING THIS COLUMN FOR THIS PASS
   50000 CONTINUE
!  NOW SAVE COLUMN
   CALL pack ( zdc( idx4+1 ), ofile, filed )
 END DO
 RETURN
END SUBROUTINE mma114

