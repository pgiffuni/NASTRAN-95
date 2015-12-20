SUBROUTINE smc2rs ( zi, zs, zil, zol, nar, lasrow, rtemp,  &
        i1, i2, i3  )
     
! ZIL    = INNER LOOP TERMS (SIZE = MAXNAC * (MAXNCOL+NEXTRA)
! ZOL    = OUTER LOOP TERMS (SIZE = (MAXNCOL+NEXTRA) * 2)
! NAR    = SAVE AREA FOR ACTIVE ROWS OF PREVIOUS COLUMN
! I1     = MAXIMUM NUMBER OF ACTIVE ROWS FOR THIS COLUMN
! I2     = NUMBER OF COLUMNS ALLOCATED FOR STORAGE OF INNER AND
!          NUMBER OF ROWS ALLOCATED FOR OUTER LOOP
! I3     = MAXIMUM NUMBER OF WORDS FOR DEFINING THE ACTIVE ROWS FOR
!          ANY COLUMN
! LASROW = LAST NON-ZERO ROW INDEX FOR A GIVEN COLUMN (SIZE = MAXNCOL
!          +NEXTRA)
 
 
 INTEGER, INTENT(IN)                      :: zi(10)
 REAL, INTENT(IN OUT)                     :: zs(10)
 REAL, INTENT(OUT)                        :: zil( i1, i2 )
 REAL, INTENT(OUT)                        :: zol( i2, 2 )
 INTEGER, INTENT(IN OUT)                  :: nar( i3 )
 INTEGER, INTENT(OUT)                     :: lasrow(i2)
 REAL, INTENT(OUT)                        :: rtemp(i3)
 INTEGER, INTENT(IN)                      :: i1
 INTEGER, INTENT(IN)                      :: i2
 INTEGER, INTENT(IN OUT)                  :: i3
 
 REAL :: zoltmp
 
 
 INCLUDE           'SMCOMX.COM'
 
! GET ROW VALUES CORRESPONDING TO THE ACTIVE ROWS OF COLUMN K FOR
! EACH COLUMN KFRCOL THROUGH KLSCOL IN ORDER TO FILL INNER LOOP AND
! OUTER LOOP AREAS.
 
 
! BEGIN TO PROCESS EACH COLUMN
! FOR COLUMN K, GET OUTER LOOP TERMS
!    A(K,J) / A(J,J)
!       K = CURRENT PIVOTAL COLUMN
!       J = RANGES FROM FIRST COLUMN DATA NEEDED FOR COLUMN K TO K-1
!    (E.G.,
!         A(5,1)/A(1,1)
!         A(5,2)/A(2,2)
!         A(5,3)/A(3,3)
!         A(5,4)/A(4,4)
! ALSO, GET INNER LOOP TERMS
!    A(I,J)
!       K = CURRENT PIVOTAL COLUMN
!       I = RANGES FROM K TO LAST ACTIVE ROW OF COLUMN K
!       J = RANGES FROM FIRST COLUMN DATA NEEDED FOR COLUMN K TO K-1
!    (E.G.,
!         A(5,1) A(6,1)  .  A(N,1)
!         A(5,2) A(6,2)  .  A(N,2)
!         A(5,3) A(6,3)  .  A(N,3)
!         A(5,4) A(6,4)  .  A(N,4)
 
!  LOOP 7000 WILL BE ON K
!  LOOP 6000 WILL BE ON J
 
 ic1     = 1
 ic2     = 2
 iilrow1 = 1
!      print *,' i1,i2,i3,maxncol,maxnac=',i1,i2,i3,maxncol,maxnac
 DO  k = 1, ncol
   kk      = MOD( k, i2 )
   IF ( kk == 0 ) kk = i2
   lasrow( kk ) = 0
!      PRINT *,' SMC2RS PROCESSING COLUMN K=',K
   kcol    = k
   kdir    = k*4 - 3
   kmidx   = zi( kdir   )
   
! SEE IF DATA IS ON IN MEMORY OR ON THE SPILL FILE
   
   IF ( kmidx /= 0 ) GO TO 500
   
! DATA IS ON THE SPILL FILE
   
   CALL smcspl ( kcol, zi )
   kmidx  = zi( kdir )
   500   CONTINUE
   kfrcolp= kfrcol
   klscolp= klscol
   kfrcol = zi( kdir+1 )
   km2    = zi( kmidx+1)
   kridxn = kmidx + 4 + km2
   klscol = k - 1
   kridx  = kmidx+4
   kridxs = kridx
   krow1  = zi( kridx   )
   krown  = krow1 + zi( kridx+1 ) - 1
   karows = 0
   DO  kk = 1, km2, 2
     karows = karows + zi( kridx+kk )
   END DO
   
! IF THE PREVIOUS COLUMN DID NOT NEED DATA FROM A COLUMN PRECEEDING IT,
! THEN MUST RELOAD THE INNER AND OUTER LOOP ARRAYS
   
   IF ( klscolp < kfrcolp ) GO TO 1350
   
! NOW MUST FIND THE ROW AND COLUMN NUMBER FOR THIS PIVOT COLUMN
! THAT IS NOT ALREADY IN THE INNER LOOP AND OUTER LOOP ARRAYS.
! FIRST CHECK THAT THE FIRST REQUIRED ROW IS STORED, IF NOT THEN WE MUST
! BEGIN AS IF NOTHING STORED.  IF SOME OF THE REQUIRED ROWS ARE PRESENT,
! THEN FIND THE NEXT POSITION AND ROW NUMBER TO BE STORED IN THE INNER
! LOOP ARRAY AND THE NEXT POSITION AND COLUMN NUMBER TO BE STORED IN THE
! OUTER LOOP ARRAY.
   
! IF THE FIRST COLUMN IS LESS THAN FIRST COLUMN OF LAST PIVOT COLUMN
! THEN WE MUST LOAD THE INNER AND OUTER LOOPS FROM THE BEGINNING
   
   IF ( kfrcol < kfrcolp ) GO TO 1350
   kr      = 1
   lrow1   = nar( 1 )
   lrown   = nar( 1 ) + nar( 2 ) - 1
   
!  LROW1 = FIRST ROW OF A STRING OF CONTIGUOUS ROWS OF LAST PIVOT
!          COLUMN PROCESSED
!  LROWN = LAST ROW OF A STRING OF CONTIGUOUS ROWS OF LAST PIVOT COLUMN
!          PROCESSED
   
! FIND FIRST ROW IN INNER LOOP THAT MATCHES THE FIRST ROW REQUIRED
! FOR THIS COLUMN
   
! IF THERE IS NO MATCH FOR THE FIRST COLUMN, THEN GO TO 1350
   
   1105  CONTINUE
   IF ( lrow1 > krow1 ) GO TO 1350
   IF ( krow1 < lrown ) GO TO 1100
   
! NO OVERLAP WITH THIS STRING, GO AND GET NEXT STRING
! ADJUST 'ILLROW1' WHICH IS THE POINTER TO THE FIRST ROW IN THE INNER
! LOOP THAT CONTAINS THE VALUE OF ROW "KROW1" OF EACH COLUMN.
   
   incr    = lrown - lrow1 + 1
   iilrow1 = iilrow1 + incr
   IF ( iilrow1 > i1 ) iilrow1 = iilrow1 - i1
   kr      = kr + 2
   lrow1   = nar( kr )
   IF ( lrow1 == 0 ) GO TO 1350
   lrown   = lrow1 + nar( kr+1 ) - 1
   GO TO 1105
   1100  CONTINUE
   
! THERE IS AN OVERLAP, SET KROWB, KROWSB, AND IILROW1 TO REFLECT
! THE PROPER ROW NUMBER IN THE INNER LOOP
   
   incr    = krow1 - lrow1
   krowb   = krow1
   krowsb  = krown - krowb + 1
   kridxs  = kridx
   iilrow1 = iilrow1 + incr
   IF ( iilrow1 > i1    ) iilrow1 = iilrow1 - i1
   lrow1   = krow1
   iilrow  = iilrow1
   1120  IF ( lrow1 /= krow1 ) GO TO 1180
   IF ( lrown == krown ) GO TO 1130
   IF ( lrown < krown ) GO TO 1140
   IF ( lrown > krown ) GO TO 1150
   
! THIS SET OF ROWS MATCHES, GO AND CHECK THE NEXT SET OF ROW NUMBERS
   
   1130  CONTINUE
   incr    = krown - krowb + 1
   iilrow  = iilrow + incr
   IF ( iilrow > i1 ) iilrow = iilrow - i1
   kridx   = kridx + 2
   IF ( kridx == kridxn ) GO TO 1170
   kr      = kr + 2
   krow1   = zi( kridx )
   krowb   = krow1
   krowsb  = zi( kridx+1 )
   krown   = krow1 + krowsb -1
   kridxs  = kridx
   lrow1   = nar( kr )
   lrown   = lrow1 + nar( kr+1 ) - 1
   IF ( lrow1 == 0 ) GO TO 1180
   GO TO 1120
   
! LAST ROW NUMBERS DO NOT MATCH, KROWN GT LROWN
   
   1140  CONTINUE
   incr    = lrown  - krowb + 1
   1145  krowb   = krowb  + incr
   krowsb  = krowsb - incr
   kridxs  = kridx
   iilrow  = iilrow + incr
   IF ( iilrow > i1 ) iilrow = iilrow - i1
   GO TO 1180
   
! LAST ROW NUMBERS DO NOT MATCH, KROWN LT LROWN
   
   1150  CONTINUE
   incr    = lrown - lrow1 + 1
   GO TO 1145
   
! ROWS MATCH FOR INNER LOOP COLUMN VALUES, NOW DETERMINE THE COLUMN INDEX
! FOR THE NEXT COLUMN TO ADD TO THE INNER AND OUTER LOOP ARRAYS.
   1170  CONTINUE
   kfrcolg = klscolp+1
   iilrow  = iilrow1
   GO TO 1400
   
! NOT ALL NEEDED ROW VALUES ARE PRESENT, MUST GET NEEDED ROWS
! FOR ALL COLUMNS REQUIRED FOR THIS PIVOT COLUMN
   
   1180  CONTINUE
   kfrcolg = kfrcol
   GO TO 1400
   
! NO MATCH FOUND, WILL START LOADING THE INNER AND OUTER LOOP ARRAYS
! FROM THE BEGINNING
   
   1350  iilrow1 = 1
   iilrow  = 1
   krowb   = krow1
   krowsb  = krown - krow1 + 1
   kfrcolg = kfrcol
   1400  CONTINUE
   kridx   = kmidx+4
   DO  j = 1, km2
     nar( j ) = zi( kridx+j-1 )
   END DO
   nar( km2+1 ) = 0
   iilrowb = iilrow
   
! KFRCOL  = FIRST COLUMN NEEDED FOR PIVOT COLUMN "K"
! KLSCOL  = LAST COLUMN NEEDED FOR PIVOT COLUMN "K"
! KFRCOLG = FIRST COLUMN TO BE PLACED IN INNER/OUTER LOOP ARRAYS
! KFRCOLP = FIRST COLUMN OF LAST PIVOT COLUMN PROCESSED
! KLSCOLP = LAST COLUMN OF LAST PIVOT COLUMN PROCESSED
   
!      PRINT *,' KFRCOL,KLSCOL,KFRCOLG,KFRCOLP,KLSCOLP,KAROWS='
!      PRINT *,  KFRCOL,KLSCOL,KFRCOLG,KFRCOLP,KLSCOLP,KAROWS
!      PRINT *,' KROWB,KROWSB,IILROW1,IILROW,kridx='
!      PRINT *,  KROWB,KROWSB,IILROW1,IILROW,kridx
   
! KLSCOL WILL BE LESS THAN KFRCOLG FOR THE FIRST COLUMN AND FOR ANY
! COLUMN THAT DOES NOT NEED A PRECEEDING COLUMN OF DATA
   
   IF ( klscol < kfrcolg ) GO TO 6000
   DO  j = kfrcolg, klscol
     iilcol = MOD ( j, i2 )
     IF ( iilcol == 0 ) iilcol = i2
     jcol   = j
     jdir   = j*4 - 3
     jmidx  = zi( jdir   )
     
! SEE IF COLUMN DATA IS IN MEMORY OR ON THE SPILL FILE
     
     IF (  jmidx /= 0 ) GO TO 1500
     
! DATA IS ON THE SPILL FILE
     
     CALL smcspl ( jcol, zi )
     IF ( zi( jdir ) == 0 ) jmidx = ispill
     IF ( zi( jdir ) /= 0 ) jmidx = zi( jdir )
     1500  CONTINUE
     jridx  = jmidx + 4
     jm2    = zi( jmidx + 1 )
     jridxn = jridx + jm2
     jrowl  = zi( jridx+jm2-2 ) + zi( jridx+jm2-1 ) - 1
     jvidx  = jridxn
     
! SAVE DIAGONAL TERM FOR COLUMN J ; (ALWAYS, THE FIRST TERM)
     
     zol( iilcol, ic2 )  = 1.0 / zs( jvidx )
     
! FOR EACH COLUMN J, GET REQUIRED ROWS; I.E, ACTIVE ROWS OF COLUMN K
     
     IF ( j > klscolp ) GO TO 1530
     
! SET VARIABLES FOR ADDING ROW TERMS TO AN EXISTING COLUMN IN THE INNER LOOP
     
     kridx  = kridxs
     krow   = krowb
     krows  = krowsb
     iilrow = iilrowb
     
! SET LASROW TO ZERO IF THIS COLUMN IS BEING RELOADED INTO ZIL AND NOT
! BEING ADDED TO FROM SOME PREVIOUS COLUMN PROCESSING.
     
     IF ( iilrowb == iilrow1 ) lasrow( j ) = 0
     GO TO 1540
     1530  CONTINUE
     
!  MUST RESET KRIDX, KROW AND KROWS FOR INSERTION OF NEW COLUMN IN INNER LOOP
     
     kridx  = kmidx+4
     krow   = zi( kridx   )
     krows  = zi( kridx+1 )
     iilrow = iilrow1
     1540  CONTINUE
     krown  = krow + krows - 1
     
! JROWL IS LAST ROW TERM IN COLUMN "J".  IF THIS IS BEFORE THE FIRST ROW
! "KROW" TERM NEEDED, THEN NO MORE TERMS ARE NEEDED FROM COLUMN "J" AND
! "LASROW" WILL INDICATE THE LAST VALUE STORED FOR COLUMN "J".
     
     IF ( jrowl < krow ) CYCLE
     2000  jrow   = zi( jridx )
     jrows  = zi( jridx+1 )
     jrown  = jrow + jrows - 1
     2010  CONTINUE
     IF ( jrown < krow  ) GO TO 2895
     IF ( jrow  > krown ) GO TO 2400
     missin = krow - jrow
     
! CHECK TO SEE IF THERE ARE MISSING TERMS, I.E., TERMS CREATED DURING
! THE DECOMPOSITION.  IF THERE ARE MISSING TERMS, THEN SET THEIR VALUES
! TO BE INITIALLY ZERO.
     
     IF ( missin >= 0 ) GO TO 2050
     nzeros = IABS( missin )
     
!  STORE "NZEROS" NUMBER OF ZEROS FOR INNER LOOP TERMS
     
     iavail = i1 - ( iilrow+nzeros-1 )
     IF ( iavail < 0 ) GO TO 2022
     DO  i = 1, nzeros
       zil( iilrow+i-1, iilcol ) = 0.0
     END DO
     iilrow = iilrow + nzeros
     GO TO 2028
     2022  ilim1 = i1 - iilrow + 1
     ilim2 = nzeros - ilim1
     DO  i = 1, ilim1
       zil( iilrow+i-1, iilcol ) = 0.0
     END DO
     DO  i = 1, ilim2
       zil( i, iilcol ) = 0.0
     END DO
     iilrow = ilim2 + 1
     2028  CONTINUE
     krow  = krow  + nzeros
     krows = krows - nzeros
     2050  CONTINUE
     IF ( missin <= 0 ) GO TO 2070
     iskip  = krow  - jrow
     jvidx  = jvidx + iskip*nvterm
     jrow   = jrow  + iskip
     2070  CONTINUE
     irown  = MIN0 ( krown, jrown )
     num    = irown - krow + 1
     
!  MOVE INNER LOOP VALUES FROM IN-MEMORY LOCATION TO
!  THE INNER LOOP AREA
     
     nrows = irown - krow + 1
     IF ( nrows > ( i1 - iilrow + 1 ) ) GO TO 2120
     DO  i = 1, nrows
       zil( iilrow+i-1, iilcol ) = zs(jvidx+i-1 )
     END DO
     iilrow = iilrow + nrows
     GO TO 2180
     2120  ilim1 = i1 - iilrow + 1
     ilim2 = nrows - ilim1
     DO  i = 1, ilim1
       zil( iilrow+i-1, iilcol ) = zs( jvidx+i-1 )
     END DO
     jvtmp = jvidx + ilim1
     DO  i = 1, ilim2
       zil( i, iilcol ) = zs( jvtmp+i-1 )
     END DO
     iilrow = ilim2 + 1
     2180  CONTINUE
     lasrow( iilcol ) = iilrow
     
! IF ALL OF THE ROWS ARE NON-ZERO, SET LASROW COUNTER TO IILROW1
     
     IF ( iilrow == iilrow1 ) lasrow( iilcol ) = iilrow1
     jvidx  = jvidx + nrows
     jrow   = jrow  + nrows
     krow   = irown + 1
     krows  = krown - irown
     
! INCREMENT EITHER KROW OR JROW DEPENDING UPON WHETHER IROWN = JROWN
! OR IROWN = KROWN
     
     IF ( irown == jrown ) GO TO 2900
     GO TO 2530
     2400  CONTINUE
     
! STORE ZEROS FOR CREATED TERMS AND INCREMENT TO THE NEXT SET OF
! OF ROWS FOR THIS PIVOTAL COLUMN.
     
     iavail = i1 - ( iilrow+krows-1 )
     IF ( iavail < 0 ) GO TO 2522
     DO  i = 1, krows
       zil( iilrow+i-1, iilcol ) = 0.0
     END DO
     iilrow = iilrow + krows
     GO TO 2528
     2522  CONTINUE
     ilim2  = krows - ( i1 - iilrow + 1 )
     DO  i = iilrow1, i1
       zil( i, iilcol ) = 0.0
     END DO
     DO  i = 1, ilim2
       zil( i, iilcol ) = 0.0
     END DO
     iilrow = ilim2 + 1
     2528  CONTINUE
     
! INCREMENT THE INDEX TO THE NEXT SET OF ROWS FOR COLUMN "K"
     
     2530  kridx  = kridx + 2
     
! IF THERE ARE NO MORE ROWS FOR THIS COLUMN, THEN COLUMN IS COMPLETE
     
     IF ( kridx >= kridxn ) CYCLE
     krow   = zi( kridx )
     krows  = zi( kridx+1)
     krown  = krow + krows - 1
     GO TO 2010
     2895  CONTINUE
     
! INCREMENT "JVIDX" TO POINT TO THE CORRESPONDING VALUE TERM FOR THE
! NEXT ROW OF COLUMN "J"
     
     jvidx  = jvidx + ( jrown - jrow + 1 )*nvterm
     
! INCREMENT THE INDEX TO THE NEXT SET OF ROWS FOR COLUMN "J"
     
     2900  jridx  = jridx + 2
     IF ( jridx >= jridxn ) CYCLE
     GO TO 2000
   END DO
   IF ( k == 1 ) GO TO 6000
   
! COMPUTE THE TERMS FOR THE CURRENT COLUMN OF DATA
   
!      do 100 k = 1,n
!         do 10  i = k,n
!         temp = 0.
!         do 5  l = 1,k-1
!            temp = temp + a(i,l)*a(k,l) / a(l,l)
!    5       continue
!         a(i,k) = a(i,k) - temp
!   10    continue
   
!  THE FOLLOWING LAST COMPUTATION TAKES PLACE IN SUBROUTINE SMCOUT.
!  THE RESULTS OF THE DIVISION ARE WRITTEN TO THE OUTPUT FILE BUT
!  THE RESULTS OF THE ABOVE (WITHOUT THE DIVISION BELOW) IS
!  MAINTAINED IN MEMORY FOR REMAINING COLUMN COMPUTATIONS.
   
!         do 11 j = k+1,n
!           a(k,j) = a(j,k) / a( k,k )
!   11      continue
!  100 continue
   
!   NROWS  = NUMBER OF ROWS STORED IN INNER LOOP
!   KCOL   = LAST COLUMN NUMBER STORED IN INNER LOOP
!   KFRCOL = FIRST COLUMN NUMBER STORED IN INNER LOOP
   
   nrows = karows
   kdir  = ( kcol-1 ) * 4 + 1
   kmidx = zi( kdir )
   kridx = kmidx + 4
   km2   = zi( kmidx+1 )
   kvidx = kridx + km2
   ilim1   = iilrow1 + nrows - 1
   ilim2   = 0
   iavail  = i1 - ilim1
   IF ( iavail >= 0 ) GO TO 4010
   ilim1   = i1
   ilim2   = nrows - ( i1 - iilrow1 + 1 )
   4010  CONTINUE
   jlim1   = MOD( kfrcol, i2 )
   jlim2   = MOD( klscol, i2 )
   IF ( jlim1 == 0 ) jlim1 = i2
   IF ( jlim2 == 0 ) jlim2 = i2
   jlim4   = 0
   IF ( kfrcol == k ) GO TO 6000
   IF ( jlim2 >= jlim1 ) GO TO 4015
   jlim4   = jlim2
   jlim2   = i2
   4015  CONTINUE
!      PRINT *,' JLIM1,JLIM2,JLIM4,IILROW1=',JLIM1,JLIM2,JLIM4,IILROW1
!      PRINT *,' ILIM1,ILIM2,JLIM1,JLIM2,JLIM4,IILROW1,NROWS'
!      PRINT *,  ILIM1,ILIM2,JLIM1,JLIM2,JLIM4,IILROW1,NROWS
   IF ( k == 1 ) GO TO 4007
   
! COMPUTE THE OUTER LOOP TERM FOR THIS COLUMN J
! I.E.,   -A(K,J) / A(J,J)
!  where K = current pivot column number; J = column being processed
   
!     KAROWS = NUMBER OF ACTIVE ROWS FOR THE CURRENT PIVOTAL COLUMN
!     JCOL   = COLUMN NUMBER OF CURRENT PIVOTAL COLUMN
!     ZOL(KBC,IC1) = FIRST ACTIVE ROW ("IILROW1") TERM OF COLUMN "KBC"
!     ZOL(KBC,IC2) = DIAGONAL TERM FOR COLUMN "KBC"
   
   DO  kbc = jlim1, jlim2
     zol( kbc, ic1 ) = zil( iilrow1, kbc ) * zol( kbc, ic2 )
   END DO
   IF ( jlim4 == 0 ) GO TO 4007
   DO  kbc = 1, jlim4
     zol( kbc, ic1 ) = zil( iilrow1, kbc ) * zol( kbc, ic2 )
   END DO
   4007  CONTINUE
!      CALL KBHELPRS( KFRCOL, KLSCOL, ZOL, ZIL, I1, I2, LASROW )
   DO  i = iilrow1, ilim1
     rtemp(i) = 0.0
   END DO
   
! PROCESS COLUMNS JLIM1 THROUGH JLIM2
   
   DO  j = jlim1, jlim2
     limit = ilim1
     itest = lasrow( j )
     IF ( itest == 0 ) CYCLE
     IF ( itest > iilrow1 ) limit = itest - 1
     
! PROCESS ROWS IILROW1 THROUGH LIMIT FOR COLUMNS JLIM1 THROUGH JLIM2
     
     zoltmp = zol( j,ic1 )
     CALL smccrs ( rtemp(iilrow1), zil( iilrow1,j ), limit-iilrow1+1  &
         , zoltmp )
!      DO 4020 I = IILROW1, LIMIT
!      RTEMP(I) = RTEMP(I) + ZIL( I, J ) * ZOLTMP
!4020  CONTINUE
   END DO
   IF ( jlim4 == 0 ) GO TO 4030
   
! PROCESS ROWS IILROW1 THROUGH LIMIT FOR COLUMNS 1 THROUGH JLIM4
   
   DO  j = 1, jlim4
     itest = lasrow( j )
     IF ( itest == 0 ) CYCLE
     limit = ilim1
     IF ( itest > iilrow1 ) limit = itest - 1
     zoltmp = zol( j,ic1 )
     CALL smccrs ( rtemp(iilrow1), zil( iilrow1,j ), limit-iilrow1+1  &
         , zoltmp )
!      DO 4023 I = IILROW1, LIMIT
!      RTEMP(I) = RTEMP(I) + ZIL( I, J ) * ZOLTMP
!4023  CONTINUE
   END DO
   4030  CONTINUE
   IF ( ilim2 == 0 ) GO TO 4060
   DO  i = 1, ilim2
     rtemp(i) = 0.0
   END DO
   
! PROCESS COLUMNS JLIM1 THROUGH JLIM2
   
   DO  j = jlim1, jlim2
     itest = lasrow( j )
     IF ( itest == 0 .OR. itest > iilrow1 ) CYCLE
     limit = ilim2
     IF ( itest <= ilim2 ) limit = itest - 1
     
! PROCESS ROWS 1 THROUGH LIMIT FOR COLUMNS JLIM1 THROUGH JLIM2
     
     zoltmp = zol( j,ic1 )
     CALL smccrs ( rtemp(1), zil( 1,j ), limit, zoltmp )
!      DO 4040 I = 1, LIMIT
!      RTEMP(I) = RTEMP(I) + ZIL( I, J ) * ZOLTMP
!4040  CONTINUE
   END DO
   IF ( jlim4 == 0 ) GO TO 4046
   
! PROCESS ROWS 1 THROUGH LIMIT FOR COLUMNS 1 THROUGH JLIM4
   
   DO  j = 1, jlim4
     itest = lasrow( j )
     IF ( itest == 0 .OR. itest > iilrow1 ) CYCLE
     limit = ilim2
     IF ( itest <= ilim2 ) limit = itest - 1
     zoltmp = zol( j,ic1 )
     CALL smccrs ( rtemp(1), zil( 1,j ), limit, zoltmp )
!      DO 4043 I = 1, LIMIT
!      RTEMP(I) = RTEMP(I) + ZIL( I, J ) * ZOLTMP
!4043  CONTINUE
   END DO
   4046  CONTINUE
   4060  CONTINUE
   
! UPDATE EACH ACTIVE ROW TERM FOR COLUMN "K" BY SUBTRACTING "RTEMP"
   
   DO  i = iilrow1, ilim1
     zs( kvidx ) = zs( kvidx ) - rtemp(i)
     kvidx = kvidx + 1
   END DO
   IF ( ilim2 == 0 ) GO TO 4070
   DO  i = 1, ilim2
     zs( kvidx ) = zs( kvidx ) - rtemp(i)
     kvidx = kvidx + 1
   END DO
   4070  CONTINUE
   
! CALL SMCOUT TO WRITE OUT THE COLUMN TO THE OUTPUT LOWER TRIANGULAR
! MATRIX FILE
   
   6000  CONTINUE
   CALL smcout ( zi, zi, zs, zol( 1,ic1 ), zol( 1,ic1 ) )
 END DO
 RETURN
END SUBROUTINE smc2rs
