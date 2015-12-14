SUBROUTINE ferswd(v1,v3,vb)
     
!  The original to this subroutine was FRSW2.  It has been modified
!  to read the matrix data from memory and after this data is exhausted
!  then to read the remaining data from the file.
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: v1(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: v3(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: vb(1)
 DOUBLE PRECISION :: ,                 xl(1)     ,xljj      ,v3j  &
     ,                 zero      ,sum       ,dcore(1)
 INTEGER :: iblk(20)
 COMMON / zzzzzz / icore(1)
 COMMON  /opinv /  mcblt(7) ,mcbsma(7)
 COMMON  /system/  ksystm(65)
 COMMON  /feerim/  nidsma, nidlt    , nidorv  , nltli  &
     ,                 nsmali, ibfsma   , ibflt  &
     ,                 ibforv, smapos(7), ltpos(7)
 EQUIVALENCE       (ksystm(02),io)
 EQUIVALENCE       ( dcore(1),icore(1), xl(1) )
 DATA              zero / 0.0D0 /
 
 nrow    = mcblt(2)
 CALL ferltd(mcbsma(1),v1(1),v3(1),vb(1))
!   FORWARD SWEEP DIRECTLY ON V3
 icrow = 1
 IF ( nidlt == 0 ) GO TO 1005
 ilrow = ltpos( 1 )
 mem   = nidlt
 DO  j = 1,nrow
   icrow = j
   IF ( icrow > ilrow ) GO TO 1000
   140 icol  = icore(mem)
   IF ( icol /= j ) GO TO 180
   ji    = mem/2+2
   ntms  = icore(mem+1)
   ntmss = ntms
   ik    = icore(mem+2+2*ntms)
   IF(ik /= j) GO TO 150
   ntms  = ntms - 1
   xljj  = dcore(ji)
   ji    = ji + 1
   ik    = ik + 1
   150 IF(ntms == 0) GO TO 170
   v3j   = v3(j)
   DO  ii = 1,ntms
     v3(ik)= v3(ik) + dcore(ji) * v3j
     ik    = ik + 1
     ji    = ji + 1
   END DO
   170 mem   = mem + ntmss*2 + 4
   GO TO 140
   180 v3(j) = v3(j) / xljj
 END DO
 GO TO 3000
 1000 CONTINUE
! POSITION FILE TO APPROPRIATE COLUMN
 CALL dsspos ( mcblt, ltpos(2), ltpos(3), ltpos(4) )
 GO TO 1008
 1005 CONTINUE
 CALL REWIND ( mcblt )
 CALL skprec ( mcblt, 1 )
 1008 CONTINUE
 iblk( 1 ) = mcblt( 1 )
 
! CONTINUE WITH FORWARD SWEEP
 
 DO  j = icrow, nrow
   iblk( 8 ) = -1
   1030 CALL getstr ( *1070, iblk )
   ik   = iblk( 4 )
   ji   = iblk( 5 )
   ntms = iblk( 6 )
   IF ( ik /= j ) GO TO 1040
   ntms = ntms - 1
   xljj = xl( ji )
   ji   = ji + 1
   ik   = ik + 1
   1040 IF ( ntms == 0 ) GO TO 1060
   v3j  = v3( j )
   IF ( v3j == zero ) GO TO 1060
   DO  ii = 1, ntms
     v3( ik ) = v3( ik ) + xl(ji)*v3j
     ik   = ik + 1
     ji   = ji + 1
   END DO
   1060 CALL endget ( iblk )
   GO TO 1030
   1070 CONTINUE
   v3( j ) = v3( j ) / xljj
 END DO
 2000 CONTINUE
 
!     BACKWARD SUBSTITUTION OMIT DIAGONAL
 
 icrow = nrow
 IF ( j == 1 ) RETURN
 IF ( ilrow == nrow .AND. nidlt /= 0 ) GO TO 3000
 j     = nrow
 2090 iblk( 8 ) = -1
 2100 CALL getstb ( *2130, iblk )
 ik    = iblk( 4 )
 ji    = iblk( 5 )
 ntms  = iblk( 6 )
 IF ( ik-ntms+1 == j ) ntms = ntms - 1
 IF ( ntms == 0 ) GO TO 2120
 sum   = zero
 DO  ii = 1, ntms
   sum   = sum + xl(ji) * v3(ik)
   ji    = ji - 1
   ik    = ik - 1
 END DO
 v3( j ) = v3( j ) + sum
 2120 CALL endgtb ( iblk )
 GO TO 2100
 2130 IF ( j == 1 ) GO TO 7000
 j = j - 1
 IF ( j <= ilrow ) GO TO 3000
 GO TO 2090
!  CONTINUE BACKWARD SUBSTITUTION USING DATA FROM MEMORY
 3000 CONTINUE
 mem = mem - ntmss*2 - 4
 3200 CONTINUE
 3210 icol = icore(mem)
 IF ( icol /= j ) GO TO 3240
 ntms  = icore(mem+1)
 ntmss = ntms
 ji    = mem/2+1+ntms
 ik    = icore(mem+2+2*ntms)+ntms-1
 IF( ik-ntms+1 == j) ntms = ntms - 1
 IF( ntms == 0 ) GO TO 3230
 v3j   = v3( j )
 DO  ii = 1,ntms
   v3j   = v3j + dcore(ji) * v3(ik)
   ji    = ji-1
   ik    = ik-1
 END DO
 v3(j) = v3j
 3230 IF ( mem == nidlt ) GO TO 3250
 ntmsnx= icore(mem-1)
 mem   = mem - ntmsnx*2 - 4
 GO TO 3210
 3240 IF ( j == 1 )  GO TO 3250
 j     = j-1
 GO TO 3200
 3250 CONTINUE
 7000 CONTINUE
 RETURN
END SUBROUTINE ferswd
