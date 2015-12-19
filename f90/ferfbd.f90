SUBROUTINE ferfbd(v1,v2,v3,vb)
     
!  FERFBD is a modification of the old FRBK2 subroutine.  It has been
!  modified to read matrix data from memory until that data is exhausted
!  and then to read the remaining data from the file.
 
 
 DOUBLE PRECISION, INTENT(IN)             :: v1(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v2(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v3(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: vb(1)
 DOUBLE PRECISION :: dcore(1)
 DOUBLE PRECISION :: xl(1)     ,xljj      ,v3j       ,v2j
 INTEGER :: iblk(20)  ,smapos
 COMMON / zzzzzz / icore(1)
 COMMON / opinv  / mcblt(7)  ,mcbsma(7)
 COMMON / system / ksystm(65)
 COMMON / feerim / nidsma    ,nidlt     ,nidorv    ,nltli  &
     ,                 nsmali    ,ibfsma    ,ibflt  &
     ,                 ibforv    ,smapos(7) ,ltpos(7)
 EQUIVALENCE       ( ksystm(02),nout)
 EQUIVALENCE       ( dcore(1)  ,icore(1), xl )
 
 nrow    = mcblt(2)
 DO  i = 1,nrow
   v2(i) = v1(i)
 END DO
 ilrow = ltpos( 1 )
 icrow = nrow
!      PRINT *,' FERFBD,ILROW,NIDLT=',ILROW,NIDLT
!      PRINT *,' LTPOS=',LTPOS
 IF ( ilrow == nrow .AND. nidlt /= 0 ) GO TO 1000
 
!     BACKWARD SUBSTITUTION
 
!     POSITION FILE TO LAST COLUMN
 
 IF ( nidlt == 0 ) GO TO 12
 11  CONTINUE
 CALL dsspos ( mcblt, ltpos(5), ltpos(6), ltpos(7) )
 GO TO 16
 12  IF ( ltpos( 5 ) /= -1 ) GO TO 11
 CALL REWIND ( mcblt )
 CALL skprec ( mcblt, nrow+1 )
 CALL dscpos ( mcblt, iblock, iclr, icbp )
 ltpos( 5 ) = iblock
 ltpos( 6 ) = iclr
 ltpos( 7 ) = icbp
 16  CONTINUE
 iblk( 1 ) = mcblt( 1 )
 j       = nrow
 15  iblk(8) = -1
 icrow   = j
 IF ( j <= ilrow ) GO TO 1000
 20  CALL getstb(*50,iblk(1))
 ntms    = iblk(6)
 ji      = iblk(5)
 ik      = iblk(4)
 IF( ik - ntms + 1 /= j) GO TO 25
 ntms    = ntms - 1
 xljj    = xl(ji-ntms)
 IF(ntms == 0) GO TO 40
 25 v2j     = v2(j)
 DO  ii= 1,ntms
   v2j     = v2j + xl(ji) * v2(ik)
   ji      = ji - 1
   ik      = ik - 1
 END DO
 v2(j)   = v2j
 40 CALL endgtb(iblk(1))
 GO TO 20
 50 v2(j)   = v2(j) / xljj
 IF(j == 1) GO TO 2000
 j       = j -1
 GO TO 15
 
!     CONTINUE BACKWARD SUBSTITUTION WITH DATA IN MEMORY
 
 1000  CONTINUE
 mem     = nltli
!      PRINT *,' AT 1000,NLTLI=',NLTLI
 ntms    = icore(mem)
!      PRINT *,' ICORE(NLTLI,-1=',ICORE(NLTLI),ICORE(NLTLI-1)
 mem     = mem - 2*ntms - 3
 j       = icrow
 1015 icol    = icore(mem)
!      PRINT *,' MEM,ICORE(MEM-1,0,+1=',MEM,ICORE(MEM-1),ICORE(MEM),
!     & ICORE(MEM+1)
!      PRINT *,' ICOL,MEM,NTMS,ICROW,J=',ICOL,MEM,NTMS,ICROW,J
 IF ( icol /= j ) GO TO 1050
 ntms    = icore(mem+1)
!      PRINT *,' FERFBD,A1015,J,NTMS,ICOL=',J,NTMS,ICOL
 ntmss   = ntms
 ji      = mem/2 + 1 + ntms
 ik      = icore( mem + 2 + 2*ntms ) + ntms - 1
!      PRINT *,' FERFBD,IK=',IK
 IF( ik-ntms+1 /= j) GO TO 1025
 ntms    = ntms - 1
 xljj    = dcore(ji-ntms)
!      PRINT *,' FERFBD,XLJJ=',XLJJ
 IF(ntms == 0) GO TO 1040
 1025 v2j     = v2(j)
 DO  ii= 1,ntms
   v2j     = v2j + dcore(ji) * v2(ik)
   ji      = ji - 1
   ik      = ik - 1
 END DO
 v2(j)   = v2j
 1040 IF ( mem == nidlt ) GO TO 1050
 ntmsnx  = icore( mem-1 )
 mem     = mem - 2*ntmsnx - 4
 GO TO 1015
 1050 v2(j)   = v2(j) / xljj
 IF(j == 1) GO TO 2000
 j       = j -1
 GO TO 1015
 2000  CONTINUE
 CALL ferltd(mcbsma(1),v2(1),v3(1),vb(1) )
 
! BEGIN FORWARD SWEEP DIRECTLY ON V3
 
 icrow = 1
 IF ( nidlt == 0 ) GO TO 3005
 mem=nidlt
 DO  j = 1, nrow
   icrow = j
   IF ( j > ilrow ) GO TO 3000
   2080 icol  = icore(mem)
   IF( icol /= j ) CYCLE
   ji    = mem/2 + 2
   ntms  = icore( mem+1 )
   ntmss = ntms
   ik    = icore(mem + 2 + 2*ntms)
   IF ( ik /= j ) GO TO 2085
   ntms  = ntms - 1
   v3(j) = v3(j) / dcore(ji)
   ji    = ji + 1
   ik    = ik + 1
   2085 IF(ntms == 0) GO TO 2100
   v3j   = v3(j)
   DO  ii = 1,ntms
     v3(ik)= v3(ik) + dcore(ji) * v3j
     ik    = ik + 1
     ji    = ji + 1
   END DO
   2100 mem   = mem + 2*ntmss + 4
   GO TO 2080
 END DO
 GO TO 7000
 3000 CONTINUE
 
!     CONTINUE FORWARD SWEEP DIRECTLY ON V3
 
!     POSITION FILE TO CONTINUE READING COLUMN DATA NOT IN MEMORY
 
 CALL dsspos ( mcblt, ltpos(2), ltpos(3), ltpos(4) )
 GO TO 3008
 3005 CALL REWIND ( mcblt )
 CALL skprec ( mcblt, 1 )
 3008 CONTINUE
 DO  j = icrow, nrow
   iblk( 8 ) = -1
   3080 CALL getstr( *3120, iblk )
!      PRINT *,' GETSTR,J,IBLK(12=',J,IBLK(12)
   ik    = iblk( 4 )
   ji    = iblk( 5 )
   ntms  = iblk( 6 )
   IF ( ik /= j) GO TO 3085
   ntms  = ntms - 1
!      PRINT *,' IK,JI,XL(JI=',IK,JI,XL(JI)
   v3(j) = v3(j) / xl(ji)
   ji    = ji + 1
   ik    = ik + 1
   3085 IF(ntms == 0) GO TO 3100
   v3j   = v3(j)
   DO  ii = 1,ntms
     v3(ik)= v3(ik) + xl(ji) * v3j
     ik    = ik + 1
     ji    = ji + 1
   END DO
   3100 CALL endget(iblk(1))
   GO TO 3080
   3120 CONTINUE
 END DO
 GO TO 7000
 7000 CONTINUE
 RETURN
END SUBROUTINE ferfbd
