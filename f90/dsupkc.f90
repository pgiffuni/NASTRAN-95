SUBROUTINE dsupkc ( itin, itout, a, b )
     
 INTEGER, INTENT(IN)                      :: itin
 INTEGER, INTENT(IN)                      :: itout
 REAL, INTENT(IN)                         :: a(4)
 REAL, INTENT(OUT)                        :: b(4)
 COMMON / system / isysbf, iwr
 REAL :: aa(4), bb(4)
 INTEGER :: nwords(4)
 REAL :: rs1, rs2
 DOUBLE PRECISION :: rd1, rd2, rdi1, rdi2
 EQUIVALENCE       (aa,rs1,rd1), (bb,rs2,rd2)
 EQUIVALENCE       ( aa(3), rdi1 ), ( bb(3), rdi2 )
 DATA              nwords / 1,2,2,4/
 
 iwrd1   = nwords( itin )
 IF ( itin /= itout ) GO TO 20
!DIR$ NEXTSCALAR
 DO  k = 1, iwrd1
   b( k )  = a( k )
 END DO
 GO TO 7777
 20      IF ( itout > 64 ) GO TO 30
 itout2 = itout
 iwrd2  = nwords( itout )
 ssign  = 1.0
 GO TO 40
 30      itout2 = itout - 64
 iwrd2  = nwords( itout2 )
 ssign  = -1.0
!DIR$ NEXTSCALAR
 40      DO  k = 1, iwrd1
   aa( k ) = a( k )
 END DO
 SELECT CASE ( itin )
   CASE (    1)
     GO TO  1000
   CASE (    2)
     GO TO  2000
   CASE (    3)
     GO TO  3000
   CASE (    4)
     GO TO  4000
 END SELECT
 1000    SELECT CASE ( itout2 )
   CASE (    1)
     GO TO  1100
   CASE (    2)
     GO TO  1200
   CASE (    3)
     GO TO  1300
   CASE (    4)
     GO TO  1400
 END SELECT
 1100    rs2 = ssign * rs1
 GO TO 7000
 1200    rd2 = ssign * rs1
 GO TO 7000
 1300    bb( 1 ) = ssign * rs1
 bb( 2 ) = 0.
 GO TO 7000
 1400    rd2 = ssign * rs1
 rdi2 = 0.
 GO TO 7000
 2000    SELECT CASE ( itout2 )
   CASE (    1)
     GO TO  2100
   CASE (    2)
     GO TO  2200
   CASE (    3)
     GO TO  2300
   CASE (    4)
     GO TO  2400
 END SELECT
 2100    rs2 = ssign * rd1
 GO TO 7000
 2200    rd2 = ssign * rd1
 GO TO 7000
 2300    bb( 1 ) = ssign * rd1
 bb( 2 ) = 0.
 GO TO 7000
 2400    rd2 = ssign * rd1
 rdi2 = 0.
 GO TO 7000
 3000    SELECT CASE ( itout2 )
   CASE (    1)
     GO TO  3100
   CASE (    2)
     GO TO  3200
   CASE (    3)
     GO TO  3300
   CASE (    4)
     GO TO  3400
 END SELECT
 3100    rs2 = ssign * aa( 1 )
 GO TO 7000
 3200    rd2 = ssign * aa( 1 )
 GO TO 7000
 3300    bb( 1 ) = ssign * aa( 1 )
 bb( 2 ) = ssign * aa( 2 )
 GO TO 7000
 3400    rd2 = ssign * aa( 1 )
 rdi2 = ssign * aa( 2 )
 GO TO 7000
 4000    SELECT CASE ( itout2 )
   CASE (    1)
     GO TO  4100
   CASE (    2)
     GO TO  4200
   CASE (    3)
     GO TO  4300
   CASE (    4)
     GO TO  4400
 END SELECT
 4100    rs2 = ssign * rd1
 GO TO 7000
 4200    rd2 = ssign * rd1
 GO TO 7000
 4300    bb( 1 ) = ssign * rd1
 bb( 2 ) = ssign * rdi1
 GO TO 7000
 4400    rd2 = ssign * rd1
 rdi2 = ssign * rdi1
 GO TO 7000
!DIR$ NEXTSCALAR
 7000    DO  k = 1, iwrd2
   b( k ) = bb( k )
 END DO
 7777    RETURN
END SUBROUTINE dsupkc
