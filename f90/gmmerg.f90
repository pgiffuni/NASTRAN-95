SUBROUTINE gmmerg(filea,file11,file21,file12,file22,rpart,cpart,  &
        nsub,mrgtyp,core,lcore)
     
!     GENERAL MATRIX MERGE ROUTINE
 
 
!                  --               --
!                  I        I        I
!                  I FILE11 I FILE12 I   --     --
!                  I        I        I   I       I
!                  I-----------------I = I FILEA I
!                  I        I        I   I       I
!                  I FILE21 I FILE22 I   --     --
!                  I        I        I
!                  --               --
 
!        WHERE
 
!             RPART - ROW PARTITIONING VECTOR
!             NSUB(1) - NUMBER OF COLUMNS IN RPART 0 SUBSET
!             NSUB(2) - NUMBER OF COLUMNS IN RPART 1 SUBSET
!             NSUB(3) - NUMBER OF ROWS IN CPART 0 SUBSET
!             NSUB(4) - NUMBER OF ROWS IN CPART 1 SUBSET
!             MRGTYP - MERGE TYPE (1 .EQ. SQUARE, 2 .EQ. RECTANGULAR)
!             CPART - COLUMN PARTITION VECTOR
 
 
 
 INTEGER, INTENT(IN)                      :: filea
 INTEGER, INTENT(IN)                      :: file11
 INTEGER, INTENT(IN)                      :: file21
 INTEGER, INTENT(IN)                      :: file12
 REAL, INTENT(IN)                         :: file22
 INTEGER, INTENT(IN)                      :: rpart
 INTEGER, INTENT(IN)                      :: cpart
 INTEGER, INTENT(IN)                      :: nsub(4)
 INTEGER, INTENT(IN)                      :: mrgtyp
 INTEGER, INTENT(OUT)                     :: core(6)
 INTEGER, INTENT(IN)                      :: lcore
 INTEGER ::  rule     , name(2),  &
             rp(7)    , cp(7)
 
 COMMON / parmeg /       ia(7)    ,ia11(7)  ,ia21(7)  ,ia12(7),  &
                         ia22(7)  ,lcr      ,rule
 
 
 
 DATA NAME / 4HGMME , 4HRG   /
 
!***********************************************************************
 
!     GET TRAILERS FOR INPUTS
 
 rp(1) = rpart
 IF(rpart /= 0) CALL rdtrl(rp)
 cp(1) = cpart
 IF(cpart /= 0) CALL rdtrl(cp)
 
 DO  i=2,7
   ia(i) = 0
   ia11(i) = 0
   ia12(i) = 0
   ia21(i) = 0
   ia22(i) = 0
 END DO
 
 ia11(1) = file11
 IF(file11 /= 0) CALL rdtrl(ia11)
 IF(ia11(1) < 0) ia11(1) = 0
 ia12(1) = file12
 IF(file12 /= 0) CALL rdtrl(ia12)
 IF(ia12(1) < 0) ia12(1) = 0
 ia21(1) = file21
 IF(file21 /= 0) CALL rdtrl(ia21)
 IF(ia21(1) < 0) ia21(1) = 0
 ia22(1) = file22
 IF(file22 /= 0) CALL rdtrl(ia22)
 IF(ia22(1) < 0) ia22(1) = 0
 
!     SET UP MATRIX CONTROL BLOCK FOR OUTPUT
 
 ia(1) = filea
 ia(4) = mrgtyp
 ia(5) = MAX0(ia11(5),ia12(5),ia21(5),ia22(5))
 
!     SET UP DUMMY PARTITION VECTOR
 
 core(1) = 0
 core(2) = 1
 core(3) = ia(2)
 core(4) = 2
 core(5) = 1
 core(6) = 0
 lcr = lcore
 rule = 0
 
 IF(rpart == 0) GO TO 30
 IF(cpart == 0) GO TO 20
 
!     FULL MERGE
 
 ia(2) = nsub(1) + nsub(2)
 ia(3) = nsub(3) + nsub(4)
 CALL merge(rp,cp,core)
 GO TO 40
 
!  *  *  MERGE COLUMNS ONLY
 
 20 ia(2) = nsub(1) + nsub(2)
 ia(3) = MAX0(ia11(3),ia12(3))
 CALL merge(rp,core,core)
 GO TO 40
 
!  *  *  MERGE ROWS ONLY
 
 30 IF(cpart == 0) GO TO 1007
 ia(2) = MAX0(ia11(2),ia21(2))
 ia(3) = nsub(3) + nsub(4)
 CALL merge(core,cp,core)
 
!     WRITE TRIALER FOR OUTPUT
 
 40 CALL wrttrl(ia)
 
 RETURN
 
!     ILLEGAL INPUT - NO PARTITION VECTOR
 
 1007 CALL mesage(-7,0,NAME)
 RETURN
END SUBROUTINE gmmerg
