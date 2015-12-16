SUBROUTINE gmprtn(filea,file11,file21,file12,file22,rpart,cpart,  &
        nsub0,nsub1,core,lcore)
     
!     GENERAL MATRIX PARTION ROUTINE
 
 
!                              --               --
!                              I        I        I
!                  --     --   I FILE11 I FILE12 I
!                  I       I   I        I        I
!                  I FILEA I = I-----------------I
!                  I       I   I        I        I
!                  --     --   I FILE21 I FILE22 I
!                              I        I        I
!                              --               --
 
!        WHERE
 
!             RPART - ROW PARTITIONING VECTOR
!             CPART - COLUMN PARTITION VECTOR
 
 
 INTEGER, INTENT(IN)                      :: filea
 INTEGER, INTENT(IN)                      :: file11
 INTEGER, INTENT(IN)                      :: file21
 INTEGER, INTENT(IN)                      :: file12
 INTEGER, INTENT(IN)                      :: rpart
 INTEGER, INTENT(IN)                      :: cpart
 INTEGER, INTENT(IN)                      :: nsub0
 INTEGER, INTENT(IN)                      :: nsub1
 INTEGER, INTENT(OUT)                     :: core(6)
 INTEGER, INTENT(IN)                      :: lcore
 INTEGER :: file22 , rule     , NAME(2)  &
     ,rp(7)    ,cp(7)
 
 COMMON / parmeg /       ia(7)    ,ia11(7)  ,ia21(7)  ,ia12(7)  &
     ,ia22(7)  ,lcr      ,rule
 
 DATA NAME / 4HGMPR , 4HTN   /
 
!***********************************************************************
 
!     GET TRAILERS FOR INPUTS
 
 rp(1) = rpart
 IF(rpart /= 0) CALL rdtrl(rp)
 cp(1) = cpart
 IF(cpart /= 0) CALL rdtrl(cp)
 ia(1) = filea
 CALL rdtrl(ia)
 
!     SET UP MATRIX CONTROL BLOCKS FOR OUTPUTS
 
 ia11(1) = file11
 ia12(1) = file12
 ia21(1) = file21
 ia22(1) = file22
 
 DO  i=2,5
   ia11(i) = ia(i)
   ia12(i) = ia(i)
   ia21(i) = ia(i)
   ia22(i) = ia(i)
 END DO
 
!     SET UP DUMMY PARTITION VECTOR
 
 core(1) = 0
 core(2) = 1
 core(3) = ia(2)
 core(4) = 2
 core(5) = 1
 core(6) = 0
 
 rule = 0
 lcr = lcore
 
 IF(rpart == 0) GO TO 30
 IF(cpart == 0) GO TO 20
 
!     FULL PARTITION
 
 ia11(3) = nsub0
 ia12(3) = nsub0
 ia21(3) = nsub1
 ia22(3) = nsub1
 CALL partn(rp,cp,core)
 GO TO 40
 
!  *  *  PARTITION COLUMNS ONLY
 
 20 CALL partn(rp,core,core)
 GO TO 40
 
!  *  *  PARTITION ROWS ONLY
 
 30 IF(cpart == 0) GO TO 1007
 ia11(3) = nsub0
 ia12(3) = nsub0
 ia21(3) = nsub1
 ia22(3) = nsub1
 CALL partn(core,cp,core)
 
!     WRITE TRAILERS FOR OUTPUTS
 
 40 IF(ia11(1) /= 0) CALL wrttrl(ia11)
 IF(ia12(1) /= 0) CALL wrttrl(ia12)
 IF(ia21(1) /= 0) CALL wrttrl(ia21)
 IF(ia22(1) /= 0) CALL wrttrl(ia22)
 
 RETURN
 
!     ILLEGAL INPUT - NO PARTITION VECTOR
 
 1007 CALL mesage(-7,0,NAME)
 RETURN
END SUBROUTINE gmprtn
