SUBROUTINE gfsmrg (filea,file11,file21,file12,file22,rpart,cpart)
     
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
!             CPART - COLUMN PARTITION VECTOR
 
 
 
 INTEGER, INTENT(IN)                      :: filea
 INTEGER, INTENT(IN)                      :: file11
 INTEGER, INTENT(IN)                      :: file21
 INTEGER, INTENT(IN)                      :: file12
 INTEGER, INTENT(IN)                      :: file22
 INTEGER, INTENT(IN)                      :: rpart
 INTEGER, INTENT(IN)                      :: cpart
 INTEGER :: rule     ,core     ,NAME(2) , rp(7)    ,cp(7)
 
!     OPEN CORE
 
 COMMON / zzzzzz /       core(1)
 
!     CALCV COMMON BLOCK
 
 COMMON / patx   /       lcore    ,nsub0    ,nsub1
 
!     PARTITION - MERGE COMMON BLOCK
 
 COMMON / parmeg /       ia(7)    ,ia11(7)  ,ia21(7)  ,ia12(7) ,  &
     ia22(7)  ,lcr      ,rule
 
 DATA     NAME   /       4HGFSM   , 4HRG    /
 
 
!     GET TRAILERS FOR INPUTS
 
 rp(1) = rpart
 IF (rpart /= 0) CALL rdtrl (rp)
 cp(1) = cpart
 IF (cpart /= 0) CALL rdtrl (cp)
 
 DO  i = 2,7
   ia(i)   = 0
   ia11(i) = 0
   ia12(i) = 0
   ia21(i) = 0
   ia22(i) = 0
 END DO
 
 ia11(1) = file11
 IF (file11  /= 0) CALL rdtrl (ia11)
 IF (ia11(1) < 0) ia11(1) = 0
 ia12(1) = file12
 IF (file12  /= 0) CALL rdtrl (ia12)
 IF (ia12(1) < 0) ia12(1) = 0
 ia21(1) = file21
 IF (file21  /= 0) CALL rdtrl (ia21)
 IF (ia21(1) < 0) ia21(1) = 0
 ia22(1) = file22
 IF (file22  /= 0) CALL rdtrl (ia22)
 IF (ia22(1) < 0) ia22(1) = 0
 
!     SET UP MATRIX CONTROL BLOCK FOR OUTPUT
 
 ia(1) = filea
 ia(4) = 2
 IF (rpart /= 0 .AND. cpart /= 0) ia(4) = 1
 ia(5) = MAX0(ia11(5),ia12(5),ia21(5),ia22(5))
 
!     SET UP DUMMY PARTITION VECTOR
 
 i = 0
 core(  1) = 0
 core(i+2) = 1
 core(i+3) = ia(2)
 core(i+4) = 2
 core(i+5) = 1
 core(i+6) = 0
 
 rule = 0
 lcr = korsz(core)
 
 IF (rpart == 0) GO TO 30
 IF (cpart == 0) GO TO 20
 
!     FULL MERGE
 
 ia(2) = nsub0 + nsub1
 ia(3) = ia(2)
 CALL merge (rp,cp,core)
 GO TO 40
 
!     ROW MERGE
 
 20 ia(2) = nsub0 + nsub1
 ia(3) = MAX0(ia11(3),ia12(3))
 CALL merge (rp,core,core)
 GO TO 40
 
!     COLUMN MERGE
 
 30 IF (cpart == 0) GO TO 50
 ia(2) = MAX0(ia11(2),ia21(2))
 ia(3) = nsub0 + nsub1
 CALL merge (core,cp,core)
 
!     WRITE TRIALER FOR OUTPUT
 
 40 CALL wrttrl (ia)
 
 RETURN
 
!     ILLEGAL INPUT - NO PARTITION VECTOR
 
 50 CALL mesage (-7,0,NAME)
 RETURN
END SUBROUTINE gfsmrg
