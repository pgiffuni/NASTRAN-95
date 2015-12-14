SUBROUTINE trht1b(iof,delta)
     
 
 
 INTEGER, INTENT(IN OUT)                  :: iof
 REAL, INTENT(IN)                         :: delta
 DOUBLE PRECISION :: qblock(6),          BLOCK(2), blk(2)
 
 INTEGER :: mcb(7),   NAME(2),  iqblk(12),iblock(11)
 
 COMMON /BLANK/    beta,     tabs,     norad,    radlin
 COMMON /trhtx /         ik(7),    ib(7),    icr1,     icr2,  &
     icr3,     icr4,     icr5,     isym, icr6,     icr7
 
 EQUIVALENCE             ( iqblk(1),         qblock(1) )
 EQUIVALENCE             ( iqblk(2),         iblock(1) )
 EQUIVALENCE             ( qblock(2),        BLOCK(1) )
 EQUIVALENCE             ( qblock(5),        blk(1) )
 
! ----------------------------------------------------------------------
 
 iblock(1) =2
 BLOCK(1)= 1.0D0/delta
 BLOCK(2)= 0.0D0
 iblock(7)= 2
 blk(1) = beta
 blk(2) = 0.0D0
 CALL ssg2c(ib,ik,icr6,1,iblock)
 mcb(1)=icr6
 CALL rdtrl(mcb(1))
 IF ( mcb(4) == 6) GO TO  10
 CALL factru(*40,icr6,icr1,icr2,icr3,icr4,icr7)
 isym = 0
 GO TO  20
 
!     SYMMETRIC DECOMP
 
 10 CALL factor( icr6, icr1, icr2, icr3, icr4, icr7 )
 isym =1
 
!     LLL  IS ON ICR1
 
!     FORM  A  MATRIX
 
 20 blk(1) = -(1.0D0 - beta)
 blk(2) = 0.0
 CALL ssg2c(ib,ik,icr6,1,iblock)
 30 RETURN
 40 CALL mesage(-5,icr6,NAME)
 GO TO 30
END SUBROUTINE trht1b
