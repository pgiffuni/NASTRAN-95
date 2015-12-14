FUNCTION maskn (l2,l1)
     
!     TO BUILD AN INTEGER MASK FOR BIT MANIPULATION
!                                                                   0 OR
!     64  60        48     36   32         <--- BIT COUNT ---          1
!      +---+---------+------+----+-------------------------------------+
!       ...   ....     ....   ...  ....00000011111111111111111111000...
!      +---+---------+------+----+-------------------------------------+
!                                            /                  /
!                                            L2                L1
 
!      BIT COUNTS FROM RIHGT (L1) TO LEFT (L2).  L1=0 IS SAME AS L1=1
 
!      E.G.      L2    L1      MASK PATTERN, WITH LEADING ZERO BITS
!               ----  ----  ------------------------------------------
!                12     0    A 12 BIT MASK, RIGHT ADJUSTED
!                24     8    A 24 BIT MASK, RIGHT ADJUSTED, WITH 8
!                            TRAILING ZERO BITS.
 
!      THIS ROUTINE IS SUITABLE FOR MACHINE WORD OF ANY BIT SIZE
!      BIT PATTERN CAN ALSO INCLUDE SIGN BIT.
!      SYSTEM MASK ROUINTE, IF IT EXISTS, IS NOT USED.
 
!      WRITTEN BY G.CHAN/UNISYS  10/1992
 
 
 INTEGER, INTENT(IN OUT)                  :: l2
 INTEGER, INTENT(IN)                      :: l1
 EXTERNAL        rshift,lshift
 INTEGER :: rshift
 
 IF (l2 < l1) CALL errtrc ('maskn   ',l2-l1)
 maskn = lshift(1,l2) - 1
 IF (l1 > 1) maskn = lshift(rshift(maskn,l1-1),l1-1)
 RETURN
END FUNCTION maskn
