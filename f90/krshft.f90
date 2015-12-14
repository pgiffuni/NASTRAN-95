FUNCTION krshft (iword,n)
     
!     CHARACTER FUNCTION KRSHFT AND KLSHFT PERFORM LEFT AND RIGHT
!     SHIFTS, BY N CHARACTERS (BYTES).
!     EMPTY BYTES ARE ZERO FILLED.
 
!     NORMALLY, KRSHFT AND KLSHFT WORK ALMOST LIKE RSHIFT AND LSHFIT
!     RESPECTIVELY, EXCEPT THEY MOVE DATA BY BYTE COUNT, NOT BY BITS.
!     HOWEVER, IF THE MACHINE STORES THE BCD WORD DATA IN REVERSE ORDER
!     (SUCH AS VAX AND SILICON GRAPHICS), KRSHFT IS EQUIVALENCED TO
!     LSHFIT, AND KLSFHT TO RSHIFT.
 
 
 INTEGER, INTENT(IN OUT)                  :: iword(1)
 INTEGER, INTENT(IN)                      :: n
 EXTERNAL        lshift,   rshift
 INTEGER :: rshift
 COMMON /machin/ mac(3),   lqro
 COMMON /system/ dummy(38),nbpc
 
 IF (MOD(lqro,10) == 1) GO TO 10
 krshft = rshift(iword(1),n*nbpc)
 RETURN
 10 krshft = lshift(iword(1),n*nbpc)
 RETURN
 
 ENTRY klshft (iword,n)
!     ======================
 
 IF (MOD(lqro,10) == 1) GO TO 20
 klshft = lshift(iword(1),n*nbpc)
 RETURN
 20 klshft = rshift(iword(1),n*nbpc)
 RETURN
END FUNCTION krshft
