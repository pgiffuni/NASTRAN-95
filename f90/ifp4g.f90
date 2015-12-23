SUBROUTINE ifp4g (ibit,FILE)
     
    !     TURNS ON BIT -IBIT- IN TRAILER FOR DATA BLOCK -FILE-

    INTEGER, INTENT(IN)                      :: ibit
    INTEGER, INTENT(IN)                      :: FILE
    EXTERNAL    orf
    INTEGER :: orf, trail(7), two
    COMMON/two/ two(32)
 
    trail(1) = FILE
    CALL rdtrl (trail)
    i1 = (ibit-1)/16 + 2
    i2 = ibit - (i1-2)*16 + 16
    trail(i1) = orf(trail(i1),two(i2))
    trail(1) = FILE
    CALL wrttrl (trail)

    RETURN
END SUBROUTINE ifp4g
