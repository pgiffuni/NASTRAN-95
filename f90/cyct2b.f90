SUBROUTINE cyct2b (INPUT,outpt,ncol,iz,mcb)
     
    !     THE PURPOSE OF THIS SUBROUTINE IS TO COPY NCOL COLUMNS FROM
    !     INPUT TO OUTPUT USING CORE AT IZ -- MCB IS THE TRAILER
 
 
    INTEGER, INTENT(IN OUT)                  :: INPUT
    INTEGER, INTENT(IN OUT)                  :: outpt
    INTEGER, INTENT(IN)                      :: ncol
    INTEGER, INTENT(OUT)                     :: iz(4)
    INTEGER, INTENT(IN OUT)                  :: mcb(7)
 
    COMMON /unpakx/ itc,iik,jjk,incr1
    COMMON /packx / ita,itb,ii,jj,incr
    EQUIVALENCE     (zero,izero)
    DATA    zero  / 0.0 /
 
 
    ita = IABS(itc)
    itb = ita
    incr= incr1
    DO  i = 1,ncol
        iik = 0
        CALL unpack (*20,INPUT,iz)
        ii  = iik
        jj  = jjk
10      CALL pack (iz,outpt,mcb)
        CYCLE
   
        !     NULL COLUMN
   
20      ii = 1
        jj = 1
        iz(1) = izero
        iz(2) = izero
        iz(3) = izero
        iz(4) = izero
        GO TO 10
    END DO
 
    RETURN
END SUBROUTINE cyct2b
