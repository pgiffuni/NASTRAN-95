SUBROUTINE bitpat (icode,ibits)
     
    !     THE PURPOSE OF THIS ROUTINE IS TO TRANSFORM THE DOF WORD INTO ITS
    !     NASTRAN DIGITAL REPRESENTATION.
 
 
    INTEGER, INTENT(IN OUT)                  :: icode
    INTEGER, INTENT(OUT)                     :: ibits(2)
    EXTERNAL        orf
    INTEGER :: list(32), orf,INT(9)
    COMMON /system/ junk(38),nbpc,nbpw
    DATA    iblank/ 4H    /
    DATA    INT   / 1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9 /
 
    ibits(1) = iblank
    ibits(2) = iblank
 
    CALL decode (icode,list,n)
    IF (n == 0) RETURN
 
    j = 1
    nbits = -nbpc
    DO  i = 1,n
        nbits = nbits + nbpc
        ia = list(i)  + 1
        k  = nbpw - nbits
        ibits(j) = klshft(krshft(ibits(j),k/nbpc),k/nbpc)
        ibits(j) = orf(ibits(j),krshft(INT(ia),nbits/nbpc))
        IF (i /= 4) CYCLE
        j = 2
        nbits = -nbpc
    END DO
    RETURN
END SUBROUTINE bitpat
