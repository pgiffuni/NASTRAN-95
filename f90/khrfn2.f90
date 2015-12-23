INTEGER FUNCTION khrfn2 (word,j,izb)
     
    !     CHARACTER FUNCTION KHRFN2 RECIEVES THE J-TH BYTE OF WORD
    !     LEFT ADJUSTED IF J IS .GE. ZERO, OR RIGHT ADJUSTED IF J .LT. ZERO
    !     ZERO FILL IF IZB IS ZERO, OTHERWISE, BLANK FILL
 
    INTEGER, INTENT(IN)                      :: word(1)
    INTEGER, INTENT(IN OUT)                  :: j
    INTEGER, INTENT(IN OUT)                  :: izb
    INTEGER :: BLANK
    COMMON /system/ dummy(40), ncpw
    DATA    BLANK / 4H      /
 
    i  = 1
    khrfn2 = BLANK
    IF (izb == 0) khrfn2 = 0
    IF (j   < 0) i = ncpw
    ij = IABS(j)
    khrfn2 = khrfn1(khrfn2,i,word(1),ij)

    RETURN
END FUNCTION khrfn2
