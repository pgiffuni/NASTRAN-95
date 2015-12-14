SUBROUTINE cpystr (inblk,outblk,flag,col)
     
    !     CPYSTR COPIES A LOGICAL RECORD WRITTEN IN STRING FORMAT
    !     FROM ONE FILE TO ANOTHER FILE.
 
    !     INBLK  = 15-WORD STRING COMMUNICATION BLOCK FOR INPUT FILE
    !     OUTBLK = 15-WORD STRING COMMUNICATION BLOCK FOR OUTPUT FILE
    !     FLAG .NE. 0 MEANS 1ST CALL GETSTR HAS BEEN MADE FOR THE RECORD
    !          .EQ. 0 MEANS 1ST CALL GETSTR HAS NOT BEEN MADE
    !     COL .EQ. 0 MEANS COLUMN NUMBER IS IN INBLK(12)
    !         .NE. 0 MEANS COL IS COLUMN NUMBER
 
 
 
    INTEGER, INTENT(OUT)                     :: inblk(15)
    INTEGER, INTENT(OUT)                     :: outblk(15)
    INTEGER, INTENT(IN OUT)                  :: flag
    INTEGER, INTENT(IN)                      :: col
    INTEGER :: prc,words,rlcmpx,TYPE, prec,rc,out
    DOUBLE PRECISION :: xnd(1)
    COMMON /zzzzzz/  xns(1)
    COMMON /TYPE  /  prc(2),words(4),rlcmpx(4)
    EQUIVALENCE      (xns(1),xnd(1))
 
    !     ON OPTION, MAKE 1ST CALL TO GETSTR AND THEN INITIALIZE
 
    IF (flag /= 0) GO TO 10
    inblk(8) = -1
    CALL getstr (*50,inblk)
10  outblk(2) = inblk(2)
    outblk(3) = inblk(3)
    outblk(4) = inblk(4)
    outblk(8) = -1
    outblk(12) = col
    IF (col == 0) outblk(12) = inblk(12)
    outblk(13)= 0
    TYPE = inblk(2)
    prec = prc(TYPE)
    rc   = rlcmpx(TYPE)
 
    !     COPY A STRING
 
12  CALL putstr (outblk)
    nprev = 0
    outblk(7) = MIN0(inblk(6),outblk(6))
14  in   = inblk(5)
    out  = outblk(5)
    nstr = out + rc*(outblk(7) - nprev) - 1
    IF (prec == 2) GO TO 18
 
    DO  jout = out,nstr
        xns(jout) = xns(in)
        in = in + 1
    END DO
    GO TO 20
 
    18 DO  jout = out,nstr
        xnd(jout) = xnd(in)
        in = in + 1
    END DO
 
    !     TEST FOR END OF INPUT STRING(S)
 
20  IF (outblk(7) == inblk(6)+nprev) GO TO 30
    outblk(13) = outblk(13) + outblk(7)
    CALL endput (outblk)
    outblk(4) = outblk(4) + outblk(7)
    inblk(6)  = inblk(6)  - (outblk(7) - nprev)
    inblk(5)  = in
    GO TO 12
 
    !     INPUT STRING HAS BEEN COPIED.  GET ANOTHER STRING.
 
30  CALL endget (inblk)
    CALL getstr (*40,inblk)
 
    !     TEST FOR STRING CONTIGUOUS WITH PREVIOUS STRING.
    !     IF SO, AND IF TERMS AVAILABLE, CONCATENATE WITH PREVIOUS STRING.
 
    IF (inblk(4) /= outblk(4)+outblk(7)) GO TO 35
    IF (outblk(7) >= outblk(6)) GO TO 35
    outblk(5) = nstr + 1
    nprev     = outblk(7)
    outblk(7) = MIN0(outblk(7)+inblk(6),outblk(6))
    GO TO 14
35  outblk(13) = outblk(13) + outblk(7)
    CALL endput (outblk)
    outblk(4) = inblk(4)
    GO TO 12
 
    !     NO MORE STRINGS -  CLOSE RECORD AND RETURN
 
40  outblk(8) = 1
    CALL endput (outblk)
    outblk(13) = (outblk(13)+outblk(7))*words(TYPE)
    RETURN
 
    !     HERE IF NO STRINGS IN RECORD - MAKE A NULL RECORD
 
50  outblk(2) = 1
    outblk(3) = 0
    outblk(8) = -1
    CALL putstr (outblk)
    outblk(8) = 1
    CALL endput (outblk)

    RETURN
END SUBROUTINE cpystr
