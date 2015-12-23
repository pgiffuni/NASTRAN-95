INTEGER FUNCTION khrfn1 (word1,i,word2,j)
     
    !     CHARACTER-FUNCTIONS 1,2,3,4, AND 5 WERE WRITTEN BY G.CHAN/UNISYS
    !     TO STANDARDIZE NASTRAN BCD-WORD BYTE PROCESSING.
 
    !     NOTE - THE INPUT WORD(S) ARE INTEGERS OR REALS, HOLDING BCD TYPE
    !            DATA. (NOT CHARACTER)
    !            BYTE COUNTS FROM LEFT TO RIGHT
 
    !     THESE FIVE CHARACTER FUNCTIONS ARE COMPLETELY MACHINE INDEPENDENT
 
    !     KHRFN1 REPLACES THE I-TH BYTE OF WORD 1 BY THE J-TH BYTE OF WORD2
    !     E.G.   WORD1=ABCD,    WORD2=1234
    !            KHRFN1(WORD1,3,WORD2,2) GIVES  AB2D
 
    !     ABSOLUTE VALUES OF I AND J ARE USED
 
    !     THE CODE BELOW WORKS WITH ALL MACHINES.  HOWEVER, SEE THE
    !     SIMPLIFIED VERSION FURTHER DOWN.
 
    !     INTEGER      WORD1(1),WORD2(1),TEMP(2)
    !     CHARACTER*8  TEMP8
 
    !     TEMP(1) = WORD1(1)
    !     TEMP(2) = WORD2(1)
    !     CALL BCDKH2 (TEMP,TEMP8)
    !     II = IABS(I)
    !     JJ = IABS(J) + 4
    !     TEMP8(II:II) = TEMP8(JJ:JJ)
    !     CALL KHRBC2 (TEMP8,TEMP)
    !     KHRFN1 = TEMP(1)
 
    !     SIMPLIFIED VERSION
 
    !     FOR MACHINES (CDC, IBM, VAX, AND GRAY) THAT ALLOW EQUIVALENCE
    !     BETWEEN CHARACTERS AND INTEGER VARIABLES, THE FOLLOWING SIMPLIFIED
    !     CODE CAN BE USED.
 
    CHARACTER (LEN=4) :: tempc1,tempc2
    !     CHARACTER*n  TEMPC1,TEMPC2
    !        (WHERE n is 10 for CDC, 8 for 64-BIT UNIX and
    !                     4 for VAX and IBM)
    INTEGER, INTENT(IN)                      :: word1(1)
    INTEGER, INTENT(IN OUT)                  :: i
    INTEGER, INTENT(IN)                      :: word2(1)
    INTEGER, INTENT(IN OUT)                  :: j
    INTEGER :: temp1,temp2

    EQUIVALENCE (temp1,tempc1), (temp2,tempc2)
 
    temp1 = word1(1)
    temp2 = word2(1)
    ii = IABS(i)
    jj = IABS(j)
    tempc1(ii:ii) = tempc2(jj:jj)
    khrfn1 = temp1
    RETURN

END FUNCTION khrfn1
