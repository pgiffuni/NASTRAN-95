INTEGER FUNCTION khrfn3 (word1,word2,move,idir)
     
    !     CHARACTER FUNCTION KHRFN3 MERGES TWO WORDS, WORD1 AND WORD2, BY
    !     BYTES
 
    !     (+)MOVE IS NO. OF BYTES INVOLVED PRELIMINARY SHIFTING
    !     (-)MOVE IS NO. OF BYTES IN MERGING, NO PRELIMINARY SHIFTING.
    !     IDIR IS LEFT OR RIGHT SHIFT OF WORD2. THE VACANT BYTES ARE THEN
    !     FILLED IN BY WORD1.  (LEFT SHIFT IF IDIR=1, RIGHT SHIFT OTHERWISE)
 
    !     NOTE - KHRFN3 HANDLES ONLY 4 BYTES OF WORD. IF MACHINE WORD HAS
    !     MORE THAN 4 BYTES PER WORD, KHRFN3 DOES NOT ZERO-FILL NOR BLANK-
    !     FILL THE REST OF THE WORD. THE CALLER SHOULD MAKE THE PROPER
    !     CHOICE BY ZERO-FILL OR BLANK-FILL THE INPUT WORDS, WORD1 ADN WORD2
 
    !     THE FOLLOWING TABLE GIVES THE RESULTS OF KHRFN3 FOR VARIOUS INPUT
    !     VALUES OF MOVE AND IDIR:
 
    !        GIVEN:    WORD1=ABCD  AND  WORD2=1234  (IN BCD)
    !                 IDIR=1   IDIR.NE.1            IDIR=1   IDIR.NE.1
    !                --------------------          --------------------
    !        MOVE= 0:   1234     1234      MOVE=-0:   1234     1234
    !        MOVE= 1:   234D     A123      MOVE=-1:   123D     A234
    !        MOVE= 2:   34CD     AB12      MOVE=-2:   12CD     AB34
    !        MOVE= 3:   4BCD     ABC1      MOVE=-3:   1BCD     ABC4
    !        MOVE= 4:   ABCD     ABCD      MOVE=-4:   ABCD     ABCD
 
    !     THIS ROUTINE WAS WRITTEN BY G.CHAN TO REPLACE THE ORIGINAL VAX
    !     ROUTINE WHICH WAS VERY VERY INEFFICIENT.
 
    INTEGER, INTENT(IN)                      :: word1(1)
    INTEGER, INTENT(IN)                      :: word2(1)
    INTEGER, INTENT(IN OUT)                  :: move
    INTEGER, INTENT(IN OUT)                  :: idir
    INTEGER :: word3
 
    ncpw  = 4
    imove = IABS(move)
    iend  = ncpw - imove
    word3 = word2(1)
    IF (move < 0) THEN
        GO TO    50
    ELSE IF (move == 0) THEN
        GO TO    90
    END IF
10  word3 = word1(1)
    IF (imove >= ncpw) GO TO 90
    IF (idir  ==    1) GO TO 30
    DO  i = 1,iend
        word3 = khrfn1(word3,i+imove,word2(1),i)
    END DO
    GO TO 90
    30   DO  i = 1,iend
        word3 = khrfn1(word3,i,word2(1),i+imove)
    END DO
    GO TO 90
50  IF (idir == 1) GO TO 70
    DO  i = 1,imove
        word3 = khrfn1(word3,i,word1(1),i)
    END DO
    GO TO 90
    70   DO  i = 1,imove
        word3 = khrfn1(word3,i+iend,word1(1),i+iend)
    END DO
90  khrfn3 = word3

    RETURN
END FUNCTION khrfn3
