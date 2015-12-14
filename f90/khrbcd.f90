SUBROUTINE khrbcd (khr80,bcd4)
     
!     MOVE ONLY THE APPROPRIATE PORTION OF THIS ROUTINE TO THE MDS GROUP
 
!     VAX, IBM, AND UNIVAC VERSION
!     ============================
 
!     THESE GROUP OF ROUTINES ARE MAINLY USED BY XREAD, RCARD2, AND
!     XRCARD IN LINK1
 
!     THESE GROUP OF ROUTINES CONVERT CHARACTER STRING (IN KHR100,
!     KHR80, KHR2), TO BCD4 ARRAY, OF 4 BYTES EACH WORD.
!     SAME OPERATION AS:
 
!           READ (KHRi,15) BCD4
!        15 FORMAT (20A4),  or (25A4),  or (2A4)
 
!     (THE READ OPERATION IS I/O BOUND, AND IS SLOW IN MOST MACHINES)
 
 
 
 CHARACTER (LEN=80), INTENT(IN OUT)       :: khr80
 INTEGER, INTENT(OUT)                     :: bcd4(2)
 INTEGER :: b4(2)
 CHARACTER (LEN=100) :: khr100,k100
 CHARACTER (LEN=80) :: k80
 CHARACTER (LEN=72) :: khr72 ,k72
 CHARACTER (LEN=8) :: khr8  ,k8
 EQUIVALENCE     (k100,k80,k72,k8,b4(1))
 
!     ROUTINE KHRBCD (KAR80,BCD4)
!     ===========================
!     A80 ---> 20A4
 
 k80=khr80
 DO  i=1,20
   bcd4(i)=b4(i)
 END DO
 RETURN
 
 
 ENTRY khrbc1 (khr100,bcd4)
!     ==========================
!     A100 ---> 25A4
 
 k100=khr100
 DO  i=1,25
   bcd4(i)=b4(i)
 END DO
 RETURN
 
 
 ENTRY khrbc2 (khr8,bcd4)
!     ========================
!     A8 ---> 2A4
 
 k8=khr8
 bcd4(1)=b4(1)
 bcd4(2)=b4(2)
 RETURN
 
 
!     THE FOLLOWING ROUTINES, BCDKHi, CONVERT FROM BCD ARRAY TO
!     CHARACTER STRING.   SAME OPERATION AS:
 
!           WRITE (KHRi,25) BCD4
!        25 FORMAT (20A4),  or (18A4), or (2A4)
 
!     WHERE KHRi IS KHR80, KHR72, OR KHR8 ACCORDINGLY
!     (THE WRITE OPERATION IS I/O BOUND, AND IS SLOW IN MOST MACHINES)
 
 
 ENTRY bcdkh8 (bcd4,khr80)
!     =========================
!     20A4 ---> A80
 
 DO  i=1,20
   b4(i)=bcd4(i)
 END DO
 khr80=k80
 RETURN
 
 
 ENTRY bcdkh7 (bcd4,khr72)
!     =========================
!     18A4 ---> A72
 
 DO  i=1,18
   b4(i)=bcd4(i)
 END DO
 khr72=k72
 RETURN
 
 
 ENTRY bcdkh2 (bcd4,khr8)
!     ========================
!     2A4 ---> A8
 
 b4(1)=bcd4(1)
 b4(2)=bcd4(2)
 khr8=k8
 RETURN
!     END
 
 
 
!     SUBROUTINE KHRBCD (KHR,BCD4)
 
!     CDC VERSION
!     ===========
 
!     THIS GROUP OF ROUTINES ARE CALLED BY XREAD, RCARD2, AND XRCARD
 
!     THESE GROUP OF ROUTINES CONVERT CHARACTER STRINGS TO BCD ARRAY,
!     4 BYTES EACH WORD, AND VISE VERSA.     SIMILARY TO -
 
!        METHOD 1:
!        --------
!           READ (KHR80,10) BCD4     and       WRITE (KHR72,20) BCD4
!        10 FORMAT (20A4)                   20 FORMAT (18A4)
 
!        METHOD 2:
!        --------
!           I2=0                     and       I2=0
!           DO 10 I=1,NWDS                     DO 20 I=1,NWDS
!           I1=I2+1                            I1=I2+1
!           I2=I2+4                            I2=I2+4
!        10 BCD4(I)(1:4)=KHR(I1:I2)         20 KHR(I1:I2)=BCD4(I)(1:4)
 
!     HOWEVER THE INTERNAL-FILE READ AND WRITE (METHOD 1) AND THE
!     CHARACTER MANIPULATION (METHOD 2) ARE EXTREMELY SLOW IN CDC
!     (METHOD 1 IS ABOUT 18 TIMES SLOWER THAN SHIFT/AND/OR OPERATIONS
!     THAT ACCOMPLISH THE SAME THING. METHOD 2 IS 2 TO 4 TIMES SLOWER)
 
!     THE CALLING ROUTINES ACTUALLY PASS THE KHR ARGUMENTS IN CHARACTER
!     STRINGS (CHARACTER*100, CHARACTER*80, CHARACTER*2), WHEREAS, THEY
!     ARE PICKED UP HERE IN THIS ROUTINE AS INTEGER-BCD ARRAYS, 10 BYTES
!     EACH WORDS. ONLY THE FIRST 4 BYTES ARE USED IN NASTRAN.
 
!     (THE FOLLOWING CODE ASSUMES NO BREAK ON THE 1ST AND 4TH BCD WORDS
!     IN A GROUP OF 5)
 
!     INPUT  - KHR  = CHARACTER STRING IN CHARACTER*80, CHARACTER*100,
!                     AND CHARACTER*2
!     OUTPUT - BCD4 = BCD ARRAYS (OF DIMENSION NWDS)
!                     EACH BCD4 WORD HOLDS ONLY 4 BYTES OF DATA
 
!     INTEGER       BCD4(1),KHR(1),BLANK,BLK90
!                                     1 2 3 4 5 6 7 8 9 10
!     DATA          M1234,M12,M34 / O"77777777000000000000",
!    1                              O"77770000000000000000",
!    2                              O"00007777000000000000"/
!     DATA          M5678,M90     / O"00000000777777770000",
!    1                              O"00000000000000007777"/
!     DATA          M3456,M7890   / O"00007777777700000000",
!    1                              O"00000000000077777777"/
!     DATA          BLANK,BLK90   / O"00000000555555555555",
!    1                              O"00000000000000005555"/
!                                     1 2 3 4 5 6 7 8 9 10
 
!     SUBROUTINE KHRBCD (KHR,BCD4)
!     ============================
!     A80 ----> 20A4
 
!     NWDS = 20
!     GO TO 40
 
 
!     ENTRY KHRBC1 (KHR,BCD4)
!     =======================
!     A100 ----> 25A4
 
!     NWDS = 25
!     GO TO 40
 
 
!     ENTRY KHRBC2 (KHR,BCD4)
!     =======================
!     A8 ----> 2A4
 
!     NWDS = 2
 
!40   I   =-5
!     II  = 0
!50   I   = I+5
!     IF (I .GE. NWDS) GO TO 80
!     II  = II+1
!     NW1 = KHR(II)
!     NN  = AND(NW1,M1234)
!     BCD4(I+1) = OR(NN,BLANK)
!     NW1 = SHIFT(NW1,24)
!     NN  = AND(NW1,M1234)
!     BCD4(I+2) = OR(NN,BLANK)
!     IF (I+2 .GE. NWDS) GO TO 80
!     NW1 = SHIFT(NW1,24)
!     II=II+1
!     NW2 = KHR(II)
!     IF (I+3 .LT. NWDS) GO TO 60
!     NW2 = SHIFT(NW2,-12)
!     GO TO 70
!60   NW2 = SHIFT(NW2,12)
!     NN  = AND(NW2,M1234)
!     BCD4(I+4) = OR(NN,BLANK)
!     NW2 = SHIFT(NW2,24)
!     NN  = AND(NW2,M1234)
!     BCD4(I+5) = OR(NN,BLANK)
!     NW2 = SHIFT(NW2,12)
!70   NW1 = AND(NW1,M12)
!     NW2 = AND(NW2,M34)
!     NN  = OR(NW1,NW2)
!     BCD4(I+3) = OR(NN,BLANK)
!     GO TO 50
!80   CONTINUE
!     GO TO 140
 
 
!     ENTRY BCDKH8 (BCD4,KHR)
!     =======================
!     20A4 ----> A80
 
!     INPUT  - BCD4 = BCD ARRAYS (OF DIMENSION NWDS). BCD DATA ARE IN
!                     A4 FORMAT
!     OUTPUT - KHR  = CHARACTER STRING IN CHARACTER*80, CHARACTER*100,
!                     AND CHARACTER*2
 
!     NWDS = 20
!     GO TO 100
 
!     ENTRY BCDKH7 (BCD4,KHR)
!     =======================
!     18A4 ----> A72
 
!     NWDS = 18
!     GO TO 100
 
!     ENTRY BCDKH2 (BCD4,KHR)
!     =======================
!     2A4 ----> A8
 
!     NWDS = 2
 
!100  I   =-5
!     II  = 0
!110  I   = I+5
!     IF (I .GE. NDWS) GO TO 140
!     II  = II+1
!     NW1 = AND(BCD4(I+1),M1234)
!     NW2 = SHIFT(BCD4(I+2),-24)
!     NW2 = AND(NW2,M5678)
!     KHR(II) = OR(NW1,NW2)
!     IF (I+2 .GE. NWDS) GO TO 120
!     NW3 = SHIFT(BCD4(I+3),12)
!     NW2 = AND(NW3,M90)
!     KHR(II) = OR(KHR(II),NW2)
!     II  = II+1
!     KHR(II) = AND(NW3,M12)
!     NW3 = SHIFT(BCD4(I+4),-12)
!     NW1 = AND(NW3,M3456)
!     NM3 = SHIFT(BCD4(I+5),-36)
!     NW2 = AND(NW3,M7890)
!     NW3 = OR(NW1,NW2)
!     KHR(II) = OR(KHR(II),NW3)
!     GO TO 110
!120  KHR(II) = OR(KHR(II),BLK90)
 
!140  RETURN
!     END
 
 
!     SUBROUTINE KHRBCD (KHR80,BCD4)
 
!     64-BIT MACHINE, UNIX VERSION
!     ============================
 
!     THIS GROUP OF ROUTINES ARE CALLED BY XREAD, RCARD2, AND XRCARD
 
!     CHARACTER*100 KHR100, KDUM
!     CHARACTER*80  KHR80 , K80
!     CHARACTER*72  KHR72 , K72
!     CHARACTER*8   KHR8  , K8,   BCD4(1)
 
!     EQUIVALENCE   (KDUM,K100,K80,K72,K8)
 
 
!     SUBROUTINE KHRBCD (KHR80,BCD4)
!     ==============================
!     A80 ----> 20A4
 
!     K80  = KHR80
!     NWDS = 20
!     GO TO 100
 
 
!     ENTRY KHRBC1 (KHR100,BCD4)
!     ==========================
!     A100 ----> 25A4
 
!     KDUM = KHR100
!     NWDS = 25
!     GO TO 100
 
 
!     ENTRY KHRBC2 (KHR8,BCD4)
!     ========================
!     A8 ----> 2A4
 
!     K8   = KHR8
!     NWDS = 2
 
! 100 I2 = 0
!     DO 200 I = 1,NWDS
!     I1 = I2 + 1
!     I2 = I2 + 4
!     BCD4(I) = KDUM(I1:I2)
! 200 CONTINUE
!     GO TO 800
 
 
!     ENTRY BCDKH8 (BCD4,KHR80)
!     =========================
!     20A4 ----> A80
 
!     NWDS  = 20
!     GO TO 300
 
 
!     ENTRY BCDKH7 (BCD4,KHR72)
!     =========================
!     18A4 ----> A72
 
!     NWDS = 18
!     GO TO 300
 
 
!     ENTRY BCDKH2 (BCD4,KHR8)
!     ========================
!     2A4 ----> A8
 
!     NWDS = 2
 
! 300 I2 = 0
!     DO 400 I = 1,NWDS
!     I1 = I2 + 1
!     I2 = I2 + 4
!     KDUM(I1:I2) = BCD4(I)
! 400 CONTINUE
!     IF (NWDS-18) 500,600,700
 
! 500 KHR8 = K8
!     GO TO 800
 
! 600 KHR72 = K72
!     GO TO 800
 
! 700 KHR80 = K80
! 800 RETURN
 
END SUBROUTINE khrbcd
