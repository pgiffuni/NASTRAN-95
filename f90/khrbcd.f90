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
