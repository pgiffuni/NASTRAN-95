SUBROUTINE na12a8 (*,a,n,b,notuse)

!     THESE ROUTTNES CONVERT N A1 BCD WORDS IN A, OR N A1 CHARACTERS IN
!     C TO AN 8-BYTE BCD WORD IN B (CDC ONLY), (OR TO TWO 4-BYTE BCD
!     WORDS IN B, ALL OTHER NON-CDC MACHINES), OR AN 8-CHARACTER WORD
!     IN D, LEFT ADJUSTED.
!     CALLING ROUTINE MUST NOT USE LOGICAL*1 FOR A-ARRAY.
!     (NO SYSTEM ENCODE/DECODE FUNCTIONS ARE USED)
 
!     ENTRY POINTS   NA1 2 A8  (BCD-BYTE  VERSION)
!                    NK1 2 K8  (CHARACTER VERSION)
 
 
!     WRITTEN BY G.CHAN/SPERRY IN AUG. 1985
!     PARTICULARLY FOR XREAD ROUTINE, IN SUPPORT OF ITS NEW FREE-FIELD
!     INPUT FORMAT.  THIS SUBROUTINE IS MACHINE INDEPENDENT
 
!     LAST REVISED  8/1988
     
 INTEGER, INTENT(IN OUT)                  :: a(1)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN OUT)                  :: b(2)
 INTEGER, INTENT(IN)                      :: notuse
 INTEGER :: cdc
 CHARACTER (LEN=1) :: c(1),     t(8)
 CHARACTER (LEN=8) :: d(1)
 CHARACTER (LEN=10) :: blnk,     temp
 COMMON /xreadx/ nout
 COMMON /machin/ mach
 EQUIVALENCE     (t(1),temp)
 DATA            blnk / '          '  /, cdc / 4  /

 IF (n > 8) GO TO 40
 temp = blnk
 CALL b2k (a,temp,n)
 IF (mach /= cdc) CALL khrbc2 (temp,b(1))
!WKBD IF (MACH .EQ. CDC) B(1) = ISWAP(TEMP)

 RETURN
 
 ENTRY nk12k8 (*,c,n,d,notuse)
!     ===============================
 
 IF (n > 8) GO TO 40
 temp = blnk
 DO  i = 1,n
   t(i) = c(i)
 END DO
 d(1) = temp
 RETURN
 
 40   WRITE  (nout,50) n
 50   FORMAT ('   N.GT.8/NA12A8',i6)
 j = notuse
 
 RETURN 1
END SUBROUTINE na12a8
