SUBROUTINE ofppnt (out,nwds,FMT)
     
!WKBD LOGICAL         DEBUG
!WKBR INTEGER         OUT(NWDS), FMT(300)
 
 INTEGER, INTENT(IN OUT)                  :: out(nwds)
 INTEGER, INTENT(IN OUT)                  :: nwds
 CHARACTER (LEN=1), INTENT(OUT)           :: FMT(1200)
 
 
!WKBI
 COMMON /machin/ machx
 COMMON /system/ sysbuf,    l
!WKBD DATA    DEBUG / .FALSE. /
 
!WKBD IF (DEBUG) WRITE (L,10) (FMT(K),K=1,32)
 10   FORMAT (' FMT=',32A4)
!WKBR 5/95     IF ( MACHX.EQ.2 .OR. MACHX.EQ.5  )
 IF ( machx == 2 .OR. machx == 5 .OR. machx == 21  )  &
     WRITE (l,FMT,IOSTAT=iosxx) (out(k),k=1,nwds)
!WKBR 5/95      IF ( MACHX.NE.2 .AND. MACHX.NE.5 )
 IF ( machx /= 2 .AND. machx /= 5 .AND. machx /= 21 )  &
     CALL forwrt (FMT, out, nwds)
 RETURN
END SUBROUTINE ofppnt
