SUBROUTINE pexit
     
 EXTERNAL        link
 INTEGER :: hh,ss,date(3)
 COMMON /output/ LE(17)
 COMMON /machin/ mach
 COMMON /msgx  / nmsg
 COMMON /resdic/ irdict
 COMMON /system/ isystm(100)
 EQUIVALENCE     (isystm( 2),nout  ), (isystm(76),nosbe),  &
     (isystm(82),icpflg), (isystm(15),date  )
 
!     SEE IF ANY MESSAGES ARE IN THE QUEUE
 
 IF (nmsg   > 0) CALL msgwrt
 IF (icpflg /= 0) WRITE (irdict,10)
10 FORMAT ('$ END OF CHECKPOINT DICTIONARY')

!     JOB DONE. PRINT LAST 4 MESSAGE LINES

CALL waltim (i)
hh = i/3600
mm = (i-hh*3600)/60
ss = i - hh*3600 - mm*60
CALL cputim (i,t,0)
IF (mach == 4) i = t
IF (LE(1) == -1 .AND. LE(2) == -1) GO TO 70
WRITE  (nout,20) LE,date,hh,mm,ss
20 FORMAT (////40X,'* * * END OF JOB * * *', /1H1, /,' JOB TITLE = ',  &
17A4, /,' DATE:',i3,1H/,i2,1H/,i2, /,' END TIME:',i3,1H:,  &
    i2,1H:,i2)

!     CDC TOTAL CPU TIME IS A BIG NUMBER. DON'T PRINT IT

IF (mach == 4 .OR. LE(1) == -1) GO TO 50
IF (mach <= 5) WRITE (nout,30) i
IF (mach > 5) WRITE (nout,40) i
30 FORMAT (' TOTAL CPU TIME',i6,' SEC.')
40 FORMAT (' TOTAL WALL CLOCK TIME',i7,' SEC.')

!     FLUSH O/P BUFFERS

50 WRITE  (nout,60)
60 FORMAT (1H )

IF (mach == 4 .AND. nosbe > 0) CALL link (-1,nosbe,1)
GO TO 90

70 j = 5
IF (LE(9) >= 0) j = 3
WRITE  (nout,80) (LE(i),i=j,8)
80 FORMAT (//1X,6A4)

90 CONTINUE
CALL dbmstf
DO  i = 1,4
  CLOSE ( i )
END DO
DO  i = 7,22
  CLOSE ( i )
END DO
!WKBR 8/94 SUN  CALL EXIT
CALL EXIT( 0 )
END SUBROUTINE pexit
