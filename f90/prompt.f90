SUBROUTINE prompt
     
!     DRIVER FOR  INTERACTIVE MODULE - PROMPT
 
!        PROMPT    //S,N,PEXIT/S,N,PLOT1/S,N,PLOT2/S,N,XYPLOT/
!                    S,N,SCAN1/S,N,SCAN2/DUM1/DUM2/DUM3/DUM4 $
 
 IMPLICIT INTEGER  (a-z)
 COMMON /BLANK / param(10)
 COMMON /system/ ksysm(100)
 EQUIVALENCE     (nout ,ksysm(2)),  (soln ,ksysm(22)),  &
     (in   ,ksysm(4)),  (intra,ksysm(86))
 DATA    p,s,c,b / 1HP, 1HS, 1HC, 1H /
 
 intra=IABS(intra)
 nout = 6
 DO  i=1,10
   param(i)=0
 END DO
 
 20   WRITE (nout,110)
 READ (in,130,ERR=20) j
 IF (j < 1 .OR. j > 6) GO TO 20
 IF (soln == 3 .AND. (j == 4 .OR. j == 6)) GO TO 70
 param(j)=-1
 IF (j == 1) RETURN
 40   WRITE (nout,120)
 i=b
 READ (in,140,END=50) i
 50   IF (i == b) RETURN
 IF (i == c) GO TO 60
 IF (i /= p .AND. i /= s) GO TO 40
 IF (i == s) intra = MOD(intra,10)
 IF (i == p) intra = MOD(intra,10) + 10
 RETURN
 
 60   param(j)=0
 GO TO 20
 70   WRITE (nout,80) j,soln
 80   FORMAT (/,' OPTION',i3,' IS NOT AVAILABLE FOR SOLUTION',i3)
 GO TO 20
 
 110  FORMAT (9X,'1. EXIT', /9X,'2. STRUCTURE PLOTS - UNDEFORMED',  &
     /9X,'3. STRUCTURE PLOTS - DEFORMED', /9X,'4. XYPLOTS',  &
     /9X,'5. SCAN OUTPUT - SORT1', /9X,'6. SCAN OUTPUT - SORT2',  &
     //9X,'SELECT ONE OPTION FROM MENU -')
 120  FORMAT (/9X,'OUTPUT TO SCREEN, OR TO PRINTFILE, OR CANCEL OPTION (  &
     S/P/C) -')
 130  FORMAT (i1)
 140  FORMAT (a1)
END SUBROUTINE prompt
