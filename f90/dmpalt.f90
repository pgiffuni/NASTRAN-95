!*DECK,DMPALT
     
 SUBROUTINE dmpalt (isize, iopen      , iptape)
!******************************************************************
!                              NOTICE                             *
!                              ------                             *
!                                                                 *
!     THIS PROGRAM BELONGS TO RPK CORPORATION.  IT IS CONSIDERED  *
!  A TRADE SECRET AND IS NOT TO BE DIVULGED OR USED BY PARTIES    *
!  WHO HAVE NOT RECEIVED WRITTEN AUTHORIZATION FROM RPK.          *
!******************************************************************
 
 
 INTEGER, INTENT(IN)                      :: isize
 INTEGER, INTENT(IN)                      :: iopen(1)
 INTEGER, INTENT(IN OUT)                  :: iptape
 INTEGER :: altfil, oldalt, xalter(2)
 
 DIMENSION ialter(2), isubr(2), icard(18)
 
 COMMON /altrxx/ altfil, newalt, nogo
 COMMON /system/ isysbf, nout
 COMMON /xrgdxx/ idum(115), nmdmap
 
 DATA oldalt /0/
 DATA isubr  /4HDMPA, 4HLT  /
 DATA xalter /4HXALT, 4HER  /
 
 ipoint = 1
 CALL WRITE (iptape, xalter, 2, 1)
 IF (newalt == 0) GO TO 200
 n2dmap = 2*nmdmap
 CALL skpfil (altfil, 1)
 CALL READ (*2000, *100, altfil, iopen, isize, 1, iend)
 ireqd = n2dmap - isize
 CALL mesage (-8, ireqd, isubr)
 100 IF (iend /= n2dmap) GO TO 2100
 ipoint = ipoint + iend
 CALL REWIND (altfil)
 
 200 CALL READ (*3000, *300, altfil, iopen(ipoint), 19, 1, iflag)
 nwords = 19
 logic = 200
 GO TO 2200
 
 300 IF (iflag /= 2) GO TO 400
 logic = 300
 ASSIGN 320 TO jgoto
 GO TO 1200
 
!     PROCESS ALTER CARDS HERE
 
 320 ialter(1) = iopen(ipoint  )
 ialter(2) = iopen(ipoint+1)
 IF (ialter(2) == 0 .OR. ialter(2) >= ialter(1)) GO TO 1000
 itemp = ialter(2)
 ialter(2) = ialter(1)
 ialter(1) = itemp
 GO TO 1000
 
 400 IF (iflag /= 4) GO TO 500
 nwords = 4
 logic = 400
 IF (newalt == 0) GO TO 2200
 logic = 410
 ASSIGN 420 TO jgoto
 GO TO 1200
 
!     PROCESS INSERT CARDS HERE
 
 420 idmap1 = iopen(ipoint  )
 idmap2 = iopen(ipoint+1)
 ioccur = iopen(ipoint+2)
 ioffst = iopen(ipoint+3)
 ASSIGN 470 TO igoto
 GO TO 1500
 470 ialter(1) = inumbr
 ialter(2) = 0
 GO TO 1000
 
 500 IF (iflag /= 5) GO TO 600
 nwords = 5
 logic = 500
 IF (newalt == 0) GO TO 2200
 logic = 510
 ASSIGN 520 TO jgoto
 GO TO 1200
 
!     PROCESS DELETE CARDS WITH ONE FIELD HERE
 
 520 idmap1 = iopen(ipoint  )
 idmap2 = iopen(ipoint+1)
 ioccur = iopen(ipoint+2)
 ioffst = iopen(ipoint+3)
 icheck = iopen(ipoint+4)
 logic = 520
 IF (icheck /= 0) GO TO 2300
 ASSIGN 570 TO igoto
 GO TO 1500
 570 ialter(1) = inumbr
 ialter(2) = inumbr
 GO TO 1000
 
 600 IF (iflag /= 9) GO TO 700
 nwords = 9
 logic = 600
 IF (newalt == 0) GO TO 2200
 logic = 610
 ASSIGN 620 TO jgoto
 GO TO 1200
 
!     PROCESS DELETE CARDS WITH TWO FIELDS HERE
 
 620 idmap1 = iopen(ipoint  )
 idmap2 = iopen(ipoint+1)
 ioccur = iopen(ipoint+2)
 ioffst = iopen(ipoint+3)
 icheck = iopen(ipoint+4)
 logic = 620
 IF (icheck /= 1) GO TO 2300
 ASSIGN 670 TO igoto
 GO TO 1500
 670 ialter(1) = inumbr
 idmap1 = iopen(ipoint+5)
 idmap2 = iopen(ipoint+6)
 ioccur = iopen(ipoint+7)
 ioffst = iopen(ipoint+8)
 ASSIGN 690 TO igoto
 GO TO 1500
 690 ialter(2) = inumbr
 IF (ialter(2) == 0 .OR. ialter(2) >= ialter(1)) GO TO 1000
 itemp = ialter(2)
 ialter(2) = ialter(1)
 ialter(1) = itemp
 GO TO 1000
 
!     PROCESS DMAP STATEMENTS HERE
 
 700 nwords = iflag
 logic = 700
 IF (iflag /= 18) GO TO 2200
 CALL WRITE (iptape, iopen(ipoint), 18, 1)
 GO TO 200
 
!     WRITE ALTER CONTROL DATA ON THE NEW PROBLEM TAPE
 
 1000 IF (ialter(1) == 0) GO TO 200
 IF (ialter(1) > oldalt) GO TO 1100
 1050 nogo = 1
 WRITE (nout, 3300) icard
 GO TO 200
 1100 IF (ialter(2) /= 0 .AND. ialter(1) > ialter(2)) GO TO 1050
 CALL WRITE (iptape, ialter, 2, 1)
 oldalt = ialter(1)
 IF (ialter(2) /= 0) oldalt = ialter(2)
 GO TO 200
 
!     INTERNAL SUBROUTINE TO READ IN AN ALTER CONTROL CARD IMAGE
 
 1200 CALL READ (*2000, *1230, altfil, icard, 19, 1, iflag)
 nwords = 19
 GO TO 2200
 1230 GO TO jgoto, (320, 420, 520, 620)
 
!     INTERNAL SUBROUTINE TO FIND THE DMAP STATEMENT NUMBER
!     FOR A DMAP STATEMENT WITH A GIVEN OCCURENCE FLAG AND
!     AN OFFSET FLAG
 
 1500 icount = 0
 DO  j = 1, iend, 2
   IF (idmap1 /= iopen(j) .OR. idmap2 /= iopen(j+1)) CYCLE
   icount = icount + 1
   IF (icount < ioccur) CYCLE
   inumbr = (j+1)/2 + ioffst
   IF (inumbr >= 1 .AND. inumbr <= nmdmap) GO TO 1700
   nogo = 1
   inumbr = 0
   WRITE (nout, 3500) icard
   GO TO 1700
 END DO
 nogo = 1
 inumbr = 0
 IF (icount > 0) WRITE (nout, 3600) icard
 IF (icount == 0) WRITE (nout, 3700) icard
 1700 GO TO igoto, (470, 570, 670, 690)
 
!     ERROR MESSAGES
 
 2000 CALL mesage (-2, altfil, isubr)
 2100 WRITE (nout, 4100) iend, n2dmap
 GO TO 2900
 2200 WRITE (nout, 4200) nwords, logic
 GO TO 2900
 2300 WRITE (nout, 4300) logic
 GO TO 2900
 2900 CALL mesage (-61, 0, 0)
 
 3000 RETURN
!***********************************************************************
 3300 FORMAT ('0*** USER FATAL MESSAGE, THE DATA ON THE ',  &
     'FOLLOWING ALTER CONTROL CARD IS NOT IN PROPER ', 'SEQUENCE OR ORDER --'//  &
     5X, 18A4)
 3500 FORMAT ('0*** USER FATAL MESSAGE, ILLEGAL OFFSET FLAG ',  &
     'SPECIFIED ON THE FOLLOWING ALTER CONTROL CARD --'// 5X, 18A4)
 3600 FORMAT ('0*** USER FATAL MESSAGE, ILLEGAL OCCURENCE FLAG ',  &
     'SPECIFIED ON THE FOLLOWING ALTER CONTROL CARD --'// 5X, 18A4)
 3700 FORMAT ('0*** USER FATAL MESSAGE, NON-EXISTENT NOMINAL ',  &
     'DMAP STATEMENT SPECIFIED ON THE FOLLOWING ', 'ALTER CONTROL CARD --'//  &
     5X, 18A4)
 4100 FORMAT ('0*** SYSTEM FATAL MESSAGE, ILLEGAL NUMBER OF ',  &
     'WORDS (', i5, ') ENCOUNTERED IN THE SECOND ',  &
     'FILE OF THE ALTER SCRATCH FILE.'/ '     EXPECTED NUMBER OF WORDS = ', i5)
 4200 FORMAT ('0*** SYSTEM FATAL MESSAGE, ILLEGAL NUMBER OF ',  &
     'WORDS (', i5, ') ENCOUNTERED WHILE READING ',  &
     'A RECORD IN THE FIRST FILE OF THE ALTER SCRATCH ', 'FILE.'/  &
     '     LOGIC ERROR NO. = ', i5)
 4300 FORMAT ('0*** SYSTEM FATAL MESSAGE, ILLEGAL CONTROL WORD ',  &
     'WHILE PROCESSING THE FOLLOWING ALTER CONTROL CARD'// 5X, 18A4)
!***********************************************************************
!                              NOTICE                             *
!                              ------                             *
!                                                                 *
!     THIS PROGRAM BELONGS TO RPK CORPORATION.  IT IS CONSIDERED  *
!  A TRADE SECRET AND IS NOT TO BE DIVULGED OR USED BY PARTIES    *
!  WHO HAVE NOT RECEIVED WRITTEN AUTHORIZATION FROM RPK.          *
!******************************************************************
END SUBROUTINE dmpalt
