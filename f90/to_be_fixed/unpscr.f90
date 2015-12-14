SUBROUTINE unpscr (in,out,z,buf1,buf2,maxz,tysign,flag)
     
!     THIS ROUTINE UNPACKS A MATRIX (IN), AND TRANSFER THE DATA FROM
!     FIRST TO LAST NON-ZERO TERMS TO A SCRATCH FILE (OUT) IN VERY LARGE
!     RECORD(S), PRECEEDED BY THE FIRST AND LAST NON-ZERO TERM POINTERS.
 
!     INPPUT  - IN, + 7 TRAILER WORDS (WORDS 4,5,6, AND 7 WILL BE
!               OVERWRITTEN)
!               Z, BUF1, BUF2, MAXZ, TYSIGN, AND FLAG
!     OUTPUT  - OUT, NO TRAILER WORD WRITTEN
!               IN(4) = 10*(NO. OF RECONDS WRITTEN, HEADER RECORD
!                       EXCLUDED) + FLAG
!               IN(5) = DATA WORD TYPE UNPACKED (= 1,2,OR 4)
!               IN(6) = TOTAL NO. OF S.P. WORDS USED FOR INPUT MATRIX
!                       IN FORWARD UNPACK PASS
!               IN(7) = OUTPUT GINO NUMBER
 
!     FLAG = 1, THE MATRIX IS UNPACKED ONCE, IN FORWARD DIRECTION, THIS
!               MATRIX CAN BE IN GENERAL FORM; NEEDS NOT BE TRIANGULAR.
!     FLAG = 2, THE MATRIX IS UNPACKED FORWARD AND BACKWARD
!     FLAG = 3, THE MATRIX IS ADVANCED TO THE END AND UNPACKED BACKWARD
!               ONCE AND THEN FORWARD
!     MAXZ = n, WHERE n IS THE UPER LIMIT OF THE RECORD SIZE TO BE
!               WRITTEN (5000 MINIMUM).
!          = 0  OR LESS, OUTPUT WILL BE WRITTEN OUT IN EITHER ONE OR TWO
!               LONG RECORDS (ONE EACH FOR FORWARD AND BACKWARD UNPACK)
!     Z    =    WORKING SPACE, MINIMUM SIZE = ROW + 2 WORDS
!     TYSIGN =  (-4,-3,...,+4), IS TYPE AND SIGN FOR INPUT MATRIX UNPACK
!               NO TYPE AND SIGN CHANGE IF TYSIGN = 0.
!     BUF1, BUF2 = TWO GINO BUFFERS
!     SUBROUTINE DEBUG CAN BE ACTIVATED BY DIAG 11 OR 16
 
!     ASSUME MATRIX IN(5x5) =  a  0  0  0  0
!                              b  e  0  0  0
!                              c  f  g  0  0
!                              d  0  h  j  0
!                              0  0  i  k  l
 
!     OUTPUT FILE OUT WILL HAVE THE FOLLOWING DATA (PRECEEDED BY HEADER
!     RECORD)
 
!     FLAG 1 -  1 4 a b c d 2 3 e f 3 5 g h i 4 5 j k 5 5 l <EOF>
!     FLAG 2 -  1 4 a b c d 2 3 e f 3 5 g h i 4 5 j k 5 5 l <EOR>
!               5 5 l 4 5 j k 3 5 g h i 2 4 e f 1 4 a b d c <EOF>
!     FLAG 3 -  5 5 l 4 5 j k 3 5 g h i 2 3 e f 1 4 a b c d <EOR>
!               1 4 a b c d 2 3 e f 3 5 g h i 4 5 j k 5 5 l <EOF>
 
!     WHERE a thru l MAY BE SP, DP, CSP, OR CDP DATA
 
!     IF INPUT MATRIX IS VERY LARGE, THERE WILL BE SEVERAL LONG RECORDS
!     FOR EACH UNPACK PASS, AND EACH RECORD WILL NOT EXCEED MAXZ IN
!     LENGTH. MINIMUM OF MAXZ IS 5000. IF MAXZ IS NOT GIVEN, EACH UNPACK
!     PASS WILL GO TO ONE VERY VERY LONG RECORD. IN THIS CASE, MAXZ IS
!     SET TO 2**31
 
!     THE PURPOSE OF THIS ROUTINE IS TO AVOID UNPACKING A MATRIX TOO
!     MANY TIMES, WHILE THE MATRIX IS BEING USED REPEATEDLY.
!     SEE FBSII (REPEATEDLY CALLED BY FBS), FRBK2 (REPEATEDLY CALLED
!     BY FNXTVC), AND FRMLTD (REPEATED CALLED BY FRBK2 AND FNXTVC) IN
!     USING THIS NEW DATA FORMAT.
 
!     WRITTEN BY G.CHAN/UNISYS   11/1991
 
!     COMMENTS FROM G.C.  3/93
!     THE PRESENT UNPSCR ASSUMES THE MATRIX IS QUIT DENSE, SUCH AS THE
!     LOWER OR UPPER TRIANGULAR FACTORS. IF MATRIX IS SPARSE, SAY 33
!     PERCENT OF LESS, WE COULD WRITE THE MATRIX OUT ANOTHER WAY AND
!     SAVE LOTS OF DISC SPACE. WE COULD WRITE THE FIRST TO LAST NON-ZERO
!     TERMS IN STRING FORMS SIMILAR TO OUTPUT4 MODULE. THIS IMPROVEMENT
!     WILL BE LEFT FOR NEXT PROJECT.
 
 
 INTEGER, INTENT(IN OUT)                  :: in(7)
 INTEGER, INTENT(IN)                      :: out
 INTEGER, INTENT(OUT)                     :: z(3)
 INTEGER, INTENT(IN OUT)                  :: buf1
 INTEGER, INTENT(IN OUT)                  :: buf2
 INTEGER, INTENT(IN)                      :: maxz
 INTEGER, INTENT(IN)                      :: tysign
 INTEGER, INTENT(IN)                      :: flag
 IMPLICIT INTEGER (a-z)
 LOGICAL :: flag23,debug
 INTEGER :: nam(2),tyiijj(4),SAVE(4)
 CHARACTER (LEN=8) :: fbwd,forwd,backwd
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ sysbuf,nout
 COMMON /unpakx/ TYPE,ii,jj,incr
 COMMON /names / rd,rdrew,wrt,wrtrew,rew
 COMMON /TYPE  / rc(2),words(4)
 EQUIVALENCE     (TYPE,tyiijj(1))
 DATA    forwd , backwd / 'FORWARD','BACKWARD'/
 DATA    nam   / 4HUNPS , 2HCR  /
 
 IF (flag < 1 .OR. flag > 3 .OR. in(1) == out) GO TO 300
 CALL sswtch (11,i)
 CALL sswtch (16,j)
 debug = .false.
 IF (i+j >= 1) debug = .true.
 MAX = maxz
 IF (MAX <= 0) MAX = 1073741824
 IF (debug) WRITE (nout,5) uim
 5 FORMAT (a29,', UNPSCR DEBUG, ACTIVATED BY DIAG 11 AND/OR 16')
 IF (MAX < 5000) GO TO 280
 flag23 = flag == 2 .OR. flag == 3
 DO  i = 1,4
   SAVE(i) = tyiijj(i)
 END DO
 TYPE = in(5)
 nl   = in(2)
 IF (tysign /= 0 .AND. IABS(tysign) <= 4) TYPE = tysign
 nwds = words(IABS(TYPE))
 IF (debug) WRITE (nout,15) in(1),out,maxz,MAX,flag,nl,TYPE,nwds
 15 FORMAT (5X,'UNPSCR/@15  IN,OUT,MAXZ,MAX,FLAG,NL,TYPE,NWDS = ',  &
     2I5,2I12,i4,i7,2I4)
 incr = 1
 FORM = in(4)
 IF (flag23 .AND. FORM /= 4 .AND. FORM /= 5) GO TO 260
!                          LOWER  AND      UPPER  TRIANGULAR FACTORS
 
 FILE = out
 CALL gopen (out,z(buf2),wrtrew)
 FILE = in(1)
 CALL OPEN (*200,in,z(buf1),rdrew)
 nrec = 0
 IF (flag == 3) GO TO 90
 20 CALL fwdrec (*210,in)
 
!     UNPACK FORWARD
 
 fbwd = forwd
 tot  = 0
 sum  = 0
 DO  i = 1,nl
   ii   = 0
   CALL unpack (*60,in,z(3))
   IF (flag23 .AND. ii /= i) GO TO 220
   30 z(1) = ii
   z(2) = jj
   ll   = (jj-ii+1)*nwds + 2
   tot  = tot + ll
   sum  = sum + ll
   IF (sum <= MAX) GO TO 50
   nrec = nrec + 1
   CALL WRITE (out,0,0,1)
   sum  = sum - ll
   IF (debug) WRITE (nout,40) nrec,sum,fbwd
   40 FORMAT (5X,'UNPSCR WROTE RECORD',i5,',  NO. OF WORDS =',i9,2X,a8)
   sum  = ll
   50 CALL WRITE (out,z(1),ll,0)
   CYCLE
   60 IF (flag23) GO TO 240
   ii   = i
   jj   = i
   DO  k = 3,6
     z(k) = 0
   END DO
   GO TO 30
 END DO
 nrec = nrec + 1
 CALL WRITE (out,0,0,1)
 IF (debug) WRITE (nout,40) nrec,sum,fbwd
 IF (flag /= 2) GO TO 150
 CALL bckrec (in)
 GO TO 100
 
 90 CALL skprec (in,nl)
 
!     UNPACK BACKWARD
 
 100 fbwd = backwd
 sum  = 0
 i    = nl
 DO  j = 1,nl
   ii   = 0
   CALL unpack (*240,in,z(3))
   IF (ii /= i) GO TO 220
   z(1) = ii
   z(2) = jj
   ll   = (jj-ii+1)*nwds + 2
   sum  = sum + ll
   IF (sum <= MAX) GO TO 110
   nrec = nrec + 1
   CALL WRITE (out,0,0,1)
   sum  = sum - ll
   IF (debug) WRITE (nout,40) nrec,sum,fbwd
   sum  = ll
   110 CALL WRITE (out,z(1),ll,0)
   CALL bckrec (in)
   CALL bckrec (in)
   i    = i - 1
 END DO
 nrec = nrec + 1
 CALL WRITE (out,0,0,1)
 IF (debug) WRITE (nout,40) nrec,sum,fbwd
 IF (flag == 3) GO TO 20
 
!     END OF UNPACKING
 
!     CHANGE LAST 4 WORDS OF THE INPUT MATRIX TRAILER. PARTICULARY, SET
!     THE 7TH WORD TO NEGATIVE. NOTE, IF FLAG IS 2 OR 3, IN(4) AND IN(6)
!     TRAILER WORDS HOLD HALF OF THE ACTUAL VALUES.
!     NOTE - SINCE WRTTRL IS NOT CALLED TO REGISTER THESE TRAILER WORD
!     CHANGES, THE TRAILER WORDS ARE INTENDED FOR THE ROUTINE TO BE
!     EXECUTE NEXT.  ALSO NOTE THAT OUTPUT FILE HAS NO TRAILER.
!     LASTLY, WE NEED TO RESTORE ORIGINAL WORDS IN /UNPAKX/ PREVIOUSLY
!     SAVED.
 
 150 CALL CLOSE (in, rew)
 CALL CLOSE (out,rew)
 in(7) =-out
 in(6) = tot
 in(5) = nwds
 i     = nrec
 IF (.NOT.flag23) GO TO 160
 i     = nrec/2
 tot   = tot*2
 160 in(4) = 10*i + flag
 DO  i = 1,4
   tyiijj(i) = SAVE(i)
 END DO
 IF (.NOT.debug) GO TO 350
 WRITE  (nout,180) uim,tot,nrec,nl,in(3)
 180 FORMAT (a29,1H,,i10,' S.P. WORDS MOVED TO SCRATCH FILE BY UNPSCR',  &
     /5X,'IN',i5,' RECORDS.', 5X,'INPUT MATRIX =',i8,3H by,i7)
 GO TO 350
 
 200 j = -1
 GO TO 330
 210 j = -2
 GO TO 330
 220 WRITE  (nout,230) sfm,i,ii,jj,fbwd,flag
 230 FORMAT (a25,',  I & II MISMATCH ',3I6,3H  /,a8,i9)
 GO TO 320
 240 WRITE  (nout,250) i,fbwd,flag
 250 FORMAT ('0*** NULL COLUMN ENCOUNTERED IN TRIANGULAR FACTOR.  ',  &
     'COLUMN',i7,3X,a8,i9)
 GO TO 320
 260 CALL fname (in(1),in(2))
 WRITE  (nout,270) in(2),in(3),FORM,flag
 270 FORMAT ('0*** INPUT MATRTIX ',2A4,' IS NOT A TRIANGULAR FACTOR.',  &
     '   FORM,FLAG =',2I4)
 CALL errtrc ('UNPSCR  ',270)
 280 WRITE  (nout,290) maxz
 290 FORMAT ('0*** MAXZ ERROR ',i9,'  (TOO SMALL)')
 CALL errtrc ('UNPSCR  ',290)
 GO TO 320
 300 WRITE  (nout,310) sfm,flag,in(1),out
 310 FORMAT (a25,',  FLAG,IN(1),OUT =',3I5)
 320 j = -37
 330 CALL mesage (j,FILE,nam)
 
 350 IF (debug) WRITE (nout,360)
 360 FORMAT (' ... UNPSCR DEBUG ENDS',/)
 RETURN
END SUBROUTINE unpscr
