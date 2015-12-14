SUBROUTINE type10 (x,y,xyd,chr,nn,opt)
!
!     (X,Y) = STARTING OR ENDING POINT OF THE LINE TO BE TYPED (ALWAYS
!             LEFT-TO-RIGHT OR TOP-TO-BOTTOM)
!     XYD   = (+/-)1 IF X = STARTING OR ENDING POINT OF THE LINE
!           = (+/-)2 IF Y = STARTING OR ENDING POINT OF THE LINE
!     CHR   = CHARACTERS TO BE TYPED
!     NN    = NUMBER OF CHARACTERS
!     OPT   = -1 TO INITIATE  THE TYPING MODE
!           = +1 TO TERMINATE THE TYPING MODE
!           =  0 TO TYPE A LINE
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 INTEGER, INTENT(IN)                      :: xyd
 INTEGER, INTENT(IN OUT)                  :: chr(1)
 INTEGER, INTENT(IN)                      :: nn
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: optx,a(6),TYPE,d,pltype
 REAL :: xy(2,2),cscale
 COMMON /pltdat/ skpplt(2),xymin(2),xymax(15),cscale,skpa(3), cntchr(6),pltype
 DATA    a(6)  , TYPE,lstchr / 0, 4, 48 /
 
 IF (pltype < 0) GO TO 175
 optx = -1
 IF (opt < 0.0) THEN
   GO TO   200
 ELSE IF (opt == 0.0) THEN
   GO TO   100
 ELSE
   GO TO   150
 END IF
 100 a(5) = IFIX(cscale+.44)
 xy(1,1) = x
 xy(2,1) = y
 xy(1,2) = x
 xy(2,2) = y
 n = 1
 IF (n <= 0) n = 1
 
!     SCREEN OUT TRAILING BLANKS
 
 DO  j = 1,nn
   IF (IABS(chr(j)) /= 48) n = j
 END DO
 IF (n == 1 .AND. IABS(chr(1)) == 48) RETURN
 d = MAX0(IABS(xyd),1)
 s = cntchr(d)
 IF (xyd == -1 .OR. xyd == 2) s = -s
 
!     TYPE THE LINE
 
 loop125:  DO  j = 1,n
   xy(d,2) = xy(d,1) + s*FLOAT(j-1)
   DO  i = 1,2
     IF (xy(i,2)+.1 < xymin(i) .OR. xy(i,2)-.1 > xymax(i)) CYCLE loop125
     a(i+2) = xy(i,2) + .1
   END DO
   
!     MAKE SURE EACH CHARACTER IS A VALID CHARACTER (UNLESS NN.LE.0)
   
   k = j
   IF (xyd < 0) k = n - j + 1
   a(2) = IABS(chr(k))
   IF (nn  <= 0) GO TO 120
   IF (a(2) == 0 .OR. a(2) > lstchr) CYCLE loop125
   IF (a(2) == 0) CYCLE loop125
   
!     TYPE THE CHARACTER
   
   120 a(1) = TYPE
   IF (optx == 0) GO TO 121
   a(1) = TYPE + 10
   optx = 0
   121 CALL wplt10 (a,0)
 END DO loop125
 GO TO 200
 
!     TERMINATE THE TYPING MODE
 
 150 CALL wplt10 (a,1)
 optx = -1
 GO TO 200
 
!     DRAW THE LINE OF CHARACTERS
 
 175 CALL drwchr (x,y,xyd,chr,nn,opt)
 
 200 RETURN
END SUBROUTINE type10
