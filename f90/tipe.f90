SUBROUTINE tipe (x,y,xyd,chr,n,opt)
     
!     (X,Y) = STARTING OR ENDING POINT OF THE LINE TO BE TYPED (ALWAYS
!             LEFT-TO-RIGHT OR TOP-TO-BOTTOM.
!     XYD   = +/-1 IF X = STARTING OR ENDING POINT OF THE LINE.
!           = +/-2 IF Y = STARTING OR ENDING POINT OF THE LINE.
!     CHR   = CHARACTERS TO BE TYPED.
!     N     = NUMBER OF CHARACTERS.
!     OPT   = -1 TO INITIATE  THE TYPING MODE.
!           = +1 TO TERMINATE THE TYPING MODE.
!           =  0 TO TYPE A LINE.
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 INTEGER, INTENT(IN)                      :: xyd
 INTEGER, INTENT(IN)                      :: chr(1)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: ploter,CHAR,c(80),BLANK,lstchr, charx,d
 REAL :: xy(2,2)
 COMMON /pltdat/ model,ploter,skpplt(18),skpa(3),cntchr(2)
 COMMON /char94/ CHAR(60)
 DATA    BLANK , lstchr / 48,47 /
 
 IF (opt /= 0)  GO TO 150
 
!     OPT = 0.
 
 d = MAX0(IABS(xyd),1)
 s = cntchr(d)
 IF (xyd == -1 .OR. xyd == 2) s = -s
 xy(1,1) = x
 xy(2,1) = y
 xy(1,2) = xy(1,1)
 xy(2,2) = xy(2,1)
 
!     PRINT A MAXIMUM OF 80 CHARACTERS AT A TIME.
 
 DO  j = 1,n,80
   IF (xyd < 0) GO TO 105
   l1 = j
   l2 = l1 + 79
   IF (l2 > n) l2 = n
   GO TO 106
   105 l2 = n  - j + 1
   l1 = l2 - 79
   IF (l1 <= 0) l1 = 1
   
   106 nc = 0
   DO  l = l1,l2
     charx = chr(l)
     DO  i = 1,lstchr
       IF (charx == CHAR(i)) GO TO 111
     END DO
     i  = BLANK
     111 nc = nc + 1
     c(nc) = i
   END DO
   
!     TYPE THE -NC- CHARACTERS JUST PROCESSED.
   
   xy(d,2) = xy(d,1) + s*FLOAT(l1-1)
   CALL type10 (xy(1,2),xy(2,2),xyd,c,nc,0)
   CYCLE
 END DO
 GO TO 200
 
!     OPT = +/-1
 
 150 CALL type10 (0,0,0,0,0,opt)
 GO TO 200
 200 RETURN
END SUBROUTINE tipe
