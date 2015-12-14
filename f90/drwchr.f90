SUBROUTINE drwchr (x,y,xyd,chr,nn,opt)
     
!     (X,Y)  = STARTING OR ENDING POINT OF THE LINE TO BE TYPED (ALWAYS
!              LEFT-TO-RIGHT OR TOP-TO-BOTTOM)
!     XYD    = (+/-)1 IF X = STARTING OR ENDING POINT OF THE LINE
!            = (+/-)2 IF Y = STARTING OR ENDING POINT OF THE LINE
!     CHR    = CHARACTERS TO BE DRAWN
!     NN     = NUMBER OF CHARACTERS
!     OPT    = -1 TO INITIATE  THE TYPING MODE
!            = +1 TO TERMINATE THE TYPING MODE
!            =  0 TO TYPE A LINE
!     CSCALE = SCALE FOR CHARACTER SIZE (REAL)
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 INTEGER, INTENT(IN)                      :: xyd
 INTEGER, INTENT(IN)                      :: chr(1)
 INTEGER, INTENT(IN)                      :: nn
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: chrind,xychr,d
 REAL :: SAVE(2,2),xy(2,2),xyc(2,2),cscale
 COMMON /pltdat/ skplt(2),reg(2,2),xymax(2),edge(11),cscale, skpa(3),cntchr(2)
 COMMON /chrdrw/ lstind,chrind(60),xychr(2,1)
 DATA    lstchr/ 48 /
 
 IF (opt == 0) GO TO 100
 CALL line (0,0,0,0,0,opt)
 GO TO 200
 
 100 n = nn
 IF (n <= 0) n = 1
 d = MAX0(IABS(xyd),1)
 s = cntchr(d)
 IF (xyd == -1 .OR. xyd == 2) s = -s
 xyc(1,1) = 3.0*cscale
 xyc(2,1) = 3.0*cscale
 xy(1,1)  = x - xyc(1,1)
 xy(2,1)  = y - xyc(2,1)
 xy(1,2)  = xy(1,1)
 xy(2,2)  = xy(2,1)
 DO  i = 1,2
   SAVE(i,1)= reg(i,1)
   reg(i,1) = AMAX1(-edge(i),reg(i,1)-xyc(i,1))
   SAVE(i,2)= reg(i,2)
   reg(i,2) = AMIN1(xymax(i)+edge(i),reg(i,2)+xyc(i,1))
 END DO
 
!     TYPE THE LINE.
 
 DO  j = 1,n
   xy(d,2)  = xy(d,1) + s*FLOAT(j-1)
   
!     MAKE SURE EACH CHARACTER IS A VALID CHARACTER.
   
   i = j
   IF (xyd < 0) i = n - j + 1
   k = chr(i)
   IF (nn /= 0 .AND. k >= lstchr) CYCLE
   IF (k > lstind) CYCLE
   
!     DRAW THE CHARACTER.
   
   120 n1 = chrind(k)
   IF (n1 > 0) GO TO 121
   k  = -n1
   GO TO 120
   121 n2 = chrind(k+1)
   IF (n2 > 0) GO TO 122
   k  = k + 1
   GO TO 121
   
   122 n2 = n2 - 1
   DO  l = n1,n2
     DO  i = 1,2
       xyc(i,1) = xyc(i,2)
       xyc(i,2) = xy(i,2) + cscale*FLOAT(IABS(xychr(i,l)))
     END DO
     IF (l == n1 .OR. xychr(1,l) < 0 .OR. xychr(2,l) < 0) CYCLE
     CALL line (xyc(1,1),xyc(2,1),xyc(1,2),xyc(2,2),1,0)
   END DO
 END DO
 
 DO  i = 1,2
   reg(i,1) = SAVE(i,1)
   reg(i,2) = SAVE(i,2)
 END DO
 
 200 RETURN
END SUBROUTINE drwchr
