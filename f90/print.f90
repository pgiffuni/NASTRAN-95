SUBROUTINE PRINT (x,y,xyd,chr,n,opt)
     
!     (X,Y) = STARTING OR ENDING POINT OF THE LINE TO BE PRINTED (ALWAYS
!             LEFT-TO-RIGHT OR TOP-TO-BOTTOM).
!     CHR   = CHARACTERS TO BE PRINTED (4 PER WORD).
!     N     = NUMBER OF 4 CHARACTER WORDS.
!     XYD   = +/-1 IF X = STARTING OR ENDING POINT OF THE LINE.
!     ...   = +/-2 .. Y = ........ .. ...... ..... .. ... .....
!     OPT   = -1 TO INITIATE  THE TYPING MODE.
!     ...   = +1 .. TERMINATE ... ...... .....
!     ...   =  0 .. PRINT A LINE.
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 INTEGER, INTENT(IN)                      :: xyd
 INTEGER, INTENT(IN OUT)                  :: chr(1)
 INTEGER, INTENT(IN)                      :: n
 INTEGER, INTENT(IN OUT)                  :: opt
 EXTERNAL        orf,krshft,klshft
 INTEGER :: orf,c(80),BLANK,blnk,charx,d
 REAL :: xy(2,2)
 COMMON /pltdat/ skpplt(20),skpa(3),cntchr(2)
 COMMON /system/ skpsys(40),ncpw
 DATA    BLANK / 1H /
 
 IF (opt /= 0) GO TO 150
 blnk = krshft(klshft(BLANK,1),1)
 d = MAX0(IABS(xyd),1)
 s = cntchr(d)
 IF (xyd == -1 .OR. xyd == 2) s = -s
 xy(1,1) = x
 xy(2,1) = y
 xy(1,2) = xy(1,1)
 xy(2,2) = xy(2,1)
 
!     SEPARATE 80 CHARACTERS AT A TIME.
 
 DO  j = 1,n,20
   IF (xyd < 0) GO TO 105
   l1 = j
   l2 = l1 + 19
   IF (l2 > n) l2 = n
   GO TO 106
   105 l2 = n - j + 1
   l1 = l2 - 19
   IF (l1 <= 0) l1 = 1
   
   106 nc = 0
   DO  l = l1,l2
     DO  i = 1,4
       charx = krshft(chr(l),ncpw-i)
       nc = nc + 1
       c(nc) = orf(klshft(charx,ncpw-1),blnk)
     END DO
   END DO
   
!     TYPE THE -NC- CHARACTERS JUST SEPARATED.
   
   xy(d,2) = xy(d,1) + s*FLOAT(l1-1)
   CALL tipe (xy(1,2),xy(2,2),xyd,c,nc,0)
 END DO
 GO TO 200
 
!     OPT = +/-1
 
 150 CALL tipe (0,0,0,0,0,opt)
 200 RETURN
END SUBROUTINE PRINT
