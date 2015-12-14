SUBROUTINE typint (x,y,xyd,num,field,opt)
     
!     (X,Y) = STARTING OR ENDING POINT OF THE NUMBER TO BE TYPED (ALWAYS
!             LEFT-TO-RIGHT OR TOP-TO-BOTTOM).
!     XYD   = NO ACTION IF OPT IS ZERO
!           = (+/-)1 IF X = STARTING OR ENDING POINT OF THE NUMBER.
!           = (+/-)2 IF Y = STARTING OR ENDING POINT OF THE NUMBER.
!     NUM   = INTEGER NUMBER TO BE TYPED (AT MOST 10 DIGITS).
!     FIELD = NO ACTION IF OPT IS ZERO
!           = 1 IF THE NUMBER IS TO BE CENTERED AT (X,Y). IF XYD=1 OR 2,
!             THE NUMBER WILL BE TYPED IN THE X OR Y DIRECTION.
!           = 0 OR -1 IF THE NUMBER IS TO BE TYPED STARTING OR ENDING AT
!             (X,Y). IF FIELD = -1, FIELD WILL BE SET TO THE NUMBER OF
!             DIGITS PRINTED.
!     OPT   =-1 TO INITIATE  THE TYPING MODE.
!           =+1 TO TERMINATE THE TYPING MODE.
!           = 0 TO TYPE A LINE.
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 INTEGER, INTENT(IN)                      :: xyd
 INTEGER, INTENT(IN OUT)                  :: num
 INTEGER, INTENT(OUT)                     :: field
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: ploter,aster,dir,d(11)
 COMMON /pltdat/ model,ploter,skpplt(18),skpa(3),cntx,cnty
 DATA    aster , minus / 41,40 /
 
 IF (opt == 0) GO TO 100
 CALL tipe (0,0,0,0,0,opt)
 GO TO 200
 
!     SEPARATE THE DIGITS OF THE NUMBER (MAXIMUM OF 10).
 
 100 nd = -1
 IF (num >= 0) GO TO 110
 nd = 0
 d(1) = minus
 110 n = IABS(num)
 DO  i = 1,10
   j = n/10**(10-i)
   IF (j == 0 .AND. nd <= 0) CYCLE
   IF (j  > 9) j  = aster - 1
   IF (nd <= 0) nd = nd + 1
   nd = nd + 1
   d(nd) = j + 1
   n = n - j*10**(10-i)
 END DO
 IF (nd > 0) GO TO 112
 nd   = 1
 d(1) = 1
 
 112 xx = x
 yy = y
 IF (field > 0 .AND. nd > 1) GO TO 120
 
!     THE TYPED NUMBER IS NOT TO BE CENTERED AT (X,Y).
 
 dir = xyd
 IF (field < 0) field = nd
 GO TO 150
 
!     THE TYPED NUMBER MUST BE CENTERED AT (X,Y).
 
 120 xy = nd/2
 IF (nd/2 == (nd+1)/2) xy = xy - .5
 dir = MAX0(IABS(xyd),1)
 IF (dir == 1) xx = x - xy*cntx
 IF (dir == 2) yy = y - xy*cnty
 
!     TYPE THE NUMBER.
 
 150 CALL type10 (xx,yy,dir,d,nd,0)
 GO TO 200
 
 200 RETURN
END SUBROUTINE typint
