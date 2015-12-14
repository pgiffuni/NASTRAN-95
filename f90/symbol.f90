SUBROUTINE symbol (x,y,symx,opt)
     
!     (X,Y) = POINT AT WHICH THE SYMBOLS ARE TO BE TYPED.
!     SYMX  = SYMBOLS TO BE TYPED.
!     OPT   = -1 TO INITIATE  THE TYPING MODE.
!           = +1 TO TERMINATE THE TYPING MODE.
!           =  0 TO TYPE THE SYMBOL.
 
 
 REAL, INTENT(IN OUT)                     :: x
 REAL, INTENT(IN OUT)                     :: y
 INTEGER, INTENT(IN)                      :: symx(2)
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: sym, ploter,symbl
 COMMON /pltdat/ model,ploter
 COMMON /symbls/ nsym, symbl(20,2)
 
 IF (opt == 0) GO TO 110
 CALL tipe (0,0,0,0,0,opt)
 GO TO 200
 
 110 DO  i = 1,2
   IF (symx(i) <= 0) CYCLE
   sym = symx(i) - nsym*((symx(i)-1)/nsym)
   sym = symbl(sym,ploter)
   CALL type10 (x,y,0,sym,0,0)
   CYCLE
 END DO
 
 200 RETURN
END SUBROUTINE symbol
