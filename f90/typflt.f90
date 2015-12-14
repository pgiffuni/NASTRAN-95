SUBROUTINE typflt (x,y,xyd,v,field,opt)
     
 
!     (X,Y) = STARTING OR ENDING POINT OF THE NUMBER TO BE TYPED (ALWAYS
!             LEFT-TO-RIGHT OR TOP-TO-BOTTOM).
!     XYD   = +/-1 IF X = STARTING OR ENDING POINT OF THE NUMBER.
!           = +/-2 IF Y = STARTING OR ENDING POINT OF THE NUMBER.
!     V     = REAL NUMBER TO BE TYPED.
!     FIELD = FIELD WIDTH OF THE NUMBER (IF POSITIVE, THE NUMBER WILL BE
!             CENTERED AT (X,Y) - IF NEGATIVE, THE NUMBER WILL BE TYPED
!             STARTING OR ENDING AT (X,Y) - IF XYD = 1 OR 2, THE NUMBER
!             WILL BE TYPED IN THE X OR Y DIRECTION).
!     OPT   = -1 TO INITIATE  THE TYPING MODE.
!           = +1 TO TERMINATE THE TYPING MODE.
!           =  0 TO TYPE THE NUMBER.
 
 
 REAL, INTENT(IN)                         :: x
 REAL, INTENT(IN)                         :: y
 INTEGER, INTENT(IN)                      :: xyd
 REAL, INTENT(IN OUT)                     :: v
 INTEGER, INTENT(IN OUT)                  :: field
 INTEGER, INTENT(IN OUT)                  :: opt
 INTEGER :: ploter,dir,EXP,d(9),c(100),aster,decpnt,plus, minus,fw,expfld,tra
 DOUBLE PRECISION :: val,z
 COMMON /pltdat/  model,ploter,skpplt(18),skpa(3),cntx,cnty
 DATA    aster ,  decpnt,plus,minus / 41,44,39,40 /
 DATA    tenm2 ,  ten7,ten8  / 1.e-2, 1.e7, 1.e8  /
 
 IF (opt == 0) GO TO 20
 CALL tipe (0,0,0,0,0,opt)
 GO TO 200
 20 val = ABS(v)
 fw  = MIN0(25,IABS(field))
 IF (fw == 0) GO TO 200
 DO  i = 1,fw
   c(i) = 1
 END DO
 EXP = 0
 IF (v /= 0.) GO TO 30
 
!     INPUT VALUE = 0.
 
 fw   = MIN0(fw,2)
 nsig = 1
 c(2) = decpnt
 GO TO 150
 
 30 expfld = 0
 IF (v < 0.) GO TO 35
 
!     SINCE -V- IS POSITIVE, THE NUMBER WILL BE UNSIGNED. IF FIELD.GT.4,
!     THE NUMBER OF SIGNIFICANT DIGITS TYPED WILL BE AT LEAST -FIELD-4-.
!     IF FIELD.LE.4, -FIELD-1-.
 
 nsig = fw - 4
 IF (nsig > 0) THEN
   GO TO   100
 ELSE
   GO TO    40
 END IF
 
!     SINCE -V- IS NEGATIVE, THE NUMBER WILL BE SIGNED.  IF FIELD.GT.5,
!     THE NUMBER OF SIGNIFICANT DIGITS TYPED WILL BE AT LEAST -FIELD-5-.
!     IF FIELD.LE.5, -FIELD-2-.
 
 35 nsig = fw - 5
 IF (nsig > 0) THEN
   GO TO   100
 END IF
 
!     THE NUMBER WILL BE TYPED WITHOUT AN EXPONENT.
 
 40 nsig = nsig + 3
 expfld = 1
 
!     THE NUMBER MUST FIRST BE MULTIPLIED BY SOME POWER OF TEN (EXP)
!     SUCH THAT THE PRODUCT IS BETWEEN 10**7 AND 10**8 SO THAT IT
!     CAN BE EXPRESSED AS AN 8-SIGNIFICANT DIGIT INTEGER.
 
 100 z = 10.d0**IABS(EXP)
 IF (EXP < 0) a = val/z
 IF (EXP >= 0) a = val*z
 IF (a >= tenm2) GO TO 105
 
!     A .LT. 10**-2
 
 EXP = EXP + 10
 GO TO 100
 
 105 IF (a >= ten7 .AND. a < ten8) GO TO 115
 IF (a < ten7) GO TO 110
 
!     A .GE. 10**8
 
 EXP = EXP - 10
 GO TO 100
 
!     A .GE. 10**-2  AND  .LT. 10**7
 
 110 EXP = EXP + 1
 GO TO 100
 
!     A .GE. 10**7  AND  .LT. 10**8  (SEPARATE THE 8 SIGNIFICANT DIGITS)
 
 115 num = a
 EXP = -EXP + 7
 DO  i = 1,8
   j   = num/10**(8-i)
   d(i)= j + 1
   num = num - j*10**(8-i)
 END DO
 IF (expfld /= 0) GO TO 130
 IF (EXP >= -4 .AND. EXP <= nsig+2) GO TO 135
 
!     USE STANDARD FORMAT (-X.XXX-XX)
 
 nsig = MIN0(nsig,8)
 ASSIGN 120 TO tra
 GO TO 180
 120 n = 0
 IF (v > 0.) GO TO 121
 c(1) = minus
 n = 1
 121 c(n+1) = d(1)
 c(n+2) = decpnt
 n = n + 2
 IF (nsig == 1) GO TO 124
 DO  i = 2,nsig
   n = n + 1
   c(n) = d(i)
 END DO
 124 IF (EXP >= 0) c(n+1) = plus
 IF (EXP < 0) c(n+1) = minus
 n = n + 1
 num = IABS(EXP)
 DO  i = 1,2
   j = num/10**(2-i)
   n = fw - (2-i)
   c(n) = j + 1
   num = num - j*10**(2-i)
 END DO
 GO TO 150
 
!     STANDARD FORMAT CANNOT BE USED.
 
 130 IF (EXP < nsig .AND. EXP >= -nsig) GO TO 136
 DO  i = 1,fw
   c(i) = aster
 END DO
 GO TO 150
 
!     THE NUMBER CAN BE EXPRESSED WITHOUT AN EXPONENT.
 
 135 nsig = MIN0(8,nsig+3)
 136 ASSIGN 137 TO tra
 GO TO 180
 137 n = 1
 IF (v > 0.) GO TO 138
 c(1) = minus
 n = 2
 138 IF (EXP >= 0) GO TO 144
 
!     NEGATIVE EXPONENT
 
 j = nsig
 141 d(j+1) = d(j)
 j = j - 1
 IF (j /= 0) GO TO 141
 d(1) = 1
 ASSIGN 142 TO tra
 IF (nsig+n >= fw) GO TO 180
 nsig = nsig + 1
 142 c(n+0) = d(1)
 c(n+1) = decpnt
 n = n + 1 + IABS(EXP)
 DO  i = 2,nsig
   c(n) = d(i)
   n = n + 1
 END DO
 GO TO 150
 
!     POSITIVE EXPONENT.
 
 144 ASSIGN 145 TO tra
 IF (nsig+n >= fw) GO TO 180
 145 j = EXP + 1
 DO  i = 1,j
   c(n) = d(i)
   n = n + 1
 END DO
 c(n) = decpnt
 j = j + 1
 IF (j > nsig) GO TO 150
 DO  i = j,nsig
   n = n + 1
   c(n) = d(i)
 END DO
 
 150 xx = x
 yy = y
 IF (field > 0 .AND. nsig > 1) GO TO 155
 
!     THE TYPED NUMBER IS NOT TO BE CENTERED AT (X,Y).
 
 dir = xyd
 GO TO 160
 
!     THE TYPED NUMBER IS TO BE CENTERED AT (X,Y).
 
 155 xy = fw/2
 IF (fw/2 == (fw+1)/2) xy = xy - .5
 dir = MAX0(1,IABS(xyd))
 IF (dir == 1) xx = x - xy*cntx
 IF (dir == 2) yy = y - xy*cnty
 
!     TYPE THE NUMBER.
 
 160 CALL type10 (xx,yy,dir,c,fw,0)
 GO TO 200
 
!     ROUND THE NUMBER.
 
 180 IF (nsig == 8) GO TO 190
 IF (d(nsig+1) <= 5) GO TO 190
 j = nsig
 181 d(j) = d(j) + 1
 IF (d(j) <= 10) GO TO 190
 d(j) = 1
 j = j - 1
 IF (j /= 0) GO TO 181
 IF (d(1) /= 1) GO TO 190
 j = nsig - 1
 182 IF (j == 0) GO TO 183
 d(j+1) = d(j)
 j = j - 1
 GO TO 182
 183 d(1) = 2
 EXP = EXP + 1
 190 GO TO tra, (120,137,142,145)
 
 200 RETURN
END SUBROUTINE typflt
