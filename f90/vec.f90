SUBROUTINE vec
     
!     THE CALL TO THIS MODULE IS
!                   VEC USET  / V / C,N,X / C,N,X0 / C,N,X1 $
!          OR       VEC USETD / V / C,N,X / C,N,X0 / C,N,X1 $
 
!     ALTERNATE FORM OF THE CALL TO THIS MODULE IS
!                   VEC USET  / V / C,N,X / C,N,X0 / C,N,COMP $
!          OR       VEC USETD / V / C,N,X / C,N,X0 / C,N,COMP $
 
!     ALTERNATE FORM OF THE CALL TO THIS MODULE IS
!                   VEC USET  / V / C,N,X / C,N,COMP / C,N,X1 $
!          OR       VEC USETD / V / C,N,X / C,N,COMP / C,N,X1 $
 
!     ALTERNATE FORM OF THE CALL TO THIS MODULE IS
!                   VEC USET  / V / C,N,BITID / C,N,* / C,N,* / C,N,I $
!          OR       VEC USET  / V / C,N,BITID / C,N,X1 $
!          OR       VEC USETD / V / C,N,BITID / C,N,* / C,N,* / C,N,I $
!          OR       VEC USETD / V / C,N,BITID / C,N,X1 $
 
!     ALTERNATE FORM OF THE CALL TO THIS MODULE IS
!                   VEC USET  / V / C,N,COLUMNS / C,N,LEFT  / C,N,* /
!                                                             C,N,I $
!          OR       VEC USETD / V / C,N,COLUMNS / C,N,LEFT  / C,N,* /
!                                                             C,N,I $
!                   ( V WILL HAVE -I- COLUMNS GENERATED FROM BIT
!                     POSITIONS 1,2,3,...,I OF USET (OR USETD) WHERE
!                     THE 32 RIGHT-MOST BITS ARE CONSIDERED, COUNTING
!                     FROM LEFT TO RIGHT. )
 
!     ALTERNATE FORM OF THE CALL TO THIS MODULE IS
!                   VEC USET  / V / C,N,COLUMNS / C,N,RIGHT / C,N,* /
!                                                             C,N,I $
!          OR       VEC USETD / V / C,N,COLUMNS / C,N,RIGHT / C,N,* /
!                                                             C,N,I $
!                   ( V WILL HAVE -I- COLUMNS GENERATED FROM BIT
!                     POSITIONS 32,31,...,33-I OF USET (OR USETD) WHERE
!                     THE 32 RIGHT-MOST BITS ARE CONSIDERED, COUNTING
!                     FROM LEFT TO RIGHT. )
 
 
!     CORE REQUIREMENTS.. ONE BUFFER PLUS USET (OR USETD).
!     FOR COLUMNS OPTION, ONE GINO BUFFER PLUS 2*USET (OR USETD) REQD.
 
 
 EXTERNAL        andf
 LOGICAL :: lz,l0,l1,cols,flag1,flag2
 INTEGER :: andf,modnam(2),fi,fo,f,nam(2),t(7),two,  &
     p(2),p1,p2,p3,p4,bn,BLANK,tyin,tyou,b(2),c(2), offset,d(2),lr(2,2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / p1(2),p2(2),p3(2),p4 /zzzzzz/ x(1)  &
     /system/ lb,nout,nerr /bitpos/ bn(32,2)  &
     /packx / tyin,tyou,ii,nn,incr /two   / two(32)
 EQUIVALENCE     (nr,t(3))
 DATA    nermax, BLANK / 10,1H /
 DATA    b,c,d / 4HBITI,4HD    ,4HCOMP,4H     ,4HCOLU,4HMNS  /
 DATA    lr    / 4HRIGH,4HT    ,4HLEFT,4H     /
 DATA    modnam/ 4HVEC ,4H     /
 DATA    fi,fo , nbn   / 101,201, 32          /
 
 
 flag1  = .false.
 flag2  = .false.
 offset = 0
 nerr   = 0
 lz     = .false.
 l0     = .false.
 l1     = .false.
 lc     = korsz(x) - lb
 IF (lc <= 0) CALL mesage (-8,lc,modnam)
 ib = lc + 1
 
!     CHECK PARAMETER VALUES AND COMPUTE MASKS.
 
 IF (p1(1) /= d(1) .OR. p1(2) /= d(2)) GO TO 5
 cols = .true.
 DO  j = 1,2
   IF (p2(1) == lr(1,j) .AND. p2(2) == lr(2,j)) GO TO 4
 END DO
 j = 2
 4 j = 2*j - 3
 GO TO 13
 5 CONTINUE
 cols = .false.
 IF (p1(1) == b(1) .AND. p1(2) == b(2)) GO TO 13
 IF (p1(2) /= BLANK) GO TO 11
 DO  i = 1,nbn
   IF (p1(1) == bn(i,2)) GO TO 19
 END DO
 11 p(1) = p1(1)
 p(2) = p1(2)
 GO TO 9904
 13 lz   = .true.
 l0   = .true.
 IF (p4 < 0 .OR. p4 > 32) GO TO 9908
 IF (cols) GO TO 50
 IF (p4 > 0) GO TO 18
 IF (p2(2) /= BLANK) GO TO 21
 DO  i = 1,nbn
   IF (p2(1) == bn(i,2)) GO TO 35
 END DO
 GO TO 21
 18 maskx1 = two(p4)
 GO TO 50
 19 i = bn(i,1)
 maskx  = two(i)
 
 IF (p2(1) == c(1) .AND. p2(2) == c(2)) GO TO 23
 IF (p2(2) /= BLANK) GO TO 21
 DO  i = 1,nbn
   IF (p2(1) == bn(i,2)) GO TO 25
 END DO
 21 p(1) = p2(1)
 p(2) = p2(2)
 GO TO 9904
 23 l0   = .true.
 GO TO 26
 25 i    = bn(i,1)
 maskx0 = two(i)
 
 26 CONTINUE
 IF (p3(1) == c(1) .AND. p3(2) == c(2)) GO TO 33
 IF (p3(2) /= BLANK) GO TO 31
 DO  i = 1,nbn
   IF (p3(1) == bn(i,2)) GO TO 35
 END DO
 31 p(1) = p3(1)
 p(2) = p3(2)
 GO TO 9904
 33 l1   = .true.
 IF (l0) GO TO 9907
 GO TO 50
 35 i    = bn(i,1)
 maskx1 = two(i)
 
!     BLAST READ USET (OR USETD) INTO CORE.
 
 50 CONTINUE
 f = fi
 CALL fname (f,nam)
 CALL gopen (f,x(ib),0)
 CALL READ (*9902,*100,f,x,lc,0,nw)
 
!     INSUFFICIENT CORE - IF DESIRED, THIS ROUTINE CAN BE WRITTEN TO
!     RUN IN SMALLER CORE.
 
 lcex = 0
 70 CALL READ (*9902,*80,f,x,lc,0,nw)
 lcex = lcex + lc
 GO TO 70
 80 lcex = lcex + nw
 IF (cols) lcex = 2*lcex
 GO TO 9903
 100 CONTINUE
 CALL CLOSE (f,1)
 IF (.NOT.cols) GO TO 150
 IF (p4 <= 0) GO TO 9908
 offset = nw
 k = 1
 l = 1
 IF (j < 0) k = 32
 maskx1 = two(k)
 IF (2*nw <= lc) GO TO 150
 lcex = 2*nw - lc
 GO TO 9903
 150 CONTINUE
 
!     PREPARE OUTPUT FILE.
 
 f = fo
 CALL gopen  (f,x(ib),1)
 CALL makmcb (t,f,0,2,1)
 tyin = 1
 tyou = 1
 ii   = 1
 incr = 1
 
!     CREATE VECTOR IN CORE OCCUPIED BY USET (OR USETD).
 
 170 nr = 0
 nz = 0
 
 DO  i = 1,nw
   IF (lz) GO TO 220
   IF (andf(x(i),maskx) == 0) GO TO 400
   220 CONTINUE
   IF (.NOT.l0) GO TO 230
   IF (andf(x(i),maskx1) == 0) GO TO 370
   GO TO 300
   230 IF (.NOT.l1) GO TO 240
   IF (andf(x(i),maskx0) == 0) GO TO 300
   GO TO 370
   240 CONTINUE
   IF (andf(x(i),maskx1) == 0) GO TO 350
   IF (andf(x(i),maskx0) == 0) GO TO 300
   nerr = nerr + 1
   IF (nerr > nermax) CYCLE
   WRITE  (nout,250) ufm,i
   250 FORMAT (a23,' 2120, MODULE VEC - BOTH SUBSET BITS ARE NON-ZERO.',  &
       3X,'I =',i10)
   CYCLE
   300 nr = nr + 1
   nz = nz + 1
   x(nr+offset) = 1.0
   CYCLE
   350 CONTINUE
   IF (andf(x(i),maskx0) /= 0) GO TO 370
   nerr = nerr + 1
   IF (nerr > nermax) CYCLE
   WRITE  (nout,360) ufm,i
   360 FORMAT (a23,' 2121, MODULE VEC - BOTH SUBSET BITS ARE ZERO.',3X,  &
       'I =',i10)
   CYCLE
   370 nr = nr + 1
   x(nr+offset) = 0.0
   CYCLE
   400 IF (l0) GO TO 450
   IF (andf(x(i),maskx0) == 0) GO TO 450
   nerr = nerr + 1
   IF (nerr > nermax) GO TO 450
   WRITE  (nout,410) ufm,i
   410 FORMAT (a23,' 2122, MODULE VEC - SET X BIT IS ZERO BUT SUBSET X0',  &
       ' BIT IS NOT.  I =',i10)
   450 IF (l1) CYCLE
   IF (andf(x(i),maskx1) == 0) CYCLE
   nerr = nerr + 1
   IF (nerr > nermax) CYCLE
   WRITE  (nout,460) ufm,i
   460 FORMAT (a23,' 2123, MODULE VEC - SET X BIT IS ZERO BUT SUBSET X1',  &
       ' BIT IS NOT.  I =',i10)
 END DO
 
 IF (nerr <= 0) GO TO 540
 IF (nerr-nermax > 0) THEN
   GO TO  9906
 ELSE
   GO TO  9995
 END IF
 540 CONTINUE
 
 IF (flag1) GO TO 600
 flag1 = .true.
 IF (nr > 0) GO TO 600
 WRITE  (nout,550) uwm
 550 FORMAT (a25,' 2124, MODULE VEC - NR=0, OUTPUT WILL BE PURGED.')
 GO TO 900
 600 IF (nz > 0) GO TO 700
 IF (flag2) GO TO 700
 flag2 = .true.
 WRITE  (nout,650) uwm
 650 FORMAT (a25,' 2125, MODULE VEC - NZ=0, ONE OR MORE COLUMNS OF ',  &
     'OUTPUT MATRIX WILL BE NULL.')
 GO TO 750
 700 CONTINUE
 
!     PACK OUT COLUMN OF OUTPUT VECTOR.
 
 750 nn = nr
 CALL pack (x(offset+1),f,t)
 IF (.NOT.cols .OR. l >= p4) GO TO 800
 l = l + 1
 k = k + j
 maskx1 = two(k)
 GO TO 170
 800 CALL wrttrl (t)
 900 CALL CLOSE  (f,1)
 
 RETURN
 
!     ERROR PROCESSING.
 
 9902 WRITE  (nout,9952) ufm,f,nam
 9952 FORMAT (a23,' 2141, MODULE VEC - EOF ENCOUNTERED WHILE READING ',  &
     'GINO FILE ',i3,', DATA BLOCK ',2A4)
 GO TO 9995
 9903 WRITE  (nout,9953) ufm,lc,lcex
 9953 FORMAT (a23,' 2142, INSUFFICIENT CORE FOR MODULE VEC.  AVAILABLE',  &
     ' CORE =',i11,' WORDS.', /5X, 'ADDITIONAL CORE NEEDED =',i11,' WORDS.')
 GO TO 9995
 9904 WRITE  (nout,9954) ufm,p
 9954 FORMAT (a23,' 2143, MODULE VEC UNABLE TO IDENTIFY SET OR SUBSET ',  &
     'DESCRIPTOR ',2A4)
 GO TO 9995
 9906 WRITE  (nout,9956) ufm,nerr,nermax
 9956 FORMAT (a23,' 2145,',i8,' FATAL MESSAGES HAVE BEEN GENERATED IN',  &
     ' SUBROUTINE VEC.', /5X, 'ONLY THE FIRST',i4,' HAVE BEEN PRINTED.')
 GO TO 9995
 9907 WRITE  (nout,9957) ufm
 9957 FORMAT (a23,' 2146, BOTH OF THE SECOND AND THIRD VEC PARAMETERS ',  &
     'REQUEST COMPLEMENT.')
 GO TO 9995
 9908 WRITE  (nout,9958) ufm,p4
 9958 FORMAT (a23,' 2150, ILLEGAL VALUE FOR FOURTH PARAMETER =',i11)
 GO TO 9995
 9995 CALL mesage (-61,0,0)
 RETURN
 
END SUBROUTINE vec
