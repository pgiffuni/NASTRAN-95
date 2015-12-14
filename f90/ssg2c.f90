SUBROUTINE ssg2c (a,b,c,op,BLOCK)
     
 
 INTEGER, INTENT(IN)                      :: a
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN)                      :: c
 INTEGER, INTENT(IN OUT)                  :: op
 INTEGER, INTENT(IN)                      :: BLOCK(11)
 INTEGER :: na(2)    ,nb(2)
 DIMENSION        ia(5)    ,ib(5)    ,ic(5)    ,it(1)    ,it1(1)
 DOUBLE PRECISION :: dit1
 INTEGER :: dt1(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm      ,uwm
 COMMON /zzzzzz/  core(1)
 COMMON /saddx /  nomat    ,lcore    ,mcbs(67)
 COMMON /system/  ksystm(65)
 EQUIVALENCE      (dt1,dit1)
 EQUIVALENCE      (ksystm(55),ipr1)  ,(mcbs(1),ia(1))  ,  &
     (mcbs(8),it(1),dit),(mcbs(13),ib(1)) ,  &
     (mcbs(20),it1(1)),(mcbs(61),ic(1)) ,  &
     (ksystm(2),nout)   ,(ia5,ia(5)) ,(ib5,ib(5))
 
!     BLOCK(6) WAS NOT USED IN ORIGINAL NASTRAN. IT IS NOW USED TO FLAG
!     THE CHECKING OF THE INPUT MATRICES COMPATABILITY IF THE CALLER
!     PRESETS BLOCK(6) TO -1
 
 ia(1) = a
 CALL rdtrl (ia)
 IF (ia(1) < 0) ia(1) = 0
 ib(1) = b
 CALL rdtrl (ib)
 IF (ib(1) > 0) GO TO 10
 ib(1) = 0
 IF (ia(1) > 0) THEN
   GO TO    30
 ELSE
   GO TO   150
 END IF
 10 DO  i = 2,4
   ic(i) = ib(i)
 END DO
 GO TO 50
 30 DO  i = 2,4
   ic(i) = ia(i)
 END DO
 
 50 nomix = 0
 IF (BLOCK(6) /= -1) GO TO 70
 IF (ia5 == 0 .OR. ib5 == 0) GO TO 70
 IF ((ia5 <= 2 .AND. ib5 <= 2) .OR. (ia5 >= 3 .AND. ib5 >= 3)) GO TO 70
 IF (MAX0(ia5,BLOCK(1)) == MAX0(ib5,BLOCK(7))) GO TO 70
 nomix = 1
 CALL fname (a,na)
 CALL fname (b,nb)
 WRITE  (nout,60) uwm,na,ia(2),ia(3),ia5,ia(4),nb,ib(2),ib(3), ib5,ib(4)
 60 FORMAT (a25,', SSG2C RECEIVES TWO MIXED FILE TYPES FOR ADDING.',  &
     /,2(5X,'FILE ',2A4,'(',i6,' X',i6,') TYPE =',i3, ', FORM =',i3))
 
!     UNSY + SYM = UNSY
 
 70 IF (ic(4) /= 6) GO TO 80
 IF (ia(1) /= 0 .AND. ia(4) /= 6) ic(4) = 1
 IF (ib(1) /= 0 .AND. ib(4) /= 6) ic(4) = 1
 80 IF (op < 0) ia(2) = -ic(2)
 DO  i = 1,5
   it(i)  = BLOCK(i  )
   it1(i) = BLOCK(i+6)
 END DO
 dt1(1) = mcbs(20)
 dt1(2) = mcbs(21)
 IF (nomix /= 0) WRITE (nout,92,ERR=95) it(1),dit,it1(1),dit1
 92 FORMAT ('  MULTIPLIERS =',i3,d12.3,i8,d12.3)
 95 ic(1)  = c
 lcore  = korsz(core)
 
!     DETERMINE TYPE OF OUTPUT
 
 irc = 0
 IF (ia(1) == 0) GO TO 100
 IF (ia5 > 2 .OR.  it(1) > 2) irc = 2
 100 IF (ib(1) == 0) GO TO 110
 IF (ib5 > 2 .OR. it1(1) > 2) irc = 2
 110 CONTINUE
 iprec = ipr1
 ic(5) = irc + iprec
 nomat = 2
 IF (nomix == 0) GO TO 130
 CALL fname (ic(1),na)
 WRITE  (nout,120) na,ic(2),ic(3),ic(5),ic(4)
 120 FORMAT (5X,'FILE ',2A4,'(',i6,' X',i6,') TYPE =',i3,', FORM =',i3,  &
     5X,'(RESULTANT)')
 130 CALL sadd (core,core)
 CALL wrttrl (ic)
 150 RETURN
END SUBROUTINE ssg2c
