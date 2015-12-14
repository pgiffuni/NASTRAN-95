SUBROUTINE dumerg
     
!     DRIVER FOR DMAP MODULE UMERGE
 
!     UMERGE   USET,PHIA,PHIO/PHIF/C,N,MAJOR/C,N,SUB0/C,N,SUB1 $
 
 INTEGER :: uset,phia,phio,phif,scr1,sub0,sub1,iabit(16), NAME(2),ib(3)
 COMMON /BLANK / major(2),sub0(2),sub1(2)
 COMMON /bitpos/ ibit(32),iabit
 DATA    NAME  / 4HUMER , 4HGE       /
 DATA    uset  , phia,phio,phif,scr1 / 101,102,103,201,301/
 
 nogo = 0
 
!     DECIDE IF CHARACTERS ARE LEGAL BIT NUMBERS
 
 ib(1) = major(1)
 ib(2) = sub0(1)
 ib(3) = sub1(1)
 
 loop30:  DO  j = 1,3
   DO  i = 1,32
     IF (ib(j) /= iabit(i)) CYCLE
     ib(j) = ibit(i)
     CYCLE loop30
   END DO
   
!     INVALID
   
   CALL mesage (59,ib(j),NAME)
   nogo = 1
 END DO loop30
 
 IF (nogo == 1) CALL mesage (-7,0,NAME)
 CALL sdr1b (scr1,phia,phio,phif,ib(1),ib(2),ib(3),uset,0,0)
 RETURN
 
END SUBROUTINE dumerg
