SUBROUTINE sadd (z,dz)
     
!     TO COMPUTE MATRIX SUM WITH MULTIPLIERS
!         ACCEPTS 1 TO 5 MATRIX BLOCKS PASSED ON VIA /SADDX/
!     COMMON BLOCK /SADDX/ NOMAT,LCORE,MCBS(60),MC(7)
!         NOMAT - NUMBER OF MATRICES INPUT
!         LCORE - LENGTH OF Z ARRAY (OPEN CORE)
!         MCBS  - MATRIX CONTROL BLOCKS AND MULTIPLIERS
!                 (12 WORDS/MATRIX)
 
!                 1 - FILE NAME            7 - NOT USED
!                 2 - NUMBER OF COLUMN     8 - TYPE OF MULTIPLIER
!                 3 - NUMBER OF ROW        9 - MULTIPLIER   *  LENGTH
!                 4 - FORM OF MATRIX      10 - MULTIPLIER   *  DEPENDS
!                 5 - TYPE OF MATRIX      11 - MULTIPLIER   *  ON THE
!                 6 - MAXIMUM NUMBER OF   12 - MULTIPLIER   *  TYPE
!                     NON-ZERO ELEMENTS
 
!         MC    - MATRIX CONTROL BLOCK OF THE OUTPUT
 
 
 REAL, INTENT(OUT)                        :: z(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dz(1)
INTEGER :: END    ,eol    ,hop    ,NAME(2),one    ,prc    ,  &
    prec   ,rc     ,sysbuf ,TYPE   ,typin  ,typout
REAL :: amcb(1),alph(1)
DOUBLE PRECISION :: da(2)  ,dalph(10)      ,dmcb(1)
COMMON /packx /  typin  ,typout ,one    ,n      ,incr
COMMON /saddx /  nomat  ,lcore  ,mcbs(60)       ,mc(7)
COMMON /system/  sysbuf ,nout
COMMON /TYPE  /  prc(2) ,nwds(4),rc(4)
COMMON /zntpkx/  a(4)   ,ii     ,eol
EQUIVALENCE      (amcb(1),mcbs(9)),(alph(1),dalph(1)),  &
    (da(1) ,a(1))  ,(dmcb(1),mcbs(9)), (ntype ,mc(5)) ,(nrow,   mc(3))
DATA    NAME  /  4HSADD ,4H    /


END   = (nomat-1)*12 + 1
prec  = -nomat*2
TYPE  = -nomat*2

!     DETERMINE PRECISION TO BE USED FOR CALCULATIONS

!     NOTE - PRC ARRAY IS DIMENSIONED ONLY TO 2
!            PRC(1) = 1, PRC(2) = 2, AND
!            PRC(3) = NWDS(1) = 1, PRC(4) = NWDS(2) = 2
!            WHERE 1 MEANS S.P., 2 D.P.
!          - RC ARRAY = 1,1,2,2, WHERE 1 MEANS REAL, 2 COMPLEX

DO  i = 1,END,12
  IF (mcbs(i) /= 0) GO TO 10
  prec  = prec + 2
  TYPE  = TYPE + 2
  CYCLE
  10 j     = mcbs(i+4)
  prec  = prec  + prc(j)
  TYPE  = TYPE  + rc(j)
  j     = mcbs(i+7)
  prec  = prec  + prc(j)
  TYPE  = TYPE  + rc(j)
END DO
typin = 1
IF (TYPE > 0) typin = 3
IF (prec > 0) typin = typin + 1
num = nrow*nwds(typin)
IF (lcore < (nomat+1)*sysbuf+num+1) CALL mesage (-8,0,NAME)

!     MOVE AND CONVERT MULTIPLIERS

IF (prec > 0) GO TO 60

!     SINGLE PRECISION

j = 1
DO  i = 1,END,12
  k = mcbs(i+7)
  IF (prc(k) == 2) GO TO 30
  alph(j  ) = amcb(i  )
  alph(j+1) = amcb(i+1)
  GO TO 40
  30 k = i/2 + 1
  alph(j  ) = dmcb(k  )
  alph(j+1) = dmcb(k+1)
  40 j = j + 1
  IF (TYPE > 0) j = j + 1
END DO
IF (TYPE <= 0) alph(j+1) = 0.0
GO TO 100

!     DOUBLE PRECISION

60 j = 1
DO  i = 1,END,12
  k = mcbs(i+7)
  IF (prc(k) == 2) GO TO 70
  dalph(j  ) = amcb(i  )
  dalph(j+1) = amcb(i+1)
  GO TO 80
  70 k = i/2 + 1
  dalph(j  ) = dmcb(k  )
  dalph(j+1) = dmcb(k+1)
  80 j = j + 1
  IF (TYPE > 0) j = j + 1
END DO
IF (TYPE <= 0) dalph(j+1) = 0.0D+0

100 SELECT CASE ( typin )
  CASE (    1)
    GO TO 110
  CASE (    2)
    GO TO 120
  CASE (    3)
    GO TO 130
  CASE (    4)
    GO TO 140
END SELECT
110 ASSIGN 300 TO hop
GO TO  150
120 ASSIGN 350 TO hop
GO TO  150
130 ASSIGN 400 TO hop
GO TO  150
140 ASSIGN 450 TO hop

!     OPEN AND ASSIGN FILES

150 ibuf = lcore
DO  i = 1,END,12
  ibuf = ibuf - sysbuf
  IF (mcbs(i) == 0) CYCLE
  CALL gopen (mcbs(i),z(ibuf),0)
END DO
ibuf = ibuf - sysbuf
CALL gopen (mc,z(ibuf),1)

!     SETUP PACK PARAMETERS

one    = 1
n      = nrow
typout = ntype
incr   = 1
ncol1  = mc(2)
mc(2)  = 0
mc(6)  = 0
mc(7)  = 0

!     ADD MATRICES

DO  i = 1,ncol1
  
!     CLEAR CORE
  
  DO  j = 1,num
    z(j) = 0.0
  END DO
  
  one = n
  n   = 1
  DO  j = 1,nomat
    k   = 12*(j-1) + 1
    IF (mcbs(k  ) == 0) CYCLE
    IF (mcbs(k+1) < i) CYCLE
    CALL intpk (*900,mcbs(k),0,typin,0)
    
!     READ IN NON ZERO ELEMENT
    
    220 CALL zntpki
    IF (ii > nrow) GO TO 500
    one = MIN0(one,ii)
    n   = MAX0(n  ,ii)
    GO TO hop, (300,350,400,450)
    300 z(ii) = z(ii) + alph(j)*a(1)
    GO TO 500
    350 dz(ii) = dz(ii) + dalph(j)*da(1)
    GO TO 500
    400 ii = ii + ii - 1
    jj = j  + j  - 1
    z(ii  ) = z(ii)  +  alph(jj)*a(1) - alph(jj+1)*a(2)
    z(ii+1) = z(ii+1)+  alph(jj)*a(2) + alph(jj+1)*a(1)
    GO TO 500
    450 ii = ii + ii - 1
    jj = j  + j  - 1
    dz(ii  ) = dz(ii  ) + dalph(jj)*da(1) - dalph(jj+1)*da(2)
    dz(ii+1) = dz(ii+1) + dalph(jj)*da(2) + dalph(jj+1)*da(1)
    500 IF (eol == 0) GO TO 220
  END DO
  
!     END OF COLUMN
  
  one = MIN0(one,n)
  ll  = (one-1)*nwds(typin) + 1
  CALL pack (z(ll),mc(1),mc)
END DO

!     DONE - CLOSE FILES AND RETURN

DO  i = 1,END,12
  IF (mcbs(i) /= 0) CALL CLOSE (mcbs(i),1)
END DO
CALL CLOSE (mc,1)
RETURN
END SUBROUTINE sadd
