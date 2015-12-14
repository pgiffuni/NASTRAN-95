SUBROUTINE ampd (qjhua,qhho,skj,gki,qhh,scr1,scr2,scr3,scr4)
     
    !     THE PURPOSE OF THIS ROUTINE IS TO COMPUTE(OR RETRIEVE) QHH
 
    !     QHH EITHER EXISTS ON QHHO (AS COLUMN QCOL) OR MUST BE COMPUTED
    !     AS FOLLOWS
 
    !     1. QKH = SKJ*QJH
    !     2. QIH = GKI(T)*QKH
    !     3. QHH = 1 QIH 1
    !              1-----1
    !              1 0   1
    !              1     1
 
 
    INTEGER, INTENT(IN)                      :: qjhua
    INTEGER, INTENT(IN)                      :: qhho
    INTEGER, INTENT(IN)                      :: skj
    INTEGER, INTENT(IN OUT)                  :: gki
    INTEGER, INTENT(IN OUT)                  :: qhh
    INTEGER, INTENT(IN OUT)                  :: scr1
    INTEGER, INTENT(IN OUT)                  :: scr2
    INTEGER, INTENT(IN OUT)                  :: scr3
    INTEGER, INTENT(IN)                      :: scr4
    INTEGER :: ajjcol,qhhcol,sysbuf,FILE,  NAME(2),mcb(7), qkh
    COMMON /ampcom/ ncol,nsub,xm,xk,ajjcol,qhhcol,ngp,ngpd(2,30),  &
        mcbqhh(7),mcbqjh(7),noh
    COMMON /zzzzzz/ iz(1)
    COMMON /system/ sysbuf,nout,skp(52),iprec
    COMMON /unpakx/ itc,ii,jj,incr
    COMMON /BLANK / noue
    DATA    NAME  / 4HAMPD,4H    /
 
    ibuf1 = korsz(iz) - sysbuf + 1
    ibuf2 = ibuf1 - sysbuf
    incr  = 1
    itc   = mcbqhh(5)
 
    !     DETERMINE IF QHH EXISTS ON QHHO
 
    IF (qhhcol == 0) GO TO 100
 
    !     COPY FROM QHHO TO QHH
 
    CALL gopen (qhh,iz(ibuf1),3)
    CALL gopen (qhho,iz(ibuf2),0)
    k = qhhcol - 1
    IF (k == 0) GO TO 20
    FILE = qhho
    DO  i = 1,k
        CALL fwdrec (*910,qhho)
    END DO
20 CONTINUE
   CALL cyct2b (qhho,qhh,noh,iz,mcbqhh)
   CALL CLOSE  (qhho,1)
   CALL CLOSE  (qhh,3)
   RETURN
 
   !     QHH MUST BE COMPUTED
 
100 CONTINUE
 
    !     COPY SKJ TO SCR4 FOR PROPER M-K PAIR
 
    CALL gopen (skj,iz(ibuf1),0)
    CALL gopen (scr4,iz(ibuf2),1)
    k = ajjcol - 1
    CALL skprec (skj,k)
    mcb(1) = qjhua
    CALL rdtrl (mcb)
    ncolj  = mcb(3)
    mcb(1) = skj
    CALL rdtrl (mcb)
    mcbqjh(3) = mcb(3)
    mcb(1) = scr4
    mcb(2) = 0
    mcb(6) = 0
    mcb(7) = 0
    itc = mcb(5)
    CALL cyct2b (skj,scr4,ncolj,iz,mcb)
    CALL CLOSE  (skj,1)
    CALL CLOSE  (scr4,1)
    CALL wrttrl (mcb)
    CALL ssg2b  (scr4,qjhua,0,scr1,0,iprec,1,scr2)
 
    !     COPY SCR1(QKH) TO OUTPUT
 
    qkh = mcbqjh(1)
    IF (qkh <= 0) GO TO 200
    itc  = mcbqjh(5)
    incr = 1
    CALL gopen  (scr1,iz(ibuf1),0)
    CALL gopen  (qkh,iz(ibuf2),3)
    CALL cyct2b (scr1,qkh,noh,iz,mcbqjh)
    CALL CLOSE  (qkh,3)
    CALL CLOSE  (scr1,1)
200 CONTINUE
    CALL ssg2b (gki,scr1,0,scr3,1,iprec,1,scr2)
 
    !     COPY TO QHH
 
    CALL gopen (qhh,iz(ibuf1),3)
    CALL gopen (scr3,iz(ibuf2),0)
    itc  = mcbqhh(5)
    incr = 1
    CALL cyct2b (scr3,qhh,noh,iz,mcbqhh)
    CALL CLOSE  (scr3,1)
    CALL CLOSE  (qhh,3)
    RETURN
 
    !     ERRORS
 
910 ip1 = -2
    CALL mesage (ip1,FILE,NAME)
    GO TO 910

END SUBROUTINE ampd
