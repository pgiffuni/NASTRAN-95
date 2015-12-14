SUBROUTINE ampf (skj,gkh,ajjl,qhjl,plan,imax,scr1,scr2,scr3,scr4,  &
                 scr5,scr6,scr7,scr8,scr9,scr10)
     
    !     THE PURPOSE OF THIS ROUTINE IS TO SOLVE FOR QHJL
 
    !     THE STEPS ARE AS FOLLOWS
 
    !       I.  FOR EACH M-K PAIR
 
    !           A. FIND SKJ FROM SKJ LIST
    !                                      T
    !           B.  COMPUTE  S(K) =  SKJ(K) *GKH
 
    !           C.  FOR EACH GROUP
    !                                                  G
    !               1. BREAK  S(K) INTO GROUPS  =  S(K)
 
    !               2. SOLVE FOR  RJH
    !                                                        -1     G
    !                      D.L. AND D.L. WITH BODIES RGH= AJJ  *S(K)
    !                                                        T      G
    !                      OTHERS                    RGH= AJJ  *S(K)
 
    !               3. MERGE RESULTS
 
    !                  1    G11
    !                  1 RJH  1
    !                  1------1  =   RJH(K)
    !                  1    G21
    !                  1 RJH  1
    !                  1      1
 
 
    !           D.  APPEND  RJH ONTO GROWING  QHJL
    !                1       1       1
    !                1RJH(K1)1RJH(K2)1  =  QHJL
    !                1       1       1
    !                1       1       1
 
 
    INTEGER, INTENT(IN)                      :: skj
    INTEGER, INTENT(IN OUT)                  :: gkh
    INTEGER, INTENT(IN)                      :: ajjl
    INTEGER, INTENT(IN)                      :: qhjl
    INTEGER, INTENT(IN OUT)                  :: plan
    INTEGER, INTENT(IN)                      :: imax
    INTEGER, INTENT(IN OUT)                  :: scr1
    INTEGER, INTENT(IN)                      :: scr2
    INTEGER, INTENT(IN OUT)                  :: scr3
    INTEGER, INTENT(IN OUT)                  :: scr4
    INTEGER, INTENT(IN OUT)                  :: scr5
    INTEGER, INTENT(IN OUT)                  :: scr6
    INTEGER, INTENT(IN OUT)                  :: scr7
    INTEGER, INTENT(IN OUT)                  :: scr8
    INTEGER, INTENT(IN)                      :: scr9
    INTEGER, INTENT(IN)                      :: scr10
    INTEGER :: mcb(7),sysbuf,NAME(2), ajjcol,  &
        rjh
    COMMON /ampcom/ ncolj,nsub,xm,xk,ajjcol,qhhcol,ngp,ngpd(2,30),  &
        mcbqhh(7),mcbqjh(7),noh,idjh
    COMMON /system/ sysbuf,nout,skp(52),iprec
    COMMON /zzzzzz/ z(1)
    COMMON /cdcmpx/ dum32(32),ib
    COMMON /unpakx/ itc,ii,jj,incr
    COMMON /packx / itc1,itc2,ii1,jj1,incr1
    DATA    NAME  / 4HAMPF,1H /
 
    !     INITIALIZE
 
    ibuf1 = korsz(z) - sysbuf + 1
    ibuf2 = ibuf1 - sysbuf
    iop   = 0
    itl   = 0
    DO  iloop = 1,imax
        CALL klock (its)
        CALL gopen (plan,z(ibuf1),iop)
        iop = 2
        CALL fread (plan,xm,4,1)
        CALL CLOSE (plan,2)
   
        !     FIND  SKJ(K) IN SKJL
   
        CALL gopen (skj,z(ibuf1),0)
        CALL gopen (scr1,z(ibuf2),1)
        k = ajjcol - 1
        CALL skprec (skj,k)
        mcb(1) = skj
        CALL rdtrl (mcb)
        CALL makmcb (mcb,scr1,mcb(3),mcb(4),mcb(5))
        incr = 1
        itc  = mcb(5)
        CALL cyct2b (skj,scr1,ncolj,z,mcb)
        CALL CLOSE  (skj,1)
        CALL CLOSE  (scr1,1)
        CALL wrttrl (mcb)
        !                     T
        !     MULTIPLY  SKJ(K) *GKH  ONTO SCR2
   
        CALL ssg2b (scr1,gkh,0,scr2,1,iprec,1,scr3)
   
        !     POSITION AJJL
   
        CALL gopen (ajjl,z(ibuf1),0)
        k = ajjcol - 1
        CALL skprec (ajjl,k)
        CALL CLOSE  (ajjl,2)
   
        !     SET UP TO LOOP ON CONSTANT THEORY
   
        ngps = 1
        nth  = ngpd(1,ngps)
        ncolth = 0
135     nclold = ncolth + 1
140     IF (ngps > ngp) GO TO 150
        IF (ngpd(1,ngps) /= nth) GO TO 150
        ncolth = ncolth + ngpd(2,ngps)
        ngps   = ngps + 1
        GO TO 140
150 CONTINUE
    ionce = 0
    IF (nclold == 1 .AND. ngps > ngp) ionce = 1
    !                                 G
    !     COPY AJJL(K) TO SCR1 (AJJ(K) )
   
    CALL gopen (ajjl,z(ibuf1),2)
    CALL gopen (scr1,z(ibuf2),1)
    mcb(1) = ajjl
    CALL rdtrl (mcb)
    CALL makmcb (mcb,scr1,ncolth,mcb(4),mcb(5))
    ii   = nclold
    jj   = ncolth
    ii1  = 1
    jj1  = ncolth - nclold + 1
    itc  = mcb(5)
    itc1 = itc
    itc2 = itc
    incr = 1
    incr1= 1
    CALL ampc1 (ajjl,scr1,ncolth,z,mcb)
    CALL CLOSE (ajjl,2)
    CALL CLOSE (scr1,1)
    CALL wrttrl (mcb)
    !                                   G
    !     COPY SKJ(K)  ONTO SCR3 (SKJ(K) )
   
    CALL gopen (scr2,z(ibuf1),0)
    CALL gopen (scr3,z(ibuf2),1)
    mcb(1) = scr2
    CALL rdtrl (mcb)
    CALL makmcb (mcb,scr3,ncolth,mcb(4),mcb(5))
    itc  = mcb(5)
    itc1 = itc
    itc2 = itc
    CALL ampc1 (scr2,scr3,noh,z,mcb)
    CALL CLOSE (scr2,1)
    CALL CLOSE (scr3,1)
    CALL wrttrl (mcb)
    rjh = scr10
    IF (ionce /= 0) rjh = scr9
   
    !     BRANCH ON THEORY
   
    SELECT CASE ( nth )
        CASE (    1)
            GO TO 1000
        CASE (    2)
            GO TO 2000
        CASE (    3)
            GO TO 3000
        CASE (    4)
            GO TO 4000
        CASE (    5)
            GO TO 5000
    END SELECT
   
!     DOUBLET LATTICE--D.L. WITH SLENDER BODIES
   
1000 CONTINUE
2000 CONTINUE
     !                     G
     !     DECOMPOSE AJJ(K)
   
     ib = 0
     CALL cfactr (scr1,scr4,scr5,scr6,scr7,scr8,iopt)
     CALL cfbsor (scr4,scr5,scr3,rjh,iopt)
     GO TO 1020
   
 !     OTHER THEORIES
   
3000 CONTINUE
4000 CONTINUE
5000 CONTINUE
     CALL ssg2b (scr1,scr3,0,rjh,1,iprec,1,scr4)
   
     !     COPY ACCUMULATIVELY ONTO RJH(K)
   
1020 IF (ionce /= 0) GO TO 8000
     CALL ampc2 (rjh,scr9,scr1)
     IF (ngps > ngp) GO TO 8000
     GO TO 135
   
 !     ALL GROUPS /THEORIES COMPLETE
   
8000 CONTINUE
   
     !     COPY ONTO  QHJL
   
     CALL gopen (scr9,z(ibuf1),0)
     CALL gopen (qhjl,z(ibuf2),3)
     mcb(1) = qhjl
     CALL rdtrl (mcb(1))
     itc  = mcb(5)
     incr = 1
     CALL cyct2b (scr9,qhjl,noh,z,mcb)
     CALL CLOSE  (qhjl,2)
     CALL CLOSE  (scr9,1)
     CALL wrttrl (mcb)
   
     !     END LOOP ON M-K PAIRS
   
     IF (iloop == imax) CYCLE
   
     !     CHECK TIME
   
     CALL klock (itf)
     CALL tmtogo (itmto)
     itl= MAX0(itf-its,1,itl)
     IF (1.1*itl >= itmto) GO TO 9010
 END DO
 RETURN
 
 !     INSUFFICIENT TIME TO COMPLETE
 
9010 CALL mesage (45,imax-iloop,NAME)

 RETURN
END SUBROUTINE ampf
