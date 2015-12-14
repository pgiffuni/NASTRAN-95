SUBROUTINE dadd
     
    !     DMAP DRIVER FOR ADD--
 
    !     ADD    A,B/C/V,N,ALPHA/V,N,BETA/V,N,DALPHA/V,N,DBETA/V,N,ECHO  $
 
    !            MATRIX C = ALPHA*MATRIX A + BETA*MATRIX B
 
    !     MATRIX C IS COMPLEX IF ANY ONE OF THE MATRIX A, MATRIX B, SCALE
    !     ALPHA, OR SCLAE BETA IS COMPLEX
 
    LOGICAL :: dblea    ,dbleb
    INTEGER :: fn(2)    ,echo     ,aa(2)     ,bb(2)
    DOUBLE PRECISION :: dalpha   ,dbeta    ,dalp(2)   ,dbta(2)  ,  &
        zero     ,one      ,xx
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg /   ufm      ,uwm      ,uim
    COMMON /system/   ibuf     ,nout
    COMMON /BLANK /   alpha(2) ,beta(2)  ,dalpha(2) ,dbeta(2) , echo
    COMMON /saddx /   nomat    ,lcore    ,ia(7)     ,ita      ,  &
        alp(4)   ,ib(7)    ,itb       ,bta(4)   , cde(3,12),ic(7)
    COMMON /zzzzzz/   core(1)
    EQUIVALENCE       (alp(1),dalp(1))   ,(bta(1),dbta(1))
    DATA              in1,in2,iout1,zero /101,102,201, 0.0D+0 /
    DATA              one,xx ,x    / 1.0D+0, 1.0D+37, 1.0E+37 /
 
 
    !     SCALE FACTORS ALPHA, DALPHA, BETA AND DBETA WERE INITIALLY SET TO
    !     (1.1+37, 1.1+37) BY XMPLDD
 
    !     IN THIS ROUTINE -
    !     IF ALPHA, DALPHA, BETA AND DBETA ARE NOT SPECIFIED BY USER, THEY
    !     WILL BE SET TO -
    !     ALPHA AND DALPHA TO (1.0, 0.0), AND
    !     BETA  AND DBETA  TO (1.0, 0.0), SAME DEFAULTS AS 88 AND EARLIER
    !                                     NASTRAN VERSIONS.
    !     NOTE - DEFAULTS WERE ALL ZEROS IN 89 NASTRAN VERSION
 
    !     NOTE - THIS ROUTINE WILL CALL SADD TO DO THE ACTUAL MATRIX MULTI-
    !     PLICATION, WHICH WILL AUTOMATICALLY ADJUST THE SCALE FACTORS
    !     WHETHER THEY ARE S.P. OR D.P. (E.G. S.P. ALPHA AND BETA CAN BE
    !     USED FOR D.P. A AND B MATRICES, AND VISE VERSA)
 
    CALL fname (iout1,fn(1))
    lcore = korsz(core)
    DO  i = 1,7
        ia(i) = 0
        ib(i) = 0
        ic(i) = 0
    END DO
    ia(1) = in1
    ib(1) = in2
    CALL rdtrl (ia)
    CALL rdtrl (ib)
    IF (ia(1) < 0) ia(1) = 0
    IF (ib(1) < 0) ib(1) = 0
    IF (ia(1)+ib(1) == 0) GO TO 100
 
    !     SET DEFAULT VALUES FOR THE SCALE FACTORS
 
    !     WHEN AN ITEM IS .LT. X OR XX, THAT ITEM HAS INPUT FROM USER
 
    dblea = .true.
    dbleb = .true.
    IF (alpha(1) < x .OR. alpha(2) < x .OR. dalpha(1) < xx .OR.  &
        dalpha(2) < xx) GO TO 20
    alp(1)   = 1.0
    alp(2)   = 0.0
    alpha(1) = 1.0
    alpha(2) = 0.0
    dblea  = .false.
20  IF (beta(1) < x .OR. beta(2) < x .OR. dbeta(1) < xx .OR.  &
        dbeta(2) < xx) GO TO 25
    bta(1)  = 1.0
    bta(2)  = 0.0
    beta(1) = 1.0
    beta(2) = 0.0
    dbleb  = .false.
    IF (.NOT.dblea) GO TO 40
 
25  IF ((alpha(1) < x .OR. alpha(2) < x) .AND. (dalpha(1) < xx .OR.  &
        dalpha(2) < xx)) GO TO 120
    IF (( beta(1) < x .OR.  beta(2) < x) .AND. ( dbeta(1) < xx .OR.  &
        dbeta(2) < xx)) GO TO 120
 
    IF (dalpha(1) > xx .AND. dalpha(2) > xx) dblea = .false.
    IF ( dbeta(1) > xx .AND.  dbeta(2) > xx) dbleb = .false.
 
    DO  i = 1,2
        IF ( alpha(i) >  x)  alpha(i) = 0.0
        IF (dalpha(i) > xx) dalpha(i) = zero
        IF (  beta(i) >  x)   beta(i) = 0.0
        IF ( dbeta(i) > xx)  dbeta(i) = zero
    END DO
 
    !     MOVE ALPHA, BETA, DALPHA AND DBETA INTO ALP AND BTA ARRAYS FOR
    !     MATRIX MULTIPLICATION TO BE PERFORMED IN SADD.
 
    DO  i = 1,2
        IF (.NOT.dblea)  alp(i) =  alpha(i)
        IF (.NOT.dbleb)  bta(i) =   beta(i)
        IF (     dblea) dalp(i) = dalpha(i)
        IF (     dbleb) dbta(i) =  dbeta(i)
    END DO
 
40  IF (echo == 0) GO TO 55
    WRITE  (nout,45) uim,fn
45  FORMAT (a29,', SCALE FACTORS FOR THE OUTOUT DATA BLOCK ',2A4,  &
        ', IN ADD MODULE ARE -')
    IF (.NOT.dblea) WRITE (nout,50) alp(1) ,alp(2)
    IF (     dblea) WRITE (nout,51) dalp(1),dalp(2)
    IF (.NOT.dbleb) WRITE (nout,52) bta(1) ,bta(2)
    IF (     dbleb) WRITE (nout,53) dbta(1),dbta(2)
50  FORMAT (5X,'1ST S.F. = (',e12.5,1H,,e12.5,1H))
51  FORMAT (5X,'3RD S.F. = (',d12.5,1H,,d12.5,1H))
52  FORMAT (1H+,48X,'2ND S.F. = (',e12.5,1H,,e12.5,1H))
53  FORMAT (1H+,48X,'4TH S.F. = (',d12.5,1H,,d12.5,1H))
 
    !     ENSURE THAT THE MATRICES BEING ADDED ARE OF THE SAME ORDER
 
55  IF (ia(1) == 0 .OR. ib(1) == 0) GO TO 70
    IF (ia(2) == ib(2) .AND. ia(3) == ib(3)) GO TO 70
    CALL fname (ia(1),aa)
    CALL fname (ib(1),bb)
    WRITE  (nout,60) ufm,aa,bb,fn,ia(2),ia(3),ib(2),ib(3)
60  FORMAT (a23,' 4149, ATTEMPT TO ADD MATRICES OF UNEQUAL ORDER IN',  &
        ' MODULE ADD, ',2A4,' TO ',2A4, /5X,'INTENDED OUTOUT DATA',  &
        ' BLOCK NAME =',2A4,i7,' BY',i6,' TO',i7,' BY',i6)
    GO TO 160
70  ic(1) = iout1
    ic(2) = ia(2)
    ic(3) = ia(3)
    IF (ia(4) == 3) ic(2) = ia(3)
    IF (ia(1) /= 0) GO TO 80
    ic(2) = ib(2)
    ic(3) = ib(3)
 
    !     DETERMINE TYPE
 
80  ita = 3
    itb = 3
    IF (alp(2) == 0.0 .AND. alp(4) == 0.0) ita = 1
    IF (bta(2) == 0.0 .AND. bta(4) == 0.0) itb = 1
    ic(5) = MAX0(ia(5),ib(5),ita,itb)
    IF (ic(5) == 3 .AND. (ia(5) == 2 .OR. ib(5) == 2)) ic(5) = 4
 
    !     DETERMINE FORM
 
    ic(4) = ia(4)
    IF (ia(1) == 0) ic(4) = ib(4)
    IF (ic(4) /= 1 .OR. ic(4) /= 6) GO TO 90
    ic(4) = 6
    IF (ia(1) /= 0 .AND. ia(4) /= 6) ic(4) = 1
    IF (ib(1) /= 0 .AND. ib(4) /= 6) ic(4) = 1
    IF (ic(2) /= ic(3)) ic(4) = 2
90  IF (ia(4) == 3 .AND. ib(1) /= 0) ic(4) = ib(4)
    IF (ia(4) == 3 .AND. ib(1) == 0) ic(4) = ia(4)
 
    nomat = 2
    CALL sadd (core,core)
    CALL wrttrl (ic)
    GO TO 170
 
100 WRITE  (nout,110) ufm,fn
110 FORMAT (a23,', INPUT MATRICES NOT SPECIFIED IN ADD MODULE.',  &
        ' INTENDED OUTPUT DATA BLOCK NAME =',2A4)
    GO TO 160
 
    120 DO  i=1,2
        IF ( alpha(i) >  x) alpha(i)  = 0.0
        IF (dalpha(i) > xx) dalpha(i) = zero
        IF (  beta(i) >  x)   beta(i) = 0.0
        IF ( dbeta(i) > xx)  dbeta(i) = zero
    END DO
    WRITE  (nout,150) ufm,fn,alpha,beta,dalpha,dbeta
150 FORMAT (a23,' IN ADD MODULE. INTENDED OUTPUT DATA BLOCK =',2A4,  &
        /5X,'SCALE FACTORS ARE ERRONEOUS =',4E9.2,2X,4D10.3)
160 CALL mesage (-61,0,0)
 
170 RETURN
END SUBROUTINE dadd
