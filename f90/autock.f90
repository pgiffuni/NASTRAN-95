SUBROUTINE autock (iadd)
     
    !     THIS ROUTINE GENERATES A CHKPT OSCAR RECORD WHEN THE PRECHK
    !     OPTION IS BEING USED, THE ADDRESS IADD IS THE STARTING
    !     LOCATION OF THE OUTPUT FILE NAMES TO BE TESTED
 
 
    INTEGER, INTENT(IN)                      :: iadd
    EXTERNAL        lshift,rshift
    INTEGER :: prenam,preflg,osprc,osbot,ospnt,list(100),xchk(2),  &
        dmpcnt,rshift,casess(2),casecc(2),casei(2), oscar(1),os(5)
    COMMON /autocm/ preflg,nnames,prenam(100)
    COMMON /zzzzzz/ core(1)
    COMMON /xgpi4 / junk(2),iseqn,dmpcnt
    COMMON /autohd/ ihead
    EQUIVALENCE     (loscar,core(1),os(1)),(osprc,os(2)),  &
        (osbot,os(3)),(ospnt,os(4)),(oscar(1),os(5))
    DATA    casess/ 4HCASE, 4HSS  /
    DATA    casecc/ 4HCASE, 4HCC  /
    DATA    casei / 4HCASE, 4HI   /
    DATA    xchk  / 4HXCHK, 4H    /
    DATA    iblank/ 0             /
 
    ihead = 0
    iop   = 0
    IF (preflg > 0.0) THEN
        GO TO     5
    ELSE
        GO TO   200
    END IF
5   preflg = IABS(preflg)
    nopf = oscar(iadd)
    nwd  = 3*nopf
    ist  = iadd + 1
    ifin = ist + nwd - 1
    nlist= 0
    incr = 3
 
    SELECT CASE ( preflg )
        CASE (    1)
            GO TO 1
        CASE (    2)
            GO TO 2
        CASE (    3)
            GO TO 3
    END SELECT
 
    !     CHECK OUTPUT FILE AGAINST LIST
 
1   n2 = 2*nnames
    DO  i = ist,ifin,incr
        DO  j = 1,n2,2
            IF (prenam(j) == casess(1) .AND.prenam(j+1) == casess(2)) CYCLE
            IF (prenam(j) == casecc(1) .AND.prenam(j+1) == casecc(2)) CYCLE
            IF (prenam(j) == casei( 1) .AND.prenam(j+1) == casei( 2)) CYCLE
            IF (prenam(j) /= oscar(i) .OR. prenam(j+1) /= oscar(i+1)) CYCLE
            nlist = nlist + 1
            list(2*nlist-1) = oscar(i  )
            list(2*nlist  ) = oscar(i+1)
        END DO
    END DO
    IF (iop   == 1) GO TO 300
    IF (nlist == 0) RETURN
    GO TO 100
 
    !     PREFLG=ALL OPTION, CHECKPOINT ALL OUTPUT DATA BLOCKS
 
    2 DO  i = ist,ifin,incr
        IF (oscar(i) == iblank .AND. oscar(i+1) == iblank) CYCLE
        IF (oscar(i) == casess(1) .AND. oscar(i+1) == casess(2)) CYCLE
        IF (oscar(i) == casecc(1) .AND. oscar(i+1) == casecc(2)) CYCLE
        IF (oscar(i) == casei( 1) .AND. oscar(i+1) == casei( 2)) CYCLE
        nlist = nlist + 1
        list(2*nlist-1) = oscar(i  )
        list(2*nlist  ) = oscar(i+1)
    END DO
    IF (iop == 1) GO TO 300
    GO TO 100
 
    !     CHECK OUTPUT FILES EXCEPT THOSE IN LIST
 
3   n2 = 2*nnames
    loop30:  DO  i = ist,ifin,incr
        DO  j = 1,n2,2
            IF (prenam(j) == oscar(i) .AND.prenam(j+1) == oscar(i+1)) CYCLE loop30
        END DO
        IF (oscar(i) == iblank .AND. oscar(i+1) == iblank ) CYCLE loop30
        IF (oscar(i) == casess(1) .AND. oscar(i+1) == casess(2)) CYCLE loop30
        IF (oscar(i) == casecc(1) .AND. oscar(i+1) == casecc(2)) CYCLE loop30
        IF (oscar(i) == casei( 1) .AND. oscar(i+1) == casei( 2)) CYCLE loop30
        nlist = nlist + 1
        list(2*nlist-1) = oscar(i  )
        list(2*nlist  ) = oscar(i+1)
    END DO loop30
    IF (iop   == 1) GO TO 300
    IF (nlist == 0) RETURN
    GO TO 100
 
    !     PURGE OR EQUIV DATA BLOCK LIST MUST BE CHECKED
 
200 nwd = oscar(ospnt)
    mi  = rshift(oscar(ospnt+2),16)
    ib  = ospnt + 6
    preflg = IABS(preflg)
    ndb = oscar(ib)
    iop = 1
    IF (mi ==  9) ist  = ib + 1
    IF (mi == 10) ist  = ib + 4
    IF (mi ==  9) ifin = ist + 2*ndb - 1
    IF (mi == 10) ifin = ist + 2*ndb - 3
    nwd  = nwd - 6
    incr = 2
    nlist= 0
    SELECT CASE ( preflg )
        CASE (    1)
            GO TO 1
        CASE (    2)
            GO TO 2
        CASE (    3)
            GO TO 3
    END SELECT
300 nwd = nwd - 2*ndb - 2
    IF (mi == 10) nwd = nwd - 1
    IF (nwd <= 0 .AND. nlist /= 0) GO TO 100
    IF (nwd <= 0 .AND. nlist == 0) GO TO 999
    ndb = oscar(ifin+2)
    IF (mi ==  9) ist = ifin + 3
    IF (mi == 10) ist = ifin + 6
    ifin = ist + 2*ndb - 1
    IF (mi == 10) ifin = ifin - 2
    SELECT CASE ( preflg )
        CASE (    1)
            GO TO 1
        CASE (    2)
            GO TO 2
        CASE (    3)
            GO TO 3
    END SELECT
 
    !     UPDATE OSCAR PARAMETERS
 
100 ihead = 1
    osprc = osbot
    osbot = oscar(osbot) + osbot
    ospnt = osbot
    iseqn = oscar(osprc+1) + 1
 
    !     LOAD HEADER
 
    oscar(ospnt  ) = 6
    oscar(ospnt+1) = iseqn
    oscar(ospnt+2) = 4 + lshift(3,16)
    oscar(ospnt+3) = xchk(1)
    oscar(ospnt+4) = xchk(2)
    oscar(ospnt+5) = dmpcnt
    IF (iop == 1) oscar(ospnt+5) = oscar(ospnt+5) - 1
    oscar(ospnt+6) = nlist
    CALL xlnkhd
    IF (nlist == 0) GO TO 110
 
    !     LOAD CHKPNT INFORMATION
 
    nlist = 2*nlist
    DO  i = 1,nlist,2
        oscar(ospnt+6+i) = list(i)
        oscar(ospnt+7+i) = list(i+1)
    END DO
110 oscar(ospnt) = oscar(ospnt) + nlist + 1
999 ihead = 0

    RETURN
END SUBROUTINE autock
