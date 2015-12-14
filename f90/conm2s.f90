SUBROUTINE conm2s
     
    ! THIS SUBROUTINE COMPUTES THE CONCENTRATED MASS ELEMENTS MASS MATRIX
    ! FOR THE M2 TYPE ELEMENT
    ! SINGLE PRECISION VERSION
 
    ! ECPTNO  NAME       TYPE  DESCRIPTION
    ! ******  ****       ****  ***********
 
    !   1     IELID      I     ELEMENT ID
    !   2     IGP        I     GRID POINT NUMBER
    !   3     ICIDT2     I     COORDINATE SYSTEM ID FOR T2
    !   4     MASS       R     LUMPED MASS
    !   5     OFFSET(1)  R
    !   6     OFFSET(2)  R     X,Y, AND Z COORDINATES OF THE
    !   7     OFFSET(3)  R     OFFSET
    !   8     MMI(1,1)   R
    !   9     MMI(2,1)   R     MASS MOMENTS OF INERTIA
    !  10     MMI(2,2)   R
    !  11     MMI(3,1)   R
    !  12     MMI(3,2)   R
    !  13     MMI(3,3)   R
    !  14     ICIDT1     I     COORDINATE SYSTEM ID FOR T1
    !  15     X          R
    !  16     Y          R
    !  17     Z          R
 
    INTEGER :: dict(11), elid, estid, iecpt(14)
    REAL :: mm(36),tt(36),t(36)
    REAL :: mb, iner(6)
 
    COMMON /emgest/ ecpt(100)
    COMMON /emgdic/ dmm(2),nlocs,elid,estid
    COMMON /emgprm/ dum(15),ismb(3),iprec,nogo
 
    EQUIVALENCE (ecpt(1),iecpt(1),ielid)
    EQUIVALENCE (dict(5),dict5), (ecpt(4),mb)
    EQUIVALENCE (ecpt(5),xof),(ecpt(6),yof),(ecpt(7),zof)
    EQUIVALENCE (iner(1),ecpt(8))
 
    !     INITIALIZE
 
    IF (ismb(2) == 0) RETURN
    dict(1) = estid
    dict(2) = 1
    dict(3) = 6
    dict(4) = 63
    dict(5) = 0
    ip = iprec
 
    ! COMPUTE NON-TRANSFORMED MASS MATRIX.  INITIALIZE TO ZERO
    ! THEN FILL IN NON-ZERO TERMS
 
    DO  i=1,36
        mm(i) = 0.
    END DO
 
    icidt2 = iecpt(3)
    IF (icidt2 >= 0) GO TO 120
    icidt2 = 0
    DO  i = 1,3
        ecpt (i+4) = ecpt(i+4) - ecpt(i+14)
    END DO
 
120 mm(1) = mb
    mm(5) = mb*zof
    mm(6) = -mb*yof
    mm(8) = mb
    mm(10) = -mm(5)
    mm(12) = mb*xof
    mm(15) = mb
    mm(16) = -mm(6)
    mm(17) = -mm(12)
    mm(20) = mm(10)
    mm(21) = mm(16)
    x2 = xof**2
    y2 = yof**2
    z2 = zof**2
    mm(22) = iner(1) + (y2 + z2)*mb
    mm(23) = -iner(2) + mm(6)*xof
    mm(24) =  -iner(4)+mm(10)*xof
    mm(25) = mm(5)
    mm(27) = mm(17)
    mm(28) = mm(23)
    mm(29) = iner(3) + (x2 + z2)*mb
    mm(30) = -iner(5) + mm(6)*zof
    mm(31) = mm(6)
    mm(32) = mm(12)
    mm(34) = mm(24)
    mm(35) = mm(30)
    mm(36) = iner(6) + (x2 + y2)*mb
 
    icidt1 = iecpt(14)
 
    ! PERFORM TRANSFORMATIONS.  IF CSIDS 1 AND 2 ARE EQUAL,
    ! T1 = T2 SO MASS MATRIX IS COMPLETE
 
    IF (icidt2 == icidt1) GO TO 240
    !                            T
    ! NOT EQUAL SO COMPUTE T = (T )(T )
    !                            1   2
    ! GET T1 AND T2 IF NEEDED
    it = 18
    IF (icidt1 == 0) GO TO 130
 
    CALL transs (ecpt(14),t(1))
    GO TO 140
    ! ONLY T2 NEEDED SO T = T2
130 it = 9
140 IF (icidt2 == 0) GO TO 150
    itemp = iecpt(14)
    iecpt(14) = icidt2
    CALL transs (ecpt(14),t(10))
    iecpt(14) = itemp
 
    IF(icidt1 == 0) GO TO 210
    CALL gmmats (t(1),3,3,2,t(10),3,3,0,t(19))
    GO TO 210
 
    ! HERE T2 IS IDENTITY AND T1 IS AT T(1) SO
    ! T = T1 (TRANSPOSE).  SO INSERT INTO T
    150 DO  i = 1,3
        DO  j = 1,3
            ij = 3*(i-1) + j
            ji = i + 3*(j-1) + 18
            t(ji) = t(ij)
        END DO
    END DO
 
    ! T = (T ) (T ) IS COMPLETE. INSERT IT IN THE 6X6 TRANSFORMATION MATRIX.
    !       1    2
 
    210 DO  i = 1,36
        tt(i) = 0.
    END DO
 
    DO  i = 1,3
        ij = i + it
        tt(i) = t(ij)
        tt(i + 6) = t(ij + 3)
        tt(i + 12) = t(ij + 6)
        tt(i + 21) = t(ij)
        tt(i + 27) = t(ij + 3)
        tt(i + 33) = t(ij + 6)
    END DO
    !           T
    ! FORM T*M*T  AND STORE IN MM
 
    CALL gmmats (tt(1),6,6,0,mm(1),6,6,0,t(1))
    CALL gmmats (t(1),6,6,0,tt(1),6,6,1,mm(1))
 
240 CALL emgout (mm,mm,36,1,dict,2,ip)

    RETURN
END SUBROUTINE conm2s
