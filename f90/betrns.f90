SUBROUTINE betrns (tbe,gg,kflag,elid)
    !     &    ENTRY BETRND (TBD,GD,KFLAG,ELID)
 
    !*****
    !     SUBROUTINE WHICH CALCULATES THE TBE TRANSFORMATION
    !     MATRIX WHICH RELATES THE ELEMENT TO THE BASIC C.S.
 
    !     GG(9) OR GD(9) IS A 9X1 ARRAY WHICH STORES THE GRID PT. COORD.
    !     X(G1),Y(G1),Z(G1),X(G2),Y(G2),Z(G2),X(G3),Y(G3),Z(G3)
    !     GG(1),GG(2),GG(3),GG(4),GG(5),GG(6),GG(7),GG(8),GG(9), OR
    !     GD(1),GD(2),GD(3),GD(4),GD(5),GD(6),GD(7),GD(8),GD(9)
 
    !     KFLAG = 0, TBE (OR TBD) IS OUTPUT WITHOUT TRANSPOSING
    !           = 1, TBE (OR TBD) IS OUTPUT AFTER IT IS TRANSPOSED
    !*****
 
 
    REAL, INTENT(OUT)                        :: tbe(9)
    REAL, INTENT(IN)                         :: gg(9)
    INTEGER, INTENT(IN OUT)                  :: kflag
    INTEGER, INTENT(IN OUT)                  :: elid
 
    REAL :: rsstr(3),rstr(3),r12(3),rv(3),LEN
    DOUBLE PRECISION :: gd(9),tbd(9),dsstr(3),dstr(3),d12(3),dv(3), dtemp,led
 
    !     SINGLE PRECISION VERSION
 
    !*****
    !     CALCULATE APPROPRIATE LENGTH QUANTITIES
    !*****
    LEN = SQRT((gg(4)-gg(1))**2 + (gg(5)-gg(2))**2 +(gg(6)-gg(3))**2)
    IF (LEN == 0.0) GO TO 40
    !*****
    !     CALCULATE APPROPRIATE VECTOR QUANTITIES
    !*****
    r12(1) = (gg(4)-gg(1))/LEN
    r12(2) = (gg(5)-gg(2))/LEN
    r12(3) = (gg(6)-gg(3))/LEN
    rv(1)  = (gg(7)-gg(1))
    rv(2)  = (gg(8)-gg(2))
    rv(3)  = (gg(9)-gg(3))
    !*****
    !     CALCULATE ENTRIES INTO THE TRANSFORMATION MATRIX
    !*****
    rstr(1) = (r12(2)*rv(3) - r12(3)*rv(2))
    rstr(2) = (r12(3)*rv(1) - r12(1)*rv(3))
    rstr(3) = (r12(1)*rv(2) - r12(2)*rv(1))
 
    LEN = SQRT (rstr(1)**2 + rstr(2)**2 + rstr(3)**2)
    IF (LEN == 0.0) GO TO 40
    DO  i=1,3
        rstr(i) = rstr(i)/LEN
    END DO
 
    rsstr(1)= (rstr(2)*r12(3) - rstr(3)*r12(2))
    rsstr(2)= (rstr(3)*r12(1) - rstr(1)*r12(3))
    rsstr(3)= (rstr(1)*r12(2) - rstr(2)*r12(1))
    tbe(1)  = r12(1)
    tbe(2)  = r12(2)
    tbe(3)  = r12(3)
    tbe(4)  = rsstr(1)
    tbe(5)  = rsstr(2)
    tbe(6)  = rsstr(3)
    tbe(7)  = rstr(1)
    tbe(8)  = rstr(2)
    tbe(9)  = rstr(3)
    IF (kflag == 0) GO TO 30
    !*****
    !     TRANSPOSE TBE(9) SINCE KFLAG.NE.ZERO
    !*****
    temp   = tbe(2)
    tbe(2) = tbe(4)
    tbe(4) = temp
    temp   = tbe(3)
    tbe(3) = tbe(7)
    tbe(7) = temp
    temp   = tbe(6)
    tbe(6) = tbe(8)
    tbe(8) = temp
    GO TO 30
 
    ENTRY betrnd (tbd,gd,kflag,elid)
    !     ================================
 
    !     DOUBLE PRECISION VERSION
 
    !*****
    !     CALCULATE APPROPRIATE LENGTH QUANTITIES
    !*****
    led = DSQRT((gd(4)-gd(1))**2 + (gd(5)-gd(2))**2 +(gd(6)-gd(3))**2)
    IF (led == 0.0D+0) GO TO 40
    d12(1) = (gd(4)-gd(1))/led
    d12(2) = (gd(5)-gd(2))/led
    d12(3) = (gd(6)-gd(3))/led
    dv(1)  = (gd(7)-gd(1))
    dv(2)  = (gd(8)-gd(2))
    dv(3)  = (gd(9)-gd(3))
    !*****
    !     CALCULATE ENTRIES INTO THE TRANSFORMATION MATRIX
    !*****
    dstr(1) = (d12(2)*dv(3) - d12(3)*dv(2))
    dstr(2) = (d12(3)*dv(1) - d12(1)*dv(3))
    dstr(3) = (d12(1)*dv(2) - d12(2)*dv(1))
 
    led = DSQRT(dstr(1)**2 + dstr(2)**2 + dstr(3)**2)
    IF (led == 0.0D+0) GO TO 40
    DO  i=1,3
        dstr(i) = dstr(i)/led
    END DO
 
    dsstr(1)= (dstr(2)*d12(3) - dstr(3)*d12(2))
    dsstr(2)= (dstr(3)*d12(1) - dstr(1)*d12(3))
    dsstr(3)= (dstr(1)*d12(2) - dstr(2)*d12(1))
    tbd(1)  = d12(1)
    tbd(2)  = d12(2)
    tbd(3)  = d12(3)
    tbd(4)  = dsstr(1)
    tbd(5)  = dsstr(2)
    tbd(6)  = dsstr(3)
    tbd(7)  = dstr(1)
    tbd(8)  = dstr(2)
    tbd(9)  = dstr(3)
    IF (kflag == 0) GO TO 30
    !*****
    !     TRANSPOSE TBD(9) SINCE KFLAG.NE.ZERO
    !*****
    dtemp  = tbd(2)
    tbd(2) = tbd(4)
    tbd(4) = dtemp
    dtemp  = tbd(3)
    tbd(3) = tbd(7)
    tbd(7) = dtemp
    dtemp  = tbd(6)
    tbd(6) = tbd(8)
    tbd(8) = dtemp
30  RETURN
    !*****
    !     ZERO LENGTH ERROR, BAD GEOMETRY
    !*****
40  CALL mesage (-30,31,elid)
    RETURN
END SUBROUTINE betrns
