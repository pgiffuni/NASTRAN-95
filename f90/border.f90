SUBROUTINE border (gplst,x,u,istore,deform,b1,opcor)
     
 
    INTEGER, INTENT(IN OUT)                  :: gplst(1)
    REAL, INTENT(IN)                         :: x(3,1)
    REAL, INTENT(IN)                         :: u(2,1)
    INTEGER, INTENT(IN OUT)                  :: istore(2)
    INTEGER, INTENT(IN OUT)                  :: deform
    INTEGER, INTENT(IN OUT)                  :: b1
    INTEGER, INTENT(IN)                      :: opcor
    INTEGER :: scr2,words(2),elid, scr4
 
    DIMENSION      pt(2,3)
    COMMON /BLANK/ skip(25),scr2,scr3,scr4
    EQUIVALENCE    (words(1),nelmt),(words(2),igdpt)
 
    lcor = opcor/5 - 1
    CALL OPEN (*150,scr2,gplst(b1),0)
    CALL line (0.,0.,0.,0.,1,-1)
9   CALL fwdrec (*100,scr2)
10  CALL READ (*100,*9,scr2,iflag,1,0,m)
    IF (iflag ==  0) GO TO 100
    IF (iflag == -1) GO TO 9
    CALL fread (scr2,words,2,0)
    ie = -1
20  ie = ie + 2
    CALL READ (*100,*30,scr2,elid,1,0,m)
    CALL fread (scr2,istore(ie),2,0)
    GO TO 20
30  ione = istore(1)
    itwo = istore(2)
    IF (nelmt == 1) GO TO 50
    ie = 2*nelmt
    ie1= ie - 1
    loop37:  DO  i = 1,ie1
        IF (istore(i) == 0) CYCLE loop37
        ip1 = i + 1
        DO  j = ip1,ie
            IF (istore(i) /= istore(j)) CYCLE
            istore(i) = 0
            istore(j) = 0
            CYCLE loop37
        END DO
    END DO loop37
    j = 0
    DO  i = 1,ie
        IF (istore(i) == 0) CYCLE
        j = j + 1
        IF (j-1 > 0) THEN
            GO TO    39
        END IF
38      ione = istore(i)
        CYCLE
39      itwo = istore(i)
    END DO
    IF (j == 0) GO TO 10
50  ig = IABS(gplst(igdpt))
    IF (deform /= 0) GO TO 57
    pt(1,3) = x(2,ig)
    pt(2,3) = x(3,ig)
    GO TO 60
57  pt(1,3) = u(1,ig)
    pt(2,3) = u(2,ig)
60  ig = ione
    DO  i = 1,2
        ig = IABS(gplst(ig))
        IF (deform /= 0) GO TO 63
        pt(1,i) = x(2,ig)
        pt(2,i) = x(3,ig)
        GO TO 64
63      pt(1,i) = u(1,ig)
        pt(2,i) = u(2,ig)
64      CALL line (pt(1,i),pt(2,i),pt(1,3),pt(2,3),1,0)
        ig = itwo
    END DO
    GO TO 10
100 CALL line (0.,0.,0.,0.,1,+1)
    CALL CLOSE (scr2,1)

150 RETURN
END SUBROUTINE border
