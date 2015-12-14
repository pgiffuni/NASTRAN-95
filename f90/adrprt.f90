SUBROUTINE adrprt (casecc,pkf,spline,sila,useta,freq,nfreq,  &
    ncore,nload)
     
    !     ADRPRT FORMATS PKF BY USER SET REQUEST FOR EACH FREQUENCY
 
 
    INTEGER, INTENT(IN)                      :: casecc
    INTEGER, INTENT(IN)                      :: pkf
    INTEGER, INTENT(IN OUT)                  :: spline
    INTEGER, INTENT(IN)                      :: sila
    INTEGER, INTENT(IN OUT)                  :: useta
    REAL, INTENT(IN OUT)                     :: freq(1)
    INTEGER, INTENT(IN)                      :: nfreq
    INTEGER, INTENT(IN)                      :: ncore
    INTEGER, INTENT(IN)                      :: nload
    EXTERNAL        andf
    INTEGER :: sysbuf,out,nam(2), lsp(2),andf,setno,all,extid,trl(7)
    REAL :: z(1),buf(12),tsave(96)
    COMMON /system/ sysbuf,out,dum1(6),nlpp
    COMMON /unpakx/ ito,ii,nn,incr
    COMMON /two   / itwo(32)
    COMMON /bitpos/ ibit(64)
    COMMON /output/ head(96)
    COMMON /zzzzzz/ iz(1)
    EQUIVALENCE     (z(1),iz(1))
    DATA    iaero / 176 /, lcs /200/
    DATA    nhfssu, nam /  4HFSSU,4HADRP,4HRT  /
    DATA    lsp   / 200,2/
 
    !     CORE LAYOUT
    !       FREQ LIST          NFREQ
    !       SPLINE TRIPLETS    3*K POINTS
    !       SILS FOR K POINTS  1 PER K
    !       USET MASKS         6*K POINTS
    !       CASECC RECORD      TRL(4) LONG
    !       LOAD VECTOR        K SIZE
    !       BUFFERS            2 * SYSBUF
 
    mask = ibit(19)
    mask = itwo(mask)
    DO  i = 1,96
        tsave(i) = head(i)
    END DO
    ibuf1 = ncore - sysbuf - 1
    ibuf2 = ibuf1 - sysbuf
    izspl = nfreq
    nr = ibuf2 - izspl
    CALL preloc (*1000,z(ibuf1),spline)
    CALL locate (*1000,z(ibuf1),lsp,dum)
    CALL READ (*1000,*10,spline,z(izspl+1),nr,0,nwr)
    GO TO 999
10  nspl  = nwr
    izsil = izspl + nwr
    CALL CLOSE (spline,1)
 
    !     FIND SMALLEST SILGA POINTER (-1+NEXTRA = NSKIP ON SILA)
 
    ismal = 1000000
    DO  i = 1,nspl,3
        ismal = MIN0(ismal,iz(izspl+i+1))
    END DO
    ismal = ismal - 1
    trl(1)= sila
    CALL rdtrl (trl)
    IF (trl(1) < 0) GO TO 1000
    nextra = trl(3)
    CALL gopen (sila,z(ibuf1),0)
    nskip = ismal + nextra
    nr = ibuf2 - izsil
    CALL READ (*1000,*1000,sila,z(izsil+1),-nskip,0,nwr)
    CALL READ (*1000,*30,sila,z(izsil+1),nr,0,nwr)
    GO TO 999
30  nsil = nwr
    CALL CLOSE (sila,1)
    izuset = izsil + nwr
    nr = ibuf2 - izuset
    nskip = iz(izsil+1) -1
    CALL OPEN (*1000,useta,z(ibuf1),0)
    CALL fwdrec (*1000,useta)
    CALL READ (*1000,*1000,useta,z(izuset+1),-nskip,0,nwr)
    CALL READ (*1000,*40,useta,z(izuset+1),nr,0,nwr)
    GO TO 999
40  icc = izuset + nwr
    CALL CLOSE (useta,1)
 
    !     ADJUST SILA AND USET POINTERS FOR SHRUNKEN LISTS
 
    DO  i = 1,nspl,3
        iz(izspl+i+1) = iz(izspl+i+1) - ismal
    END DO
    DO  i = 1,nsil
        iz(izsil+i) = iz(izsil+i) - nskip
    END DO
    CALL bug (nhfssu,60,z,icc)
    trl(1) = casecc
    CALL rdtrl (trl)
    lcc = trl(4) + 1
    izvect = icc + lcc
    trl(1) = pkf
    CALL rdtrl (trl)
    ito = 3
    ii  = 1
    nn  = trl(3)
    incr  = 1
    nvect = trl(3)*2
    iend  = izvect + nvect
    IF (iend > ibuf2) GO TO 999
    CALL OPEN (*1000,casecc,z(ibuf1),0)
    CALL fwdrec (*1000,casecc)
    CALL OPEN (*1000,pkf,z(ibuf2),0)
    CALL fwdrec (*1000,pkf)
 
    !     LOOP OVER NLOAD (CASECC RECORDS)
    !     THEN LOOP OVER NFREQ  (PKF COLUMNS)
    !     OUTPUT K POINTS FOR SET LIST
 
    DO  k = 1,nload
        CALL READ (*1000,*65,casecc,z(icc+1),lcc,1,nwr)
65      setno = iz(icc+iaero)
        all = 0
        DO  i = 1,96
            head(i) = z(icc+i+38)
        END DO
        IF(setno < 0.0) THEN
            GO TO    70
        ELSE IF (setno == 0.0) THEN
            GO TO   250
        ELSE
            GO TO    80
        END IF
70      all = 1
        GO TO 100
80      isetno = lcs + iz(icc+lcs) + 1 + icc
90      iset = isetno + 2
        nset = iz(isetno+1) + iset - 1
        IF (iz(isetno) == setno) GO TO 100
        isetno = nset +1
        IF (isetno < izvect) GO TO 90
        all = 1
        100 DO  j = 1,nfreq
            nlppp = nlpp
            CALL unpack (*110,pkf,z(izvect+1))
            GO TO 120
110         CALL zeroc (z(izvect+1),nvect)
     
            !     PRINT LOOP
     
120         IF (all == 0) GO TO 150
            ASSIGN 140 TO iret
            l = 1
            GO TO 181
140         l = l + 3
            IF (l >= nspl) CYCLE
            GO TO 181
150         i = iset
155         IF (i == nset) GO TO 170
            IF (iz(i+1) > 0) GO TO 170
            id = iz(i  )
            n  =-iz(i+1)
            i  = i+1
            ASSIGN 160 TO iret1
            GO TO 180
160         id = id + 1
            IF (id <= n) GO TO 180
            GO TO 175
170         id = iz(i)
            ASSIGN 175 TO iret1
            GO TO 180
175         i = i + 1
            IF (i <= nset) GO TO 155
            CYCLE
     
            !     LOCATE ELEMENT THEN  PRINT DATA
     
180         ASSIGN 190 TO iret
            CALL bisloc (*190,id,iz(izspl+1),3,nspl/3,l)
181         extid = iz(izspl+l)
            ipsil = iz(izspl+l+1)
            irow  = iz(izspl+l+2) *2 - 1 + izvect
            ipuset= iz(izsil+ipsil) + izuset - 1
            GO TO 200
190         GO TO iret1, (160,175)
     
            !     PRINT
     
200         IF (nlppp < nlpp) GO TO 210
            CALL page1
            WRITE  (out,201) j,freq(j)
201         FORMAT (44X,42HAERODYNAMIC loads  (UNIT dynamic pressure),  /  &
                30X,7HVECTOR ,i8,10X,12HFREQUENCY = ,1P,e14.6,7H  hertz,  /,  &
                11H box OR    ,12X,7HT1 / r1,23X,7HT2 / r2,23X,7HT3 / r3, /,  &
                11H body elmt., 3(4X,4HREAL,10X,12HIMAGINARY   ))
            nlppp = 1
            210 DO  m = 1,6
                mm = m*2 - 1
                buf(mm  ) = 0.0
                buf(mm+1) = 0.0
                IF (andf(iz(ipuset+m),mask) == 0) CYCLE
                buf(mm  ) = z(irow  )
                buf(mm+1) = z(irow+1)
                irow = irow + 2
            END DO
            WRITE  (out,221) extid,buf
221         FORMAT (1H0,i10,6(1P,e15.6), /11X,6(1P,e15.6))
            nlppp = nlppp + 3
            GO TO iret, (140,190)
        END DO
250     IF (k == nload) CYCLE
        CALL REWIND (pkf)
        CALL skprec (pkf,1)
    END DO
 
    !     CLOSE UP AND RETURN
 
1000 CALL CLOSE (casecc,1)
    CALL CLOSE (pkf,1)
    CALL CLOSE (sila,1)
    CALL CLOSE (spline,1)
    DO  i = 1,96
        head(i) = tsave(i)
    END DO
    CALL page2 (1)
    RETURN
 
    !     ERROR MESSAGES
 
999 CALL mesage (8,0,nam)
    GO TO 1000
END SUBROUTINE adrprt
