SUBROUTINE algap (ifname,ifnm)
     
    !     THIS ROUTINE IS A MODIFIED VERSION OF SUBROUTINE TABPCH. IT WILL
    !     ONLY PUNCH ONE TABLE INTO DTI CHARDS.
 
    !     CONTINUATION CARD CHARACTERS ARE - AL.
 
    !     SINGLE FIELD CARDS WILL BE MADE UNLESS REAL NUMBERS ARE TO BE MADE
    !     ALL REAL NUMBERS ARE ASSUMED TO BE SINGLE PRECISION.
 
    !  $MIXED_FORMATS
 
 
    INTEGER, INTENT(IN OUT)                  :: ifname
    INTEGER, INTENT(IN)                      :: ifnm
    INTEGER :: sysbuf   ,iz(10)   ,NAME(2)  ,INT(2)   ,ireal(2) ,  &
               mcb(7)   ,FILE     ,tabnm(2) ,dti(2)   ,dtis(2)  ,  &
               idata(20),endrec(2),out      ,iform(20),BLANK    ,  &
               ibcd(2)  ,intd(2)  ,ibcdd(2) ,pform(30),ll(4)    ,  &
               form(30,2)         ,forms(30,2)

    REAL :: rdata(20)

    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg / ufm      ,uwm      ,uim
    COMMON /system/ ksystm(100)
    COMMON /zzzzzz/ z(1)

    EQUIVALENCE     (ksystm( 1),sysbuf),(ksystm(2),out ),  &
                    (ksystm(91),lpunch),(iz(1)    ,z(1)),  &
                    (idata(1),rdata(1))

    DATA    BLANK / 1H             /
    DATA    dti   / 4HDTI , 1H     /
    DATA    dtis  / 4HDTI*, 1H     /
    DATA    endrec/ 4HENDR, 4HEC   /
    DATA    forms / 4H(2A4, 26*4H    ,4H,1H+ ,4HA2,i,4H5)  ,  &
                    4H(a1,, 4HA2,i   ,4H5    ,24*4H        ,  &
                    4H,1H+  , 4HA2,i, 4H5)   /
    DATA    ibcd  / 4H,2A4, 1H     /
    DATA    ibcdd / 4H,2A4, 4H,8X  /
    DATA    INT   / 4H,i8 , 1H     /
    DATA    intd  / 4H,i16, 1H     /
    DATA    iplus / 1H+            /
    DATA    ireal / 4H,e16, 4H.9   /
    DATA    istar / 1H*            /
    DATA    NAME  / 4HALGA, 4HP    /
    DATA    n1    / 2HAL           /
    DATA    ll    / 3, 1, 3, 2     /

 
    nz   = korsz(z)
    ibuf = nz - sysbuf + 1
    nz   = ibuf - 1
    IF (nz <= 10) CALL mesage (-8,0,NAME)
    nread = nz/2  - 2
    nlist = nread + 3
    DO  j = 1,2
        DO  i = 1,30
            FORM(i,j) = forms(i,j)
        END DO
    END DO
 
    !     FOR EACH  TABLE DEFINED
 
    mcb(1) = ifnm
    CALL rdtrl(mcb)
    IF (mcb(1) <= 0) GO TO 310
 
    !     TABLE EXISTS SET IT UP
 
    FILE = ifnm
    CALL OPEN (*320,FILE,iz(ibuf),0)
    CALL READ (*340,*350,FILE,iz(1),-2,0,ilen)
    CALL fname (ifname,tabnm)
    irecno = 0
    ichr   = n1
    iz(3)  = 0
 
    !     SET UP FIRST RECORD
 
    iz(1) = tabnm(1)
    iz(2) = tabnm(2)
    iz(4) = mcb(2)
    iz(5) = mcb(3)
    iz(6) = mcb(4)
    iz(7) = mcb(5)
    iz(8) = mcb(6)
    iz(9) = mcb(7)
    CALL READ (*290,*10,FILE,iz(10),nread,0,ilen)
    CALL mesage (-8,0,NAME)
10  ilen = ilen + 11
11  iz(ilen-1) = endrec(1)
    iz(ilen  ) = endrec(2)
    GO TO 40
 
    !     BRING IN NEXT RECORD
 
20  CALL READ (*290,*30,FILE,iz(4),nread,0,ilen)
    CALL mesage (-8,0,NAME)
30  iz(3) = iz(3) + 1
    IF (ilen == 0) GO TO 20
    ilen = ilen + 5
    GO TO 11
 
    !     BUILD FORMAT VECTOR  1= INTEGER, 2 =BCD, 3=REAL
 
    40 DO  k = 1,ilen
        m = nlist + k - 1
        j = numtyp(iz(k))
        iz(m) = ll(j+1)
    END DO
 
    !     MOVE DATA/FORMAT TO DATA AREA 8 FIELDS AT A TIME--SET D.F. FLAG
 
    id = 1
    IF = nlist
    ifrs = 1
 
    !     HERE FOR EIGHT MORE WORDS
 
60  idf = 0
    idt = 1
    ift = 1
    nf  = 1
 
    !     HERE  FOR EACH FIELD
 
70  idata(idt) = iz(id)
    iform(ift) = iz(IF)
    IF (iform(ift) == 3) idf = 1
    IF (iform(ift) /= 2) GO TO 80
 
    !     BCD IS TWO WORDS
 
    idata(idt+1) = iz(id+1)
 
    !     MAY BE FALSE BCD, CHECK FORMAT OF SECOND WORD ALSO
    !     ( SOME REAL NUMBER BIT PATTERNS LOOK LIKE BCD ).
 
    IF (iz(IF+1) == 2) GO TO 100
 
    !     SECOND WORD IS NOT BCD, ASSUME FIRST WORD IS REAL.
 
    idf = 1
    iform(ift) = 3
    GO TO 80
100 idt = idt + 2
    ift = ift + 1
    id  = id  + 2
    IF  = IF  + 2
    GO TO 90
 
    !     REAL OR INTEGER
 
80  idt = idt + 1
    ift = ift + 1
    id  = id  + 1
    IF  = IF  + 1
 
    !     BUMP FIELD COUNTER
 
90  nf = nf + 1
    IF (nf > 8) GO TO 110
    IF (id < ilen) GO TO 70
 
    !     FILL  WITH BLANKS
 
    idata(idt  ) = BLANK
    idata(idt+1) = BLANK
    iform(ift  ) = 2
    GO TO 100
 
    !     PUNCH OUT 8 FIELDS OF DATA
 
110 idt = 0
    IF (idf /= 0) GO TO 200
 
    !     SINGLE FIELD CARD
 
    nf = 1
120 m = 2*nf + 2
    IF (iform(nf)-2 < 0) THEN
        GO TO   130
    ELSE IF (iform(nf)-2 == 0) THEN
        GO TO   150
    ELSE
        GO TO   160
    END IF
 
    !     INTEGER
 
130 FORM(m  ,ifrs) = INT(1)
    FORM(m+1,ifrs) = INT(2)
 
    !     GET NEXT ITEM
 
    idt = idt + 1
140 nf  = nf  + 1
    IF (nf <= 8) GO TO 120
    GO TO 170
 
    !     BCD
 
150 FORM(m  ,ifrs) = ibcd(1)
    FORM(m+1,ifrs) = ibcd(2)
    idt = idt + 2
    GO TO 140
 
    !     REAL NOT LEGAL
 
160 CALL mesage (-61,0,NAME)
    RETURN
 
    !     PUNCH OUT SINGLE CARD
 
170 IF (ifrs /= 1) GO TO 190
    DO  j = 1,30
        pform(j) = FORM(j,1)
    END DO
    WRITE (lpunch,pform,ERR=173) dti,(rdata(m),m=1,idt),ichr,irecno
173 irecno = irecno + 1
    ifrs = 2
    DO  j = 1,30
        FORM(j,1) = forms(j,1)
    END DO
180 IF (id >= ilen) GO TO 20
    GO TO 60
 
    !     CONTINUATION CARD
 
190 ircnm1 = irecno - 1
    DO  j = 1,30
        pform(j) = FORM(j,2)
    END DO
    WRITE (lpunch,pform,ERR=193) iplus,ichr,ircnm1,(rdata(m),m=1,idt),ichr,irecno
193 irecno = irecno + 1
    DO  j = 1,30
        FORM(j,2) = forms(j,2)
    END DO
    GO TO 180
 
    !     DOUBLE FIELD CARDS
 
200 nf = 1
    is = 1
    it = 4
    idt= 0
    m  = 2
210 m  = m + 2
    IF (iform(nf)-2 < 0) THEN
        GO TO   211
    ELSE IF (iform(nf)-2 == 0) THEN
        GO TO   240
    ELSE
        GO TO   250
    END IF
 
    !     INTEGER
 
211 FORM(m  ,ifrs) = intd(1)
    FORM(m+1,ifrs) = intd(2)
220 idt = idt + 1
230 nf  = nf  + 1
    IF (m <= 8) GO TO 210
    GO TO 260
 
    !     BCD
 
240 FORM(m  ,ifrs) = ibcdd(1)
    FORM(m+1,ifrs) = ibcdd(2)
    idt = idt + 2
    GO TO 230
 
    !     REAL
 
250 FORM(m  ,ifrs) = ireal(1)
    FORM(m+1,ifrs) = ireal(2)
    GO TO 220
 
    !     PUNCH OUT DOUBLE FIELD CARD
 
260 IF (ifrs /= 1) GO TO 280
    DO  j = 1,30
        pform(j) = FORM(j,1)
    END DO
    WRITE (lpunch,pform,ERR=263) dtis,(rdata(m),m=is,idt),ichr,irecno
263 irecno = irecno + 1
    DO  j = 1,30
        FORM(j,1) = forms(j,1)
    END DO
    ifrs = 2
270 it = 8
    m  = 2
    is = idt + 1
    GO TO 210
 
    !     CONTINUATION CARD
 
280 ircnm1 = irecno - 1
    DO  j = 1,30
        pform(j) = FORM(j,2)
    END DO
    WRITE (lpunch,pform,ERR=283)  &
        istar,ichr,ircnm1,(rdata(m),m=is,idt),ichr,irecno
283 irecno = irecno + 1
    DO  j = 1,30
        FORM(j,2) = forms(j,2)
    END DO
    IF (it == 4) GO TO 270
    GO TO 180
 
    !     CLOSE OFF FILES
 
290 CALL CLOSE (FILE,1)
    WRITE  (out,300) uim,tabnm,irecno
300 FORMAT (a29,' 4015.', /5X,'TABLE NAMED ',2A4,' PUNCHED ONTO',i9,  &
        ' CARDS.')
310 CONTINUE
    WRITE  (lpunch,311)
311 FORMAT (1H , /,1H , /,1H )
    RETURN
 
    !     ERROR MESAGES
 
320 ip1 = -1
330 CALL mesage (ip1,FILE,NAME)
    CALL mesage (-61,0,NAME)
340 ip1 =-2
    GO TO 330
350 ip1 =-3
    GO TO 330

END SUBROUTINE algap
