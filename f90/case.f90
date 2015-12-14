SUBROUTINE case
     
    !     CASE READS THE CASE CONTROL DATA BLOCK AND WRITES A NEW
    !     DATA BLOCK WHICH CONTAINS ONLY THOSE RECORDS WHICH DESCRIBE THE
    !     CURRENT CASE IN THE LOOP. ADDITIONALLY, THE LOOP CONTROL PARAMETER
    !     IS SET.

 
    INTEGER :: app    ,count  ,sysbuf,casecc,casexx,FILE  ,z    ,  &
               buf1   ,buf2   ,rfmts ,branch,buf   ,error(2)
    INTEGER :: buf3   ,psdl
    DIMENSION       nam(2) ,buf(20),mcb(7),rfmts(40)
    COMMON /BLANK / app(2) ,count  ,loop
    COMMON /system/ sysbuf
    COMMON /names / rd     ,rdrew  ,wrt   ,wrtrew,clsrew
    COMMON /zzzzzz/ z(1)
 
    !     DATA DESCRIBING DATA BLOCK FILE NAMES AND POSITION
    !     OF PARAMETERS IN THE CASE CONTROL RECORD.
 
    DATA   casecc / 101/ ,casexx /201/ ,ik2pp  /139/ ,im2pp /141/ ,  &
           ib2pp  / 143/ ,itfl   / 15/ ,psdl   /102/ ,irand /163/
    DATA   error  / 4HPSDL,4HCASE/
    DATA   ifreq  / 14/  ,imeth  /  5/
 
    !     DATA DEFINING RIGID FORMATS.
 
    DATA   nrigds / 10   /, rfmts  /  &
        4HSTAT,4HICS , 4HREIG,4HEN  , 4HDS0 ,4H    ,  &
        4HDS1 ,4H    , 4HFREQ,4H    , 4HTRAN,4HSNT ,  &
        4HBKL0,4H    , 4HBKL1,4H    , 4HCEIG,4HEN  , 4HPLA ,4H    , 20*0  /
 
    !     MISC DATA
 
    DATA   nam    / 4HCASE,4H    /, mcb    / 7*0  /
 
    !     PERFORM BUFFER ALLOCATION.
 
    buf1  = korsz(z) - sysbuf + 1
    buf3  = buf1 - sysbuf
    buf2  = buf3 - sysbuf
    iry   = 0
    m8    = -8
    IF (count <= 0) count = 1
    loop  = 1
    iocnt = count
 
    !     SET PARAMETER FOR APPROACH.
 
    n = 2*nrigds - 1
    DO  i = 1,n,2
        IF (rfmts(i) == app(1)) GO TO 30
    END DO
    CALL mesage (30,75,app)
    i = 19
30  branch = (i+1)/2
 
    !     OPEN CASECC. SKIP RECORDS ALREADY PROCESSED. OPEN CASEXX.
    !     WRITE HEADER RECORD. THEN BRANCH ON APPROACH.
 
    FILE = casecc
    CALL OPEN (*130,casecc,z(buf1),rdrew)
    DO  i = 1,count
        CALL fwdrec (*140,casecc)
    END DO
    FILE = casexx
    CALL OPEN  (*130,casexx,z(buf2),wrtrew)
    CALL fname (casexx,buf)
    CALL WRITE (casexx,buf,2,1)
    SELECT CASE ( branch )
        CASE (    1)
            GO TO 120
        CASE (    2)
            GO TO 50
        CASE (    3)
            GO TO 120
        CASE (    4)
            GO TO 120
        CASE (    5)
            GO TO 50
        CASE (    6)
            GO TO 100
        CASE (    7)
            GO TO 120
        CASE (    8)
            GO TO 120
        CASE (    9)
            GO TO 50
        CASE (   10)
            GO TO 120
    END SELECT
 
    !     COMPLEX EIGENVALUES OR FREQUENCY RESPONSE.
 
50  CALL READ (*140,*60,casecc,z,buf2,1,ncc)
    CALL mesage (m8,0,nam)
60  buf(1) = z(ik2pp  )
    buf(2) = z(ik2pp+1)
    buf(3) = z(im2pp  )
    buf(4) = z(im2pp+1)
    buf(5) = z(ib2pp  )
    buf(6) = z(ib2pp+1)
    buf(7) = z(itfl)
    irset  = z(irand)
    ifrqst = z(ifreq)
    imrqst = z(imeth)
    IF (branch == 5 .AND. irset /= 0) iry = 1
    IF (iry == 0) GO TO 70
 
    !     BUILD LIST OF UNIQUE LOAD ID-S
 
    FILE = psdl
    CALL OPEN (*68,psdl,z(buf3),rdrew)
    CALL fwdrec (*90,psdl)
    ils  = buf2
    ilf  = buf2 - 1
61  CALL READ (*90,*66,psdl,z(ncc+1),6,0,j)
    IF (z(ncc+1) /= irset) GO TO 61
    j = 1
    iload = z(ncc+2)
    IF (ils == ilf+1) GO TO 63
    65 DO  i = ils,ilf
        IF (z(i) == iload) GO TO 64
    END DO
 
    !     NEW LOAD ID
 
63  ils = ils - 1
    z(ils) = iload
64  IF (j == 0) GO TO 61
    j = 0
    iload = z(ncc+3)
    GO TO 65
 
    !     END OF PSDL RECORD
 
66  CALL CLOSE (psdl,clsrew)
    IF (ils == ilf+1) CALL mesage (-31,irset,error(1))
    buf2 = ils - 1
    GO TO 70
 
    !     NO PSDL IS EQUIVALENT TO NO RANDOM
 
68  iry = 0
70  CALL WRITE (casexx,z,ncc,1)
    count = count + 1
    IF (iry == 0) GO TO 71
 
    !     CHECK  SUBCASE ID-S
 
    DO  i = ils,ilf
        IF (z(1) == z(i)) GO TO 74
    END DO
    GO TO 71
 
    !     MARK USED
 
74  z(i) = -z(i)
71 CONTINUE
   CALL READ (*90,*80,casecc,z,buf2,1,ncc)
   CALL mesage (m8,0,nam)
80 IF (z(ik2pp) /= buf(1) .OR. z(ik2pp+1) /= buf(2) .OR.  &
       z(im2pp) /= buf(3) .OR. z(im2pp+1) /= buf(4) .OR.  &
       z(ib2pp) /= buf(5) .OR. z(ib2pp+1) /= buf(6)) GO TO 120
   IF (z(itfl) /= buf(7)) GO TO 120
   IF (z(imeth) /= 0 .AND. z(imeth) /= imrqst) GO TO 120
 
   !     TEST FOR CHANGED FREQUENCY SET
 
   IF (z(ifreq) /= ifrqst .AND. branch == 5) GO TO 120
   GO TO 70
90 count = -1
   GO TO 120
 
   !     TRANSIENT RESPONSE.
 
100 CALL READ (*140,*110,casecc,z,buf2,1,ncc)
   CALL mesage (m8,0,nam)
110 CALL WRITE (casexx,z,ncc,1)
   count = count + 1
   CALL READ (*90,*120,casecc,z,buf2,1,ncc)
   GO TO 120
 
   !     CLOSE FILES. WRITE TRAILER. RETURN.
 
120 CALL CLOSE (casecc,clsrew)
   CALL CLOSE (casexx,clsrew)
   mcb(1) = casexx
   mcb(2) = count
   CALL wrttrl (mcb)
   IF (count <= 1 .AND. iocnt == 1) loop = -1
 
   !     CHECK ALL PSDL ACCOUNTED FOR
 
   IF (iry == 0) GO TO 125
   nogo = 0
   DO   i = ils,ilf
       IF (z(i) < 0) CYCLE
       nogo = -1
       CALL mesage (33,z(i),nam)
   END DO
   IF (nogo < 0) CALL mesage (-7,0,nam)

125 RETURN
 
   !     FATAL FILE ERRORS.
 
130 n = -1
   GO TO 150
140 n = -2
   FILE = casecc
150 CALL mesage (n,FILE,nam)
   GO TO 150
 
END SUBROUTINE case
