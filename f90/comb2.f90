SUBROUTINE comb2
     
    !     COMB2 PERFORMS THE TRANSFORMATION AND ADDITION OF STIFFNESS, MASS,
    !     OR LOAD MATRICES FOR THE PHASE 2 SUBSTRUCTURE COMBINE OPERATION

    !     NOVEMBER 1973
 
 
    LOGICAL :: addflg
    INTEGER :: tflag      ,signab     ,signc      ,prec      ,  &
        scr1     ,scr2       ,scr3       ,scr4       ,scr5      ,  &
        rule     ,typin      ,typout     ,acomb      ,amcb(7,7) ,  &
        hmcb(6,7),sof1       ,sof2       ,sof3       ,dry       ,  &
        buf1     ,rdsof      ,rc         ,use        ,horg      ,  &
        pvec     ,rfiles     ,iz(1)      ,rect       ,rsp       ,  &
        name(2)  ,rsofar     ,kmp(5)     ,TYPE       ,kmpitm(5) ,  &
        blank    ,sysbuf     ,cpv(7)     ,rpv(7)     ,scr6      ,  &
        xxxx     ,scr7       ,pora       ,papp
 
    DOUBLE PRECISION :: dbz(1)
    DIMENSION         mcbtrl(7)
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg /   ufm        ,uwm        ,uim        ,sfm
    COMMON /BLANK /   dry        ,TYPE(2)    ,pora(2)    ,namess(2,7),  &
        acomb      ,use(14)    ,rfiles(3)  ,kk        , kn         ,jn
    COMMON /system/   sysbuf     ,nout
    COMMON /names /   rd         ,rdrew      ,wrt        ,wrtrew    ,  &
        rew        ,norew      ,eofnrw     ,rsp       ,  &
        rdp        ,csp        ,cdp        ,square    , rect
    COMMON /packx /   typin      ,typout     ,irow       ,nrow      , incr
    COMMON /parmeg/   mcbp(7)    ,mcbp11(7)  ,mcbp21(7)  ,mcbp12(7) ,  &
        mcbp22(7)  ,mrgz       ,rule
    COMMON /mpy3tl/   mcba2(7)   ,mcbb2(7)   ,mcbc2(7)   ,mcbd2(7)  ,  &
        scr5       ,scr6       ,scr7       ,lkore     ,  &
        icode      ,iprec      ,dummy(13)
    COMMON /mpyadx/   mcba(7)    ,mcbb(7)    ,mcbc(7)    ,mcbd(7)   ,  &
        lcore      ,tflag      ,signab     ,signc     , prec       ,mscr       ,dumm
    COMMON /zzzzzz/   z(1)

    EQUIVALENCE      (dbz(1),z(1),iz(1)),(pvec,kmpitm(3))

    DATA   NAME   /  4HCOMB,4H2   /
    DATA   horg   /  4HHORG       /
    DATA   BLANK  /  4H           /
    DATA   xxxx   /  4HXXXX       /
    DATA   papp   /  4HPAPP       /
    DATA   kmp    /  4HK   , 4HM   , 4HP   , 4HB   , 4HK4   /
    DATA   kmpitm /  4HKMTX, 4HMMTX, 4HPVEC, 4HBMTX, 4HK4MX /
 
    !     INITIALIZE
 
    DO  i = 1,14
        IF (namess(i,1) == xxxx .OR. namess(i,1) == 0) namess(i,1) = BLANK
    END DO
    acomb = 201
    scr1  = 301
    scr2  = 302
    scr3  = 303
    scr4  = 304
    scr5  = 305
    scr6  = 306
    scr7  = 307
    signab= 1
    signc = 1
    prec  = 0
    mscr  = scr5
    icode = 0
    rule  = 0
    rdsof = 1
    nogo  = 0
    rfiles(1) = acomb
    nsize = 0
    mcbp21(1) = 0
    mcbp22(1) = 0
    rsofar = 0
    kn = 1
    jn =-1
    DO  i = 1,5
        IF (TYPE(1) == kmp(i)) GO TO 20
    END DO
    WRITE (nout,6302) sfm,TYPE
    IF (dry < 0) RETURN
 
    dry  =-2
    item = 0
    GO TO 30
20  item = kmpitm(i)
    IF (item == pvec) item = pora(1)
    IF (dry < 0) RETURN
 
30  lcore = korsz(z) - 1
    lkore = lcore
    buf1  = lcore - sysbuf + 1
    sof1  = buf1  - sysbuf
    sof2  = sof1  - sysbuf - 1
    sof3  = sof2  - sysbuf
    IF (sof3 > 0) GO TO 40
    CALL mesage (8,0,NAME)
    dry =-2
    RETURN
 
40  CALL sofopn (z(sof1),z(sof2),z(sof3))
 
    !     GRAB THE MATRIX CONTROL BLOCKS
 
    nmat = 0
    DO  i = 1,7
        IF (namess(1,i) == BLANK) CYCLE
        amcb(1,i) = 100 + i
        CALL rdtrl (amcb(1,i))
        IF (amcb(1,i) > 0) GO TO 135
   
        !     NO GINO FILE.  CHECK SOF
   
        CALL softrl (namess(1,i),item,amcb(1,i))
        rc = amcb(1,i)
        SELECT CASE ( rc )
            CASE (    1)
                GO TO 130
            CASE (    2)
                GO TO 110
            CASE (    3)
                GO TO 115
            CASE (    4)
                GO TO 120
            CASE (    5)
                GO TO 120
        END SELECT
110     nogo = 1
        WRITE (nout,6301) sfm,namess(1,i),namess(2,i),item
        CYCLE
115     IF (TYPE(1) == kmp(3)) CYCLE
120     nogo = 1
        CALL smsg (rc-2,item,namess(1,i))
        CYCLE
   
    !     MATRIX FOUND ON SOF
   
130 CONTINUE
    amcb(1,i) = 0
   
    !     GRAB THE MCB OF THE TRANSFORMATION MATRIX
   
135 CALL softrl (namess(1,i),horg,mcbtrl)
    rc = mcbtrl(1)
    SELECT CASE ( rc )
        CASE (    1)
            GO TO 160
        CASE (    2)
            GO TO 140
        CASE (    3)
            GO TO 150
        CASE (    4)
            GO TO 155
        CASE (    5)
            GO TO 155
    END SELECT
140 nogo = 1
    CALL smsg (1,horg,namess(1,i))
    CYCLE
150 nogo = 1
    WRITE (nout,6303) sfm,namess(1,i),namess(2,i)
    CYCLE
155 nogo = 1
    CALL smsg (rc-2,horg,namess(1,i))
    CYCLE
    160 DO  it = 1,6
        hmcb(it,i) = mcbtrl(it+1)
    END DO
    nmat = nmat + 1
    use(2*nmat-1) = i
    den = FLOAT(amcb(7,i))/10000.
    use(2*nmat) = amcb(2,i)*amcb(3,i)*den
   
    !     CHECK COMPATIBILITY OF DIMENSIONS
   
    IF (nsize == 0) nsize = hmcb(1,i)
    IF (hmcb(1,i) == nsize .AND. hmcb(2,i) == amcb(2,i) .AND.  &
        hmcb(2,i) == amcb(3,i)) CYCLE
    IF (item == pvec .OR. item == papp .AND. hmcb(1,i) == nsize .AND.  &
        hmcb(2,i) == amcb(3,i)) CYCLE
    nogo = 1
    WRITE (nout,6304) sfm,i,namess(1,i),namess(2,i)
END DO
IF (nogo == 0) GO TO 175
 
174 dry =-2
WRITE  (nout,177) amcb,hmcb
177 FORMAT ('0*** COMB2 MATRIX TRAILER DUMP', //7(4X,7I10/), /7(11X,6I10/))
GO TO 9999
 
175 IF (nmat == 0) GO TO 9999
 
!     DETERMINE PRECISION FOR FINAL MATRIX
 
iprc = 1
ityp = 0
DO  i = 1,nmat
    IF (amcb(5,i) == 2 .OR. amcb(5,i) == 4) iprc = 2
    IF (amcb(5,i) >= 3) ityp = 2
END DO
iprec = ityp + iprc
 
IF (item == pvec .OR. item == papp) GO TO 300
!                                               ******
!                                                         *
!     PROCESS STIFFNESS, MASS OR DAMPING MATRICES         *
!                                                         *
!                                               ******
 
!     IF NMAT IS ODD, PUT FIRST RESULT ON ACOMB.  IF EVEN, PUT IT ON
!     SCR4.  FINAL RESULT WILL THEN BE ON ACOMB.
 
CALL sort (0,0,2,2,use,2*nmat)
irf = 1
IF ((nmat/2)*2 == nmat) irf = 2
iform = 6
rfiles(2) = scr4
addflg =.false.
 
DO  kk = 1,nmat
    j  = 2*kk - 1
    jn = jn + 2
    inuse = use(jn)
   
    !     MOVE TRANSFORMATION MATRIX TO SCR2
   
    CALL mtrxi (scr2,namess(1,inuse),horg,z(buf1),rc)
   
    !     IF INPUT MATRIX IS ON SOF, MOVE IT TO SCR1
   
    mcbb2(1) = 100 + inuse
    IF (amcb(1,inuse) > 0) GO TO 180
    mcbb2(1) = scr1
    CALL mtrxi (scr1,namess(1,inuse),item,z(buf1),rc)
   
    !     PERFORM TRIPLE MULTIPLY  H(T)*INPUT*H
   
180 CALL sofcls
    mcba2(1) = scr2
    mcbc2(1) = 0
    IF (addflg) mcbc2(1) = rfiles(3-irf)
    addflg = .true.
    DO  j = 2,7
        mcba2(j) = hmcb(j-1,inuse)
        mcbb2(j) = amcb(j,inuse)
    END DO
    IF (mcbb2(4) <= 2) iform = mcbb2(4)
    CALL makmcb (mcbd2,rfiles(irf),hmcb(1,inuse),iform,iprec)
   
    CALL mpy3dr (z)
   
    CALL wrttrl (mcbd2)
    DO  j = 2,7
        mcbc2(j) = mcbd2(j)
    END DO
    irf = 3 - irf
    CALL sofopn (z(sof1),z(sof2),z(sof3))
END DO
GO TO 9999
!                                          ******
!                                                    *
!     PROCESS LOAD MATRICES                          *
!                                                    *
!****                                           ******
300 mcbc(1) = 0
mcba(1) = scr2
tflag   = 1
mrgz    = lcore
prec    = 0
 
!     SELECT FIRST RESULT FILE SO THAT FINAL RESULT WILL WIND UP ON
!     ACOMB
 
rfiles(2) = scr3
rfiles(3) = scr4
irf = 1
IF (nmat == 2 .OR. nmat == 5) irf = 2
IF (nmat == 3 .OR. nmat == 6) irf = 3
IF (nmat == 1) mcbp(1) = acomb
 
!     CREATE COLUMN PARTITIONING VECTOR FOR ALL MERGES
!     (VECTOR IS ALWAYS NULL)
 
CALL makmcb (cpv,scr6,nsize,rect,rsp)
typin  = rsp
typout = rsp
irow = 1
nrow = 1
incr = 1
CALL gopen (scr6,z(buf1),wrtrew)
CALL pack  (0,scr6,cpv)
CALL CLOSE (scr6,rew)
addflg =.true.
 
DO  kk = 1,nmat
    inuse = use(2*kk-1)
   
    !     COPY TRANSFORMATION MATRIX TO SCR2
   
    CALL mtrxi (scr2,namess(1,inuse),horg,z(buf1),rc)
   
    !     IF LOAD MATRIX IS ON SOF, COPY IT TO SCR1
   
    mcbb(1) = 100 + inuse
    IF (amcb(1,inuse) > 0) GO TO 330
    mcbb(1) = scr1
    CALL mtrxi (scr1,namess(1,inuse),item,z(buf1),rc)
   
    !     MULTIPLY (HT * A) AND STORE RESULT ON RFILES(IRF)
    !     (ACOMB,SCR3, OR SCR4)
   
330 CALL sofcls
    DO  j = 2,7
        mcba(j) = hmcb(j-1,inuse)
        mcbb(j) = amcb(j,inuse)
    END DO
    IF (mcbb(6) == 0 .OR. mcba(3) == mcbb(3)) GO TO 350
    i = kk
    nogo = 1
    WRITE (nout,6304) sfm,i,namess(1,i),namess(2,i)
    GO TO 174
350 CALL makmcb (mcbd,rfiles(irf),hmcb(1,inuse),rect,iprec)
   
    CALL mpyad (z,z,z)
   
    IF (addflg) GO TO 390
   
    !     COMPUTE ROW PARTITIONING VECTOR TO MERGE RESULT OF THIS MULTIPLY
    !     WITH ALL PREVIOUS RESULTS
   
    k = amcb(2,inuse)
    CALL makmcb (rpv,scr5,rsofar+k,rect,rsp)
    IF (k > lcore) GO TO 9008
    DO  j = 1,k
        z(j)  = 1.0E0
    END DO
    typin = rsp
    typout= rsp
    irow  = rsofar + 1
    nrow  = rsofar + k
    incr  = 1
    CALL gopen (scr5,z(buf1),wrtrew)
    CALL pack  (z,scr5,rpv)
    CALL CLOSE (scr5,rew)
   
    !     MERGE MATRICES   STORE RESULT ON NEXT AVAILABLE RFILE
   
    j = MOD(irf,3) + 1
    CALL makmcb (mcbp,rfiles(j),nsize,rect,iprec)
    mcbp(2) = rpv(3)
    j = MOD(j,3) + 1
    mcbp11(1) = rfiles(j)
    mcbp12(1) = rfiles(irf)
    DO  j = 2,7
        mcbp11(j) = mcbp(j)
        mcbp12(j) = mcbd(j)
    END DO
   
    CALL merge (rpv,cpv,z)
   
    irf = MOD(irf,3) + 1
    GO TO 395
    390 DO  j = 2,7
        mcbp(j) = mcbd(j)
    END DO
395 rsofar = rsofar + amcb(2,inuse)
    addflg =.false.
    irf = MOD(irf,3) + 1
    CALL sofopn (z(sof1),z(sof2),z(sof3))
END DO
CALL wrttrl (mcbp)
GO TO 9999
 
!     DIAGNOSTICS
 
6301 FORMAT (a25,' 6301, DATA MISSING IN GO MODE FOR SUBSTRUCTURE ',  &
    2A4,' ITEM ',a4)
6302 FORMAT (a25,' 6302, ',2A4,' IS ILLEGAL MATRIX TYPE FOR MODULE ',  &
    'COMB2')
6303 FORMAT (a25,' 6303, H OR G TRANSFORMATION MATRIX FOR SUBSTRUCTURE'  &
    ,       1X,2A4,' CANNOT BE FOUND ON SOF')
6304 FORMAT (a25,' 6304, MODULE COMB2 INPUT MATRIX NUMBER ',i2,  &
    ' FOR SUBSTRUCTURE ,2A4,28H HAS INCOMPATIBLE DIMENSIONS')
9008 CALL mesage (8,0,NAME)
 
!     NORMAL COMPLETION
 
9999 CALL sofcls

RETURN
END SUBROUTINE comb2
