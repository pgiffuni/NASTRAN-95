SUBROUTINE xsfa (x)
     
    !     ENTRY SIZE NUMBERS,  1=FIAT, 2=SOS, 3=MD, 4=DPD
 
    !     REVISED  8/89,  SEE XSFABD
 
    IMPLICIT INTEGER (a-z)
    EXTERNAL        lshift,rshift,andf,orf,complf
    DIMENSION       iprt(23),nsfa(3),ddbn(1),dfnu(1),fcum(1),  &
                    fcus( 1),fdbn(1),fequ(1),FILE(1),fknd(1),  &
                    fmat( 1),fntu(1),fpun(1),fon (1),ford(1),  &
                    minp( 1),mlsn(1),mout(1),mscr(1),sal (1),  &
                    sdbn( 1),sntu(1),sord(1),pfil(2,3)
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm

    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /blank / ibnk(1)
    COMMON /machin/ mch
    COMMON /xfiat / fiat(7)
    COMMON /xfist / fist
    COMMON /xdpl  / dpd(6)
    COMMON /zzzzzz/ buf1
    COMMON /system/ ibufsz,outtap,dum(17),pltflg,dum1,thislk,dum2,  &
                    icfiat,dmm(14),nbpc,nbpw,ncpw
    COMMON /ixsfa / lmt3,bff,pad,idefr1,idefr2
    COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(5)

    EQUIVALENCE             (dpd  (1),dnaf    ),(dpd  (2),dmxlg   ),  &
        (dpd  (3),dculg   ),(dpd  (4),ddbn (1)),(dpd  (6),dfnu (1)),  &
        (fiat (1),funlg   ),(fiat (2),fmxlg   ),(fiat (3),fculg   ),  &
        (fiat (4),fequ (1)),(fiat (4),FILE (1)),(fiat (4),ford (1)),  &
        (fiat (5),fdbn (1)),(fiat (7),fmat (1)),(md   (1),mlgn    ),  &
        (md   (2),mlsn (1)),(md   (3),minp (1)),(md   (4),mout (1)),  &
        (md   (5),mscr (1)),(sos  (1),slgn    ),(sos  (2),sdbn (1)),  &
        (sos  (4),sal  (1)),(sos  (4),sntu (1)),(sos  (4),sord (1)),  &
        (xf1at(1),fntu (1)),(xf1at(1),fon  (1)),(xf1at(2),fpun (1)),  &
        (xf1at(3),fcum (1)),(xf1at(4),fcus (1)),(xf1at(5),fknd (1))
    EQUIVALENCE             (comm (1),almsk   ),(comm (2),apndmk  ),  &
        (comm (3),cursno  ),(comm (4),entn1   ),(comm (5),entn2   ),  &
        (comm (6),entn3   ),(comm (7),entn4   ),(comm (8),flag    ),  &
        (comm (9),fnx     ),(comm(10),lmsk    ),(comm(11),lxmsk   ),  &
        (comm(12),macsft  ),(comm(13),rmsk    ),(comm(14),rxmsk   ),  &
        (comm(15),s       ),(comm(16),scornt  ),(comm(17),tapmsk  ),  &
        (comm(18),thcrmk  ),(comm(19),zap     )
    
    DATA  oscar1, oscar2/ 4HXOSC, 4HAR  /, pool  / 4HPOOL  /
    DATA  nsfa  / 4HXSFA, 4H    , 4H    /, ns14  / 4HNS14  /
    DATA  ibegn,  iend  / 4HBEGN, 4HEND /
    DATA  plus  / 1H+   /
    DATA  pfil  / 4HPLTP, 4HAR  , 4HGPSE, 4HTS  , 4HELSE, 4HTS   /

    CALL xsfadd
    nsfa(3) = ibegn
    CALL conmsg (nsfa,3,0)

    !     ALMSK  = O 377777777777     Z 7FFFFFFF
    almsk  = rshift(complf(0),1)

    !     THCRMK = O 777777000000     Z FFFFFF00
    thcrmk = lshift(almsk,nbpw-(3*nbpc))

    !     S      = O 400000000000     Z 80000000
    s      = lshift(1,nbpw-1)

    !     MACSFT = SHIFT COUNT TO PLACE INTEGER IN 4TH FROM LEFT CHARACTER
    macsft = (ncpw-4)*nbpc

    entn1  = icfiat
    cursno = x

    !     GET OSCAR FILE POSITION AND SAVE IN FNOS
    !     ALSO SAVE RECORD POSITION IN RNOS

    CALL xpolck (oscar1,oscar2,fnos,nx)
    IF (fnos == 0) GO TO 920
    fnx  = fnos
    rnos = cursno
    CALL xsosgn
    IF (mlgn == 0) GO TO 930
    CALL xclean

    !     INITIALIZE PRIOR TO FIRST MODULE ALLOCATION

    ASSIGN 670 TO itest

    lmt1  = mlgn *entn3
    lmt8  = funlg*entn1
    lmt8p1= lmt8 + 1
    DO  i = 1,lmt8,entn1
        IF (andf(tapmsk,FILE(i)) /= 0) GO TO 120
    END DO
    tapmsk = 0

    !     LOOP THRU ALL MODULES IN SOS

120 i = 1
125 totio = minp(i)+ mout(i)
    totf  = totio  + mscr(i)
    alcnt = 0
    lmt2  = lmt3 + 1
    lmt4  = lmt3 + minp(i)*entn2
    lmt5  = lmt4 + mout(i)*entn2
    lmt3  = lmt3 + totf   *entn2
    lmt9  = fculg*entn1
    nfculg= lmt9 + 1
    itiord= lshift(mlsn(i),16)
    DO  j = 1,lmt9,entn1
        fcum(j) = 0
    END DO

    !     SEQUENCE THRU SOS (ONE MODULE) LOOK FOR NAME MATCH + LTU COMPARE

150 flag = 0
    DO  k = lmt2,lmt3,entn2
        IF (sal(k) < 0)  CYCLE
        itpflg = andf(tapmsk,sntu(k))
  
        !     SEQUENCE THRU FIAT (NAME MATCH)
  
        DO  f1 = 1,lmt9,entn1
            IF (sdbn(k) /= fdbn(f1) .OR. sdbn(k+1) /= fdbn(f1+1)) CYCLE
            IF (fpun(f1) < 0) GO TO 680
            fntu(f1) = orf(andf(s,fon(f1)),sntu(k))
            fcum(f1) = -1
            fcus(f1) = -1
            IF (fknd(f1) == 0) fknd(f1) = 1
            GO TO 230
        END DO
        IF (mlsn(i) < 0) CYCLE
        IF (k    <= lmt4) CYCLE
        IF (andf(apndmk,sord(k)) == apndmk) CYCLE
  
        !     SEQUENCE THRU FIAT (LTU COMPARE)
  
        loop220:  DO  f1 = 1,lmt9,entn1
            IF (itiord <= andf(lmsk,ford(f1))) CYCLE loop220
            IF (fon (f1) < 0) CYCLE loop220
            IF (fcum(f1) < 0) CYCLE loop220
            IF (fdbn(f1) == 0) CYCLE loop220
            IF (andf(rmsk,FILE(f1)) == rmsk) CYCLE loop220
            IF (andf(lmsk,ford(f1)) == lmsk) CYCLE loop220
            IF (itpflg /= 0 .AND. andf(tapmsk,FILE(f1)) == 0) CYCLE loop220
            IF (fequ(f1) >= 0) GO TO 210
            fil = andf(rmsk,FILE(f1))
            DO  l = 1,lmt9,entn1
                IF (fequ(l) >= 0) CYCLE
                IF (f1      == l) CYCLE
                IF (fil    /= andf(rmsk,FILE(l))) CYCLE
                IF (itiord <= andf(lmsk,ford(l))) CYCLE loop220
                IF (fon(l)  < 0) CYCLE loop220
                IF (fcum(l) < 0) CYCLE loop220
            END DO
210         IF (fculg+pad >= fmxlg) GO TO 680
            fon(f1) = orf(s,fon(f1))
            fdbn(nfculg  ) = sdbn(k  )
            fdbn(nfculg+1) = sdbn(k+1)
            ford(nfculg  ) = orf(andf(lxmsk,sord(k)),andf(rxmsk,FILE(f1)))
            fntu(nfculg  ) = sntu(k)
            fcum(nfculg  ) = -1
            fcus(nfculg  ) = -1
            fknd(nfculg  ) = 2
            nfculg = nfculg+ entn1
            fculg  = fculg + 1
            GO TO 230
        END DO loop220
        CYCLE
230     sal(k) = orf(s,sal(k))
        alcnt  = alcnt + 1
    END DO
    IF (alcnt == totf) GO TO 600

    !     SEQUENCE THRU SOS (ONE MODULE) LOOK FOR BLANK FILES + GREATER NTU

    DO  k = lmt2,lmt3,entn2
        IF (sal(k) < 0) CYCLE
        IF (flag /= 0 .AND. k > lmt4 .AND. k <= lmt5) GO TO 150
        iapflg = 0
        iunpfg = 0
        IF (andf(apndmk,sord(k)) == apndmk) iapflg = -1
        itpflg = andf(tapmsk,sntu(k))
  
        !     SEQUENCE THRU FIAT-UNIQUE (BLANK FILES)
  
        IF (bff < 0) GO TO 390
        DO  f1 = 1,lmt8,entn1
            IF (fdbn(f1) /= 0) CYCLE
            IF (itpflg /= 0 .AND. andf(tapmsk,FILE(f1)) == 0) CYCLE
            IF (k > lmt4 .AND. iapflg == 0) GO TO 310
            CALL xpolck (sdbn(k),sdbn(k+1),fn,nx)
            IF (iapflg /= 0 .AND. fn == 0) GO TO 310
            IF (fn     /= 0) GO TO 300
            IF (pltflg /= 0) GO TO 280
            DO  ip = 1,3
                IF (sdbn(k) == pfil(1,ip) .AND. sdbn(k+1) == pfil(2,ip)) GO TO 300
            END DO
280         IF (thislk /= ns14) CALL mesage (22,0,sdbn(k))
            GO TO 320
300         fpun(f1) = fn
            iunpfg   = f1
310         fdbn(f1  ) = sdbn(k  )
            fdbn(f1+1) = sdbn(k+1)
            ford(f1) = orf(andf(lxmsk,sord(k)),FILE(f1))
            fntu(f1) = sntu(k)
            fcum(f1) =-1
            fcus(f1) =-1
            fknd(f1) = 3
320         sal(k) = orf(s,sal(k))
            alcnt  = alcnt + 1
            GO TO 540
        END DO
        IF (itpflg == 0) bff = -1
  
        !     SEQUENCE THRU FIAT (GREATEST NTU) FOR POOLING
  
390     IF (mlsn(i) < 0) GO TO 680
  
        !     BEFORE PERMITTING POOLING CHECK IF AT LEAST ONE MODULE IS ALLOCATE
  
        IF (i /= 1) GO TO 620
400     mxntu  = cursno
        mxntui = 0
        DO  f1 = 1,lmt8,entn1
            IF (fcus(f1) < 0) CYCLE
            IF (idefr2   < 0) GO TO 420
            IF (fmat(f1) /= 0 .OR. fmat(f1+1) /= 0 .OR. fmat(f1+2) /= 0) GO TO 410
            IF (entn1 == 11 .AND. (fmat(f1+5) /= 0 .OR. fmat(f1+6) /= 0 .OR.  &
                fmat(f1+7) /= 0)) GO TO 410
            GO TO 420
410         idefr1 = -1
            CYCLE
420         IF (fknd(f1) < 0) CYCLE
            IF (fpun(f1) /= 0) CYCLE
            IF (itpflg /= 0 .AND. andf(tapmsk,FILE(f1)) == 0) CYCLE
            trial = andf(fntu(f1),rmsk)
            IF (trial <= mxntu) CYCLE
            mxntu  = trial
            mxntui = f1
        END DO
        IF (mxntui /= 0) GO TO 463
  
        !     FILE NOT FOUND - HAS A PASS BEEN DEFERRED
  
        IF (idefr1 == 0) GO TO 680
  
        !     PASS HAS BEEN DEFERRED - TRY IT NOW
  
        idefr1 = 0
        idefr2 =-1
        DO  ix = 1,lmt8,entn1
            fknd(ix) = IABS(fknd(ix))
        END DO
        GO TO 400
  
        !     A GREATER NTU FILE EXISTS
  
463     n = 1
  
        !     SEARCH FOR EQUIV OR STACKED MATCH
  
        fil = andf(rmsk,FILE(mxntui))
        DO  j = lmt8p1,lmt9,entn1
            IF (fil /= andf(rmsk,FILE(j))) CYCLE
    
            !     A MATCH IS FOUND, IS MATCHED FILE USED IN CURRENT SEG
    
            IF (fcus(j) < 0) GO TO 490
    
            !     IF MATCHED FILE HAS NTU LESS - TEST AND SET DEFER FLAG
    
            IF (idefr2 < 0) GO TO 465
            IF (fmat(j) /= 0 .OR.  fmat(j+1) /= 0 .OR. fmat(j+2) /= 0) GO TO 464
            IF (entn1 == 11 .AND. (fmat(j+5) /= 0 .OR. fmat(j+6) /= 0 .OR.  &
                fmat(j+7) /= 0)) GO TO 464
            IF (andf(rmsk,fntu(1)) >= andf(rmsk,fntu(mxntui))) GO TO 465
464         idefr1 = -1
            GO TO 490
    
            !     MATCHED FILE IS O.K. - IS IT EQUIV OR STACKED
    
465         IF (fequ(j) >= 0) GO TO 467
            fknd(j) = 7
            n = n + 1
            CYCLE
    
            !     STACKED - WIPE OUT MATCH (IF EMPTY)
    
467         IF (fmat(j) /= 0 .OR.  fmat(j+1) /= 0 .OR. fmat(j+2) /= 0) GO TO 490
            IF (entn1 == 11 .AND. (fmat(j+5) /= 0 .OR. fmat(j+6) /= 0 .OR.  &
                fmat(j+7) /= 0)) GO TO 490
            FILE(j  ) = 0
            fdbn(j  ) = 0
            fdbn(j+1) = 0
        END DO
        fpun(mxntui) = orf(s,n)
        IF (k > lmt4 .AND. iapflg == 0) GO TO 520
        CALL xpolck (sdbn(k),sdbn(k+1),fn,nx)
        IF (iapflg /= 0 .AND. fn == 0) GO TO 520
        IF (fn /= 0) GO TO 500
        IF (thislk /= 14) CALL mesage (22,0,sdbn(k))
        GO TO 530
490     IF (fknd(mxntui) == 0) fknd(mxntui) = 9
        fknd(mxntui) = -IABS(fknd(mxntui))
        GO TO 400
500     fpun(nfculg) = fn
        iunpfg = nfculg
520     IF (fculg+pad >= fmxlg) GO TO 680
        fon(mxntui ) = orf(s,fon(mxntui))
        ford(nfculg) = orf(andf(rxmsk,FILE(mxntui)),andf(lxmsk,sord(k)))
        fknd(nfculg) = orf(fknd(nfculg),5)
        fdbn(nfculg  ) = sdbn(k  )
        fdbn(nfculg+1) = sdbn(k+1)
        fntu(nfculg) = sntu(k)
        fcum(nfculg) = -1
        fcus(nfculg) = -1
        nfculg = nfculg+ entn1
        fculg  = fculg + 1
530     sal(k) = orf(s,sal(k))
        alcnt  = alcnt+ 1
540     IF (iunpfg   == 0) CYCLE
        IF (dfnu(nx) >= 0) CYCLE
        CALL xpleqk (nx,iunpfg)
        lmt9  = fculg*entn1
        nfculg= lmt9 + 1
    END DO

    !     MODULE ALLOCATION COMPLETE

600 cursno = andf(rmsk,mlsn(i)) + 1

    !     END OF I MODULE PSEUDO LOOP

    i = i + entn3
    IF (i <= lmt1) GO TO 125

620 CALL xpunp
    CALL xdph

    !     REPOSITION OSCAR FOR SEM

    CALL xpolck (oscar1,oscar2,fnos,nx)
    IF (fnos == 0) GO TO 920
630 CALL OPEN (*940,pool,buf1,0)
    IF (fnos /= 1) CALL skpfil (pool,fnos-1)
    DO  j = 1,rnos
        CALL fwdrec (*950,pool)
    END DO
    CALL CLOSE (pool,2)

655 CONTINUE

    !     DUMP FIAT IF SENSE SWITCH 2 IS ON

    CALL sswtch (2,ix)
    IF (ix /= 1) GO TO itest, (670,715)
    CALL page1
    CALL page2 (-4)
    WRITE  (outtap,660) fiat(1),fiat(2),fiat(3),x,cursno
660 FORMAT (15H0FIAT after sfa,3I4,12H  oscar str ,i4,6H, stp ,i4, //,  &
        ' EQ AP  LTU  TP  UNIT  NTU  OF SG KN TR DATA-BLK      *',  &
        6X,'*   TRAILER   *      *      *  PRI BLKS   SEC FLS/BLKS',  &
        3X,'TER FLS/BLKS')
    ii = fiat(3)*entn1
    DO  ix = 1,ii,entn1
        iprt( 1) = rshift(fequ(ix),nbpw-1)
        iprt( 2) = rshift(andf(apndmk,ford(ix)),30)
        iprt( 3) = rshift(andf(lmsk  ,ford(ix)),16)
        iprt( 4) = rshift(andf(tapmsk,FILE(ix)),15)
        iprt( 5) = andf(rmsk,FILE(ix))
        iprt( 6) = andf(rmsk,fntu(ix))
        iprt( 7) = rshift(fon(ix),nbpw-1)
        iprt( 8) = fcus(ix)
        iprt( 9) = fknd(ix)
        iprt(10) = rshift(andf(tapmsk,fntu(ix)),15)
        iprt(11) = fdbn(ix  )
        iprt(12) = fdbn(ix+1)
        IF (iprt(11) /= 0) GO TO 661
        iprt(11) = nsfa(2)
        iprt(12) = nsfa(2)
661     IF (entn1 == 11) GO TO 662
        iprt(13) = rshift(fmat(ix),16)
        iprt(14) = andf(rxmsk,fmat(ix))
        iprt(15) = rshift(fmat(ix+1),16)
        iprt(16) = andf(rxmsk,fmat(ix+1))
        iprt(17) = rshift(fmat(ix+2),16)
        iprt(18) = andf(rxmsk,fmat(ix+2))
        GO TO 663
662     iprt(13) = fmat(ix  )
        iprt(14) = fmat(ix+1)
        iprt(15) = fmat(ix+2)
        iprt(16) = fmat(ix+5)
        iprt(17) = fmat(ix+6)
        iprt(18) = fmat(ix+7)
663     iprt(19) = rshift(fmat(ix+3),16)
        itemp    = andf(fmat(ix+3),rxmsk)
        iprt(20) = rshift(itemp,8)
        iprt(21) = rshift(fmat(ix+4),16)
        iprt(22) = itemp - iprt(20)*2**8
        iprt(23) = andf(rxmsk,fmat(ix+4))
        CALL page2 (-1)
        WRITE  (outtap,664) (iprt(iy),iy=1,23)
664     FORMAT (1H ,2(i2,1X),i5,1X,i2,2(1X,i5),4(1X,i2),1X,2A4,6I7,  &
            4X,i5,1X,2(7X,i2,1H/,i5))
    END DO
    CALL xflszd (0,blksiz,0)
    CALL page2 (-2)
    WRITE  (outtap,628) blksiz
628 FORMAT (30X,20H each BLOCK contains,i5,7H words.)
    WRITE  (outtap,666)
666 FORMAT (52H pool FILE contents   EQ    size   FILE   DATA BLOCK)
    ii = dpd(3)*3
    DO  ix = 1,ii,3
        iprt(1) = rshift(dfnu(ix),nbpw-1)
        iprt(2) = rshift(dfnu(ix),16)
        iprt(3) = andf(rxmsk,dfnu(ix))
        iprt(4) = ddbn(ix  )
        iprt(5) = ddbn(ix+1)
        CALL page2 (-1)
        WRITE  (outtap,667) (iprt(iy),iy=1,5)
667     FORMAT (22X,i2,i7,i7,3X,2A4)
    END DO
    CALL dbmdia
    CALL dbmstf

    GO TO itest, (670,715)

670 j = mch
    IF (IABS(ibnk(entn1*5))/1000 /= j .AND. j > 6) comm(4) = j
    x = cursno
    nsfa(3) = iend
    CALL conmsg (nsfa,3,0)
    RETURN

    !     MODULE ALLOCATION INCOMPLETE

680 IF (i      /= 1) GO TO 620
    IF (itpflg == 0) GO TO 700

    !     LOOKING FOR A TAPE + AT LEAST ONE TAPE EXISTS

    noaval = 0
    DO  m = 1,lmt8,entn1
        IF (andf(tapmsk,FILE(m)) == 0) CYCLE
        IF (andf(tapmsk,fntu(m)) == 0) GO TO 710
        noaval = 1
    END DO
    IF (noaval == 0) GO TO 700
    tapmsk = 0
    GO TO 790
700 cursno = 0
    GO TO 630

    !     A TAPE FILE EXIST CONTAINING A D.B. NOT REQUIRING A TAPE  -
    !     FREE THAT TAPE***  CHECK FOR EQUIV AND LTU D.B. ON SAME UNIT

710 n = 1

    ASSIGN 715 TO itest
    GO TO 655
715 CONTINUE
    ASSIGN 670 TO itest

    trial = andf(rmsk,FILE(m))
    lmt   = lmt8 + 1
    DO  j = lmt,lmt9,entn1
        IF (trial /= andf(rmsk,FILE(j))) CYCLE
        inam1 = fdbn(j  )
        inam2 = fdbn(j+1)
        IF (fequ(m) < 0 .AND. fequ(j) < 0) GO TO 720
        fdbn(j) = almsk
        GO TO 725
720     n = n + 1
        725 DO  l = lmt2,lmt3,entn2
            IF (inam1 == sdbn(l) .AND. inam2 == sdbn(l+1)) GO TO 740
        END DO
        CYCLE
  
        !     TURN OFF ALLOC FLAG
  
740     sal(l) = orf(almsk,sal(l))
        alcnt  = alcnt - 1
    END DO
    inam1 = fdbn(m  )
    inam2 = fdbn(m+1)
    DO  l = lmt2,lmt3,entn2
        IF (inam1 == sdbn(l) .AND. inam2 == sdbn(l+1)) GO TO 770
    END DO
    GO TO 780
770 sal(l) = orf(almsk,sal(l))
    alcnt  = alcnt - 1
780 fpun(m)= orf(s,n)
    CALL xpunp
    fdbn(m  ) = sdbn(k  )
    fdbn(m+1) = sdbn(k+1)
    ford(m  ) = orf(andf(lxmsk,sord(k)),andf(rxmsk,FILE(m)))
    fknd(m  ) = 8

    CALL sswtch (2,ix)
    IF (ix /= 1) GO TO 790
    CALL page2 (-2)
    WRITE  (outtap,785)
785 FORMAT (38H0* xsfa repeats TO use freed tape FILE)

790 bff = 0
    GO TO 150

920 WRITE  (outtap,921) sfm
921 FORMAT (a25,' 1001, OSCAR NOT FOUND IN DPL')
    GO TO  1000
930 WRITE  (outtap,931) sfm
931 FORMAT (a25,' 1002, OSCAR CONTAINS NO MODULES')
    GO TO  1000
940 WRITE  (outtap,941) sfm
941 FORMAT (a25,' 1003, POOL COULD NOT BE OPENED')
    GO TO  1000
950 WRITE  (outtap,951) sfm
951 FORMAT (a25,' 1004, ILLEGAL EOF ON POOL')
1000 CALL mesage (-37,0,nsfa)

    RETURN
END SUBROUTINE xsfa
