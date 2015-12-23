SUBROUTINE ifs2p (*,*,*)
     
    IMPLICIT INTEGER (a-z)
    EXTERNAL        lshift,orf
    LOGICAL :: abort,flush,flshal,ec,INT,dmiflg,fphys,fphys1,  &
        baddat,badfor,iax,iaxf,fphys2,lharm,secd,fthru
    INTEGER :: nam(2),onm(2),nm(2),t(7),ihill(2),ihoff(2),  &
        itsai(2),istrs(2),istrn(2),iall(2),isym(2), imem(2),isymm(2)
    REAL :: xm(100),z(100),xl,xl1,x1,x2,zseq,zseq1,oldxm3
    DOUBLE PRECISION :: da(2)
    CHARACTER (LEN=25) :: sfm
    CHARACTER (LEN=29) :: uim
    CHARACTER (LEN=25) :: uwm
    CHARACTER (LEN=23) :: ufm
    COMMON /machin/ mach
    COMMON /xmssg / ufm,uwm,uim,sfm
    COMMON /system/ ksystm(77)
    COMMON /xpfist/ ipfist
    COMMON /xfist / ifist(1)
    COMMON /xfiat / ifiat(2)
    COMMON /ifpx1 / nt1,t1(2,1)
    COMMON /ifpdta/ id,n,k,kx,ky,i(100),m(100),mf(100),m1(100),  &
        m1f(100),kn,baddat,badfor,nopen,nparam,iax,nax,  &
        iaxf,naxf,lharm,knt,slotdf(5),gc(7),ll(6)
    COMMON /zblpkx/ a(4),i0
    COMMON /zzzzzz/ ibuf(1)
    COMMON /cifs2p/ fphys,fphys1,km,dmiflg,ibcds,fthru,fphys2
    COMMON /xdpl  / p(3)
    COMMON /l15_l8/ l15,l8
 
    !     P(1) = NEXT AVAILABLE FILE ON POOL
    !     P(2) = TOTAL NUMBER OF POSSIBLE ENTRYS
    !     P(3) = CURRENT NUMBER OF ENTRYS PRESENT
    !     P(4) - P(3*P(2)+3) = THREE WORDS FOR EACH ENTRY AS FOLLOWS...
    !            1.  NAME(1)
    !            2.  NAME(2)
    !            3.  EQUIV FLAG, SIZE/1000, FILE NO. ON POOL
 
    EQUIVALENCE   (ksystm( 1) , nbuf  )  , (ksystm(24) , icfiat) ,  &
        (ksystm( 2) , nout  )  , (ksystm(55) , kprec ) ,  &
        (ksystm( 3) , abort )  , (ksystm(77) , bandit) ,  &
        (nrows,t(3)),(ifo,t(4)),(ty2,t(5)),(z(1),i(1)) , (xm(1),m(1)),(da(1),a(1))
 
    DATA   nam    / 4HISF2,   4HP    /
    DATA   endt   / 4HENDT /, skip   / 4HSKIP /, pool   / 4HPOOL /
    DATA   bcdblk / 4H     /, bcddet / 4HDET  /, bcdsdt / 4HSDET /,  &
        bcdudt / 4HUDET /, bcdinv / 4HINV  /, bcdsin / 4HSINV /,  &
        bcduin / 4HUINV /, bcdgiv / 4HGIV  /, bcdmgv / 4HMGIV /,  &
        bcdhes / 4HHESS /, bcdfer / 4HFEER /, bcdmas / 4HMASS /,  &
        bcdmax / 4HMAX  /, bcdpoi / 4HPOIN /,  &
        bcdq   / 4H-q   /, bcdt   / 4HT    /, bcdz   / 4H-x   /,  &
        bcdll  / 4HLL   /, bcdsl  / 4HSL   /, bcdls  / 4HLS   /
    DATA   thru   / 4HTHRU /, eigr   / 4HEIGR /, eigb   / 4HEIGB /
    DATA   dmi    / 4H dmi /, dti    / 4H dti /, dmig   / 4HDMIG /
    DATA   endrc1,  endrc2 /  4HENDR , 4HEC   /
    DATA   iscr1  / 301    /, icomp  / 1      /
    DATA   ihill  , ihoff   , itsai  , istrs   , istrn  /  &
        4HHILL , 4H      , 4HHOFF , 4H      , 4HTSAI , 4H      ,  &
        4HSTRE , 4HSS    , 4HSTRA , 4HIN   /
    DATA   iall   , isym    , imem   , isymm  /  &
        4HALL  , 4H      , 4HSYM  , 4H      , 4HMEM  , 4H      , 4HSYMM , 4HEM   /
    DATA   iyes,    ino    /  4HYES  , 4HNO   /
 
    !     =======================================================
    !     DMI AND DMIG MUST ACCOMODATE ALL KINDS OF SPECIAL FORMS
    !     E.G., IDENTITY MATRIX
    !     =======================================================
 
    IF (k > 100) GO TO 81
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5, 850, 850, 870,   5, 890,   5,  &
        5,   5, 920, 920, 920, 960, 920,   5,   5,   5 ), k
81  IF (kx > 100) GO TO 82
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,1190,1200,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5, 920, 920,   5,   5,   5,   5,   5, 920,  &
        960,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5, 920,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,1000,   5,   5,  &
        920,2900,   5,   5,   5,   5,   5,   5,   5,2000 ), kx
82  IF (ky > 100) GO TO 83
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        3100,3300,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,4100,  &
        4300,4500,4700,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5 ), ky
83  kz = k - 300
    IF (kz > 60) GO TO 5
    GO TO (     5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,   5,   5,   5,   5,   5,   5,   5,  &
        5,   5,   5,3200,   5,   5, 920,   5,   5,   5 ), kz
5   CALL page2 (2)
    WRITE  (nout,6) sfm
6   FORMAT (a25,' 322, ILLEGAL ENTRY TO IFS2P.')
    abort  =.true.
    RETURN 1
7   badfor =.true.
    RETURN 1
8   baddat =.true.
    RETURN 1
    3 DO  l = 1,n
        i(l) = m(l)
    END DO
2   RETURN
9   RETURN 3
 
    !*******       85-EIGR, 86-EIGB      ***********************************
 
850 IF (m(1) <= 0) GO TO 8
    IF (m(3) /= bcdblk .AND. m(2) /= bcdfer) GO TO 8
    IF (m(2) /= bcddet .AND. m(2) /= bcdsdt .AND. m(2) /= bcdudt .AND.  &
        m(2) /= bcdinv .AND. m(2) /= bcdsin .AND. m(2) /= bcduin .AND.  &
        m(2) /= bcdgiv .AND. m(2) /= bcdmgv .AND. m(2) /= bcdfer) GO TO 8
    IF (m(2) == bcdfer .AND. (m(3) /= bcdblk .AND. m(3) /= bcdq .AND.  &
        m(3) /= bcdz)) GO TO 8
    IF (m(10)+m(11) == 0) GO TO 852
    IF (m(10) /= bcdblk .OR. m(11) /= bcdblk) GO TO 860
852 nm(1) = eigr
    nm(2) = bcdmas
    IF (k == 85) GO TO 855
    nm(1) = eigb
    nm(2) = bcdmax
855 m(10) = nm(2)
    m(12) = 0
    m(13) = 0
    CALL mesage (30,222,nm)
    GO TO 865
860 IF ((m(10) /= bcdmas .OR. m(11) /= bcdblk) .AND.  &
        (m(10) /= bcdmax .OR. m(11) /= bcdblk) .AND.  &
        (m(10) /= bcdpoi .OR. m(11) /= bcdt )) GO TO 8
    IF (m(10) /= bcdpoi .AND. (m(12) /= 0 .OR. m(13) /= 0)) GO TO 8
    IF (m(10) == bcdpoi .AND. (m(12) <= 0 .OR. m(13) < 0)) GO TO 8
865 IF (m(6) == 0 .AND. m(2) /= bcdgiv .AND. m(2) /= bcdmgv .AND.  &
        m(2) /= bcdfer) GO TO 8
    IF (k == 86 .AND. (m(2) == bcdgiv .OR. m(2) == bcdmgv)) GO TO 8
    IF ((m(2) == bcddet .OR. m(2) == bcdsdt) .AND. xm(4) < 0.0) GO TO 8
    IF (m(2) == bcdudt .AND. xm(4) < 0.0) GO TO 8
    IF (k == 85 .AND. m(2) /= bcdgiv .AND. m(2) /= bcdmgv .AND.  &
        xm(4) < 0.0) GO TO 8
    IF (m(2) /= bcdgiv .AND. m(2) /= bcdmgv .AND. m(2) /= bcdfer .AND.  &
        xm(5) <= 0.0) GO TO 8
    IF (m(2) /= bcdgiv .AND. m(2) /= bcdmgv .AND. m(2) /= bcdfer .AND.  &
        xm(4) >= xm(5)) GO TO 8
    n = 18
    GO TO 3
 
    !*****         87-EIGC          **************************************
 
870 IF (km /= 0) GO TO 872
    IF (mf(1) /= 1 .OR. mf(2) /= 3 .OR. mf(3) /= 3 .OR. mf(4) /= 1  &
        .AND. mf(4) /= 0 .OR. mf(5) /= 1 .AND. mf(5) /= 0 .OR.  &
        mf(6) /= 2 .AND. mf(6) /= 0 .OR. mf(7) /= 0 .OR. mf(8) /= 0) GO TO 875
    IF (m(1) <= 0 .OR. m(2) /= bcddet .AND. m(2) /= bcdinv .AND.  &
        m(2) /= bcdhes .AND. m(2) /= bcdfer .OR. m(4) /= bcdmax .AND.  &
        (m(4) /= bcdpoi .OR. m(5) /= bcdt) .OR. xm(8) < 0.) GO TO 875
    IF (m(4) == bcdmax .AND. (m(6) /= 0 .OR. m(7) /= 0)) GO TO 875
    IF (m(4) == bcdpoi .AND. (m(6) <= 0 .OR. m(7) < 0)) GO TO 875
    IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 875
    n = 10
    874 DO  l = 1,n
        i(l) = m(l)
    END DO
    GO TO 876
    872 DO  l = 1,5
        IF (mf(l) /= 2 .AND. mf(l) /= 0) GO TO 875
    END DO
    IF (mf(6) /= 1 .AND. mf(6) /= 0 .OR. mf(7) /= 1 .AND. mf(7) /= 0  &
        .OR. mf(8) /= 0) GO TO 875
    IF (xm(5) <= 0.) xm(5) = 1.0
    IF (m(6) < 0 .OR. m(7) < 0) GO TO 875
    n = 7
    GO TO 874
875 baddat = .true.
876 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 877
    DO  l = 1,7
        n = n + 1
        i(n) =-1
    END DO
    km = 0
    kn = 0
    GO TO 9
877 kn = 1
    km = 1
    GO TO 9
 
    !*******       -BLANK CARD-        *************************************
 
890 IF (ibcds /= 0)RETURN 2
    ibcds = 1
    CALL page2 (2)
    WRITE  (nout,891) uwm
891 FORMAT (a25,' 324, BLANK CARD(S) IGNORED.')
    RETURN 2
 
    !*******       93- TABLEM1, 94-TABLEM2, 95-TABLEM3  ********************
    !              133-TABLED1,134-TABLED2,140-TABLED3
    !              162-TABDMP1, 97-TABLES1,191-TABRND1
    !              357-TABLEM5
    !                 (TABLEM5 IS DESIGNED FOR THERMAL COEFFICIENT WHICH IS
    !                          FUNCTION OF TIME
    !                          THIS PROJECT TEMPORARY HALTS HERE  6/90)
 
920 IF (km /= 0) GO TO 933
    i2 = m(1)
    items = 0
    n = 8
    IF (m(1) <= 0) baddat = .true.
    IF (mf(1) /= 1) badfor = .true.
    i(1) = i2
    DO  l = 2,7
        IF (mf(l) /= 0 .AND. mf(l) /= 2) badfor = .true.
        i(l) = m(l)
    END DO
 
    !     LOGARITHMIC SCALE
    !     I(8) = 0, LINEAR-LINEAR SCALE (ALL TABLES)
    !          = 1, LOG-LOG SCALE (TABLE-1 ONLY)
    !          = 2, LINEAR-LOG SCALE (TABLE-1, TABLE-2 AND TABLE-3)
    !          = 3, LOG-LINEAR SCALE (TABLE-1 ONLY)
    !     TABLE-1 INCLUDES TABLED1, TABLEM1, TABLES1, TABDMP1 AND TABRND1
    !     TABLE-2 INCLUDES TABLED2 AND TABLEM2
 
    i(8) = 0
    IF (mf(8) /= 3) GO TO 930
    IF (m(8) == bcdll) i(8) = 1
    IF (m(8) == bcdsl) i(8) = 2
    IF (m(8) == bcdls) i(8) = 3
    IF (m(8) /= bcdsl .AND. (k == 94 .OR. k == 95 .OR. k == 134 .OR.  &
        k == 140)) baddat = .true.
 
930 IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 932
    kn = 1
    km = 1
    GO TO 966
932 baddat = .true.
    kn = 0
    km = 0
    GO TO 966
933 l1 = 0
    DO  l = 1,7,2
        IF (mf(l) == 3 .OR. mf(l+1) == 3) GO TO 937
        IF (mf(l) /= 0 .AND. mf(l) /= 2 .OR. mf(l+1) /= 0 .AND.  &
            mf(l+1) /= 2) GO TO 943
        items = items + 1
        n  = n + 2
        l1 = l1 + 2
        i(n-1) = m(l1-1)
        i(n  ) = m(l1  )
        IF (items > 2) GO TO 935
        IF (items > 1) GO TO 934
        x1 = z(n-1)
        xl = x1
        GO TO 936
934     x2  = z(n-1)
        xl1 = xl
        xl  = x2
        zseq= SIGN(1.0,x2-x1)
        IF (x2 == x1) baddat = .true.
        GO TO 936
935     xl1 = xl
        xl  = z(n-1)
        zseq1 = SIGN(1.0,xl-xl1)
        IF (zseq1 /= zseq .AND. xl /= xl1) baddat = .true.
936     CYCLE
937     IF (mf(l) == 3) GO TO 938
        l1  = l1 + 1
        lp1 = l1
        kword1 = 0
        GO TO 939
938     l1  = l1 + 2
        lp1 = l1 - 1
        kword1 = m(lp1)
939     IF (mf(l+1) == 3) GO TO 941
        l1  = l1 + 1
        lp2 = l1
        kword2 = 0
        GO TO 942
941     l1  = l1 + 2
        lp2 = l1 - 1
        kword2 = m(lp2)
942     IF (kword1 == endt .OR. kword2 == endt) GO TO 961
        IF (kword1 == skip .OR. kword2 == skip) CYCLE
        baddat = .true.
        GO TO 956
    END DO
    GO TO 956
943 badfor = .true.
    GO TO 956
956 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 966
    kn = 0
    km = 0
    baddat = .true.
    GO TO 966
961 n = n + 2
    i(n-1) = -1
    i(n  ) = -1
    IF (xl  == xl1) baddat = .true.
    IF (items < 2) baddat = .true.
958 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 965
    kn = 0
    km = 0
    GO TO 966
965 kn = 1
    km = 1
    baddat = .true.
966 IF (baddat .OR. badfor) GO TO 968
    GO TO 2
968 m(1) = i2
    GO TO 8
 
    !*******      96-TABLEM4, 141-TABLED4     ******************************
 
960 IF (km /= 0) GO TO 964
    items = 0
    i2 = m(1)
    n  = 8
    IF (m(1)  <= 0) baddat = .true.
    IF (mf(1) /= 1) badfor = .true.
    i(1) = i2
    IF (m(3) == 0) baddat = .true.
    DO  l = 2,8
        IF (mf(l) /= 0 .AND. mf(l) /= 2 .OR. l >= 6 .AND. mf(l) /= 0)  &
            badfor = .true.
        i(l) = m(l)
    END DO
    i(8) = 0
    IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 963
    kn = 1
    km = 1
    GO TO 966
963 baddat = .true.
    kn = 0
    km = 0
    GO TO 966
964 l1 = 0
    DO  l = 1,8
        kword1 = 0
        IF (mf(l) == 3) GO TO 969
        IF (mf(l) /= 0 .AND. mf(l) /= 2) GO TO 943
        n = n + 1
        items = items + 1
        l1   = l1 + 1
        i(n) = m(l1)
        CYCLE
969     l1 = l1 + 2
        kword1 = m(l1-1)
        IF (kword1 == endt) GO TO 959
        baddat = .true.
        GO TO 956
    END DO
    GO TO 956
959 n = n + 1
    i(n) = -1
    IF (items < 1) baddat = .true.
    GO TO 958
 
    !*****     188-TABRNDG       **************************************
 
1000 IF (m(1) < 0) GO TO 8
    IF (m(2) < 1 .OR. m(2) > 2) GO TO 8
    i(1) = m(1)
    i(2) = m(2)
    i(3) = m(3)
    i(4) = m(4)
    i(5) = 0
    i(6) = 0
    i(7) = 0
    i(8) = 0
    i(9) =-1
    i(10)=-1
    n    = 10
    GO TO 2
 
    !******         119-DMI          ************************************
 
1190 IF (km /= 0) GO TO 8150
    IF (fphys) GO TO 8005
    IF (m(1) == nm(1) .AND. m(2) == nm(2)) GO TO 8100
    ASSIGN 8010 TO r
    GO TO 8973
8005 IF (p(1) > 1) dmiflg = .true.
    fphys = .false.
    nm(1) = 0
    nm(2) = 0
    IF (bandit /= -1 .AND. bandit /= -2) bandit = +9
8010 flush = .false.
    flshal= .false.
    ec    = .true.
    secd  = .false.
    t(1)  = iscr1
    DO  l = 2,7
        t(l) = 0
    END DO
    IF (m(3) /= 0) flush = .true.
    onm(1) = nm(1)
    onm(2) = nm(2)
    IF (mf(1) /= 3 .OR. m(1) == onm(1) .AND. m(2) == onm(2)) flush = .true.
    nm(1)  = m(1)
    nm(2)  = m(2)
    iprint = 0
    j0     = 0
    IF (p(1) <= p(2)) GO TO 8020
    flush  = .true.
    flshal = .true.
8020 ASSIGN 8025 TO r1
    ASSIGN 8030 TO r
    GO TO 200
8025 flush = .true.
8030 IF (flush) GO TO 8960
    ifo = m(4)
    ty1 = m(5)
    ty2 = m(6)
    IF (ty2 == 0 .AND. MOD(ty1,2) == 1) ty2 = ty1 + kprec - 1
    IF (ty2 == 0 .AND. MOD(ty1,2) == 0) ty2 = ty1
    IF (mach /= 12) GO TO 8033
    IF (ty2 == 2 .OR. ty2 == 4) ty2 = ty2 - 1
8033 CONTINUE
     IF (ty1 < 1 .OR. ty1 > 4 .OR. ty2 < 1 .OR. ty2 > 4) GO TO 8950
     IF (ty1 >= 3 .AND. ty2 <= 2) WRITE (nout,8035) uwm,dmi,nam(1), nam(2),knt
8035 FORMAT (a25,' 327A, ',a4,' CARD ',2A4,', SORTED CARD COUNT =',i7,  &
         ' SPECIFYING COMPLEX DATA INPUT', /5X,  &
         'AND REAL MATRIX OUTPUT MAY NOT MAKE SENSE',/)
     nrows = m(8)
     ncols = m(9)
     IF (ifo   > 8) GO TO 8950
     IF (mf(6) /= 0) GO TO 8950
     IF (nrows <= 0 .OR. ncols <= 0) GO TO 8950
     IF ((ifo == 1 .OR. ifo == 6 .OR. ifo == 8) .AND. (nrows /= ncols))  &
         GO TO 8950
     nbuf2 = 2*nbuf
     CALL OPEN (*9997,iscr1,ibuf(nbuf2+1),1)
     CALL WRITE (iscr1,nm,2,1)
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 8950
     IF (ifo == 8) GO TO 8040
     IF (m1(1) /= t1(1,k) .OR. m1(2) /= t1(2,k)) GO TO 8950
     GO TO 8960
8040 IF (m1(1) == t1(1,k) .AND. m1(2) == t1(2,k) .AND.  &
         m1(3) == nm(1)   .AND. m1(4) == nm(2) ) GO TO 8950
     GO TO 8960
8100 IF (.NOT.ec) GO TO 8950
     IF (flush) GO TO 8960
     ec = .false.
     IF (m(3) <= j0) GO TO 8950
8130 j0 = j0 + 1
     IF (m(3) == j0) GO TO 8140
     CALL bldpk (ty1,ty2,iscr1,0,0)
     CALL bldpkn (iscr1,0,t)
     GO TO 8130
8140 i0 = 1
     l1 = 4
     l1f=-1
     l2 = 9
     IF (ty1 == 2 .OR. ty1 == 4) l2 = 14
     IF (mf(3) /= 1 .OR. m(4) < i0) GO TO 8950
     i0  = m(4) - 1
     INT = .false.
     CALL bldpk (ty1,ty2,iscr1,0,0)
     GO TO 8155
8150 IF (j0 <= 0 .OR. j0 > ncols) GO TO 8950
     l1  = 1
     l1f = 0
     l2  = 8
     IF (ty1 == 2 .OR. ty1 == 4) l2 = 16
8155 l  = l1
8156 lf = l + l1f
     IF (fthru) GO TO 8192
     IF (mf(lf) == 0) GO TO 8300
     IF (mf(lf) == 2 .OR. mf(lf) == 4) GO TO 8180
     IF (mf(lf) == -32767) GO TO 8291
     IF (INT) GO TO 8950
     IF (mf(lf) == 3) GO TO 8191
     IF (mf(lf) /= 1 .OR. m(l) < i0 .OR. m(l) > nrows) GO TO 8950
     i0  = m(l)
     INT = .true.
     GO TO 8290
8180 SELECT CASE ( ty1 )
         CASE (    1)
             GO TO 8181
         CASE (    2)
             GO TO 8182
         CASE (    3)
             GO TO 8183
         CASE (    4)
             GO TO 8184
     END SELECT
     !   . REAL SINGLE PRECISION
8181 IF (mf(lf) == 4) GO TO 8950
     IF (flush .OR. m(l) == 0) GO TO 8190
     a(1) = m(l)
     GO TO 8185
     !   . REAL DOUBLE PRECISION
8182 IF (mf(lf) == 2) GO TO 8950
     a(1) = m(l  )
     a(2) = m(l+1)
     l    = l + 1
     l1f  = l1f - 1
     IF (flush .OR. da(1) == 0.0D0) GO TO 8190
     GO TO 8185
     !   . COMPLEX SINGLE PRECISION
8183 IF (mf(lf) == 4) GO TO 8950
     IF (secd) GO TO 8186
     a(1) = m(l)
     secd = .true.
     GO TO 8290
8186 a(2) = m(l)
     secd = .false.
     IF (a(1) == 0 .AND. a(2) == 0 .OR. flush) GO TO 8190
     GO TO 8185
     !   . COMPLEX DOUBLE PRECISION
8184 IF (mf(lf) == 2) GO TO 8950
     IF (secd) GO TO 8187
     a(1) = m(l  )
     a(2) = m(l+1)
     l    = l + 1
     l1f  = l1f - 1
     secd = .true.
     GO TO 8290
8187 a(3) = m(l  )
     a(4) = m(l+1)
     l    = l + 1
     l1f  = l1f - 1
     secd = .false.
     IF (flush .OR. da(1) == 0.0D0 .AND. da(2) == 0.0D0) GO TO 8190
 
     !     PACK AN ELEMENT
 
8185 CALL zblpki
8190 INT = .false.
     i0  = i0 + 1
     IF (i0 > nrows) GO TO 8300
     IF (l+1   > l2) GO TO 8290
     IF (mf(lf+1) /= 3) GO TO 8290
     l = l + 1
8191 IF (m(l) /= thru) GO TO 8950
     fthru = .true.
     l1f = l1f - 1
     l2  = l2  + 1
     l   = l   + 1
     IF (l >= l2) GO TO 8291
     l  = l + 1
     lf = l + l1f
8192 IF (mf(lf) /= 1 .OR. m(l) < i0 .OR. m(l) > nrows) GO TO 8950
8193 CALL zblpki
     i0 = i0 + 1
     IF (i0 <= m(l)) GO TO 8193
     fthru = .false.
     IF (i0 > nrows) GO TO 8300
8290 l = l + 1
     IF (l <= l2) GO TO 8156
8291 IF (m1(1) == 0 .AND. m1(2) == 0 .OR. INT) GO TO 8960
     GO TO 8315
8300 IF (l == l2) GO TO 8310
     lf = lf + 1
     DO  lx = lf,8
         IF (mf(lx) == -32767) EXIT
         IF (mf(lx) /= 0) GO TO 8950
     END DO
8310 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 8950
8315 IF (flush) GO TO 8320
     IF (secd ) GO TO 8950
     IF (fthru) GO TO 8950
     CALL bldpkn (iscr1,0,t)
8320 ec = .true.
     GO TO 8960
8950 flush = .true.
8960 IF (m1(1) == 0 .AND. m1(2) == 0 .OR. m1(1) == t1(1,k) .AND.  &
         m1(2) == t1(2,k)) GO TO 8970
     ASSIGN 8970 TO r
     GO TO 8973
8970 n = 0
     IF (.NOT.flush .OR. iprint /= 0) GO TO 1226
     CALL page2 (2)
     WRITE  (nout,8971) ufm,nm(1),nm(2),knt
8971 FORMAT (a23,' 325, BAD DATA OR FORMAT OR NON-UNIQUE NAME. DMI ',  &
         2A4,10X,' SORTED CARD COUNT =',i7)
     iprint = 1
     GO TO 1226
8973 IF (flshal) GO TO 8993
     IF (flush ) GO TO 8987
     IF (ifo == 8) GO TO 8995
8975 j0 = j0 + 1
     IF (j0 > ncols) GO TO 8977
     CALL bldpk  (ty1,ty2,iscr1,0,0)
     CALL bldpkn (iscr1,0,t)
     GO TO 8975
8977 IF (ncols == t(2)) GO TO 8978
     flush = .true.
     GO TO 8987
8978 CONTINUE
     CALL CLOSE (iscr1,1)
     CALL wrttrl (t)
     CALL rdtrl  (t)
     IF (icfiat == 11) GO TO 8982
     DO  lx = 1,3
         t(lx+1) = orf(lshift(t(2*lx),16),t(2*lx+1))
     END DO
     j = 3
     GO TO 8985
8982 j = 6
8985 CALL WRITE (pool,nm,2,0)
     CALL WRITE (pool,t(2),j,1)
     IF (l8 /= 0) WRITE (nout,8986) nm,dmi,(t(ip+1),ip=1,j)
8986 FORMAT ('0*** DIAG  8 MESSAGE -- TRAILER FOR DATA BLOCK ',2A4,  &
         ' (VIA ',a4,' CARDS) = ',5I7,i9)
     CALL gopen (iscr1,ibuf(2*nbuf+1),2)
     CALL cpyfil (iscr1,pool,ibuf(3*nbuf+1),nopen,nwords)
     CALL CLOSE (iscr1,1)
     CALL eof (pool)
     dmiflg = .true.
     p(1) = p(1) + 1
8987 ip   = 3*p(3) + 4
     p(ip  ) = nm(1)
     p(ip+1) = nm(2)
     IF (flush) nwords = 0
     p(ip+2) = orf(lshift(nwords/1000,16),p(1)-1)
     p(3   ) = p(3) + 1
     IF (.NOT.flush) GO TO 8992
     CALL CLOSE (iscr1,1)
     CALL eof (pool)
     p(1) = p(1) + 1
     CALL skpfil (pool,-1)
     IF (dmiflg) CALL eof (pool)
8990 abort = .true.
8992 GO TO r, (8010,8970)
8993 WRITE  (nout,8994) sfm,nm(1),nm(2)
8994 FORMAT (a25,' 326, NO ROOM IN /XDPL/ FOR DMI ',2A4)
     CALL page2 (2)
     GO TO 8990
8995 t(2) = ncols
     GO TO 8978
9997 CALL mesage (-1,iscr1,nm)
 
     !******          120-DMIG          ********************************
 
1200 IF (.NOT.fphys1) GO TO 1202
     fphys1 = .false.
     nm(1)  = 0
     nm(2)  = 0
1202 ierr   = 0
     IF (km   /= 0) GO TO 1208
     IF (m(3) == 0) GO TO 1206
     IF (m(1) /= nm(1) .OR. m(2) /= nm(2)) GO TO 1218
     IF (mf(2) /= 1 .OR. mf(3) /= 1 .AND. mf(3) /= 0) GO TO 1218
     IF (m(3) <= 0 .OR. m(4) < 0 .OR. m(4) > 6) GO TO 1218
     IF (mf(4) /= 0) GO TO 1218
     IF (mf(5) /= 1 .OR. mf(6) /= 1 .AND. mf(6) /= 0) GO TO 1218
     IF (m(6) <= 0 .OR. m(7) < 0 .OR. m(7) > 6) GO TO 1218
     IF (mf(7)+ity1 /= 4 .AND. mf(7) /= 0) GO TO 1218
     IF ((ty1 == 1 .OR.  ty1 == 2)  .AND. mf(8) /= 0 .OR.  &
         ty1 == 3 .AND. mf(8) /= 2 .AND. mf(8) /= 0 .OR.  &
         ty1 == 4 .AND. mf(8) /= 4 .AND. mf(8) /= 0) GO TO 1218
     n = 5
     i(n-4) = m(3)
     i(n-3) = m(4)
     i(n-2) = m(6)
     i(n-1) = m(7)
     i(n  ) = m(8)
     IF (ty1 == 1) GO TO 1204
     n = 6
     i(n) = m(9)
     IF (ty1 /= 4) GO TO 1204
     n = 8
     i(n-1) = m(10)
     i(n  ) = m(11)
1204 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1230
     n = n + 2
     i(n-1) = -1
     i(n  ) = -1
     GO TO 1216
1206 IF (mf(1) /= 3 .OR. m(1) == nm(1) .AND. m(2) == nm(2)) GO TO 1218
     ifo = m(4)
     ty1 = m(5)
     ity1= 2*MOD(ty1,2)
     ty2 = m(6)
     IF (ty2 == 0 .AND. MOD(ty1,2) == 1) ty2 = ty1 + kprec - 1
     IF (ty2 == 0 .AND. MOD(ty1,2) == 0) ty2 = ty1
     IF (mach /= 12) GO TO 1207
     IF (ty2 == 2.OR.ty2 == 4) ty2 = ty2 - 1
1207 CONTINUE
     IF (ty1 <= 0 .OR. ty1 > 4 .OR. ty2 <= 0 .OR. ty2 > 4) GO TO 1218
     IF (ty1 >= 3 .AND. ty2 <= 2) WRITE (nout,8035) uwm,dmig,nm(1), nm(2),knt
     IF (ifo /= 1 .AND. ifo /= 2 .AND. ifo /= 6) GO TO 1218
     IF (ty2 == 1 .AND. ty1 == 3) GO TO 1218
     nm(1) = m(1)
     nm(2) = m(2)
     IF (mf(6) /= 0 .OR. mf(7) /= 0 .OR. mf(8) /= 0) GO TO 1220
     IF (m1f(2) /= 3 .OR. m1(3) /= nm(1) .OR. m1(4) /= nm(2)) GO TO 1220
     m(6) = ty2
     n = 9
     GO TO 3
1208 lf = 1
     l  = 1
1210 IF (m(l) /= 0 .OR. m(l+1) /= 0 .OR. m(l+2) /= 0 .OR. m(l+3) /= 0)  &
         GO TO 1212
     lf = lf + 4
     l  = l  + 4
     GO TO 1214
1212 IF (m(l) <= 0 .OR. m(l+1) < 0 .OR. m(l+1) > 6) GO TO 1220
     IF (mf(lf) /= 1 .OR. mf(lf+1) /= 1 .AND. mf(lf+1) /= 0) GO TO 1220
     ierr = 1
     IF (mf(lf+2)+ity1 /= 4 .AND. mf(lf+2) /= 0) GO TO 1220
     IF (mf(lf+3) /= 0 .AND. ty1 /= 3 .AND. ty1 /= 4) GO TO 1220
     n = n + 3
     i(n-2) = m(l  )
     i(n-1) = m(l+1)
     i(n  ) = m(l+2)
     lf = lf + 4
     l  = l  + 4
     IF (ty1 == 1) GO TO 1214
     n = n + 1
     i(n) = m(l-1)
     IF (ty1 == 2) l = l + 1
     IF (ty1 /= 4) GO TO 1214
     n = n + 2
     i(n-1) = m(l  )
     i(n  ) = m(l+1)
     l  = l + 2
1214 IF (lf <= 7) GO TO 1210
     IF (n  <= 0) GO TO 1220
     IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1230
     n = n + 2
     i(n-1) = -1
     i(n  ) = -1
1216 IF (m1(1) == t1(1,k) .AND. m1(2) == t1(2,k) .AND.  &
         m1(3) == nm(1  ) .AND. m1(4) == nm(2 )) GO TO 1228
     n = n + 2
     i(n-1) = -1
     i(n  ) = -1
     GO TO 1228
1218 nm(1) = m(1)
     nm(2) = m(2)
1220 abort = .true.
     CALL page2 (2)
     WRITE  (nout,1222) ufm,nm(1),nm(2),knt
1222 FORMAT (a23,' 327, BAD DATA OR FORMAT OR NON-UNIQUE NAME. DMIG ',  &
         2A4,10X,' SORTED CARD COUNT =',i7)
     IF (ierr == 1) WRITE (nout,1224)
1224 FORMAT (5X,'INPUT MATRIX TYPE (TIN) AND INPUT DATA (XIJ OR YIJ) ',  &
         'ARE NOT CONSISTANT')
1226 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 1230
1228 km = 0
     kn = 0
     GO TO 2
1230 kn = 1
     km = 1
     GO TO 2
 
     !*******        200 - DTI       ****************************************
 
2000 IF (km /= 0) GO TO 2120
     IF (fphys2) GO TO 2010
     IF (m(1) == nm(1) .AND. m(2) == nm(2)) GO TO 2100
     ASSIGN 2020 TO r
     GO TO 2300
2010 IF (p(1) > 1) dmiflg = .true.
     fphys2 = .false.
     nm(1)  = 0
     nm(2)  = 0
2020 flush  = .false.
     flshal = .false.
     IF (m(3) /= 0) flush = .true.
     onm(1) = nm(1)
     onm(2) = nm(2)
     IF (mf(1) /= 3 .OR. m(1) == onm(1) .AND. m(2) == onm(2)) flush = .true.
     nm(1)  = m(1)
     nm(2)  = m(2)
     iprint = 0
     nwords = 2
     j0 = 0
     IF (p(1) <= p(2)) GO TO 2050
     flush  = .true.
     flshal = .true.
2050 ASSIGN 2055 TO r1
     ASSIGN 2056 TO r
     GO TO 200
2055 flush = .true.
2056 IF (flush) GO TO 2195
     itrlt = 0
     DO  l = 2,7
         itrlt = itrlt+m(l+2)
         IF (icfiat == 8  .AND. (m(l+2) < 0 .OR. m(l+2) > 65535)) flush = .true.
         !     2147483647 = 2**31-1
         IF (icfiat == 11 .AND. (m(l+2) < 0 .OR. m(l+2) > 2147483647))  &
             flush = .true.
         t(l) = m(l+2)
     END DO
     IF (itrlt /= 0) GO TO 2080
     DO  l = 2,7
         t(l) = 32767
     END DO
2080 CONTINUE
     CALL WRITE (pool,nm,2,0)
     IF (icfiat == 11) GO TO 2087
     DO  lx = 1,3
         t(lx+1) = orf(lshift(t(2*lx),16),t(2*lx+1))
     END DO
     l = 3
     GO TO 2090
2087 l = 6
2090 CALL WRITE (pool,t(2),l,1)
     CALL WRITE (pool,nm,2,0)
     IF (l8 /= 0) WRITE (nout,8986) nm,dti,(t(ip+1),ip=1,j)
     IF (m1(1) == t1(1,k) .AND. m1(2) == t1(2,k)) CALL WRITE (pool,nm,0,1)
     GO TO 2200
2100 j0 = j0 + 1
     IF (m(3) /= j0) GO TO 2190
     l1  = 4
     l1f =-1
     GO TO 2150
2120 l1 = 1
     l1f= 0
2150 l  = l1
     lf = l + l1f
2160 IF (mf(lf) == 3 .AND. m(l) == endrc1 .AND. m(l+1) == endrc2) GO TO 2180
     IF (mf(lf) > 2) l = l + 1
     l  = l  + 1
     lf = lf + 1
     IF (mf(lf) >= 0) GO TO 2160
     CALL WRITE (pool,m(l1),l-l1,0)
     nwords = nwords + l - l1
     IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 2190
     GO TO 2200
2180 CALL WRITE (pool,m(l1),l-l1,1)
     nwords = nwords + l - l1
     IF (m1(1) /= 0 .OR. m1(2) /= 0) GO TO 2195
2190 flush = .true.
2195 IF (m1(1) == 0 .AND. m1(2) == 0 .OR.  &
         m1(1) == t1(1,k) .AND. m1(2) == t1(2,k)) GO TO 2200
     ASSIGN 2200 TO r
     GO TO 2300
2200 n = 0
     IF (.NOT.flush .OR. iprint /= 0) GO TO 1226
     CALL page2 (2)
     WRITE (nout,2350) ufm,nm(1),nm(2),knt
     iprint = 1
     GO TO 1226
2300 IF (flshal) GO TO 2370
     IF (flush ) GO TO 2330
     CALL eof (pool)
     dmiflg = .true.
     p(1)   = p(1) + 1
2330 ip = 3*p(3) + 4
     p(ip  ) = nm(1)
     p(ip+1) = nm(2)
     IF (flush) nwords = 0
     p(ip+2) = orf(lshift(nwords/1000,16),p(1)-1)
     p(3) = p(3) + 1
     IF (.NOT.flush) GO TO 2365
     CALL page2 (2)
     WRITE  (nout,2350) ufm,nm(1),nm(2),knt
2350 FORMAT (a23,' 317, BAD DATA OR FORMAT OR NON-UNIQUE NAME FOR DTI '  &
         ,       2A4,10X,'SORTED CARD COUNT =',i7)
     CALL eof (pool)
     p(1) = p(1) + 1
     CALL skpfil (pool,-1)
     IF (dmiflg) CALL skpfil (pool,+1)
2360 abort = .true.
2365 GO TO r, (2020,2200)
2370 WRITE  (nout,2380) sfm,nm(1),nm(2)
2380 FORMAT (a25,' 318, NO ROOM IN /XDPL/ FOR DTI ',2A4)
     CALL page2 (2)
     GO TO 2360
 
 !     ******************************************************************
 
 !     CHECK NAME FOR UNIQUENESS AMONG DMI CARDS, DTI CARDS, ETC. AND
 !     RESERVED NAMES
 
200 CONTINUE
 
    !     CHECK  FIST, FIAT, DPL FOR A NAME MATCH
 
    DO  ii = 1,ipfist
        IF (nm(1) == ifist(2*ii+1) .AND. nm(2) == bcdblk) GO TO 250
    END DO
    nfiat = icfiat*ifiat(2) - 2
    DO  ii = 4,nfiat,icfiat
        IF (nm(1) == ifiat(ii) .AND. nm(2) == ifiat(ii+1)) GO TO 250
    END DO
    ndpl = p(3)*3 + 1
    DO  ii = 4,ndpl,3
        IF (nm(1) == p(ii) .AND. nm(2) == p(ii+1)) GO TO 250
    END DO
    GO TO r,  (8030,2056)
250 GO TO r1, (8025,2055)
 
    !*******      192-PLOAD4     ****************************************
 
2900 IF (km == 1) GO TO 2940
    km = 1
    kn = 1
    IF (mf(1) /= 1) badfor = .true.
    IF (.NOT.(mf(2) == 2 .AND. mf(3) == 1 .AND. mf(4) == 0 .AND.  &
        mf(5) == 0 .AND. mf(6) == 0)) GO TO 2905
 
    !     SPECIAL - ALLOWING PLOAD4 TO TAKE ON PLOAD2 FORMAT
    !     (PLOAD4,SID,P1,E1,blank,blank,blank,"THRU",E2) FOR QUICK INPUT
    !     DATA SWITCHING.  INTERCHAGNE 2ND AND 3RD FIELDS
 
    mf(2) = 1
    mf(3) = 2
    l     = m(2)
    m(2)  = m(3)
    m(3)  = l
2905 IF (mf(2) /= 1) badfor = .true.
    DO  l = 3,6
        IF (mf(l) /= 2 .AND. mf(l) /= 0) badfor = .true.
    END DO
    IF (mf(7) /= 3 .AND. mf(7) /= 0 .AND.  &
        .NOT.(mf(7) == 1 .AND. m(7) == 0)) badfor = .true.
    IF (mf(8) /= 1 .AND. mf(8) /= 0) badfor = .true.
    IF (mf(7) == 0 .AND. mf(8) /= 0) badfor = .true.
    IF (mf(7) == 3 .AND. mf(8) /= 1) badfor = .true.
    IF (m(1) <= 0) baddat = .true.
    IF (m(2) <= 0) baddat = .true.
    IF (mf(7) == 3 .AND. m(7) /= thru) baddat = .true.
    IF (mf(7) == 3 .AND. m(9) <=   0) baddat = .true.
    IF (mf(7) == 3 .AND. m(9) <= m(2)) baddat = .true.
    l1 = 0
    IF (mf(7) == 3) l1 = 1
    DO  l = 1,6
        i(l) = m(l)
    END DO
    i(7) = -1
    IF (l1 == 1) i(7) = 0
    i(8) = m(l1+8)
    n = 8
    IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
    DO  l = 9,12
        i(l) = 0
    END DO
    n  = 12
    km = 0
    kn = 0
    GO TO 9
 
2940 IF (mf(1) > 1) badfor = .true.
    DO  l = 2,4
        IF (mf(l) /= 2 .AND. mf(l) /= 0) badfor = .true.
    END DO
    IF (mf(1) == 0) m(1) = 0
    IF (m(1)  < 0) baddat = .true.
    DO  l = 1,4
        i(l) = m(l)
    END DO
    n  = 4
    km = 0
    kn = 0
    GO TO 9
 
    !*******       261-CQUAD4    ****************************************
 
3100 IF (mf(2) == 0) m(2) = m(1)
    i(1) = m(1)
    DO  l = 2,6
        IF (mf(l) /= 1) badfor = .true.
        IF (m(l)  <= 0) baddat = .true.
        i(l) = m(l)
    END DO
    l1 = 6
    DO  l = 11,14
        l1 = l1 + 1
        i(l1) = m(l)
    END DO
    IF (mf(7) /= 1 .AND. mf(7) /= 2 .AND. mf(7) /= 0) badfor = .true.
    IF (mf(7) == 1 .AND. (m(7) < 0 .OR. m(7) >= 1000000)) baddat = .true.
    i(11) = m(7)
    i(12) = 0
    IF (mf(7) == 1) i(12) = 1
    i(13) = m(8)
    n = 13
    GO TO 9
 
    !*******       354-CTRIA3      **************************************
 
3200 IF (mf(2) == 0) m(2) = m(1)
    i(1) = m(1)
    DO  l = 2,5
        IF (mf(l) /= 1) badfor = .true.
        IF (m(l)  <= 0) baddat = .true.
        i(l) = m(l)
    END DO
    IF (mf(6) /= 1 .AND. mf(6) /= 2 .AND. mf(6) /= 0) badfor = .true.
    IF (mf(6) == 1 .AND. (m(6) < 0 .OR. m(6) >= 1000000)) baddat = .true.
    i( 6) = m(11)
    i( 7) = m(12)
    i( 8) = m(13)
    i( 9) = m(6)
    i(10) = 0
    i(11) = m(7)
    IF (mf(6) == 1) i(10) = 1
    n = 11
    GO TO 9
 
    !*******        262-MAT8      ****************************************
 
3300 IF (mf(2) == 0 .OR. mf(3) == 0 .OR. mf(5) == 0) GO TO 7
    IF (m(1) <= 0) GO TO 8
    IF (xm(2) == 0.0 .OR. xm(3) == 0.0) GO TO 8
    IF (xm(5) <= 0.0) GO TO 8
    IF (mf(12) == 2 .AND. xm(12) <= 0.0) GO TO 8
    IF (mf(14) == 2 .AND. xm(14) <= 0.0) GO TO 8
    IF (mf(16) == 2 .AND. xm(16) <= 0.0) GO TO 8
    IF (mf(13) == 0) xm(13) = xm(12)
    IF (mf(15) == 0) xm(15) = xm(14)
    n = 18
    GO TO 3
 
    !*******        280-PCOMP     ****************************************
 
4100 kn = 1
    IF (icomp > 1) GO TO 4140
    icomp = 2
    IF (mf(1) /= 1) badfor = .true.
    IF (mf(2) /= 2 .AND. mf(2) /= 0) badfor = .true.
    IF (mf(3) /= 2 .AND. mf(3) /= 0) badfor = .true.
    IF (mf(4) /= 2 .AND. mf(4) /= 0) badfor = .true.
    IF (mf(5) /= 3 .AND. mf(5) /= 0) badfor = .true.
    l = 0
    IF (mf(5) == 3) l = 1
    IF (mf(6) /= 0) badfor = .true.
    IF (mf(7) /= 0) badfor = .true.
    IF (mf(8) /= 3 .AND. mf(8) /= 0   ) badfor = .true.
    IF (m(1) <= 0 .OR. m(1) >= 1000000) baddat = .true.
    IF (mf(5) == 3 .AND. xm(4) <= 0.0 ) baddat = .true.
    failur = -1
    IF (mf(5)  == 0) failur = 0
    IF (failur == 0) GO TO 4120
    IF (m(5) == ihill(1) .AND. m(6) == ihill(2)) failur = 1
    IF (m(5) == ihoff(1) .AND. m(6) == ihoff(2)) failur = 2
    IF (m(5) == itsai(1) .AND. m(6) == itsai(2)) failur = 3
    IF (m(5) == istrs(1) .AND. m(6) == istrs(2)) failur = 4
    IF (m(5) == istrn(1) .AND. m(6) == istrn(2)) failur = 5
    IF (failur == -1) baddat = .true.
4120 lamopt = -1
    IF (mf(8)  == 0) lamopt = 0
    IF (lamopt == 0) GO TO 4130
    IF (m(l+8) == iall (1) .AND. m(l+9) == iall (2)) lamopt = 0
    IF (m(l+8) == isym (1) .AND. m(l+9) == isym (2)) lamopt = 1
    IF (m(l+8) == imem (1) .AND. m(l+9) == imem (2)) lamopt = 2
    IF (m(l+8) == isymm(1) .AND. m(l+9) == isymm(2)) lamopt = 3
    IF (lamopt == -1) baddat = .true.
4130 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4135
    badfor = .true.
    kn = 0
    icomp = 1
    n = 0
    GO TO 9
4135 i(1) = m(1)
    i(2) = m(2)
    i(3) = m(3)
    i(4) = m(4)
    i(5) = failur
    i(6) = 0
    i(7) = 0
    i(8) = lamopt
    n    = 8
    GO TO 9
 
4140 n = 0
    DO  l = 1,2
        l1 = 4*(l-1)
        l2 = l1
        DO  l3 = 1,4
            IF (mf(l1+l3) /= 0) GO TO 4160
        END DO
        IF (l == 1) badfor = .true.
        GO TO 4195
4160    IF (l == 2 .AND. mf(4) == 3) l2 = l2 + 1
        IF (icomp == 3) GO TO 4170
        icomp = 3
        IF (mf(1) /= 1) badfor = .true.
        IF (mf(2) /= 2) badfor = .true.
        IF (mf(3) /= 2 .AND. mf(3) /= 0) badfor = .true.
        IF (m(1)  <=  0) baddat = .true.
        IF (xm(2) <= 0.0) baddat = .true.
        GO TO 4180
4170    IF (mf(l1+1) /= 1 .AND. mf(l1+1) /= 0) badfor = .true.
        IF (mf(l1+2) /= 2 .AND. mf(l1+2) /= 0) badfor = .true.
        IF (mf(l1+3) /= 2 .AND. mf(l1+3) /= 0) badfor = .true.
        IF (mf(l1+1) == 1 .AND. m (l2+1) <= 0) baddat = .true.
        IF (mf(l1+1) == 0) m(l2+1) = iold1
        IF (mf(l1+2) == 2 .AND. xm(l2+2) <= 0.0) baddat = .true.
        IF (mf(l1+2) == 0) m(l2+2) = iold2
        IF (mf(l1+3) == 0) m(l2+3) = iold3
4180    IF (mf(l1+4) /= 3 .AND. mf(l1+4) /= 0) badfor = .true.
        IF (mf(l1+4) == 3 .AND. (m(l2+4) /= iyes .AND. m(l2+4) /= ino))  &
            baddat = .true.
        iout = 0
        IF (m(l2+4) == iyes) iout = 1
        i(n+1) = m(l2+1)
        i(n+2) = m(l2+2)
        i(n+3) = m(l2+3)
        i(n+4) = iout
        iold1  = m(l2+1)
        iold2  = m(l2+2)
        iold3  = m(l2+3)
        n = n + 4
    END DO
    IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
4195 kn = 0
    icomp  = 1
    i(n+1) =-1
    n = n + 1
    GO TO 9
 
    !*******        281-PCOMP1    ****************************************
 
4300 kn = 1
    IF (icomp > 1) GO TO 4340
    icomp = 2
    IF (mf(1) /= 1) badfor = .true.
    IF (mf(2) /= 2 .AND. mf(2) /= 0) badfor = .true.
    IF (mf(3) /= 2 .AND. mf(3) /= 0) badfor = .true.
    IF (mf(4) /= 2 .AND. mf(4) /= 0) badfor = .true.
    IF (mf(5) /= 3 .AND. mf(5) /= 0) badfor = .true.
    l = 0
    IF (mf(5) == 3) l = 1
    IF (mf(6) /= 1) badfor = .true.
    IF (mf(7) /= 2) badfor = .true.
    IF (mf(8) /= 3 .AND. mf(8) /= 0   ) badfor = .true.
    IF (m(1) <= 0 .OR. m(1) >= 1000000) baddat = .true.
    IF (mf(5) == 3 .AND. xm(4) <= 0.0 ) baddat = .true.
    failur = -1
    IF (mf(5)  == 0) failur = 0
    IF (failur == 0) GO TO 4320
    IF (m(5) == ihill(1) .AND. m(6) == ihill(2)) failur = 1
    IF (m(5) == ihoff(1) .AND. m(6) == ihoff(2)) failur = 2
    IF (m(5) == itsai(1) .AND. m(6) == itsai(2)) failur = 3
    IF (m(5) == istrs(1) .AND. m(6) == istrs(2)) failur = 4
    IF (m(5) == istrn(1) .AND. m(6) == istrn(2)) failur = 5
    IF (failur == -1) baddat = .true.
4320 IF (m(l+6) <=  0) baddat = .true.
    IF (xm(l+7) <= 0.0) baddat = .true.
    lamopt = -1
    IF (mf(8)  == 0) lamopt = 0
    IF (lamopt == 0) GO TO 4330
    IF (m(l+8) == iall (1) .AND. m(l+9) == iall (2)) lamopt = 0
    IF (m(l+8) == isym (1) .AND. m(l+9) == isym (2)) lamopt = 1
    IF (m(l+8) == imem (1) .AND. m(l+9) == imem (2)) lamopt = 2
    IF (m(l+8) == isymm(1) .AND. m(l+9) == isymm(2)) lamopt = 3
    IF (lamopt == -1) baddat = .true.
4330 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4335
    badfor = .true.
    kn = 0
    icomp = 1
    n = 0
    GO TO 9
4335 i(1) = m(1)
    i(2) = m(2)
    i(3) = m(3)
    i(4) = m(4)
    i(5) = failur
    i(6) = m(l+6)
    i(7) = m(l+7)
    i(8) = lamopt
    n    = 8
    GO TO 9
 
4340 n = 0
    DO  l = 1,8
        IF (mf(l) /= 0) GO TO 4360
        IF (l == 1) badfor = .true.
        GO TO 4395
4360    IF (icomp == 3) GO TO 4370
        icomp = 3
        IF (mf(1) /= 2) badfor = .true.
        GO TO 4380
4370    IF (mf(l) /= 2 .AND. mf(l) /= 0) badfor = .true.
        IF (mf(l) == 0) m(l) = iold1
4380    i(n+1) = m(l)
        iold1  = m(l)
        n = n + 1
    END DO
    IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
4395 kn = 0
    icomp  = 1
    i(n+1) =-1
    n = n + 1
    GO TO 9
 
    !*******        282-PCOMP2    ****************************************
 
4500 kn = 1
    IF (icomp > 1) GO TO 4540
    icomp = 2
    IF (mf(1) /= 1) badfor = .true.
    IF (mf(2) /= 2 .AND. mf(2) /= 0) badfor = .true.
    IF (mf(3) /= 2 .AND. mf(3) /= 0) badfor = .true.
    IF (mf(4) /= 2 .AND. mf(4) /= 0) badfor = .true.
    IF (mf(5) /= 3 .AND. mf(5) /= 0) badfor = .true.
    l = 0
    IF (mf(5) == 3) l = 1
    IF (mf(6) /= 1) badfor = .true.
    IF (mf(7) /= 0) badfor = .true.
    IF (mf(8) /= 3 .AND. mf(8) /= 0   ) badfor = .true.
    IF (m(1) <= 0 .OR. m(1) >= 1000000) baddat = .true.
    IF (mf(5) == 3 .AND. xm(4) <= 0.0 ) baddat = .true.
    failur = -1
    IF (mf(5)  == 0) failur = 0
    IF (failur == 0) GO TO 4520
    IF (m(5) == ihill(1) .AND. m(6) == ihill(2)) failur = 1
    IF (m(5) == ihoff(1) .AND. m(6) == ihoff(2)) failur = 2
    IF (m(5) == itsai(1) .AND. m(6) == itsai(2)) failur = 3
    IF (m(5) == istrs(1) .AND. m(6) == istrs(2)) failur = 4
    IF (m(5) == istrn(1) .AND. m(6) == istrn(2)) failur = 5
    IF (failur == -1) baddat = .true.
4520 IF (m(l+6) <=  0) baddat = .true.
    lamopt = -1
    IF (mf(8)  == 0) lamopt = 0
    IF (lamopt == 0) GO TO 4530
    IF (m(l+8) == iall (1) .AND. m(l+9) == iall (2)) lamopt = 0
    IF (m(l+8) == isym (1) .AND. m(l+9) == isym (2)) lamopt = 1
    IF (m(l+8) == imem (1) .AND. m(l+9) == imem (2)) lamopt = 2
    IF (m(l+8) == isymm(1) .AND. m(l+9) == isymm(2)) lamopt = 3
    IF (lamopt == -1) baddat = .true.
4530 IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 4535
    badfor = .true.
    kn = 0
    icomp = 1
    n = 0
    GO TO 9
4535 i(1) = m(1)
    i(2) = m(2)
    i(3) = m(3)
    i(4) = m(4)
    i(5) = failur
    i(6) = m(l+6)
    i(7) = 0
    i(8) = lamopt
    n    = 8
    GO TO 9
 
4540 n = 0
    DO  l = 1,4
        l1 = 2*(l-1)
        DO  l3 = 1,2
            IF (mf(l1+l3) /= 0) GO TO 4560
        END DO
        IF (l == 1) badfor = .true.
        GO TO 4595
4560    IF (icomp == 3) GO TO 4570
        icomp = 3
        IF (mf(1)  /= 2) badfor = .true.
        IF (mf(2)  /= 2) badfor = .true.
        IF (xm(1) <= 0.0) baddat = .true.
        GO TO 4580
4570    IF (mf(l1+1) /= 2 .AND. mf(l1+1) /= 0) badfor = .true.
        IF (mf(l1+2) /= 2 .AND. mf(l1+2) /= 0) badfor = .true.
        IF (mf(l1+1) == 2 .AND. xm(l1+1) <= .0) baddat = .true.
        IF (mf(l1+1) == 0) m(l1+1) = iold1
        IF (mf(l1+2) == 0) m(l1+2) = iold2
4580    i(n+1) = m(l1+1)
        i(n+2) = m(l1+2)
        iold1  = m(l1+1)
        iold2  = m(l1+2)
        n = n + 2
    END DO
    IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
4595 kn = 0
    icomp  = 1
    i(n+1) =-1
    n = n + 1
    GO TO 9
 
    !*******        283-PSHELL    ****************************************
 
4700 IF (km == 1) GO TO 4740
    km = 1
    kn = 1
    IF (mf( 1) /= 1                  ) badfor = .true.
    IF (mf( 2) /= 1 .AND. mf( 2) /= 0) badfor = .true.
    IF (mf( 3) /= 2 .AND. mf( 3) /= 0) badfor = .true.
    IF (mf( 4) /= 1 .AND. mf( 4) /= 0) badfor = .true.
    IF (mf( 5) /= 2 .AND. mf( 5) /= 0) badfor = .true.
    IF (mf( 6) /= 1 .AND. mf( 6) /= 0) badfor = .true.
    IF (mf( 7) /= 2 .AND. mf( 7) /= 0) badfor = .true.
    IF (mf( 8) /= 2 .AND. mf( 8) /= 0) badfor = .true.
    IF (m(1) <= 0) baddat = .true.
    IF (mf(2) == 1 .AND.  m(2) <= 0) baddat = .true.
    IF (mf(4) == 1 .AND.  m(4) <= 0) baddat = .true.
    IF (mf(4) /= 0 .AND. mf(5) == 0) xm(5) = 1.0
    IF (mf(6) == 1 .AND.  m(6) <= 0) baddat = .true.
    IF (mf(6) /= 0 .AND. mf(4) == 0) baddat = .true.
    IF (mf(6) /= 0 .AND. mf(7) == 0) xm(7) = 0.833333
    DO  l = 2,6,2
        IF (m(l) == 0 .AND. xm(l+1) > 0.0) baddat = .true.
    END DO
    DO  l = 1,8
        i(l) = m(l)
    END DO
    iolmf2 = mf(2)
    iolmf4 = mf(4)
    ioldm2 =  m(2)
    ioldm4 =  m(4)
    ioldm6 =  m(6)
    oldxm3 = xm(3)
    n = 8
    IF (m1(1) == 0 .AND. m1(2) == 0) GO TO 9
    z( 9) = -0.5*oldxm3
    z(10) =  0.5*oldxm3
    DO  l = 11,17
        i(l) = 0
    END DO
    n  = 17
    km = 0
    kn = 0
    GO TO 9
 
4740 IF (mf(1) /= 2 .AND. mf(1) /= 0) badfor = .true.
    IF (mf(2) /= 2 .AND. mf(2) /= 0) badfor = .true.
    IF (mf(3) /= 1 .AND. mf(3) /= 0) badfor = .true.
    IF (mf(4) /= 1 .AND. mf(4) /= 2 .AND. mf(4) /= 0) badfor = .true.
    IF (mf(5) /= 1 .AND. mf(5) /= 2 .AND. mf(5) /= 0) badfor = .true.
    IF (mf(6) /= 2 .AND. mf(6) /= 0) badfor = .true.
    IF (mf(1) == 0) xm(1) = -0.5*oldxm3
    IF (mf(2) == 0) xm(2) =  0.5*oldxm3
    IF (mf(3) == 1 .AND. m(3) <= 0) baddat = .true.
    IF (mf(3) /= 0 .AND. (iolmf2 == 0 .OR. iolmf4 == 0)) baddat = .true.
    IF (mf(3) /= 0 .AND. (m(3) == ioldm2 .OR. m(3) == ioldm4)) baddat = .true.
    IF (mf(4) == 1 .AND. m(4) < 0) baddat = .true.
    IF (mf(5) == 1 .AND. m(5) < 0) baddat = .true.
    IF (ioldm2 == 0 .AND. ioldm4 == 0 .AND.  &
        ioldm6 == 0 .AND. m(3) == 0) baddat = .true.
    DO  l = 1,4
        i(l) = m(l)
    END DO
    i(5) = 0
    IF (mf(4) == 1) i(5) = 1
 
    !     I(6) IS THE INTEGRATION ORDER (SET TO 0)
 
    !     NOTE
    !     ----
 
    !     THE INTEGRATION ORDER IS NOT USED IN THE PROGRAM,
    !     BUT THIS WORD IS REQUIRED BECAUSE OF THE DESIGN
    !     OF THE EST DATA FOR THE CQUAD4 ELEMENT.
 
    i(6) = 0
    i(7) = m(5)
    i(8) = 0
    IF (mf(5) == 1) i(8) = 1
    i(9) = m(6)
    n  = 9
    km = 0
    kn = 0
    GO TO 9
 
END SUBROUTINE ifs2p
