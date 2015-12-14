SUBROUTINE ascm06 (NAME,iphase,isol,nogo)
     
!     RECOVER, MRECOVER COMMAND DMAP DATA
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,2),xtra(15),subnam(2),rdmap(18,17),  &
     rdmap1(18,9),rdmap2(18,8),oct(3,1),oct1(3,1),  &
     ptbs(7,31),ptbs1(7,18),ptbs2(7,13)
 COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(543)
 EQUIVALENCE    (rdmap1(1,1),rdmap(1, 1)),(ptbs1(1,1),ptbs(1,1)),  &
     (rdmap2(1,1),rdmap(1,10)),(ptbs2(1,1),ptbs(1,19)), (oct1(1,1),oct(1,1))
 DATA comnd   / 4HRECO    , 17    , 15    ,  1    , 31    ,  0  ,  &
     4HMREC    , 17    , 15    ,  1    , 31    ,  0  /
 DATA slash   / 1H/       /
 DATA rdmap 1 /  &
     4HFILE,4H    ,4H  u1,4H=app,4HEND/,4HU2=a,4HPPEN,4HD/u3,4H=app,  &
4HEND/,4HU4=a,4HPPEN,4HD/u5,4H=app,4HEND ,4H$   ,4H    ,4H    ,  &
    4HPARA,4HM   ,4H  //,4H*add,4H*/il,4HOOP/,4H0/0 ,4H$   ,4H    ,  &
    4H    , 8*4H    , 4HLABE,4HL   ,4H  lb,4HSTP ,4H$   ,13*4H    ,  &
    4HRCOV,4HR   ,4H  ca,4HSESS,4H,geo,4HM4,k,4HGG,m,4HGG,p,4HGG,u,  &
    4HGV ,,4HDIT,,4HDLT,,4HBGG,,4HK4GG,4H,ppf,4H/oug,4HV1 ,,4H    ,  &
    4H    ,4H    ,4H  op,4HG1,o,4HQG1,,4HU1,u,4H2,u3,4H,u4,,4HU5/s,  &
    4H,n,d,4HRY/s,4H,n,i,4HLOOP,4H/stp,4H!*NA,4HMEFS,4HS */,4H    ,  &
    4H    ,4H    ,4H  ns,4HOL/n,4HEIGV,4H/s,n,4H,lui,4H/s,n,4H,u1n,  &
    4H/s,n,4H,u2n,4H/s,n,4H,u3n,4H/s,n,4H,u4n,4H/s,n,4H,u5n,4H/   ,  &
    4H    ,4H    ,4H  s,,4HN,no,4HSORT,4H2/v,,4HY,ut,4HHRES,4HH/v,,  &
    4HY,pt,4HHRES,4HH/v,,4HY,qt,4HHRES,4HH $ , 3*4H    ,  &
    4HEQUI,4HV   ,4H  ou,4HGV1 ,4H,oug,4HV /n,4HOSOR,4HT2/o,4HQG1,,  &
    4HOQG/,4HNOSO,4HRT2 ,4H$   , 5*4H    ,  &
    4HEQUI,4HV   ,4H  op,4HG1,o,4HPG/n,4HOSOR,4HT2 $,4H    ,4H    ,  &
    4H    , 8*4H        /
DATA rdmap 2 /  &
    4HCOND,4H    ,4H  ns,4HT2ST,4HP,no,4HSORT,4H2 $ ,4H    ,4H    ,  &
    4H    , 8*4H    ,  &
    4HSDR3,4H    ,4H  ou,4HGV1 ,4H,opg,4H1,oq,4HG1,,,4H,/ou,4HGV ,,  &
    4HOPG,,4HOQG,,4H,, $, 6*4H    ,  &
    4HLABE,4HL   ,4H  ns,4HT2ST,4HP $ ,13*4H    ,  &
    4HOFP ,4H    ,4H  ou,4HGV ,,4HOPG,,4HOQG,,4H,,//,4HS,n,,4HCARD,  &
    4HNO $, 8*4H    , 4HCOND,4H    ,4H  lb,4HBSTP,4H,ilo,4HOP $,12*4H    ,  &
    4HREPT,4H    ,4H  lb,4HSTP,,4H100 ,4H$   ,12*4H    ,  &
    4HLABE,4HL   ,4H  lb,4HBSTP,4H $  ,13*4H    ,  &
    4HSOFO,4H    ,4H  ,u,4H1,u2,4H,u3,,4HU4,u,4H5//-,4H1!*X,4HXXXX,  &
    4HXXX*,4H $  , 7*4H        /
DATA xtra    / 4HPRIN,4HSAVE,4HDISP,4HOLOA,4HSPCF,4HMODE,4HRANG ,  &
    4HSUBC,4HSORT,4HBASI,4HVELO,4HACCE,4HENER,4HUIMP , 4HSTEP       /
DATA oct  1  / 9    ,    262144    ,         0     /
DATA ptbs 1  / 3  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
    4  , 11  , 15  ,  2  ,4HCASE  ,         0  ,  0  ,  &
    4  , 18  , 18  ,  5  ,4HGORL  ,         0  ,  0  ,  &
    4  , 24  , 27  ,  0  ,4HNAME  ,         1  ,  0  ,  &
    4  , 28  , 31  ,  0  ,4HNAME  ,         2  ,  0  ,  &
    4  , 32  , 32  ,  3  ,4HPVEC  ,         0  ,  0  ,  &
    4  , 36  , 36  ,  4  ,4HUVEC  ,         0  ,  0  ,  &
    4  , 41  , 44  ,  0  ,4HNAME  ,    458752  ,  0  ,  &
    4  , 45  , 48  ,  0  ,4HNAME  ,    458752  ,  0  ,  &
    4  , 49  , 52  ,  0  ,4HNAME  ,    458768  ,  0  ,  &
    4  , 53  , 57  ,  0  ,4HNAME  ,    458784  ,  0  ,  &
    4  , 58  , 58  ,  3  ,4HPFTL  ,         0  ,  0  ,  &
    4  , 62  , 62  ,  6  ,4HOVEC  ,         0  ,  0  ,  &
    5  , 11  , 15  ,  0  ,4HNAME  ,    262144  ,  0  ,  &
    5  , 54  , 54  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
    5  , 57  , 59  ,  8  ,4HNAME  ,         0  ,  0  ,  &
    6  , 11  , 11  ,  4  ,4HSOL   ,         0  ,  0  ,  &
    6  , 16  , 21  ,  0  ,4HNAME  ,   1769472  ,  0  /
DATA ptbs 2  / 8  , 11  , 11  ,  6  ,4HOVEC  ,         0  ,  0  ,  &
    8  , 18  , 18  ,  5  ,4HOVC2  ,         0  ,  0  ,  &
    10  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
    11  , 11  , 11  ,  6  ,4HOVEC  ,         0  ,  0  ,  &
    11  , 18  , 22  ,  0  ,4HNAME  ,    262144  ,  0  ,  &
    11  , 31  , 31  ,  5  ,4HOVC2  ,         0  ,  0  ,  &
    11  , 37  , 40  ,  0  ,4HNAME  ,    262144  ,  0  ,  &
    12  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
    13  , 11  , 11  ,  5  ,4HOVC2  ,         0  ,  0  ,  &
    13  , 17  , 20  ,  0  ,4HNAME  ,    262144  ,  0  ,  &
    14  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
    15  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
    16  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  /
DATA subnam  / 4HASCM,2H06  /

!     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION))

rdmap(15,5) = khrfn1(rdmap(15,5),1,slash,1)
rdmap(8,17) = khrfn1(rdmap(8,17),2,slash,1)

!     VALIDATE COMMAND AND SET POINTERS

DO  i = 1,2
  IF (NAME == comnd(1,i)) GO TO 20
END DO
GO TO 70
20 icomnd = i
irdm   = 1
nrdm   = comnd(2,icomnd)
ixtra  = irdm  + 18*nrdm
nxtra  = comnd(3,icomnd)
ioct   = ixtra + nxtra
noct   = comnd(4,icomnd)
iptbs  = ioct  + 3*noct
nptbs  = comnd(5,icomnd)
iph    = iptbs + 7*nptbs
nph    = comnd(6,icomnd)

!     MOVE RDMAP DATA

k = 0
IF (nrdm == 0) GO TO 35
DO  j = 1,nrdm
  DO  i = 1,18
    k = k + 1
    idat(k) = rdmap(i,j)
  END DO
END DO
35 CONTINUE

!     MOVE XTRA DATA

IF (nxtra == 0) GO TO 45
DO  i = 1,nxtra
  k = k + 1
  idat(k) = xtra(i)
END DO
45 CONTINUE

!     MOVE OCT DATA

IF (noct == 0) GO TO 55
DO  j = 1,noct
  DO  i = 1,3
    k = k + 1
    idat(k) = oct(i,j)
  END DO
END DO
55 CONTINUE

!     MOVE PTBS DATA

IF (nptbs == 0) GO TO 65
DO  j = 1,nptbs
  DO  i = 1,7
    k = k + 1
    idat(k) = ptbs(i,j)
  END DO
END DO
65 CONTINUE

RETURN

!     INPUT ERROR

70 CALL mesage (7,0,subnam)
nogo = 1
RETURN

END SUBROUTINE ascm06
