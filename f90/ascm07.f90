SUBROUTINE ascm07 (NAME,iphase,isol,nogo)
     
!     BRECOVER COMMAND DMAP DATA
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,1),subnam(2),rdmap(18,21),rdmap1(18,9),  &
     rdmap2(18,9),rdmap3(18,3),oct(3,13),oct1(3,13),  &
     ptbs(7,26),ptbs1(7,18),ptbs2(7,8)
 COMMON /phas37/ ipas37(6)
 COMMON /asdbd / irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs,  &
     iph,nph,idat(605)
 EQUIVALENCE     (rdmap1(1,1),rdmap(1, 1)),(oct1(1,1),oct(1,1)),  &
     (rdmap2(1,1),rdmap(1,10)),(ptbs1(1,1),ptbs(1,1)),  &
     (rdmap3(1,1),rdmap(1,19)),(ptbs2(1,1),ptbs(1,19))
 DATA comnd    / 4HBREC  , 21    ,  0    , 13    , 26    ,  6  /
 DATA slash    / 1H/     /
 DATA rdmap1  / 4HALTE,4HR   ,4H  (s,4HOLVE,4H) $ ,13*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*nop,4H*/al,4HWAYS,4H=-1 ,4H$   ,4H    ,  &
     4H    , 8*4H    ,  &
     4HSSG1,4H    ,4H  sl,4HT,bg,4HPDT,,4HCSTM,4H,sil,4H,est,4H,mpt,  &
     4H,gpt,4HT,ed,4HT,mg,4HG,ca,4HSECC,4H,dit,4H,/pg,4H,,,,,4H/   ,  &
     4H    ,4H    ,4H  lu,4HSET/,4HNSKI,4HP $ ,12*4H    ,  &
     4HSSG2,4H    ,4H  us,4HET,g,4HM,ys,4H,kfs,4H,GO,,4H,pg/,4HQR,p,  &
     4HO,ps,4H,pl ,4H$   , 6*4H    ,  &
     4HRCOV,4HR3  ,4H  ,p,4HG,ps,4H,po,,4HYS/u,4HAS ,,4HQAS,,4HPGS,,  &
     4HPSS,,4HPOS,,4HYSS,,4HLAMA,4H/sol,4HN!*N,4HAME ,4H   *,4H/   ,  &
     4H    ,4H    ,4H  no,4HUE $,14*4H    ,  &
     4HEQUI,4HV   ,4H  pg,4HS,pg,4H/alw,4HAYS ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  ps,4HS,ps,4H/alw,4HAYS ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H        /
 DATA rdmap2 /  &
     4HEQUI,4HV   ,4H  po,4HS,po,4H/alw,4HAYS ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  ys,4HS,ys,4H/alw,4HAYS ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H    , 4HCOND,4H    ,4H  lb,4HSSTP,4H,omi,4HT $ ,12*4H    ,  &
     4HFBS ,4H    ,4H  lo,4HO,,p,4HOS/u,4HOOV/,4H1/1/,4HPREC,4H/0 $,  &
     4H    , 8*4H    , 4HLABE,4HL   ,4H  lb,4HSSTP,4H $  ,13*4H    ,  &
     4HOFP ,4H    ,4H  la,4HMA,,,4H,,,/,4H/car,4HDNO ,4H$   ,4H    ,  &
     4H    , 8*4H    , 4HALTE,4HR   ,4H  (s,4HDR1),4H $  ,13*4H    ,  &
     4HUMER,4HGE  ,4H  us,4HET,q,4HAS,/,4HQGS/,4H*g*/,4H*a*/,4H*o* ,  &
     4H$   , 8*4H    ,  &
     4HADD ,4H    ,4H  qg,4H ,qg,4HS/qg,4HT/  ,4H(1.0,4H,0.0,4H)/(1,  &
     4H.0,0,4H.0) ,4H$   ,6*4H      /
 DATA rdmap3 /  &
     4HEQUI,4HV   ,4H  qg,4HT,qg,4H /al,4HWAYS,4H $  ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  ca,4HSECC,4H,cas,4HEXX/,4HALWA,4HYS $,4H    ,  &
     4H    , 8*4H    , 4HALTE,4HR   ,4H  (r,4HEPT),4H $  ,13*4H        /
 DATA oct1   / 3    ,    983040    ,        12  ,  &
     4    ,    983040    ,        12  , 5    ,    524288    ,        12  ,  &
     8    ,   1835008    ,        12  , 9    ,   1835008    ,        12  ,  &
     10    ,   1835008    ,        12  , 11    ,   1835008    ,        12  ,  &
     12    ,   1835008    ,        12  , 13    ,   1835008    ,        12  ,  &
     14    ,   1835008    ,        12  , 15    ,   1769472    ,         0  ,  &
     20    ,    458752    ,         0  , 21    ,    458752    ,         0  /
 DATA ptbs1  / 1  , 11  , 11  ,  0  ,     1  ,         0  ,  0  ,  &
     5  ,  1  ,  1  ,  0  ,4HNAME  ,         0  ,  0  ,  &
     5  , 19  , 21  ,  0  ,4HNAME  ,   1048576  ,  0  ,  &
     5  , 33  , 35  ,  0  ,4HNAME  ,         0  ,  0  ,  &
     5  , 36  , 38  ,  0  ,4HNAME  ,         0  ,  0  ,  &
     5  , 42  , 44  ,  0  ,4HNAME  ,         0  ,  0  ,  &
     6  , 12  , 14  ,  0  ,4HNAME  ,    524300  ,  0  ,  &
     6  , 15  , 17  ,  0  ,4HNAME  ,    524300  ,  0  ,  &
     6  , 18  , 20  ,  0  ,4HNAME  ,   1572876  ,  0  ,  &
     6  , 21  , 23  ,  0  ,4HNAME  ,   1572876  ,  0  ,  &
     6  , 24  , 24  ,  4  ,4HUAPH  ,         0  ,  0  ,  &
     6  , 33  , 33  ,  3  ,4HPGVC  ,         0  ,  0  ,  &
     6  , 37  , 37  ,  3  ,4HPSVC  ,         0  ,  0  ,  &
     6  , 41  , 44  ,  0  ,4HNAME  ,   1572876  ,  0  ,  &
     6  , 45  , 48  ,  0  ,4HNAME  ,   1572876  ,  0  ,  &
     6  , 49  , 49  ,  4  ,4HDYNT  ,    196608  ,  0  ,  &
     6  , 54  , 54  ,  4  ,4HSOL   ,         0  ,  0  ,  &
     6  , 60  , 60  ,  8  ,4HNAME  ,         0  ,  0  /
 DATA ptbs2  / 7  , 11  , 15  ,  0  ,4HNAME  ,    458752  ,  0  ,  &
     12  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     13  , 29  , 29  ,  4  ,4HPREC  ,         0  ,  0  ,  &
     14  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     16  , 11  , 11  ,  0  ,     2  ,         0  ,  0  ,  &
     18  , 11  , 11  ,  3  ,4HQVEC  ,         0  ,  0  ,  &
     19  , 15  , 15  ,  3  ,4HQVEC  ,         0  ,  0  ,  &
     21  , 11  , 11  ,  0  ,     3  ,         0  ,  0  /
 DATA subnam  / 4HASCM,2H07  /
 
!     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION))
 
 rdmap(15,6) = khrfn1(rdmap(15,6),2,slash,1)
 
!     VALIDATE COMMAND AND SET POINTERS
 
 IF (NAME /= comnd(1,1)) GO TO 1000
 icomnd = 1
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
 
!     MOVE PHASE 3 DATA
 
 IF (iphase /= 3 .OR. nph == 0) GO TO 200
 DO  i = 1,nph
   k = k + 1
   idat(k) = ipas37(i)
 END DO
 
 200 RETURN
 
!     INPUT ERROR
 
 1000 CALL mesage (7,0,subnam)
 nogo = 1
 RETURN
 
END SUBROUTINE ascm07
