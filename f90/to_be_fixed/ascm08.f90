SUBROUTINE ascm08 (NAME,iphase,isol,nogo)
     
!     SOLVE COMMAND DMAP DATA FOR DYNAMIC ANALYSIS
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,1),subnam(2),isave(21),rdmap(18,55),  &
     rdmap1(18,9),rdmap2(18,9),rdmap3(18,9),  &
     rdmap4(18,9),rdmap5(18,9),rdmap6(18,9),  &
     rdmap7(18,1),oct(3,23),oct1(3,18),oct2(3,5),  &
     ptbs(7,25),ptbs1(7,18),ptbs2(7,7)
 COMMON /phas28/ ipas28(14)
 COMMON /asdbd / irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs,  &
     iph,nph,idat(1248)
 EQUIVALENCE     (rdmap1(1,1),rdmap(1, 1)),(oct1(1,1),oct(1,1)),  &
     (rdmap2(1,1),rdmap(1,10)),(oct2(1,1),oct(1,19)),  &
     (rdmap3(1,1),rdmap(1,19)),(ptbs1(1,1),ptbs(1,1)),  &
     (rdmap4(1,1),rdmap(1,28)),(ptbs2(1,1),ptbs(1,19)),  &
     (rdmap5(1,1),rdmap(1,37)), (rdmap6(1,1),rdmap(1,46)),  &
     (rdmap7(1,1),rdmap(1,55))
 DATA comnd    / 4HSOLV    , 55    ,  0    , 23    , 25    , 14 /
 DATA slash    / 1H/       /
 DATA isave    / 4,11,3, 13,10,1, 13,14,3, 13,16,2, 54,8,2, 54,9,2, 54,10,2 /
 DATA rdmap 1  / 4HALTE,4HR   ,4H  (g,4HP1) ,4H$   ,13*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*nop,4H*/al,4HWAYS,4H=-1 ,4H$   ,4H    ,  &
     4H    ,8*4H    ,  &
     4HSGEN,4H    ,4H  ca,4HSECC,4H,geo,4HM3,g,4HEOM4,4H,dyn,4HAMIC,  &
     4HS/ca,4HSESS,4H,cas,4HEI,g,4HPL,e,4HQEXI,4HN,gp,4HDT, ,4H    ,  &
     4H    ,4H    ,4H  bg,4HPDT,,4HSIL,,4HGE3S,4H,ge4,4HS,dy,4HNS/s,  &
     4H,n,d,4HRY!*,4HNAME,4HSOLS,4H*/S,,4HN,LU,4HSET/,4H    ,4H    ,  &
     4H    ,4H    ,4H  s,,4HN,no,4HGPDT,4H $  ,12*4H    ,  &
     4HPURG,4HE   ,4H  cs,4HTM $,14*4H    ,  &
     4HEQUI,4HV   ,4H  GE,4H3S,g,4HEOM3,4H/alw,4HAYS/,4HGE4S,4H,geo,  &
     4HM4/a,4HLWAY,4HS/ca,4HSEI,,4HCASE,4HCC/a,4HLWAY,4HS/  ,4H    ,  &
     4H    ,4H    ,4H  dy,4HNS,d,4HYNAM,4HICS/,4HALWA,4HYS $,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4HSTP,,4HDRY ,4H$   ,12*4H     /
 DATA rdmap 2  / 4HALTE,4HR   ,4H  (p,4HLOT),4H $  ,13*4H    ,  &
     4HALTE,4HR   ,4H  (c,4HOND),4H $  ,13*4H    ,  &
     4HALTE,4HR   ,4H  (g,4HPWG),4H $  ,13*4H    ,  &
     4HSOFI,4H    ,4H  /k,4HNOS,,4HMNOS,4H,bno,4HS,k4,4HNOS,,4H/dry,  &
     4H!*NA,4HMESO,4HLS*/,4H*KMT,4HX*!*,4HMMTX,4H*!*B,4HMTX*,4H/   ,  &
     4H    ,4H    ,4H  *k,4H4MX*,4H $  ,13*4H    ,  &
     4HEQUI,4HV   ,4H  kn,4HOS,k,4HGG/n,4HOKGG,4HX $ ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4H2K,n,4HOKGG,4HX $ ,12*4H    ,  &
     4HADD ,4H    ,4H  kg,4HGX,k,4HNOS/,4HKGG/,4H(1.0,4H,0.0,4H)/(1,  &
     4H.0,0,4H.0) ,4H$   ,6*4H    , 4HLABE,4HL   ,4H  lb,4H2K $,14*4H       /
 DATA rdmap 3  /  &
     4HEQUI,4HV   ,4H  mn,4HOS,m,4HGG/n,4HOMGG,4H $  ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4H2M,n,4HOMGG,4H $  ,12*4H    ,  &
     4HADD ,4H    ,4H  mg,4HG,mn,4HOS/m,4HGGX/,4H(1.0,4H,0.0,4H)/(1,  &
     4H.0,0,4H.0) ,4H$   ,6*4H    ,  &
     4HEQUI,4HV   ,4H  mg,4HGX,m,4HGG/a,4HLWAY,4HS $ ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HLABE,4HL   ,4H  lb,4H2M $,14*4H    ,  &
     4HEQUI,4HV   ,4H  bn,4HOS,b,4HGG/n,4HOBGG,4H $  ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4H2B,n,4HOBGG,4H $  ,12*4H    ,  &
     4HADD ,4H    ,4H  bg,4HG,bn,4HOS/b,4HGGX/,4H(1.0,4H,0.0,4H)/(1,  &
     4H.0,0,4H.0) ,4H$   ,6*4H    ,  &
     4HEQUI,4HV   ,4H  bg,4HGX,b,4HGG/a,4HLWAY,4HS $ ,4H    ,4H    ,  &
     4H    ,8*4H         /
 DATA rdmap 4  / 4HLABE,4HL   ,4H  lb,4H2B $,14*4H    ,  &
     4HEQUI,4HV   ,4H  k4,4HNOS,,4HK4GG,4H/nok,4H4GG ,4H$   ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4H2K4,,4HNOK4,4HGG $,12*4H    ,  &
     4HADD ,4H    ,4H  k4,4HGG,k,4H4NOS,4H/k4g,4HGX/ ,4H(1.0,4H,0.0,  &
     4H)/(1,4H.0,0,4H.0) ,4H$   ,5*4H    ,  &
     4HEQUI,4HV   ,4H  k4,4HGGX,,4HK4GG,4H/alw,4HAYS ,4H$   ,4H    ,  &
     4H    ,8*4H    , 4HLABE,4HL   ,4H  lb,4H2K4 ,4H$   ,13*4H    ,  &
     4HLABE,4HL   ,4H  lb,4HSTP ,4H$   ,13*4H    ,  &
     4HCHKP,4HNT  ,4H  mg,4HG,bg,4HG,k4,4HGG $,12*4H    ,  &
     4HALTE,4HR   ,4H  (p,4HARAM,4H) $ ,13*4H           /
 DATA rdmap 5  /  &
     4HPARA,4HM   ,4H  //,4H*AND,4H*/md,4HEMA/,4HNOUE,4H/nom,4H2PP ,  &
     4H$   ,8*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*add,4H*/kd,4HEK2/,4H1/0 ,4H$   ,4H    ,  &
     4H    ,8*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*add,4H*/no,4HMGG/,4H1/0 ,4H$   ,4H    ,  &
     4H    ,8*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*add,4H*/no,4HBGG/,4H1/0 ,4H$   ,4H    ,  &
     4H    ,8*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*add,4H*/no,4HK4GG,4H/1/0,4H $  ,4H    ,  &
     4H    ,8*4H    , 4HALTE,4HR   ,4H  (e,4HQUIV,4H) $ ,13*4H    ,  &
     4HEQUI,4HV   ,4H  k2,4HDD,k,4HDD/k,4HDEK2,4H $  ,4H    ,4H    ,  &
     4H    ,8*4H    ,  &
     4HEQUI,4HV   ,4H  m2,4HDD,m,4HDD/n,4HOMGG,4H $  ,4H    ,4H    ,  &
     4H    ,8*4H    ,  &
     4HEQUI,4HV   ,4H  b2,4HDD,b,4HDD/n,4HOBGG,4H $  ,4H    ,4H    ,  &
     4H    ,8*4H         /
 DATA rdmap 6  / 4HALTE,4HR   ,4H  (s,4HDR2),4H $  ,13*4H    ,  &
     4HEQUI,4HV   ,4H  up,4HVF,u,4HPVC/,4HNOA ,4H$   ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4HL19,,4HNOA ,4H$   ,12*4H    ,  &
     4HSDR1,4H    ,4H  us,4HETD,,4H,udv,4HF,,,,4HGOD,,4HGMD,,4H,,,/,  &
     4HUPVC,4H,,/1,4H/dyn,4HAMIC,4HS $ , 4*4H    ,  &
     4HLABE,4HL   ,4H  lb,4HL19 ,4H$   ,13*4H    ,  &
     4HCHKP,4HNT  ,4H  up,4HVC $,14*4H    ,  &
     4HEQUI,4HV   ,4H  up,4HVC,u,4HGV/n,4HOUE ,4H$   ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4HUE,n,4HOUE ,4H$   ,12*4H    ,  &
     4HUPAR,4HTN  ,4H  us,4HET,u,4HPVC/,4HUGV,,4HUEV,,4H,!*P,4H*!*G,  &
     4H*!*E,4H* $ ,7*4H       /
 DATA rdmap 7  / 4HLABE,4HL   ,4H  lb,4HUE $,14*4H      /
 DATA oct 1    / 15    ,         0    ,         1  ,  &
     16    ,         0    ,         1  , 17    ,         0    ,         1  ,  &
     18    ,         0    ,         1  , 19    ,         0    ,         2  ,  &
     20    ,         0    ,         2  , 21    ,         0    ,         2  ,  &
     22    ,         0    ,         2  , 23    ,         0    ,         2  ,  &
     24    ,         0    ,        16  , 25    ,         0    ,        16  ,  &
     26    ,         0    ,        16  , 27    ,         0    ,        16  ,  &
     28    ,         0    ,        16  , 29    ,         0    ,        32  ,  &
     30    ,         0    ,        32  , 31    ,         0    ,        32  ,  &
     32    ,         0    ,        32  /
 DATA oct 2   / 33    ,         0    ,        32  ,  &
     38    ,         0    ,         1  , 39    ,         0    ,         2  ,  &
     40    ,         0    ,        16  , 41    ,         0    ,        32  /
 DATA ptbs 1  / 1  , 11  , 11  ,  5  ,     1  ,         0  ,  0  ,  &
     4  , 43  , 45  ,  8  ,4HNAME  ,         0  ,  0  ,  &
     9  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     10  , 11  , 11  ,  6  ,     2  ,         0  ,  0  ,  &
     11  , 11  , 11  ,  6  ,     3  ,         0  ,  0  ,  &
     12  , 11  , 11  ,  6  ,     4  ,         0  ,  0  ,  &
     13  , 12  , 13  ,  3  ,4HNANO  ,         1  , -1  ,  &
     13  , 17  , 18  ,  3  ,4HNANO  ,         2  , -1  ,  &
     13  , 22  , 23  ,  3  ,4HNANO  ,        16  , -1  ,  &
     13  , 27  , 29  ,  3  ,4HNANO  ,        32  , -1  ,  &
     13  , 37  , 39  ,  8  ,4HNAME  ,         0  ,  0  ,  &
     15  , 11  , 12  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     17  , 16  , 17  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     19  , 11  , 12  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     21  , 15  , 16  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     24  , 11  , 12  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     26  , 15  , 16  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     29  , 11  , 13  ,  3  ,4HNANO  ,         0  ,  0  /
 DATA ptbs 2  / 31  , 16  , 18  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     34  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     36  , 11  , 11  ,  7  ,     5  ,         0  ,  0  ,  &
     42  , 11  , 11  ,  7  ,     6  ,         0  ,  0  ,  &
     46  , 11  , 11  ,  6  ,     7  ,         0  ,  0  ,  &
     47  , 11  , 11  ,  4  ,4HDVEC  ,         0  ,  0  ,  &
     49  , 18  , 18  ,  4  ,4HDVEC  ,         0  ,  0  /
 DATA subnam  / 4HASCM,2H08  /
 
!     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION))
 
 DO  l = 1,21,3
   i = isave(l+1)
   j = isave(l  )
   k = isave(l+2)
   rdmap(i,j) = khrfn1(rdmap(i,j),k,slash,1)
 END DO
 
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
 
!     MOVE PHASE 2 DATA
 
 IF (iphase /= 2 .OR. nph == 0) GO TO 100
 DO  i = 1,nph
   k = k + 1
   idat(k) = ipas28(i)
 END DO
 GO TO 200
 100 CONTINUE
 
 200 RETURN
 
!     INPUT ERROR
 
 1000 CALL mesage (7,0,subnam)
 nogo = 1
 RETURN
 
END SUBROUTINE ascm08
