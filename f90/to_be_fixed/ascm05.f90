SUBROUTINE ascm05 (NAME,iphase,isol,nogo)
     
!     SOLVE COMMAND DMAP DATA
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,1),subnam(2),rdmap(18,28),rdmap1(18,9),  &
     rdmap2(18,9),rdmap3(18,9),rdmap4(18,1),oct(3,5),  &
     oct1(3,5),ptbs(7,20),ptbs1(7,18),ptbs2(7,2)
 COMMON /phas25/ ipas25(14)
 COMMON /asdbd / irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs,  &
     iph,nph,idat(673)
 EQUIVALENCE     (rdmap1(1,1),rdmap(1, 1)),(oct1(1,1),oct(1,1)),  &
     (rdmap2(1,1),rdmap(1,10)),(ptbs1(1,1),ptbs(1,1)),  &
     (rdmap3(1,1),rdmap(1,19)),(ptbs 2(1,1),ptbs(1,19))  &
     ,               (rdmap4(1,1),rdmap(1,28))
 DATA comnd   / 4HSOLV    , 28    ,  0    ,  5    , 20    , 14  /
 DATA slash   /  1H/       /
 DATA rdmap 1 / 4HALTE,4HR   ,4H  (g,4HP1) ,4H$   ,13*4H    ,  &
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
 DATA rdmap 2 / 4HALTE,4HR   ,4H  (p,4HLOT),4H $  ,13*4H    ,  &
     4HALTE,4HR   ,4H  (c,4HOND),4H $  ,13*4H    ,  &
     4HCOND,4H    ,4H  lb,4HSOL,,4HNOSI,4HMP $,12*4H    ,  &
     4HALTE,4HR   ,4H  (o,4HPTP),4H $  ,13*4H    ,  &
     4HCOND,4H    ,4H  lb,4HSOL,,4HNOMG,4HG $ ,12*4H    ,  &
     4HALTE,4HR   ,4H  (s,4HMA3),4H $  ,13*4H    ,  &
     4HLABE,4HL   ,4H  lb,4HSOL ,4H$   ,13*4H    ,  &
     4HSOFI,4H    ,4H  /k,4HNOS,,4HMNOS,4H,,,/,4HDRY/,4H*nam,4HESOL,  &
     4HS*!*,4HKMTX,4H*!*M,4HMTX*,4H $  , 4*4H    ,  &
     4HEQUI,4HV   ,4H  kn,4HOS,k,4HGG/n,4HOSIM,4HP $ ,4H    ,4H    ,  &
     4H    ,8*4H      /
 DATA rdmap 3 /  &
     4HEQUI,4HV   ,4H  mn,4HOS,m,4HGG/n,4HOSIM,4HP $ ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HCOND,4H    ,4H  lb,4HSTP,,4HNOSI,4HMP $,12*4H    ,  &
     4HADD ,4H    ,4H  kg,4HGX,k,4HNOS/,4HKGG/,4H(1.0,4H,0.0,4H)/(1,  &
     4H.0,0,4H.0) ,4H$   ,6*4H    ,  &
     4HADD ,4H    ,4H  mg,4HG,mn,4HOS/m,4HGGX/,4H(1.0,4H,0.0,4H)/(1,  &
     4H.0,0,4H.0) ,4H$   ,6*4H    ,  &
     4HEQUI,4HV   ,4H  mg,4HGX,m,4HGG/a,4HLWAY,4HS $ ,4H    ,4H    ,  &
     4H    ,8*4H    , 4HLABE,4HL   ,4H  lb,4HSTP ,4H$   ,13*4H    ,  &
     4HCHKP,4HNT  ,4H  mg,4HG $ ,14*4H    ,  &
     4HALTE,4HR   ,4H  (g,4HP4) ,4H$   ,13*4H    ,  &
     4HCOND,4H    ,4H  lb,4HSEND,4H,dry,4H $  ,12*4H     /
 DATA rdmap 4 / 4HALTE,4HR   ,4H  (s,4HDR2),4H $  ,13*4H            /
 DATA oct 1   / 18    ,         0    ,         1  ,  &
     19    ,         0    ,         2  , 21    ,         0    ,         1  ,  &
     22    ,         0    ,         2  , 23    ,         0    ,         2  /
 DATA ptbs 1  / 1  , 11  , 11  ,  5  ,     1  ,         0  ,  0  ,  &
     4  , 43  , 45  ,  8  ,4HNAME  ,         0  ,  0  ,  &
     9  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     10  , 11  , 11  ,  6  ,     2  ,         0  ,  0  ,  &
     11  , 11  , 11  ,  6  ,     6  ,         0  ,  0  ,  &
     12  , 50  , 50  ,  0  ,4HSOL   ,         0  ,  0  ,  &
     13  , 11  , 11  ,  6  ,     7  ,         0  ,  0  ,  &
     14  , 50  , 50  ,  0  ,4HMSKP  ,         0  ,  0  ,  &
     15  , 11  , 11  ,  6  ,     3  ,         0  ,  0  ,  &
     17  , 12  , 13  ,  3  ,4HNANO  ,         1  , -1  ,  &
     17  , 17  , 18  ,  3  ,4HNANO  ,         2  , -1  ,  &
     17  , 28  , 30  ,  8  ,4HNAME  ,         0  ,  0  ,  &
     18  , 11  , 12  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     19  , 11  , 12  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     20  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     21  , 16  , 17  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     22  , 15  , 16  ,  3  ,4HNANO  ,         0  ,  0  ,  &
     24  , 13  , 13  ,  3  ,4HSTEP  ,         0  ,  0  /
 DATA ptbs 2  / 26  , 11  , 11  ,  5  ,     4  ,         0  ,  0  ,  &
     28  , 11  , 11  ,  6  ,     5  ,         0  ,  0  /
 DATA subnam  / 4HASCM,2H05  /
 
!     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION))
 
 rdmap(11, 4) = khrfn1(rdmap(11, 4),3,slash,1)
 rdmap(10,17) = khrfn1(rdmap(10,17),3,slash,1)
 rdmap(12,17) = khrfn1(rdmap(12,17),2,slash,1)
 
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
 DO  i = 1,4
   k = k + 1
   idat(k) = ipas25(i)
 END DO
 DO  i = 9,14
   k = k + 1
   idat(k) = ipas25(i)
 END DO
 DO  i = 5,8
   k = k + 1
   idat(k) = ipas25(i)
 END DO
 GO TO 200
 100 CONTINUE
 
 200 RETURN
 
!     INPUT ERROR
 
 1000 CALL mesage (7,0,subnam)
 nogo = 1
 RETURN
 
END SUBROUTINE ascm05
