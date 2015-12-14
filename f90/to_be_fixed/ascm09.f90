SUBROUTINE ascm09 (NAME,iphase,isol,nogo)
     
!     MREDUCE COMMAND DMAP DATA
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,1),xtra(13),subnam(2),isave(30),  &
     rdmap(18,25),rdmap1(18,9),rdmap2(18,9),  &
     rdmap3(18,7),oct(3,16),oct1(3,16),ptbs(7,53),  &
     ptbs1(7,18),ptbs2(7,18),ptbs3(7,17)
 COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(882)
 EQUIVALENCE    (rdmap1(1,1),rdmap(1, 1)),(oct1(1,1),oct(1,1)),  &
     (rdmap2(1,1),rdmap(1,10)),(ptbs1(1,1),ptbs(1,1)),  &
     (rdmap3(1,1),rdmap(1,19)),(ptbs2(1,1),ptbs(1,19)), (ptbs3(1,1),ptbs(1,37))
 DATA comnd   / 4HMRED    , 25    , 13    , 16    , 53    ,  0  /
 DATA slash   /  1H/       /
 DATA isave   /  &
     1,15,1,  2,11,2,  4,12,1,  4,16,3,  5, 5,1, 19, 7,3, 19, 8,3,  &
     19, 9,3, 22,15,2, 24, 6,2    /
 DATA rdmap 1 /  &
     4HMRED,4H1   ,4H  ca,4HSECC,4H,geo,4HM4,d,4HYNAM,4HICS,,4HCSTM,  &
     4H/use,4HTR,e,4HEDR,,4HEQST,4H,dmr,4H!*NA,4HMEA ,4H  */,4H    ,  &
     4H    ,4H    ,4H  s,,4HN,dr,4HY/st,4HP/s,,4HN,no,4HFIX/,4HS,n,,  &
     4HSKIP,4HM!*R,4HEAL*,4H $  , 5*4H    ,  &
     4HCOND,4H    ,4H  lb,4HM3ST,4HP,dr,4HY $ ,12*4H    ,  &
     4HSOFI,4H    ,4H  /k,4HNOA,,4HMNOA,4H,pno,4HA,bn,4HOA,k,4H4NOA,  &
     4H/s,n,4H,dry,4H!*NA,4HMEA ,4H  */,4H*KMT,4HX*!*,4HMMTX,4H*/  ,  &
     4H    ,4H    ,4H  *p,4HVEC*,4H!*BM,4HTX*/,4H*K4M,4HX* $,4H    ,  &
     4H    , 8*4H    ,  &
     4HCOND,4H    ,4H  lb,4HM2ST,4HP,sk,4HIPM ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  kn,4HOA,k,4HFFX/,4HNOFI,4HX $ ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  mn,4HOA,m,4HFFX/,4HNOFI,4HX $ ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  bn,4HOA,b,4HFFX/,4HNOFI,4HX $ ,4H    ,4H    ,  &
     4H    , 8*4H        /
 DATA rdmap 2 /  &
     4HEQUI,4HV   ,4H  k4,4HNOA,,4HK4FF,4HX/no,4HFIX ,4H$   ,4H    ,  &
     4H    , 8*4H    ,  &
     4HCOND,4H    ,4H  lb,4HM1ST,4HP,no,4HFIX ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HSCE1,4H    ,4H  us,4HETR,,4HKNOA,4H,mno,4HA,bn,4HOA,k,4H4NOA,  &
     4H/kff,4HX,kf,4HSX,k,4HSSX,,4HMFFX,4H,bff,4HX,k4,4HFFX ,4H$   ,  &
     4HLABE,4HL   ,4H  lb,4HM1ST,4HP $ ,13*4H    ,  &
     4HREAD,4H    ,4H  kf,4HFX,m,4HFFX,,4HBFFX,4H,k4f,4HFX,e,4HEDR,,  &
     4HUSET,4HR,/l,4HAMAR,4H,phi,4HR,mi,4HR,oe,4HIGR/,4H*MOD,4HES*/,  &
     4H    ,4H    ,4H  NE,4HIGVS,4H $  ,13*4H    ,  &
     4HOFP ,4H    ,4H  la,4HMAR,,4HOEIG,4HR,,,,4H,// ,4H$   ,4H    ,  &
     4H    , 8*4H    ,  &
     4HEQUI,4HV   ,4H  ph,4HIR,p,4HHIS/,4HNOFI,4HX $ ,4H    ,4H    ,  &
     4H    , 8*4H    ,  &
     4HCOND,4H    ,4H  lb,4HM2ST,4HP,no,4HFIX ,4H$   ,4H    ,4H    ,  &
     4H    , 8*4H        /
 DATA rdmap 3 /  &
     4HUMER,4HGE  ,4H  us,4HETR,,4HPHIR,4H,/ph,4HIS!*,4HN*!*,4HF*!*,  &
     4HS* $, 8*4H    , 4HLABE,4HL   ,4H  lb,4HM2ST,4HP $ ,13*4H    ,  &
     4HMRED,4H2   ,4H  ca,4HSECC,4H,lam,4HAR,p,4HHIS,,4HEQST,4H,use,  &
     4HTR,k,4HNOA,,4HMNOA,4H,bno,4HA,k4,4HNOA,,4HPNOA,4H,dmr,4H,   ,  &
     4H    ,4H    ,4H  qs,4HM/kn,4HOB,m,4HNOB,,4HBNOB,4H,k4n,4HOB,p,  &
     4HNOB,,4HPONO,4HB/st,4HP/s,,4HN,dr,4HY!*P,4HVEC*,4H $  ,4H    ,  &
     4HLABE,4HL   ,4H  lb,4HM3ST,4HP $ ,13*4H    ,  &
     4HLODA,4HPP  ,4H  pn,4HOB,p,4HONOB,4H/!*N,4HAMEB,4H   *,4H/S,N,  &
     4H,dry,4H $  , 7*4H    ,  &
     4HCOND,4H    ,4H  fi,4HNIS,,4HDRY ,4H$   ,12*4H        /
 DATA xtra   / 4HNAME,4HBOUN,4HFIXE,4HMETH,4HRANG,4HNMAX,4HOLDM ,  &
     4HOLDB,4HUSER,4HOUTP,4HRGRI,4HRNAM,4HRSAV        /
 DATA oct 1  / 6    ,         8    ,         0  ,  &
     7    ,         8    ,         1  , 8    ,         8    ,         2  ,  &
     9    ,         8    ,        16  , 10    ,         8    ,        32  ,  &
     11    ,         8    ,         0  , 12    ,         8    ,         0  ,  &
     13    ,         8    ,         0  , 14    ,         8    ,         0  ,  &
     15    ,         8    ,         0  , 16    ,         8    ,         0  ,  &
     17    ,         8    ,         0  , 18    ,         8    ,         0  ,  &
     19    ,         8    ,         0  , 20    ,         8    ,         0  ,  &
     24    ,         0    ,         8  /
 DATA ptbs 1  / 1  , 59  , 59  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     2  , 19  , 19  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     3  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     4  , 12  , 13  ,  3  ,4HNONA  ,         1  , -1  ,  &
     4  , 17  , 18  ,  3  ,4HNONA  ,         2  , -1  ,  &
     4  , 22  , 23  ,  3  ,4HNONA  ,        12  , -1  ,  &
     4  , 27  , 28  ,  3  ,4HNONA  ,        16  , -1  ,  &
     4  , 32  , 34  ,  3  ,4HNONA  ,        32  , -1  ,  &
     4  , 47  , 47  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     5  , 12  , 12  ,  4  ,4HPITM  ,         0  ,  0  ,  &
     6  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     7  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     8  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     9  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     10  , 11  , 13  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     11  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     12  , 17  , 18  ,  3  ,4HNONA  ,         1  ,  0  ,  &
     12  , 22  , 23  ,  3  ,4HNONA  ,         2  ,  0  /
 DATA ptbs 2  / 12  , 27  , 28  ,  3  ,4HNONA  ,        16  ,  0  ,  &
     12  , 32  , 34  ,  3  ,4HNONA  ,        32  ,  0  ,  &
     12  , 38  , 42  ,  0  ,4HNAMA  ,         1  ,  0  ,  &
     12  , 43  , 47  ,  0  ,4HNAMA  ,         1  ,  0  ,  &
     12  , 48  , 52  ,  0  ,4HNAMA  ,         1  ,  0  ,  &
     12  , 53  , 57  ,  0  ,4HNAMA  ,         2  ,  0  ,  &
     12  , 58  , 62  ,  0  ,4HNAMA  ,        16  ,  0  ,  &
     12  , 63  , 68  ,  0  ,4HNAMA  ,        32  ,  0  ,  &
     13  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     14  , 11  , 15  ,  0  ,4HNAMA  ,         1  ,  0  ,  &
     14  , 16  , 20  ,  0  ,4HNAMA  ,         2  ,  0  ,  &
     14  , 21  , 25  ,  0  ,4HNAMA  ,        16  ,  0  ,  &
     14  , 26  , 31  ,  0  ,4HNAMA  ,        32  ,  0  ,  &
     18  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     20  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     21  , 18  , 23  ,  0  ,4HNAMA  ,        55  ,  0  ,  &
     21  , 24  , 28  ,  0  ,4HNAMA  ,        55  ,  0  ,  &
     21  , 40  , 41  ,  3  ,4HNONA  ,         1  ,  0  /
 DATA ptbs 3  / 21  , 45  , 46  ,  3  ,4HNONA  ,         2  ,  0  ,  &
     21  , 50  , 51  ,  3  ,4HNONA  ,        16  ,  0  ,  &
     21  , 55  , 57  ,  3  ,4HNONA  ,        32  ,  0  ,  &
     21  , 61  , 62  ,  3  ,4HNONA  ,        12  ,  0  ,  &
     22  , 11  , 14  ,  0  ,4HNAMA  ,    131072  ,  0  ,  &
     22  , 15  , 16  ,  3  ,4HNONB  ,         1  , -1  ,  &
     22  , 20  , 21  ,  3  ,4HNONB  ,         2  , -1  ,  &
     22  , 25  , 26  ,  3  ,4HNONB  ,        16  , -1  ,  &
     22  , 30  , 32  ,  3  ,4HNONB  ,        32  , -1  ,  &
     22  , 36  , 37  ,  3  ,4HNONB  ,        12  , -1  ,  &
     22  , 41  , 43  ,  3  ,4HNONB  ,        12  , -1  ,  &
     22  , 47  , 47  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     22  , 60  , 60  ,  4  ,4HPITM  ,        12  ,  0  ,  &
     23  , 15  , 15  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     24  , 11  , 12  ,  3  ,4HNONB  ,         0  ,  0  ,  &
     24  , 16  , 18  ,  3  ,4HNONB  ,         0  ,  0  ,  &
     24  , 24  , 24  ,  8  ,4HNAMB  ,         0  ,  0  /
 DATA subnam  / 4HASCM,2H09  /
 
!     RESTORE TO ORIGINAL DATA BY REPLACEING ! BY / IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION))
 
 DO  l = 1,30,3
   i = isave(l+1)
   j = isave(l  )
   k = isave(l+2)
   rdmap(i,j) = khrfn1(rdmap(i,j),k,slash,1)
 END DO
 
!     VALIDATE COMMAND AND SET POINTERS
 
 IF (NAME /= comnd(1,1)) GO TO 70
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
 
END SUBROUTINE ascm09
