SUBROUTINE ascm04 (NAME,iphase,isol,nogo)
     
!     REDUCE COMMAND DMAP DATA
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN OUT)                  :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,1),xtra(4),subnam(2),isave(54),  &
     rdmap(18,23),rdmap1(18,9),rdmap2(18,9),  &
     rdmap3(18,5),oct(3,16),oct1(3,16),ptbs(7,67),  &
     ptbs1(7,18),ptbs2(7,18),ptbs3(7,18),ptbs4(7,13)
 COMMON /asdbd/ irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs, iph,nph,idat(935)
 EQUIVALENCE    (rdmap1(1,1),rdmap(1,1)),(rdmap2(1,1),rdmap(1,10)),  &
     (rdmap3(1,1),rdmap(1,19)),(oct1(1,1),oct(1,1)),  &
     (ptbs 1(1,1),ptbs(1, 1)),(ptbs2(1,1),ptbs(1,19)),  &
     (ptbs 3(1,1),ptbs(1,37)),(ptbs4(1,1),ptbs(1,55))
 DATA comnd   / 4HREDU    , 23    ,  4    , 16    , 67    ,  0  /
 DATA slash   / 1H/       /
 DATA isave   /  &
     1,14,1,  3,12,1,  3,16,3,  4, 5,1,  8, 8,2,  8,11,1,  9, 7,1,  &
     10,10,3, 12, 7,1, 14, 7,1, 16, 7,2, 16,10,1, 19, 7,2, 19,10,1,  &
     20, 7,2, 20,10,1, 21, 7,1, 23, 6,2    /
 DATA rdmap 1 /  &
     4HREDU,4HCE  ,4H  ca,4HSECC,4H,geo,4HM4/p,4HVNOA,4H,uss,4HTP,i,  &
     4HNSTP,4H/stp,4H/s,n,4H,dry,4H!*PV,4HEC* ,4H$   ,4H    ,4H    ,  &
     4HCOND,4H    ,4H  lb,4HRSTP,4H,dry,4H $  ,12*4H    ,  &
     4HSOFI,4H    ,4H  /k,4HNOA,,4HMNOA,4H,pno,4HA,bn,4HOA,k,4H4NOA,  &
     4H/s,n,4H,dry,4H!*NA,4HME00,4H0A*/,4H*KMT,4HX*!*,4HMMTX,4H*/  ,  &
     4H    ,4H    ,4H  *p,4HVEC*,4H!*BM,4HTX*/,4H*K4M,4HX* $,4H    ,  &
     4H    , 8*4H    , 4HCOND,4H    ,4H  lb,4HRSTP,4H,dry,4H $  ,12*4H    ,  &
     4HSMP1,4H    ,4H  us,4HSTP,,4HKNOA,4H,,,/,4HGONO,4HA,kn,4HOB,k,  &
     4HONOA,4H,lon,4HOA,,,4H,,, ,4H$   , 4*4H    ,  &
     4HMERG,4HE   ,4H  GO,4HNOA,,4HINST,4HP,,,,4H,pvn,4HOA/g,4HNOA/,  &
     4H1/ty,4HP/2 ,4H$   , 6*4H    ,  &
     4HSOFO,4H    ,4H  ,g,4HNOA,,4HLONO,4HA,,,,4H//dr,4HY!*N,4HAME0,  &
     4H00A*,4H!*HO,4HRG*/,4H*LMT,4HX* $, 4*4H    ,  &
     4HSOFO,4H    ,4H  ,k,4HNOB,,4H,,,/,4H/dry,4H!*NA,4HME00,4H0B*/,  &
     4H*kmt,4HX* $, 7*4H       /
 DATA rdmap 2 /  &
     4HSOFI,4H    ,4H  /g,4HNOA,,4H,,,/,4HS,n,,4HDRY/,4H*nam,4HE000,  &
     4HA*!*,4HHORG,4H* $ , 6*4H    ,  &
     4HMPY3,4H    ,4H  gn,4HOA,m,4HNOA,,4H/mno,4HB/0/,4H0 $ ,4H    ,  &
     4H    , 8*4H    ,  &
     4HSOFO,4H    ,4H  ,m,4HNOB,,4H,,,/,4H/dry,4H!*NA,4HME00,4H0B*/,  &
     4H*mmt,4HX* $, 7*4H    ,  &
     4HMPY3,4H    ,4H  gn,4HOA,b,4HNOA,,4H/bno,4HB/0/,4H0 $ ,4H    ,  &
     4H    , 8*4H    ,  &
     4HSOFO,4H    ,4H  ,b,4HNOB,,4H,,,/,4H/dry,4H!*NA,4HME00,4H0B*/,  &
     4H*bmt,4HX* $, 7*4H    ,  &
     4HMPY3,4H    ,4H  gn,4HOA,k,4H4NOA,4H,/k4,4HNOB/,4H0/0 ,4H$   ,  &
     4H    , 8*4H    ,  &
     4HSOFO,4H    ,4H  ,k,4H4NOB,4H,,,,,4H//dr,4HY!*N,4HAME0,4H00B*,  &
     4H!*K4,4HMX* ,4H$   , 6*4H    ,  &
     4HPART,4HN   ,4H  pn,4HOA,,,4HPVNO,4HA/po,4HNOA,,4H,,/1,4H/1/2,  &
     4H $  , 8*4H    ,  &
     4HMPYA,4HD   ,4H  gn,4HOA,p,4HNOA,,4H/pno,4HB/1/,4H1/0/,4H1 $ ,  &
     4H    , 8*4H        /
 DATA rdmap 3 /  &
     4HSOFO,4H    ,4H  ,p,4HONOA,4H,,,,,4H//dr,4HY!*N,4HAME0,4H00A*,  &
     4H!*PO,4HVE* ,4H$   , 6*4H    ,  &
     4HSOFO,4H    ,4H  ,p,4HVNOA,4H,,,,,4H//dr,4HY!*N,4HAME0,4H00A*,  &
     4H!*UP,4HRT* ,4H$   , 6*4H    ,  &
     4HSOFO,4H    ,4H  ,p,4HNOB,,4H,,,/,4H/dry,4H!*NA,4HME00,4H0B*/,  &
     4H*pve,4HC* $, 7*4H    , 4HLABE,4HL   ,4H  lb,4HRSTP,4H $  ,13*4H    ,  &
     4HLODA,4HPP  ,4H  pn,4HOB,p,4HONOA,4H/!*N,4HAME0,4H00B*,4H/S,N,  &
     4H,dry,4H $  , 7*4H      /
 DATA xtra    / 4HOUTP,4HNAME,4HBOUN,4HRSAV        /
 DATA oct 1   / 6    ,         0    ,         1  ,  &
     7    ,         0    ,         1  , 8    ,         0    ,         1  ,  &
     9    ,         0    ,         1  , 10    ,         1    ,        62  ,  &
     11    ,         0    ,         2  , 12    ,         0    ,         2  ,  &
     13    ,         0    ,        16  , 14    ,         0    ,        16  ,  &
     15    ,         0    ,        32  , 16    ,         0    ,        32  ,  &
     17    ,         0    ,        12  , 18    ,         0    ,        12  ,  &
     19    ,         0    ,        12  , 21    ,         0    ,        12  ,  &
     23    ,         0    ,         8  /
 DATA ptbs 1  / 1  , 24  , 26  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     1  , 32  , 32  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     1  , 38  , 38  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     1  , 41  , 42  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     1  , 53  , 55  ,  4  ,4HPITM  ,         0  ,  0  ,  &
     2  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     3  , 12  , 13  ,  3  ,4HNONA  ,         1  , -1  ,  &
     3  , 17  , 18  ,  3  ,4HNONA  ,         2  , -1  ,  &
     3  , 22  , 23  ,  3  ,4HNONA  ,        12  , -1  ,  &
     3  , 27  , 28  ,  3  ,4HNONA  ,        16  , -1  ,  &
     3  , 32  , 34  ,  3  ,4HNONA  ,        32  , -1  ,  &
     3  , 45  , 47  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     4  , 11  , 12  ,  4  ,4HPITM  ,         0  ,  0  ,  &
     5  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     6  , 11  , 13  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     6  , 17  , 18  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     6  , 25  , 27  ,  3  ,4HNONA  ,         0  , -1  ,  &
     6  , 31  , 32  ,  3  ,4HNONB  ,         0  , -1  /
 DATA ptbs 2  / 6  , 36  , 38  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     6  , 42  , 44  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     7  , 11  , 13  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     7  , 17  , 19  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     7  , 26  , 28  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     7  , 32  , 33  ,  3  ,4HNONA  ,         0  , -1  ,  &
     7  , 38  , 39  ,  3  ,4HPREC  ,         0  ,  0  ,  &
     8  , 12  , 13  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     8  , 17  , 19  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     8  , 17  , 21  ,  0  ,4HRSAV  ,         0  ,  0  ,  &
     8  , 30  , 32  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     8  , 48  , 54  ,  0  ,4HRSAV  ,         0  ,  0  ,  &
     9  , 12  , 13  ,  3  ,4HNONB  ,         0  ,  0  ,  &
     9  , 25  , 27  ,  8  ,4HNAMB  ,         0  ,  0  ,  &
     10  , 12  , 13  ,  3  ,4HNONA  ,         0  , -1  ,  &
     10  , 28  , 30  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     11  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     11  , 16  , 17  ,  3  ,4HNONA  ,         0  ,  0  /
 DATA ptbs 3  / 11  , 22  , 23  ,  3  ,4HNONB  ,         0  , -1  ,  &
     12  , 12  , 13  ,  3  ,4HNONB  ,         0  ,  0  ,  &
     12  , 25  , 27  ,  8  ,4HNAMB  ,         0  ,  0  ,  &
     13  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     13  , 16  , 17  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     13  , 22  , 23  ,  3  ,4HNONB  ,         0  , -1  ,  &
     14  , 12  , 13  ,  3  ,4HNONB  ,         0  ,  0  ,  &
     14  , 25  , 27  ,  8  ,4HNAMB  ,         0  ,  0  ,  &
     15  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     15  , 16  , 18  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     15  , 23  , 25  ,  3  ,4HNONB  ,         0  , -1  ,  &
     16  , 12  , 14  ,  3  ,4HNONB  ,         0  ,  0  ,  &
     16  , 26  , 28  ,  8  ,4HNAMB  ,         0  ,  0  ,  &
     17  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     17  , 17  , 19  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     17  , 23  , 25  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     18  , 11  , 12  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     18  , 16  , 17  ,  3  ,4HNONA  ,         0  ,  0  /
 DATA ptbs 4  / 18  , 22  , 23  ,  3  ,4HNONB  ,         0  , -1  ,  &
     19  , 12  , 14  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     19  , 26  , 28  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     19  , 37  , 39  ,  4  ,4HPOIT  ,         0  ,  0  ,  &
     20  , 12  , 14  ,  3  ,4HNONA  ,         0  ,  0  ,  &
     20  , 26  , 28  ,  8  ,4HNAMA  ,         0  ,  0  ,  &
     21  , 12  , 13  ,  3  ,4HNONB  ,         0  ,  1  ,  &
     21  , 25  , 27  ,  8  ,4HNAMB  ,         0  ,  0  ,  &
     21  , 36  , 38  ,  4  ,4HPITM  ,         0  ,  0  ,  &
     22  , 14  , 14  ,  3  ,4HSTEP  ,         0  ,  0  ,  &
     23  , 11  , 12  ,  3  ,4HNONB  ,         0  ,  1  ,  &
     23  , 16  , 18  ,  3  ,4HNONA  ,         0  ,  1  ,  &
     23  , 22  , 24  ,  8  ,4HNAMB  ,         0  ,  0  /
 
 DATA subnam  / 4HASCM,2H04  /
 
!     RESTORE TO ORIGINAL DATA BY REPLACEING !* BY /* IN RDMAP ARRAY
!     (SEE ASCM01 FOR EXPLANATION)
 
 DO  l = 1,51,3
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
 
END SUBROUTINE ascm04
