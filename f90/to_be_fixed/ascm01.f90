SUBROUTINE ascm01 (NAME,iphase,isol,nogo)
     
!     SUBSTRUCTURE COMMAND DMAP DATA
 
!     COMMENTS FROM G.CHAN/UNISYS   8/1991
!     IN SOME UNIX MACHINES, SUCH AS SiliconGraphics, THE FORTRAN
!     COMPILER IS A SUBSET OF THE C COMPILER. THE SYMBOL /* IS A COMMENT
!     MARKER FOR C, AND ANYTHING AFTER /* IS NOT PASS OVER TO THE
!     FORTRAN COMPILER. THEREFORE, ALL /* SYMBOLS IN RDMAP ARRAY ARE
!     REPLACED BY
!     THE ! WILL BE CHANGED BACK TO / IN THE EXECUTABLE CODE.
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME
 INTEGER, INTENT(IN)                      :: iphase
 INTEGER, INTENT(IN OUT)                  :: isol
 INTEGER, INTENT(OUT)                     :: nogo
 INTEGER :: comnd(6,3),xtra(3),subnam(2),isave(21),  &
     rdmap(18,29),rdmap1(18,9),rdmap2(18,9),  &
     rdmap3(18,9),rdmap4(18,2),oct(3,13),oct1(3,13), ptbs(7,16),ptbs1(7,16)
 COMMON /phas11/ ipas11(8)
 COMMON /phas31/ ipas31(2)
 COMMON /asdbd / irdm,nrdm,ixtra,nxtra,ioct,noct,iptbs,nptbs,  &
     iph,nph,idat(684)
 EQUIVALENCE     (rdmap1(1,1),rdmap(1, 1)),(oct1(1,1),oct(1,1)),  &
     (rdmap2(1,1),rdmap(1,10)),(ptbs1(1,1),ptbs(1,1)),  &
     (rdmap3(1,1),rdmap(1,19)), (rdmap4(1,1),rdmap(1,28))
 DATA comnd    / 4HSUBS    , 29    ,  3    , 13    , 16    ,  8  ,  &
     4HSUBS    ,  8    ,  1    ,  0    ,  3    ,  0  ,  &
     4HSUBS    ,  8    ,  1    ,  0    ,  3    ,  2  /
 DATA slash    / 1H/       /
 DATA isave    / 3,13,1, 19,8,2, 26,13,3, 26,15,2, 26,17,1, 27,5,1, 28,4,3  /
 DATA rdmap 1  / 4HALTE,4HR   ,4H  (b,4HEGIN,4H) $ ,13*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*nop,4H*/al,4HLWAY,4HS=-1,4H $  ,4H    ,  &
     4H    ,8*4H    ,  &
     4HSGEN,4H    ,4H  ca,4HSECC,4H,,,/,4HCASE,4HSS,c,4HASEI,4H,,,,,  &
     4H,,,,,4H/s,n,4H,dry,4H!*XX,4HXXXX,4HXX*/,4HS,N,,4HLUSE,4HT/  ,  &
     4H    ,4H    ,4H  s,,4HN,no,4HGPDT,4H $  ,12*4H    ,  &
     4HEQUI,4HV   ,4H  ca,4HSEI,,4HCASE,4HCC/a,4HLLWA,4HYS $,4H    ,  &
     4H    ,8*4H    , 4HALTE,4HR   ,4H  (a,4HFTGP,4H4) $,13*4H    ,  &
     4HPARA,4HM   ,4H  //,4H*add,4H*/dr,4HY/-1,4H /0 ,4H$   ,4H    ,  &
     4H    ,8*4H    , 4HLABE,4HL   ,4H  lb,4HSBEG,4H $  ,13*4H    ,  &
     4HCOND,4H    ,4H  lb,4HLIS,,4HDRY ,4H$   ,12*4H    /
 DATA rdmap 2  /  &
     4HSSG1,4H    ,4H  sl,4HT,bg,4HPDT,,4HCSTM,4H,sil,4H,est,4H,mpt,  &
     4H,gpt,4HT,ed,4HT,mg,4HG,ca,4HSECC,4H,dit,4H,/pg,4H,,,,,4H/   ,  &
     4H    ,4H    ,4H  lu,4HSET/,4HNSKI,4HP $ ,12*4H    ,  &
     4HCHKP,4HNT  ,4H  pg,4H $  ,14*4H    ,  &
     4HALTE,4HR   ,4H  (s,4HOLVE,4H) $ ,13*4H    ,  &
     4HSSG2,4H    ,4H  us,4HET,g,4HM,,k,4HFS,g,4HO,,p,4HG/qr,4H,po,,  &
     4HPS,p,4HL $ ,7*4H    , 4HCHKP,4HNT  ,4H  po,4H,ps,,4HPL $,13*4H    ,  &
     4HLABE,4HL   ,4H  lb,4HLIS ,4H$   ,13*4H    ,  &
     4HALTE,4HR   ,4H  (s,4HDR) ,4H$   ,13*4H    ,  &
     4HSUBP,4HH1  ,4H  ca,4HSECC,4H,eqe,4HXIN,,4HUSET,4H,bgp,4HDT,c,  &
     4HSTM,,4HGPSE,4HTS,e,4HLSET,4HS//s,4H,n,d,4HRY/ ,4H    ,4H    /
 DATA rdmap 3  /  &
     4H    ,4H    ,4H  *n,4HAME ,4H   *,4H/xpl,4HOTID,4H !*P,4HVEC*,  &
     4H $  ,8*4H    , 4HCOND,4H    ,4H  lb,4HSEND,4H,dry,4H $  ,12*4H    ,  &
     4HEQUI,4HV   ,4H  pg,4H,pl/,4HNOSE,4HT $ ,12*4H    ,  &
     4HCOND,4H    ,4H  lb,4HL10,,4HNOSE,4HT $ ,12*4H    ,  &
     4HSSG2,4H    ,4H  us,4HET,g,4HM,ys,4H,kfs,4H,GO,,4H,pg/,4HQR,p,  &
     4HO,ps,4H,pl ,4H$   ,6*4H    ,  &
     4HCHKP,4HNT  ,4H  po,4H,ps,,4HPL $,13*4H    ,  &
     4HLABE,4HL   ,4H  lb,4HL10 ,4H$   ,13*4H    ,  &
     4HSOFO,4H    ,4H  ,k,4HAA,m,4HAA,p,4HL,ba,4HA,k4,4HAA//,4HS,n,,  &
     4HDRY/,4H*xxx,4HXXXX,4HX*!*,4HKMTX,4H*!*M,4HMTX*,4H!*PV,4HEC*/,  &
     4H    ,4H    ,4H  *b,4HMTX*,4H!*K4,4HMX* ,4H$   ,4H    ,4H    ,  &
     4H    ,8*4H    /
 DATA rdmap 4  /  &
     4HLODA,4HPP  ,4H  pl,4H,/!*,4HNAME,4H    ,4H*/S,,4HN,DR,4HY $ ,  &
     4H    ,8*4H    ,  &
     4HEQUI,4HV   ,4H  ca,4HSESS,4H,cas,4HECC/,4HALWA,4HYS $,4H    ,  &
     4H    ,8*4H    /
 DATA xtra     / 4HSAVE,4HNAME,4HRUN      /
 DATA oct 1    / 9    ,    524288    ,         0     ,  &
     10    ,    983040    ,        12     ,  &
     11    ,    983040    ,        12     ,  &
     12    ,    983040    ,        12     ,  &
     14    ,    983040    ,        12     ,  &
     15    ,    983040    ,        12     ,  &
     16    ,    524288    ,         0     ,  &
     21    ,   1572864    ,        12     ,  &
     22    ,   1572864    ,        12     ,  &
     23    ,   1572864    ,        12     ,  &
     24    ,   1572864    ,        12     ,  &
     25    ,   1572864    ,        12     , 28    ,    524288    ,         8     /
 DATA ptbs 1   / 1  , 11  , 11  ,  7  ,     4  ,         0  ,  0  ,  &
     6  , 11  , 11  ,  8  ,     1  ,         0  ,  0  ,  &
     7  , 22  , 23  ,  3  ,4HRUN   ,         0  ,  0  ,  &
     13  , 11  , 11  ,  7  ,     2  ,         0  ,  0  ,  &
     17  , 11  , 11  ,  5  ,     3  ,         0  ,  0  ,  &
     19  , 11  , 12  ,  8  ,4HNAME  ,         0  ,  0  ,  &
     19  , 21  , 22  ,  8  ,4HSAVE  ,         0  ,  0  ,  &
     19  , 30  , 32  ,  4  ,4HPITM  ,    524300  ,  0  ,  &
     26  , 12  , 15  ,  0  ,4HNAME  ,         1  ,  0  ,  &
     26  , 16  , 19  ,  0  ,4HNAME  ,         2  ,  0  ,  &
     26  , 20  , 22  ,  0  ,4HNAME  ,        12  ,  0  ,  &
     26  , 23  , 26  ,  0  ,4HNAME  ,        16  ,  0  ,  &
     26  , 27  , 31  ,  0  ,4HNAME  ,        32  ,  0  ,  &
     26  , 40  , 42  ,  8  ,4HNAME  ,         0  ,  0  ,  &
     26  , 65  , 67  ,  4  ,4HPITM  ,    524288  ,  0  ,  &
     28  , 15  , 17  ,  8  ,4HNAME  ,         0  ,  0  /
 DATA subnam   / 4HASCM,2H01  /
 
!     RESTORE ORIGINAL DATA BY REPLACING ! BY / IN RDMAP
 
 DO  l = 1,21,3
   i = isave(l+1)
   j = isave(l  )
   k = isave(l+2)
   rdmap(i,j) = khrfn1(rdmap(i,j),k,slash,1)
 END DO
 
!     VALIDATE COMMAND AND SET POINTERS
 
 IF (NAME /= comnd(1,iphase)) GO TO 1000
 icomnd = iphase
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
 
!      MOVE XTRA DATA
 
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
 
!     MOVE PHASE 1 DATA
 
 IF (iphase /= 1 .OR. nph == 0) GO TO 80
 DO  i = 3,8
   k = k + 1
   idat(k) = ipas11(i)
 END DO
 DO  i = 1,2
   k = k + 1
   idat(k) = ipas11(i)
 END DO
 GO TO 200
 80 CONTINUE
 
!     MOVE PHASE 3 DATA
 
 IF (iphase /= 3 .OR. nph == 0) GO TO 200
 DO  i = 1,nph
   k = k + 1
   idat(k) = ipas31(i)
 END DO
 
 200 RETURN
 
!     INPUT ERROR
 
 1000 CALL mesage (7,0,subnam)
 nogo = 1
 RETURN
 
END SUBROUTINE ascm01
