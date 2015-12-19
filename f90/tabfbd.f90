BLOCK DATA tabfbd
!TABFBD
! TABFTX - BLOCK DATA PROGRAM FOR MODULE TABPRT
 
 INTEGER :: hx, re
 
 INTEGER :: hx01(32)
 INTEGER :: hx02(32)
 INTEGER :: hx03(32)
 INTEGER :: hx04(32)
 INTEGER :: hx05(32)
 INTEGER :: hx06(32)
 INTEGER :: hx07(32)
 INTEGER :: hx08(32)
 INTEGER :: hx09(32)
 INTEGER :: hx10(32)
 INTEGER :: hx11(32)
 INTEGER :: hx12(32)
 INTEGER :: hx13(32)
 INTEGER :: hx14(32)
 INTEGER :: hx15(32)
 INTEGER :: hx16(32)
 INTEGER :: hx17(32)
 INTEGER :: hx18(32)
 INTEGER :: hx19(32)
 INTEGER :: hx20(32)
 INTEGER :: hx21(32)
 INTEGER :: hx22(32)
 
 COMMON /tabftx/ la,na(2,21)  ,  hx(32,40)  , re(21)
 
 EQUIVALENCE (hx01(1),hx(1, 1))
 EQUIVALENCE (hx02(1),hx(1, 2))
 EQUIVALENCE (hx03(1),hx(1, 3))
 EQUIVALENCE (hx04(1),hx(1, 4))
 EQUIVALENCE (hx05(1),hx(1, 5))
 EQUIVALENCE (hx06(1),hx(1, 6))
 EQUIVALENCE (hx07(1),hx(1, 7))
 EQUIVALENCE (hx08(1),hx(1, 8))
 EQUIVALENCE (hx09(1),hx(1, 9))
 EQUIVALENCE (hx10(1),hx(1,10))
 EQUIVALENCE (hx11(1),hx(1,11))
 EQUIVALENCE (hx12(1),hx(1,12))
 EQUIVALENCE (hx13(1),hx(1,13))
 EQUIVALENCE (hx14(1),hx(1,14))
 EQUIVALENCE (hx15(1),hx(1,15))
 EQUIVALENCE (hx16(1),hx(1,16))
 EQUIVALENCE (hx17(1),hx(1,17))
 EQUIVALENCE (hx18(1),hx(1,18))
 EQUIVALENCE (hx19(1),hx(1,19))
 EQUIVALENCE (hx20(1),hx(1,20))
 EQUIVALENCE (hx21(1),hx(1,21))
 EQUIVALENCE (hx22(1),hx(1,22))
 
!-----------------------------------------------------------------------
 
 DATA la / 9 /
 DATA re /1,1,1,1,1,1,1,0,1,1 ,1,1,1,1,1,1,1,1,1,1  &
     ,1                    /
 
 DATA na / 4HBGPD,4HT     ,  4HGPL ,4H      ,  4HCSTM,4H    &
     , 4HGPLD,4H      ,  4HEQEX,4HIN    ,  4HEQDY,4HN     &
     , 4HGPDT,4H      ,  4HGPTT,4H      ,  4HGPCT,4H      &
     , 4H*10*,4H****  ,  4H*11*,4H****  ,  4H*12*,4H****  &
     , 4H*13*,4H****  ,  4H*14*,4H****  ,  4H*15*,4H****  &
     , 4H*16*,4H****  ,  4H*17*,4H****  ,  4H*18*,4H****  &
     , 4H*19*,4H****  ,  4H*20*,4H****  ,  4H*21*,4H**** /
 
 
 DATA hx01/4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx02/4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4HFORM,4HATTE,4HD li,4HST o,4HF ta,4HBLE ,4HDATA  &
     ,4H blo,4HCK  ,4H****,4H****,4H    ,4H( re,4HCORD,4H****  &
     ,4H  ) ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx03/4H    ,4H    ,4H    ,4H    ,4HINTE,4HRNAL,4H    ,4H coo  &
     ,4HRDIN,4HATE ,4H    ,4H    ,4H    ,4H coo,4HRDIN,4HATES  &
     ,4H in ,4HBASI,4HC co,4HORDI,4HNATE,4H sys,4HTEM ,4H      &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx04/4H    ,4H    ,4H    ,4H    ,4H   i,4HD   ,4H    ,4H sys  &
     ,4HTEM ,4HID  ,4H    ,4H    ,4H   x,4H    ,4H    ,4H    &
     ,4H    ,4H   y,4H    ,4H    ,4H    ,4H    ,4H   z,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx05/4H    ,4H  in,4HTERN,4HAL  ,4H    ,4H ext,4HERNA,4HL gr  &
     ,4HID  ,4H    ,4H ext,4HERNA,4HL gr,4HID  ,4H    ,4H ext  &
     ,4HERNA,4HL gr,4HID  ,4H    ,4H ext,4HERNA,4HL gr,4HID  &
     ,4H    ,4H ext,4HERNA,4HL gr,4HID  ,4H    ,4H    ,4H    /
 
 DATA hx06/4H    ,4H    ,4H id ,4H    ,4H    ,4H OR ,4HSCAL,4HAR i  &
     ,4HD   ,4H    ,4H OR ,4HSCAL,4HAR i,4HD   ,4H    ,4H OR  &
     ,4HSCAL,4HAR i,4HD   ,4H    ,4H OR ,4HSCAL,4HAR i,4HD    &
     ,4H    ,4H OR ,4HSCAL,4HAR i,4HD   ,4H    ,4H    ,4H    /
 
 DATA hx07/4H    ,4H  in,4HTERN,4HAL  ,4H    ,4H   e,4HXTER,4HNAL  &
     ,4HGRID,4H   s,4HEQUE,4HNCE ,4H    ,4H    ,4HEXTE,4HRNAL  &
     ,4H gri,4HD   ,4HSEQU,4HENCE,4H    ,4H    ,4H ext,4HERNA  &
     ,4HL gr,4HID  ,4H seq,4HUENC,4HE   ,4H    ,4H    ,4H    /
 
 DATA hx08/4H    ,4H    ,4H id ,4H    ,4H    ,4H   o,4HR sc,4HALAR  &
     ,4H id ,4H    ,4HNUMB,4HER  ,4H    ,4H    ,4HOR s,4HCALA  &
     ,4HR id,4H    ,4H num,4HBER ,4H    ,4H    ,4H OR ,4HSCAL  &
     ,4HAR i,4HD   ,4H  nu,4HMBER,4H    ,4H    ,4H    ,4H    /
 
 DATA hx09/4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx10/4H    ,4H    ,4H    ,4H  n ,4H    ,4H   i,4HD   ,4H    &
     ,4HTYPE,4H    ,4H    ,4H r(i,4H,1) ,4H    ,4H    ,4H    &
     ,4H r(i,4H,2) ,4H    ,4H    ,4H    ,4H r(i,4H,3) ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4HT(i),4H    ,4H    /
 
 DATA hx11/4H   e,4HXTER,4HNAL ,4H    ,4HEXTE,4HRNAL,4H gri,4HD   &
     ,4HINTE,4HRNAL,4H    ,4H ext,4HERNA,4HL gr,4HID  ,4H INT  &
     ,4HERNA,4HL   ,4H  ex,4HTERN,4HAL g,4HRID ,4H  in,4HTERN  &
     ,4HAL  ,4H   e,4HXTER,4HNAL ,4HGRID,4H   i,4HNTER,4HNAL /
 
 DATA hx12/4H   s,4HORT ,4HID  ,4H    ,4HOR s,4HCALA,4HR id,4H    &
     ,4H num,4HBER ,4H    ,4H OR ,4HSCAL,4HAR i,4HD   ,4H  nu  &
     ,4HMBER,4H    ,4H  OR,4H sca,4HLAR ,4HID  ,4H   n,4HUMBE  &
     ,4HR   ,4H   o,4HR sc,4HALAR,4H id ,4H    ,4HNUMB,4HER  /
 
 DATA hx13/4H   i,4HNTER,4HNAL ,4H    ,4H    ,4HCOOR,4HDINA,4HTE  &
     ,4H    ,4H    ,4H coo,4HRDIN,4HATES,4H in ,4HDEFI,4HNING  &
     ,4H coo,4HRDIN,4HATE ,4HSYST,4HEM  ,4H    ,4H  di,4HSPLA  &
     ,4HCEME,4HNT c,4HOOR-,4H    ,4HCONS,4HTRAI,4HNT  ,4H    /
 
 DATA hx14/4H    ,4H  id,4H    ,4H    ,4H    ,4HSYST,4HEM  ,4H    &
     ,4H    ,4H    ,4H  x ,4H    ,4H    ,4H    ,4H    ,4HY     &
     ,4H    ,4H    ,4H    ,4H z  ,4H    ,4H    ,4H  di,4HNATE  &
     ,4H sys,4HTEM ,4HID  ,4H    ,4H   c,4HODE ,4H    ,4H    /
 
 DATA hx15/4H   i,4HTERN,4HAL  ,4H    ,4H   t,4HEMPE,4HRATU,4HRE  &
     ,4H    ,4H    ,4HDEFA,4HULT ,4HTEMP,4HERAT,4HURE ,4H    &
     ,4H    ,4H   r,4HECOR,4HD nu,4HMBER,4H for,4H    ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx16/4H    ,4HINDE,4HX   ,4H    ,4H    ,4H set,4H id ,4H    &
     ,4H    ,4H    ,4H    ,4H  OR,4H fla,4HG   ,4H    ,4H    &
     ,4H    ,4H add,4HITIO,4HNAL ,4HTEMP,4H. da,4HTA  ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx17/4H  su,4HBSEQ,4HUENT,4H REC,4HORDS,4H of ,4H g p,4H t t  &
     ,4H  te,4HMPER,4HATUR,4HE da,4HTA a,4HRE l,4HISTE,4HD un  &
     ,4HDER ,4HSET ,4HID a,4HND e,4HLEME,4HNT t,4HYPE ,4HBY e  &
     ,4HLEME,4HNT i,4HD   ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx18/4H   r,4HECOR,4HD nu,4HMBER,4H   t,4HEMPE,4HRATU,4HRE s  &
     ,4HET i,4HD   ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx19/4H    ,4H    ,4H    ,4H  f ,4HO r ,4HM a ,4HT t ,4HE d  &
     ,4H  l ,4HI s ,4HT   ,4HO f ,4H  t ,4HA b ,4HL e ,4H  d  &
     ,4HA t ,4HA   ,4HB l ,4HO c ,4HK   ,4HG p ,4HC t ,4H    &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx20/4H  re,4HCORD,4H    ,4HPIVO,4HT  c,4HONNE,4HCTIN,4HG   &
     ,4H    ,4H    ,4H    ,4H    ,4H sor,4HTED ,4HLIST,4H of   &
     ,4H s i,4H l  ,4HNUMB,4HERS ,4HOF c,4HONNE,4HCTED,4H poi  &
     ,4HNTS ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    ,4H    /
 
 DATA hx21/4H  nu,4HMBER,4H    ,4HS i ,4HL   ,4H num,4HBER ,4H    &
     ,4H( 1 ,4H)   ,4H  ( ,4H2 ) ,4H    ,4H( 3 ,4H)   ,4H  (  &
     ,4H4 ) ,4H    ,4H( 5 ,4H)   ,4H  ( ,4H6 ) ,4H    ,4H( 7  &
     ,4H)   ,4H  ( ,4H8 ) ,4H    ,4H( 9 ,4H)   ,4H ( 1,4H0 ) /
 
 DATA hx22/4H   s,4HORT ,4HID  ,4H    ,4HOR s,4HCALA,4HR id,4H   c  &
     ,4HODED,4H sil,4H    ,4H OR ,4HSCAL,4HAR i,4HD   ,4HCODE  &
     ,4HD si,4HL   ,4H  OR,4H sca,4HLAR ,4HID  ,4H cod,4HED s  &
     ,4HIL  ,4H   o,4HR sc,4HALAR,4H id ,4H  co,4HDED ,4HSIL /
 
END BLOCK DATA
