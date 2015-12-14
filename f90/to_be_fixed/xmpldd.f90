SUBROUTINE xmpldd
     
 
!     MPL    = MODULE PROPERTIES TABLE
!     LMPL   = LENGTH OF MPL TABLE
!     MPLPNT = POINTER TO AN MPL ENTRY
 
!     DESCRIPTION OF VARIABLES EQUIVALENCED TO /XGPI2/ ENTRIES
!     EQUIVALENCE (LMPL,LORDNL),(MPLPNT,IORBOT),(MPL,IORDNL)
 
!     IORDNL = TABLE USED TO COMPUTE FILE ORDNALS AND NUT VALUES
!     LORDNL = LENGTH OF IORDNL TABLE
!     IORBOT = POINTER TO LAST ENTRY MADE IN IORDNL TABLE
 
!     ==================================================================
 
!     NOTE DATA ITEMS HAVE BLANK WORDS TO FACILITATE ADDITIONS
!     CHANGE DIMENSIONS ONLY IF THE BLANKS ARE DEPLETED BY ADDITIONS
 
!     ANY FOLLOWING DATA LINE ENDS WITH BCD BLANKS SHOULD BE FOLLOWED BY
!     COMMA, SO THAT A STRIPPING ROUTINE WOULD NOT STRIP OFF THOSE BLANK
 
!     ==================================================================
 
!                        LOAD /XGPI2/
!                   MODULE PROPTERIES LIST (MPL)
 
 REAL :: x(2,20)
 DOUBLE PRECISION :: xx(20)    , xxx(20)
 DIMENSION        mpl01( 68), mpl02(161), mpl03(135), mpl04(152),  &
     mpl05(138), mpl06(162), mpl07(200), mpl08(137),  &
     mpl09(173), mpl10( 93), mpl11(116), mpl12(135),  &
     mpl13(150), mpl14(151), mpl15(135), mpl16( 53),  &
     mpl17(144), mpl18(169), mpl19(193), mpl20(186),  &
     mpl21(196), mpl22(119), mpl( 3166)
 COMMON /xgpi2 /  lmpl, mplpnt    , imp( 3166)
 COMMON /xgpi2x/  xxx
 EQUIVALENCE      (xx(1),x(1,1))
 EQUIVALENCE      (mpl(   1),mpl01(1)) ,(mpl(  69),mpl02(1)) ,  &
     (mpl( 230),mpl03(1)) ,(mpl( 365),mpl04(1)) ,  &
     (mpl( 517),mpl05(1)) ,(mpl( 655),mpl06(1)) ,  &
     (mpl( 817),mpl07(1)) ,(mpl(1017),mpl08(1)) ,  &
     (mpl(1154),mpl09(1)) ,(mpl(1327),mpl10(1)) ,  &
     (mpl(1420),mpl11(1)) ,(mpl(1536),mpl12(1)) ,  &
     (mpl(1671),mpl13(1)) ,(mpl(1821),mpl14(1)) ,  &
     (mpl(1972),mpl15(1)) ,(mpl(2107),mpl16(1)) ,  &
     (mpl(2160),mpl17(1)) ,(mpl(2304),mpl18(1)) ,  &
     (mpl(2473),mpl19(1)) ,(mpl(2666),mpl20(1)) ,  &
     (mpl(2852),mpl21(1)) ,(mpl(3048),mpl22(1))
 
 DATA    lmplx          /    3166      /
 
 DATA    x(1,1)         /    -1.0      /
 DATA    xx(2)          /    -1.0D+0   /
 DATA    x(1,3),x(2,3)  /  2*-1.0      /
 DATA    xx(4),xx(5)    /  2*-1.0D+0   /
 DATA    x(1,6)         /     1.0      /
 DATA    x(1,7),x(2,7)  / 1.0,0.0      /
 DATA    xx(8)          /     0.0D+0   /
 DATA    x(1,9),x(2,9)  /     2*0.0    /
 DATA    x(1,10)        /     30.0     /
 DATA    x(1,11)        /     0.001    /
 DATA    x(1,12)        /     0.55     /
 DATA    x(1,13)        /     0.01     /
 DATA    x(1,14)        /     0.00001  /
 DATA    x(1,15)        /     1.01     /
 DATA    x(1,16)        /     0.80     /
 DATA    xx(17)         /     1.1D+37  /
 DATA    x(1,18),x(2,18)/   2*1.1E+37  /
 
 
 DATA    mpl01          / 4,   4HFILE,4H    , 5  &
     ,  4,   4HBEGI,4HN   , 5 ,  4,   4HCHKP,4HNT  , 4  &
     ,  4,   4HLABE,4HL   , 5 ,  4,   4HREPT,4H    , 3  &
     ,  4,   4HJUMP,4H    , 3 ,  4,   4HCOND,4H    , 3  &
     ,  4,   4HSAVE,4H    , 4 ,  4,   4HPURG,4HE   , 4  &
     ,  4,   4HEQUI,4HV   , 4 ,  4,   4HEND ,4H    , 3  &
     ,  4,   4HEXIT,4H    , 3 , 20,   19*0  &
     /
 
!      IN NEXT 12 LINES, '1' MAY MEAN NUMERIC ONE, OR A VERTICAL BAR
 
!      NO. OF WORDS        1  I  O  S   ---PARAMETERS
!      OF THIS DMAP           N  U  C   1  NEGATIVE FOR NO DEFAULT
!      LINE               OR  P  T  R   1  POSITIVE INDICATES DEFAULT TO
!        1                    U  P  A   1   1 = INTEGER    NEXT VALUE(S)
!        1   DAMP NAME     2  T  U  T   1   2 = RSP
!        1   1                1  T  C   1   3 = BCD    (NEXT WORD(S)
!        1   1     1. I/O DB  1  1  H   1   4 = RDP    AFTER 2,4,5,6 ARE
!        1   1     2. NO OUT- 1  1  1   1   5 = CSP    POINTER TO DEF-
!        1   1        PUT DB  1  1  1   1   6 = CDP    AULT VALUE(S) IN
!        1   1             1  1  1  1   1              X OR XX ARRAYS)
!        1   1             1  1  1  1   1   REF. PROG. MAN. SEC 2.4.2.2
 DATA  mpl02 / 25, 4HADD ,4H    , 1, 2, 1, 0,  5,  2*18, 5, 2*18, 6,  4*17  &
     ,  6,  4*17, 1, 0  &
     , 22, 4HADD5,4H    , 1, 5, 1, 0,  5,  7, 7, 5, 7, 7, 5,7,7, 5,7,7  &
     ,  5,  7, 7 , 10, 4HAMG ,4H    , 1, 2, 4, 5,  3* -1  &
     , 11, 4HAMP ,4H    , 1,10, 3,14,  2* -1, 1,-1  &
     , 12, 4HAPD ,4H    , 1, 8,12, 5,  3* -1, 2, 9  &
     , 10, 4HBMG ,4H    , 1, 4, 1, 1, -1, -1,-5  &
     , 12, 4HCASE,4H    , 1, 2, 1, 0, -3,  1, 1, 1,-1  &
     , 15, 4HCYCT,4H1   , 1, 1, 2, 3, -3, -3,-1,-1, 1, 1,1, 1  &
     , 16, 4HCYCT,4H2   , 1, 6, 5, 6, -3, -1,-1, 1,-1, 1,1, 1,1  &
     , 10, 4HCEAD,4H    , 1, 5, 4,12, -1,  1, 1  &
     , 11, 4HCURV,4H    , 1, 6, 2, 5,  1, -1, 1, 0 ,  7, 6*0  &
     /
 
 DATA  mpl03 / 10, 4HDDR ,4H    , 1, 1, 1, 0, -3,-3,-3  &
     ,  7, 4HDDR1,4H    , 1, 2, 1, 1  &
     , 14, 4HDDR2,4H    , 1, 9, 3, 6, -3, 1,-1, 1,-1,   1,-1  &
     ,  7, 4HDDRM,4HM   , 1,11, 5, 7  &
     , 21, 4HDECO,4HMP  , 1, 1, 2, 4,  1, 0, 1, 0, 4,   8, 8, 5, 9, 9  &
     ,  1, 0, 1, 0 , 12, 4HDIAG,4HONAL, 1, 1, 1, 0,  3,4HCOLU,4HMN  , 2, 6  &
     , 19, 4HDPD ,4H    , 1, 4,11, 4,  9*-1, 1, 1,-1  &
     , 16, 4HDSCH,4HK   , 1, 3, 0, 3,  2*-2, 7*-1  &
     ,  8, 4HDSMG,4H1   , 1,10, 1, 1, -1  &
     , 11, 4HDSMG,4H2   , 1,11, 7, 0,  1, 0,-1,-1 , 10, 9*0  &
     /
 
 DATA  mpl04 /  &
     33, 4HDUMM,4HOD1 , 1, 1, 2, 3,  1,-1, 1,-1, 1,-1, 1,-1, 2,1, 2,1  &
     ,  3, 4HABCD , 4HEFGH , 4, 2,2, 5,3 ,  3, 6, 4, 4, 5, 5  &
     , 33, 4HDUMM,4HOD2 , 1, 8, 8,10,  1,-1, 1,-1, 1,-1, 1,-1, 2,1, 2,1  &
     ,  3, 4HABCD , 4HEFGH,  4, 2,2, 5,3 ,  3, 6, 4, 4, 5, 5  &
     , 33, 4HDUMM,4HOD3 , 1, 8, 8,10,  1,-1, 1,-1, 1,-1, 1,-1, 2,1, 2,1  &
     ,  3, 4HABCD , 4HEFGH , 4, 2,2, 5,3 ,  3, 6, 4, 4, 5, 5  &
     , 33, 4HDUMM,4HOD4 , 1, 8, 8,10,  1,-1, 1,-1, 1,-1, 1,-1, 2,1, 2,1  &
     ,  3, 4HABCD , 4HEFGH,  4, 2,2, 5,3 ,  3, 6, 4, 4, 5, 5  &
     , 20, 19*0 /
 
 DATA  mpl05 / 11, 4HEMA1,4H    , 1, 5, 1, 2,  1,-1, 2, 6  &
     , 43, 4HEMG ,4H    , 1, 6, 7, 4,  1,-1, 1,-1, 1,-1, 1,-1, 1,-1  &
     ,  1,-1, 1,-1, 1,-1, 1,-1, 1,-1 ,  1,-1, 1,-1, 1,-1, 1,-1, 1,-1  &
     ,  1,-1, 2, 9, 2, 9 , 11, 4HFA1 ,4H    , 1, 6, 4, 6,  2*-1, 1, 0  &
     , 12, 4HFA2 ,4H    , 1, 3, 4, 0, -1,-2, 3, 4HYES   ,4H     ,  &
     15, 4HFBS ,4H    , 1, 3, 1, 1,  1, 0, 1, 1, 1, 0, 1, 0  &
     , 13, 4HFRLG,4H    , 1, 8, 5, 4, -3, 1,-1, 3, 4HFREQ,4H    ,  &
     17, 4HFRRD,4H    , 1,11, 4, 8, -3,-3,-1,-1,-1,-1,-1,-1 ,  1, 1  &
     , 16, 15*0 /
 
 DATA  mpl06 / 9, 4HGI  ,4H    , 1, 8, 1, 6,  2*-1  &
     , 24, 4HGKAD,4H    , 1,10, 8, 6, -3,-3,-3,-2,-2,-2, 11*-1  &
     , 21, 4HGKAM,4H    , 1, 9, 4, 4, -1,-1, 2, 9, 2, 1,  4*-1 , 1, 1, 1,-1  &
     , 11, 4HGP1 ,4H    , 1, 3, 6, 2, -1,-1, 1, 1  &
     ,  7, 4HGP2 ,4H    , 1, 2, 1, 4  &
     , 12, 4HGP3 ,4H    , 1, 3, 2, 2, -1, 1, 1, 1, 1  &
     , 24, 4HGP4 ,4H    , 1, 7, 5, 2,  9*-1, 1, 1, 1,-1, 1,0,  1,0  &
     , 10, 4HGPCY,4HC   , 1, 3, 1, 2, -3, 1, 1  &
     ,  8, 4HGPFD,4HR   , 1, 9, 2, 4, -3  &
     , 20, 4HDUMM,4HOD5 , 1, 5, 5, 0, -1, 1, 0, 1, 0, 1, 0,1,0,1,0,1,0  &
     , 11, 4HGPWG,4H    , 1, 4, 1, 4,  1,-1, 2, 6 ,  5, 4*0  &
     /
 
 DATA  mpl07 / 13, 4HINPU,4HT   , 1, 5, 5, 0,  1,-1, 1, 0, 1, 0  &
     , 17, 4HINPU,4HTT1 , 1, 0, 5, 0,  1, 0, 1, 0, 3, 4HXXXX,4HXXXX  &
     ,  3, 4H     , 4H       ,  &
     21, 4HINPU,4HTT2 , 1, 0, 5, 0,  1, 0, 1,14, 3, 4HXXXX,4HXXXX  &
     ,  1, 0, 1, 0, 3, 4H    ,4H    ,  &
     13, 4HINPU,4HTT3 , 1, 5, 5, 0,  1,-11, 1, 0, 1, 0  &
     , 16, 4HINPU,4HTT4 , 1, 0, 5, 0,  1, 1, 1,14, 3, 4HXXXX,4HXXXX ,  1, 0  &
     , 29, 4HMATG,4HEN  , 1, 1, 1, 0,  1, 0, 1, 0, 1, 0, 1, 0, 1, 0  &
     , 1,  0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0  &
     , 18, 4HMATG,4HPR  , 2, 4, 0, 0, -3, 3,  4H    , 4H     , 3  &
     ,  4HALL ,4H    , 2, 9, 1, 0  &
     , 19, 4HMATP,4HRN  , 2, 5, 0, 0,  1, 0, 1, 0, 1, 0, 1, 0, 1, 0 ,  1, 0  &
     , 11, 4HMATP,4HRT  , 2, 1, 0, 0,  1, 0, 1, 0  &
     ,  7, 4HMCE1,4H    , 1, 2, 1, 7 ,  7, 4HMCE2,4H    , 1, 6, 4, 6  &
     , 17, 4HMERG,4HE   , 1, 6, 1, 0,  1,-1, 1, 0, 1, 0, 1, 0, 1, 0 , 12, 11*0  /
 
 DATA  mpl08 / 20, 4HMODA,4H    , 1, 0, 4, 0,  5*-2, 5* -1, -2, -1, -1  &
     , 10, 4HMODA,4HCC  , 1, 6, 5, 0,  3,4HTRAN,4H        ,  &
     18, 4HMODB,4H    , 1, 3, 4, 0,  4*-2, 3* -1, -2,  3* -1  &
     ,  8, 4HMODC,4H    , 1, 2, 0, 0, -1  &
     , 14, 4HMPYA,4HD   , 1, 3, 1, 1, -1, 1,  1, 1,  1,  1,  0  &
     , 14, 4HMTRX,4HIN  , 1, 5, 3, 7, -1, 1, -1, 1, -1,  1, -1  &
     , 11, 4HOFP ,4H    , 2, 6, 0, 0,  1, 0,  1,-1  &
     , 10, 4HOPTP,4HR1  , 1, 5, 1, 1,  3*-1  &
     , 12, 4HOPTP,4HR2  , 1, 3, 2, 0,  3*-1,  1, 0 , 20, 19*0  &
     /
 
 DATA  mpl09 / 9, 4HOUTP,4HUT  , 2, 1, 0, 0,  1,-1  &
     , 14, 4HOUTP,4HUT1 , 2, 5, 0, 0,  1, 0, 1, 0,3,4HXXXX,4HXXXX  &
     , 21, 4HOUTP,4HUT2 , 2, 5, 0, 0,  1, 0, 1,14,3,4HXXXX,4HXXXX,  &
     1, 0, 1, 0,3,4H    ,4H    ,  &
     22, 4HOUTP,4HUT3 , 2, 5, 0, 0,  1, 0,-3,   3,3HXXX,1H ,3,3HXXX  &
     ,  1H  , 3,3HXXX,1H ,3,3HXXX,1H ,  &
     13, 4HOUTP,4HUT4 , 2, 5, 0, 0,  1,-1, 1,14,1, 1  &
     , 14, 4HPARA,4HM   , 1, 0, 0, 0, -3, 1, 1, 1,1, 1,1  &
     , 30, 4HPARA,4HML  , 2, 1, 0, 0, -3, 1, 1, 1,1, 2,9, 1,0,   4,8,8  &
     ,  3, 4H(voi ,4HD)  , 5,9,9, 6,4*8  &
     , 25, 4HPARA,4HMR  , 2, 0, 0, 0, -3, 2, 9, 2,9, 2,9, 5,9,9, 5,9,9  &
     ,  5, 9, 9, 1,0  &
     , 19, 4HPART,4HN   , 1, 3, 4, 0,  1,-1, 1, 0,1, 0,1, 0,1,0, 1,0 , 6, 5*0  &
     /
 
 DATA  mpl10 / 17, 4HMRED,4H1   , 1, 4, 4, 1, -3,-1,-1,-1, -1,-3, 1,0, 2,9  &
     , 14, 4HMRED,4H2   , 1,12, 6,11, -1,-1, 3,4H    ,4H     , 1,0  &
     , 12, 4HCMRE,4HD2  , 1,11, 6,11, -1,-1, 3,4H    ,4H     ,  &
     13, 4HPLA1,4H    , 1, 7, 4, 0,  5*-1,-5  &
     ,  8, 4HPLA2,4H    , 1, 3, 3, 0, -1 ,  9, 4HPLA3,4H    , 1, 6, 2, 1, -1,-1  &
     , 10, 4HPLA4,4H    , 1, 6, 2, 1, -1,-1,-5 , 10, 9*0  &
     /
 
 DATA  mpl11 / 14, 4HPLOT,4H    , 1,13, 1, 4,  3*-1, 1, 1, 1, 0  &
     , 10, 4HPLTS,4HET  , 1, 4, 4, 2, -1, 1,-1  &
     , 11, 4HPLTT,4HRAN , 1, 2, 2, 0,  1, 0, 1, 0  &
     ,  7, 4HPRTM,4HSG  , 2, 1, 0, 0  &
     , 13, 4HPRTP,4HARM , 2, 0, 0, 0, -1, 3, 4HXXXX ,4HXXXX, 1,0  &
     ,  9, 4HRAND,4HOM  , 1, 9, 2, 0,  1,-1 ,  7, 4HRBMG,4H1   , 1, 3, 6, 1  &
     , 11, 4HRBMG,4H2   , 1, 1, 1, 4,  1, 1, 2, 6  &
     ,  7, 4HRBMG,4H3   , 1, 3, 1, 2 ,  7, 4HRBMG,4H4   , 1, 4, 1, 3  &
     , 20, 19*0 /
 
 DATA  mpl12 / 13, 4HREAD,4H    , 1, 7, 4,10, -3,-1, 1, 1, 2,  6  &
     , 14, 4HRMG ,4H    , 1, 4, 3, 6,  2, 9, 2, 9, 1, -1, -1  &
     , 24, 4HSCAL,4HAR  , 2, 1, 0, 0,  1, 1, 1, 1, 2,  9,  4,8,8  &
     ,  5, 9, 9, 6, 4*8 ,  7, 4HSCE1,4H    , 1, 5, 6, 1  &
     ,  9, 4HSDR1,4H    , 1,11, 3, 6, -1,-3  &
     , 16, 4HSDR2,4H    , 1,16, 8, 3, -3, 1, 1, 1,-1,  1,-1, 1,1  &
     ,  7, 4HSDR3,4H    , 1, 6, 6, 8  &
     , 11, 4HSDRH,4HT   , 1,10, 1, 3,  2, 9, 1,-1  &
     , 23, 4HSEEM,4HAT  , 2, 5, 0, 0,  3,4HPRIN,4HT   , 1,0, 1,100  &
     ,  3,4HM   ,4H    , 1,1, 2,9, 2,9 , 11, 10*0  &
     /
 
 DATA  mpl13 /  &
     26, 4HSETV,4HAL  , 2, 0, 0, 0, -1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1  &
     ,  1,-1, 1,-1, 1,-1, 1,-1 , 11, 4HSMA1,4H    , 1, 5, 3, 2, -1,-1, 1,-1  &
     , 32, 4HSMA2,4H    , 1, 5, 2, 2, -2,-1,-1, 1,-1, 1,-1, 1,-1  &
     ,  1,-1, 1,-1, 1,-1, 1,-1, 1,-1 ,  1,-1, 1,-1, 1,-1  &
     , 10, 4HSMA3,4H    , 1, 2, 1, 7, -1,-1,-1 ,  7, 4HSMP1,4H    , 1, 5, 9, 7  &
     ,  7, 4HSMP2,4H    , 1, 3, 1, 6  &
     , 22, 4HSMPY,4HAD  , 1, 6, 1, 2, -1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0  &
     ,  1, 0, 1, 0 , 15, 4HSOLV,4HE   , 1, 2, 1, 5,  1, 0, 1, 1, 1, 0, 1, 0  &
     , 20, 19*0 /
 
 DATA  mpl14 / 13, 4HSSG1,4H    , 1,12, 5, 7, -1,-1, 1, 0, 2,14  &
     ,  7, 4HSSG2,4H    , 1, 7, 4, 4  &
     , 13, 4HSSG3,4H    , 1, 6, 4, 2, -1,-1, 1, 1, 1, 1  &
     ,  8, 4HSSG4,4H    , 1,11, 2, 5, -1  &
     , 23, 4HSSGH,4HT   , 1,17, 3, 5,  1,-1, 1,-1, 2,11, 2, 9, 1, 4  &
     ,  1,-1, 1, 0, 1, 0  &
     , 15, 4HTA1 ,4H    , 1, 8, 8, 4, -1,-1, 1, 1,-1,-1, 1, 1  &
     , 22, 4HTABP,4HCH  , 2, 5, 0, 0,  3,     4HAA  ,4H      ,3  &
     ,  4HAB  ,4H    ,3,4HAC  ,4H    ,3 ,  4HAD  ,4H    ,3,4HAE  ,4H    ,  &
     13, 12*0 , 12, 4HTABP,4HRT  , 2, 1, 0, 0, -3, 1, 0, 1, 0  &
     ,  7, 4HTABP,4HT   , 2, 5, 0, 0 , 18, 17*0  &
     /
 
 DATA  mpl15 /  &
     17, 4HTIME,4HTEST, 1, 0, 0, 2,  1,50, 1,50, 1, 2, 1, 1, 1, 511  &
     , 13, 4HTRD ,4H    , 1, 8, 3, 9, -3, 3*-1, 1,-1  &
     , 17, 4HTRHT,4H    , 1,10, 2, 7,  2,12, 2, 9, 1,-1, 1,-1, 2, 9  &
     , 11, 4HTRLG,4H    , 1,15, 6, 9,  1,-1, 1, 0  &
     ,  9, 4HTRNS,4HP   , 1, 1, 1, 8,  1, 0  &
     , 10, 4HUMER,4HGE  , 1, 3, 1, 1, -3,-3,-3  &
     , 10, 4HUPAR,4HTN  , 1, 2, 4, 1, -3,-3,-3  &
     , 14, 4HVDR ,4H    , 1, 7, 2, 2, -3,-3,-1, 1, 0,-1,-1  &
     , 16, 4HVEC ,4H    , 1, 1, 1, 0, -3, 3, 4HCOMP ,1H   , 3  &
     ,        4HCOMP, 1H   , 1, 0 , 18, 17*0  &
     /
 
 DATA  mpl16 / 7, 4HXYPL,4HOT  , 2, 1, 0, 2  &
     ,  7, 4HXYPR,4HNPLT, 2, 1, 0, 0  &
     , 19, 4HXYTR,4HAN  , 1, 6, 1, 5,  3,4HTRAN,4HS   ,3,4HSOL ,4H    ,  &
     1, 0, 1,0, 1,1 , 20, 19*0  &
     /
 
 DATA  mpl17 /  &
     13, 4HCOMB,4H1    , 1, 2, 1,10,  1,     0,-1,  3,4H    ,4H    ,  &
     35, 4HCOMB,4H2    , 1, 7, 1, 7, -1,    -3,     3,4H    ,4H    ,3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,3,4H    ,4H    ,3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,1,0  &
     , 36, 4HEXIO,4H     , 2, 0, 0, 2,  2*-1,  5*-3,  3  &
     ,  4HALL ,4H    ,3,4HWHOL,4HESOF,3 ,  4HXXXX,4HXXXX,3,4HXXXX,4HXXXX,3  &
     ,  4HXXXX,4HXXXX,3,4HXXXX,4HXXXX ,  1, 0,  1, 0  &
     , 39, 4HRCOV,4HR    , 1,11, 8, 9,  3*-1, -3,-1,  1, 0,    1,0,   3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,3,4H    ,4H    ,3  &
     ,  4H    ,4H    ,1,-1,2,9,2,9, 2,9
 
!     THE BCD PARAMETER   IN THE NEXT MODULE IS A DUMMY SINCE WE NEED
!     11 WORDS IN THIS SPACE
  , 11, 4HEMFL,4HD    , 1,10, 1, 1, -1,   3,4H    ,4H    ,  &
     10, 9*0 /
 
 DATA  mpl18 / 11, 4HRCOV,4HR3   , 1, 4, 7, 3, -1, -3, 1,-1  &
     , 16, 4HREDU,4HCE   , 1, 2, 3, 2,  1,  0, 1, 0,  3,4H    ,4H    ,1 ,  0  &
     , 11, 4HSGEN,4H     , 1, 4,10, 0, -1, -3,-1,-1  &
     , 25, 4HSOFI,4H     , 1, 0, 5, 0,  1, -1,-3,     3,4H    ,4H    ,3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,3,4H    ,4H    ,  &
     25, 4HSOFO,4H     , 2, 5, 0, 0,  1, -1,-3,     3,4H    ,4H    ,3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,3,4H    ,4H    ,  &
     37, 4HSOFU,4HT    , 2, 0, 0, 1,  1, -1, 2*-3,  1,0     ,3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,3,4H    ,4H    ,3  &
     ,  4H    ,4H    ,3,4H    ,4H    ,3 ,  4H    ,4H    ,3,4H    ,4H    ,  &
     15, 4HSUBP,4HH1   , 2, 7, 0, 1,  1,0,-3,1,0,   3,4H    ,4H    ,  &
     11, 4HPLTM,4HRG   , 1, 2, 6, 1, -3,3*-1, 18, 17*0  &
     /
 
 DATA  mpl19 / 9, 4HCOPY,4H    , 1, 1, 1, 0,  1, -1  &
     ,  9, 4HSWIT,4HCH  , 2, 2, 0, 0,  1, -1  &
     , 11, 4HMPY3,4H    , 1, 3, 1, 3,  1,  0, 1, 0  &
     , 33, 4HSDCM,4HPS  , 1, 4, 2, 6,  1,  0, 1, 0, 1,20, 1,0,   1,0  &
     ,  3, 1HL,1H  , 1, 0, 5,9,9, 4,8,8 ,  1,  0, 3,4H non,  1HE  &
     ,  9, 4HLODA,4HPP  , 2, 2, 0, 8, -3, -1 ,  7, 4HGPST,4HGEN , 1, 2, 1, 0  &
     , 15, 4HEQMC,4HK   , 1,12, 1, 7,  1,  0, 1,-1, -1, 3,4H non,1HE  &
     , 11, 4HADR ,4H    , 1, 7, 1, 5, -2,  2, 9,-3  &
     , 12, 4HFRRD,4H2   , 1, 6, 1, 9, -2,  2, 9, 2, 9  &
     , 14, 4HGUST,4H    , 1,10, 1, 7, -1,  2, 9, 2, 9, 2, 9  &
     ,  9, 4HIFT ,4H    , 1, 4, 2, 0,  1,  1  &
     ,  9, 4HLAMX,4H    , 1, 2, 1, 0,  1,  0  &
     , 11, 4HEMA ,4H    , 1, 3, 1, 2,  1, -1, 2, 6  &
     ,  9, 4HANIS,4HOP  , 1, 5, 1, 0,  1,  1 , 25, 24*0  &
     /
 DATA  mpl20 / 11, 4HGENC,4HOS  , 1, 2, 1, 0, -1, -1,-1,-1  &
     ,  8, 4HDDAM,4HAT  , 1, 2, 1, 0, -2  &
     ,  9, 4HDDAM,4HPG  , 1, 2, 1, 0, -1, -1  &
     , 11, 4HNRLS,4HUM  , 1, 2, 2, 3, -1, -1,-1,-1  &
     ,  9, 4HGENP,4HART , 1, 1, 4, 0, -1, -1  &
     , 10, 4HCASE,4HGEN , 1, 1, 1, 0, -1, -1,-1  &
     , 21, 4HDESV,4HEL  , 1, 2, 5, 0,  14*-2
 
!     3 DUMMY PARAMETERS IN PROLATE SO THAT AXLOOP CAN HAVE A PARAMETER
!     IN THE SAME POSITION IN BOTH SSG1 AND PROLATE
  , 15, 4HPROL,4HATE , 1,10, 1, 2,  1, -1, 1,-1, 1,-1, 2,14  &
     ,  8, 4HMAGB,4HDY  , 1, 2, 1, 0, -1  &
     ,  9, 4HCOMB,4HUGV , 1, 1, 5, 0, -1, -1  &
     , 14, 4HFLBM,4HG   , 1, 9, 4, 7, -1, -1, 5, 7, 7, 1, 0  &
     , 17, 4HGFSM,4HA   , 1,14, 5, 8, -1, -1, 2, 6, 1,-1, 1,-1, 1,-1  &
     , 10, 4HTRAI,4HLER , 2, 1, 0, 0, -3, -1,-1
!RLBR 12/29/93  SPR 93010 & 93011
!    4, 24, 4HSCAN,4H    , 1, 3, 1, 1,  3,4H    , 4H     , 1, 0, 1,20  &
 , 24, 4HSCAN,4H    , 1, 5, 2, 1,  3,4H    , 4H     , 1, 0, 1,20  &
     ,  2,  9, 2, 9, 1, 0, 1, 0, 1, 0 , 10, 9*0  &
     /
 DATA  mpl21 / 9, 4HPLTH,4HBDY , 1, 6, 4, 3, -1, -3  &
     ,  9, 4HVARI,4HAN  , 1, 5, 5, 3, -3, -2  &
     , 21, 4HFVRS,4HTR1 , 1, 8, 8,10,  13*-1,-2  &
     , 15, 4HFVRS,4HTR2 , 1, 8, 8,10,  8* -1  &
     , 29, 4HALG ,4H    , 1, 7, 2, 4,  1, -1, 1, -1, 1, -1, 1, -1, 1, 0  &
     ,  1,  0, 2,  6, 2,  9, 2,  6, 2, 6 ,  2,  6  &
     , 20, 4HAPDB,4H    , 1, 7, 5, 5, -1, -1, 2, 15, 2, 16, 1, -1  &
     ,  3,  4HCOSI,4HNE    ,-1, -1  &
     , 27, 4HPROM,4HPT1 , 2, 0, 0, 0,  1,  0, 1,  0, 1,  0, 1,  0, 1, 0  &
     ,  1,  0, 1,  0, 1,  0, 1,  0, 1, 0 ,  7, 4HSITE,4HPLOT, 2, 0, 0, 0  &
     , 16, 4HINPU,4HTT5 , 1, 0, 5, 0,  1,  0, 1, 11, 3,4HXXXX,4HXXXX ,  1,  0  &
     , 36, 4HOUTP,4HUT5 , 2, 5, 0, 0,  1,  0, 1, 11, 3,4HXXXX,4HXXXX  &
     ,  1,  0, 1,  0, 1,  0, 1,  0, 1, 0 ,  1,  0, 1,  0, 1,  0, 1,  0, 1, 0  &
     ,  1,  0 ,  7, 6*0  /
 DATA  mpl22 /  &
     34, 4HPARA,4HMD  , 2, 0, 0, 0, -3,  4,8,8, 4,8,8, 4,8,8, 6,4*8  &
     ,  6,  4*8,   6,4*8, 1,0  &
     , 12, 4HGINO,4HFILE, 1, 0, 1, 1, -1,  1,0,   1,999999  &
     , 13, 4HDATA,4HBASE, 2, 7, 0, 1,  1, 11,  1, 0,  1, 0  &
     , 16, 4HNORM,4H    , 1, 1, 1, 0,  1,  0,  1, 0,  2, 9, 3  &
     ,  4HMAX , 4H      , 13, 4HVECG,4HRB  , 1, 3, 1, 0,  1,  0,  1, 0,  1, 0  &
     , 21, 4HAUTO,4HASET, 1, 6, 2, 1,  1, -1,  1, -1, 1, -1, 1, -1  &
     ,  1, -1,  1, -1, 1, -1 , 10, 9*0  /
 
!     INITIALIZE /XGPI2/
 
 lmpl = lmplx
 DO  i = 1,lmpl
   imp(i) = mpl(i)
 END DO
 
!     INITIALIZE /XGPI2X/
 
 DO  i = 1,20
   xxx(i) = xx(i)
 END DO
 
 RETURN
END SUBROUTINE xmpldd
