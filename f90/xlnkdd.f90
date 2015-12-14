SUBROUTINE xlnkdd
     
!     LINK SPECIFICATION TABLE
 
!     A LINK TABLE ENTRY CONTAINS AN EXECUTABLE DMAP INSTRUCTION NAMES,
!     THEIR CORRESPONDING SUBROUTINE ENTRY POINT NAMES, AND THE LINK OR
!     LINKS WHERE THEY RESIDE
!     EACH BIT IN THE LINK FLAG SPECIFIES A LINK NUMBER. BIT 1 (RIGHT
!     MOST) SPECIFIES LINK ONE, BIT 2 SPECIFIES LINK 2, ETC.
!     BIT ON SPECIFIES MODULE IS IN THAT LINK, BIT OFF MEANS IT IS NOT.
!     EXAMPLE - SUPPOSE MODULE X IS IN LINKS 2,4 AND 5. ITS LINK
!               FLAG=32(8).
 
!     LLINK = LENGTH OF LINK TABLE.
 
!     SET SENSE SWITCH 28 TO GENERATE ALL FORTRAN CODE BELOW.
 
 DIMENSION       link (960), link01(90), link02(90), link03(90),link04(90),  &
     link05(90), link06(90), link07(90),link08(90),  &
     link09(90), link10(90), link11(70)
 COMMON /xlkspc/ llink     , klink(970)
 EQUIVALENCE     (link(  1), link01(1)), (link( 91),link02(1)),  &
     (link(181), link03(1)), (link(271),link04(1)),  &
     (link(361), link05(1)), (link(451),link06(1)),  &
     (link(541), link07(1)), (link(631),link08(1)),  &
     (link(721), link09(1)), (link(811),link10(1)), (link(901), link11(1))
 DATA    llinkx / 970 /
 DATA    link01 / 4HCHKP,4HNT  , 4HXCHK,4H    , 32767 ,  &
     4HREPT,4H    , 4HXCEI,4H    , 32767 ,  &
     4HJUMP,4H    , 4HXCEI,4H    , 32767 ,  &
     4HCOND,4H    , 4HXCEI,4H    , 32767 ,  &
     4HSAVE,4H    , 4HXSAV,4HE   , 32766 ,  &
     4HPURG,4HE   , 4HXPUR,4HGE  , 32767 ,  &
     4HEQUI,4HV   , 4HXEQU,4HIV  , 32767 ,  &
4HEND ,4H    , 4HXCEI,4H    , 32767 ,  &
    4HEXIT,4H    , 4HXCEI,4H    , 32767 , 4HADD ,4H    , 4HDADD,4H    , 72    ,  &
    4HADD5,4H    , 4HDADD,4H5   , 64    , 4HAMG ,4H    , 4HAMG ,4H    , 256   ,  &
    4HAMP ,4H    , 4HAMP ,4H    , 256   , 4HAPD ,4H    , 4HAPD ,4H    , 256   ,  &
    4HBMG ,4H    , 4HBMG ,4H    , 512   , 4HCASE,4H    , 4HCASE,4H    , 512   ,  &
    4HCEAD,4H    , 4HCEAD,4H    , 1024  , 4HCYCT,4H1   , 4HCYCT,4H1   , 64    /
DATA    link02 / 4HCYCT,4H2   , 4HCYCT,4H2   , 64    ,  &
    4HDDR ,4H    , 4HDDR ,4H    , 128   , 4HDDR1,4H    , 4HDDR1,4H    , 2048  ,  &
    4HDDR2,4H    , 4HDDR2,4H    , 2048  , 4HDDRM,4HM   , 4HDDRM,4HM   , 2048  ,  &
    4HDECO,4HMP  , 4HDDCO,4HMP  , 64    , 4HDIAG,4HONAL, 4HDIAG,4HON  , 16384 ,  &
    4HDPD ,4H    , 4HDPD ,4H    , 32    , 4HDSCH,4HK   , 4HDSCH,4HK   , 64    ,  &
    4HDSMG,4H1   , 4HDSMG,4H1   , 4096  , 4HDSMG,4H2   , 4HDSMG,4H2   , 8     ,  &
    4HDUMM,4HOD1 , 4HDUMO,4HD1  , 4     , 4HDUMM,4HOD2 , 4HDUMO,4HD2  , 64    ,  &
    4HDUMM,4HOD3 , 4HDUMO,4HD3  , 64    , 4HDUMM,4HOD4 , 4HDUMO,4HD4  , 64    ,  &
    4HEMA1,4H    , 4HEMA1,4H    , 128   , 4HEMG ,4H    , 4HEMG ,4H    , 128   ,  &
    4HFA1 ,4H    , 4HFA1 ,4H    , 1024  /
DATA    link03 / 4HFA2 ,4H    , 4HFA2 ,4H    , 1024  ,  &
    4HFBS ,4H    , 4HDFBS,4H    , 64    , 4HFRLG,4H    , 4HFRLG,4H    , 512   ,  &
    4HFRRD,4H    , 4HFRRD,4H    , 512   , 4HGI  ,4H    , 4HGI  ,4H    , 256   ,  &
    4HGKAD,4H    , 4HGKAD,4H    , 512   , 4HGKAM,4H    , 4HGKAM,4H    , 512   ,  &
    4HGP1 ,4H    , 4HGP1 ,4H    , 2     , 4HGP2 ,4H    , 4HGP2 ,4H    , 2     ,  &
    4HGP3 ,4H    , 4HGP3 ,4H    , 2     , 4HGP4 ,4H    , 4HGP4 ,4H    , 8     ,  &
    4HGPCY,4HC   , 4HGPCY,4HC   , 64    , 4HGPFD,4HR   , 4HGPFD,4HR   , 4096  ,  &
    4HDUMM,4HOD5 , 4HDUMO,4HD5  , 64    , 4HGPWG,4H    , 4HGPWG,4H    , 8     ,  &
    4HINPU,4HT   , 4HINPU,4HT   , 2     , 4HINPU,4HTT1 , 4HINPT,4HT1  , 2     ,  &
    4HINPU,4HTT2 , 4HINPT,4HT2  , 2     /
DATA    link04 / 4HINPU,4HTT3 , 4HINPT,4HT3  , 2     ,  &
    4HINPU,4HTT4 , 4HINPT,4HT4  , 2     , 4HMATG,4HEN  , 4HMATG,4HEN  , 64    ,  &
    4HMATG,4HPR  , 4HMATG,4HPR  , 16    , 4HMATP,4HRN  , 4HMATP,4HRN  , 32766 ,  &
    4HMATP,4HRT  , 4HPRTI,4HNT  , 32766 , 4HMCE1,4H    , 4HMCE1,4H    , 8     ,  &
    4HMCE2,4H    , 4HMCE2,4H    , 8     , 4HMERG,4HE   , 4HMERG,4HE1  , 64    ,  &
    4HMODA,4H    , 4HMODA,4H    , 64    , 4HMODA,4HCC  , 4HMODA,4HCC  , 2048  ,  &
    4HMODB,4H    , 4HMODB,4H    , 64    , 4HMODC,4H    , 4HMODC,4H    , 64    ,  &
    4HMPYA,4HD   , 4HDMPY,4HAD  , 64    , 4HMTRX,4HIN  , 4HMTRX,4HIN  , 512   ,  &
    4HOFP ,4H    , 4HOFP ,4H    , 8192  , 4HOPTP,4HR1  , 4HOPTP,4HR1  , 2     ,  &
    4HOPTP,4HR2  , 4HOPTP,4HR2  , 128   /
DATA    link05 / 4HOUTP,4HUT  , 4HOUTP,4HT   , 8192  ,  &
    4HOUTP,4HUT1 , 4HOUTP,4HT1  , 8192  , 4HOUTP,4HUT2 , 4HOUTP,4HT2  , 8192  ,  &
    4HOUTP,4HUT3 , 4HOUTP,4HT3  , 8192  , 4HOUTP,4HUT4 , 4HOUTP,4HT4  , 8192  ,  &
    4HPARA,4HM   , 4HQPAR,4HAM  , 32766 , 4HPARA,4HML  , 4HPARA,4HML  , 32766 ,  &
    4HPARA,4HMR  , 4HQPAR,4HMR  , 32766 , 4HPART,4HN   , 4HPART,4HN1  , 64    ,  &
    4HMRED,4H1   , 4HMRED,4H1   , 16384 , 4HMRED,4H2   , 4HMRED,4H2   , 16384 ,  &
    4HCMRE,4HD2  , 4HCMRD,4H2   , 16384 , 4HPLA1,4H    , 4HPLA1,4H    , 4     ,  &
    4HPLA2,4H    , 4HPLA2,4H    , 4096  , 4HPLA3,4H    , 4HPLA3,4H    , 4096  ,  &
    4HPLA4,4H    , 4HPLA4,4H    , 4096  , 4HPLOT,4H    , 4HDPLO,4HT   , 2     ,  &
    4HPLTS,4HET  , 4HDPLT,4HST  , 2     /
DATA    link06 / 4HPLTT,4HRAN , 4HPLTT,4HRA  , 2     ,  &
    4HPRTM,4HSG  , 4HPRTM,4HSG  , 2     , 4HPRTP,4HARM , 4HPRTP,4HRM  , 128   ,  &
    4HRAND,4HOM  , 4HRAND,4HOM  , 8192  , 4HRMG ,4H    , 4HRMG ,4H    , 16    ,  &
    4HRBMG,4H1   , 4HRBMG,4H1   , 8     , 4HRBMG,4H2   , 4HRBMG,4H2   , 8     ,  &
    4HRBMG,4H3   , 4HRBMG,4H3   , 8     , 4HRBMG,4H4   , 4HRBMG,4H4   , 8     ,  &
    4HREAD,4H    , 4HREIG,4H    , 32    , 4HSCAL,4HAR  , 4HSCAL,4HAR  , 16384 ,  &
    4HSCE1,4H    , 4HSCE1,4H    , 8     , 4HSDR1,4H    , 4HSDR1,4H    , 2048  ,  &
    4HSDR2,4H    , 4HSDR2,4H    , 4096  , 4HSDR3,4H    , 4HSDR3,4H    , 8192  ,  &
    4HSDRH,4HT   , 4HSDRH,4HT   , 4096  , 4HSEEM,4HAT  , 4HSEEM,4HAT  , 2     ,  &
    4HSETV,4HAL  , 4HSETV,4HAL  , 32766 /
DATA    link07 / 4HSMA1,4H    , 4HSMA1,4H    , 4     ,  &
    4HSMA2,4H    , 4HSMA2,4H    , 4     , 4HSMA3,4H    , 4HSMA3,4H    , 8     ,  &
    4HSMP1,4H    , 4HSMP1,4H    , 8     , 4HSMP2,4H    , 4HSMP2,4H    , 8     ,  &
    4HSMPY,4HAD  , 4HSMPY,4HAD  , 64    , 4HSOLV,4HE   , 4HSOLV,4HE   , 64    ,  &
    4HSSG1,4H    , 4HSSG1,4H    , 16    , 4HSSG2,4H    , 4HSSG2,4H    , 16    ,  &
    4HSSG3,4H    , 4HSSG3,4H    , 16    , 4HSSG4,4H    , 4HSSG4,4H    , 16    ,  &
    4HSSGH,4HT   , 4HSSGH,4HT   , 16    , 4HTA1 ,4H    , 4HTA1 ,4H    , 2     ,  &
    4HCURV,4H    , 4HCURV,4H    , 4096  , 4HTABP,4HCH  , 4HTABP,4HCH  , 32766 ,  &
    4HTABP,4HRT  , 4HTABF,4HMT  , 32766 , 4HTABP,4HT   , 4HTABP,4HT   , 32766 ,  &
    4HTIME,4HTEST, 4HTIMT,4HST  , 256   /
DATA    link08 / 4HTRD ,4H    , 4HTRD ,4H    , 1024  ,  &
    4HTRHT,4H    , 4HTRHT,4H    , 1024  , 4HTRLG,4H    , 4HTRLG,4H    , 16    ,  &
    4HTRNS,4HP   , 4HDTRA,4HNP  , 64    , 4HUMER,4HGE  , 4HDUME,4HRG  , 64    ,  &
    4HUPAR,4HTN  , 4HDUPA,4HRT  , 64    , 4HVDR ,4H    , 4HVDR ,4H    , 2048  ,  &
    4HVEC ,4H    , 4HVEC ,4H    , 64    , 4HXYPL,4HOT  , 4HXYPL,4HOT  , 2     ,  &
    4HXYPR,4HNPLT, 4HXYPR,4HPT  , 8192  , 4HXYTR,4HAN  , 4HXYTR,4HAN  , 2     ,  &
    4HCOMB,4H1   , 4HCOMB,4H1   , 16384 , 4HCOMB,4H2   , 4HCOMB,4H2   , 16384 ,  &
    4HEXIO,4H    , 4HEXIO,4H    , 16384 , 4HRCOV,4HR   , 4HRCOV,4HR   , 16384 ,  &
    4HRCOV,4HR3  , 4HRCOV,4HR3  , 16384 , 4HREDU,4HCE  , 4HREDU,4HCE  , 16384 ,  &
    4HSGEN,4H    , 4HSGEN,4H    , 16384 /
DATA    link09 / 4HSOFI,4H    , 4HSOFI,4H    , 16384 ,  &
    4HSOFO,4H    , 4HSOFO,4H    , 16384 , 4HSOFU,4HT   , 4HSOFU,4HT   , 16384 ,  &
    4HSUBP,4HH1  , 4HSUBP,4HH1  , 16384 , 4HPLTM,4HRG  , 4HPLTM,4HRG  , 16384 ,  &
    4HCOPY,4H    , 4HCOPY,4H    , 64    , 4HSWIT,4HCH  , 4HSWIT,4HCH  , 64    ,  &
    4HMPY3,4H    , 4HMPY3,4H    , 64    , 4HSDCM,4HPS  , 4HDDCM,4HPS  , 64    ,  &
    4HLODA,4HPP  , 4HLODA,4HPP  , 16384 , 4HGPST,4HGEN , 4HGPST,4HGN  , 8     ,  &
    4HEQMC,4HK   , 4HEQMC,4HK   , 2048  , 4HADR ,4H    , 4HADR ,4H    , 512   ,  &
    4HFRRD,4H2   , 4HFRRD,4H2   , 512   , 4HGUST,4H    , 4HGUST,4H    , 512   ,  &
    4HIFT ,4H    , 4HIFT ,4H    , 512   , 4HLAMX,4H    , 4HLAMX,4H    , 256   ,  &
    4HEMA ,4H    , 4HEMA ,4H    , 128   /
DATA    link10 / 4HANIS,4HOP  , 4HANIS,4HOP  , 2     ,  &
    4HEMFL,4HD   , 4HEMFL,4HD   , 4096  , 4HGENC,4HOS  , 4HGENC,4HOS  , 4096  ,  &
    4HDDAM,4HAT  , 4HDDAM,4HAT  , 4096  , 4HDDAM,4HPG  , 4HDDAM,4HPG  , 4096  ,  &
    4HNRLS,4HUM  , 4HNRLS,4HUM  , 4096  , 4HGENP,4HART , 4HGENP,4HAR  , 4096  ,  &
    4HCASE,4HGEN , 4HCASE,4HGE  , 4096  , 4HDESV,4HEL  , 4HDESV,4HEL  , 4096  ,  &
    4HPROL,4HATE , 4HPROL,4HAT  , 4096  , 4HMAGB,4HDY  , 4HMAGB,4HDY  , 16    ,  &
    4HCOMB,4HUGV , 4HCOMU,4HGV  , 4096  , 4HFLBM,4HG   , 4HFLBM,4HG   , 8     ,  &
    4HGFSM,4HA   , 4HGFSM,4HA   , 8     , 4HTRAI,4HLER , 4HTRAI,4HL   , 8     ,  &
    4HSCAN,4H    , 4HSCAN,4H    , 8192  , 4HPLTH,4HBDY , 4HPTHB,4HDY  , 2     ,  &
    4HVARI,4HAN  , 4HVARI,4HAN  , 8192  /
DATA    link11 / 4HFVRS,4HTR1 , 4HFVRS,4HT1  , 64    ,  &
    4HFVRS,4HTR2 , 4HFVRS,4HT2  , 64    , 4HALG ,4H    , 4HALG ,4H    , 32    ,  &
    4HAPDB,4H    , 4HAPDB,4H    , 256   , 4HPROM,4HPT1 , 4HPROM,4HPT  , 8194  ,  &
    4HSITE,4HPLOT, 4HOLPL,4HOT  , 2     , 4HINPU,4HTT5 , 4HINPT,4HT5  , 2     ,  &
    4HOUTP,4HUT5 , 4HOUTP,4HT5  , 8192  , 4HPARA,4HMD  , 4HQPAR,4HMD  , 32766 ,  &
    4HGINO,4HFILE, 4HGINO,4HFL  , 32766 , 4HDATA,4HBASE, 4HDBAS,4HE   , 8202  ,  &
    4HNORM,4H    , 4HNORM,4HAL  , 16    , 4HVECG,4HRB  , 4HGRBV,4HEC  , 64    ,  &
    4HAUTO,4HASET, 4HAASE,4HT   , 8     /

!     INITIALIZE /XLKSPC/

llink = llinkx
DO  i = 1,llink
  klink(i) = link(i)
END DO
RETURN
END SUBROUTINE xlnkdd
