BLOCK DATA sdr2bd
!SDR2BD
 IMPLICIT INTEGER (a-z)
 
 INTEGER :: rfmts(40)
 
 COMMON /sdr2x1/ ieigen,ieldef,itload,isymfl,iloads,idispl,istr  ,  &
     ielf  ,iacc  ,ivel  ,ispcf ,ittl  ,ilsym ,ifrout, isload,idload,isorc
 
 COMMON /sdr2x2/ casecc,cstm  ,mpt   ,dit   ,eqexin,sil   ,gptt  ,  &
     edt   ,bgpdt ,pg    ,qg    ,ugv   ,est   ,phig  ,  &
     eigr  ,opg1  ,oqg1  ,ougv1 ,oes1  ,oef1  ,pugv1 ,  &
     oeigr ,ophig ,pphig ,esta  ,gptta ,harms ,xycdb , scr3  ,pcomps,oes1l ,oef1l
 
COMMON /sdr2x4/ nam(2),END   ,mset  ,icb(7),ocb(7),mcb(7),dtype(8)  &
    ,               icstm ,ncstm ,ivec  ,ivecn ,temp  ,deform,FILE  ,  &
    buf1  ,buf2  ,buf3  ,buf4  ,buf5  ,any   ,all   ,  &
    tloads,eldef ,symflg,branch,ktype ,loads ,spcf  ,  &
    displ ,vel   ,acc   ,stress,force ,kwdest,kwdedt,  &
    kwdgpt,kwdcc ,nrigds,sta(2),rei(2),ds0(2),ds1(2),  &
    frq(2),trn(2),bk0(2),bk1(2),cei(2),pla(22)      ,  &
    nrings,nharms,axic  ,knset ,isopl ,strspt,ddrmm , isopl8

EQUIVALENCE     (sta(1),rfmts(1))

!*****
!     DATA DEFINING POSITIONS OF PARAMETERS IN A CASE CONTROL RECORD.
!*****
DATA  ieigen/  5/,ieldef/  6/,itload/  7/,isymfl/ 16/,iloads/ 17/,  &
    idispl/ 20/,istr  / 23/,ielf  / 26/,iacc  / 29/,ivel  / 32/,  &
    ispcf / 35/,ittl  / 39/,ilsym /200/,ifrout/145/,isload/  4/,  &
    idload/ 13/,isorc /136/
!*****
!     DATA DEFINING DATA BLOCK FILE NUMBERS.
!*****
DATA  casecc/101/,cstm  /102/,mpt   /103/,dit   /104/,eqexin/105/,  &
    sil   /106/,gptt  /107/,edt   /108/,bgpdt /109/,pg    /110/,  &
    qg    /111/,ugv   /112/,est   /113/,phig  /112/,eigr  /110/,  &
    opg1  /201/,oqg1  /202/,ougv1 /203/,oes1  /204/,oef1  /205/,  &
    pugv1 /206/,oeigr /201/,ophig /203/,pphig /206/,esta  /301/,  &
    gptta /302/,harms /137/,xycdb /114/,scr3  /303/,pcomps/116/,  &
    oes1l /207/,oef1l /208/
!*****
!     DATA DEFINING RIGID FORMATS.
!*****
DATA  nrigds/ 10   /, rfmts / 4HSTAT,4HICS ,  &
    4HREIG,4HEN  , 4HDS0 ,4H    ,  &
    4HDS1 ,4H    , 4HFREQ,4H    ,  &
    4HTRAN,4HSNT , 4HBKL0,4H    ,  &
    4HBKL1,4H    , 4HCEIG,4HEN  ,  &
    4HPLA ,4H    , 20*0         /
!*****
!     MISC. DATA.
!*****
DATA   nam  / 4HSDR2,4H    /, END/4HEND /, dtype/2,3,1,5,4,6,7,8/,  &
    mset / 1001/ ,isopl8/  0  /

END BLOCK DATA
