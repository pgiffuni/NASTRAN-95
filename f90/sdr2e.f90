SUBROUTINE sdr2e (*,ieqex,neqex)
     
!     THIS ROUTINE WHICH IS CALLED ONLY FROM SDR2D WILL PROCESS THE ESTA
!     FILE ONCE AND OUTPUT FORCE AND OR STRESS RESULTS ON OEF1 AND OR
!     OES1 WHICH ARE OPENED IN SDR2D.
 
 INTEGER, INTENT(IN OUT)                  :: ieqex
 INTEGER, INTENT(IN OUT)                  :: neqex
 INTEGER    ::  opte, a, z, branch
 EXTERNAL        andf
 LOGICAL :: eorflg,endid ,record,acstic,axic  ,again ,idstrs,  &
     idforc,eofcc ,idlyst,idlyfr,ok2wrt,heat  ,ddrmm , strain,ilogic(4)
 INTEGER :: buf(50)      ,platit(12)   ,complx(478)  ,  &
     isavef(75)   ,isaves(75)
 REAL :: zz(1) ,bufr(1)      ,tgrid(33)    ,diff1 ,diff  ,  &
     deform,frtmei,temp  ,two_to_p,fnchk
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm
 COMMON /lhpwx / lh(6) ,mtisa
 COMMON /BLANK / app(2),sort2 ,istrn ,idum1 ,comps ,idum4(4)     , strain
 COMMON /sdr2c1/ ipcmp ,npcmp ,ipcmp1,npcmp1,ipcmp2,npcmp2,nstrop
 COMMON /sdr2x1/ ieigen,ieldef,itload,isymfl,iloads,idispl,istr  ,  &
     ielf  ,iacc  ,ivel  ,ispcf ,ittl  ,ilsym ,ifrout, isload,idload,isorc
 COMMON /sdr2x2/ casecc,cstm  ,mpt   ,dit   ,eqexin,sil   ,gptt  ,  &
     edt   ,bgpdt ,pg    ,qg    ,ugv   ,est   ,phig  ,  &
     eigr  ,opg1  ,oqg1  ,ougv1 ,oes1  ,oef1  ,pugv1 ,  &
     oeigr ,ophig ,pphig ,esta  ,gptta ,harms ,idum3(3)  &
     ,               oes1l ,oef1l
 COMMON /gpta1 / nelem ,last  ,incr  ,elem(1)
COMMON /sdr2x4/ nam(2),END   ,mset  ,icb(7),ocb(7),mcb(7),dtype(8)  &
    ,               icstm ,ncstm ,ivec  ,ivecn ,temp  ,deform,FILE  ,  &
    buf1  ,buf2  ,buf3  ,buf4  ,buf5  ,any   ,all   ,  &
    tloads,eldef ,symflg,branch,ktype ,loads ,spcf  ,  &
    displ ,vel   ,acc   ,stress,force ,kwdest,kwdedt,  &
    kwdgpt,kwdcc ,nrigds,sta(2),rei(2),ds0(2),ds1(2),  &
    frq(2),trn(2),bk0(2),bk1(2),cei(2),pla(22)      ,  &
    nrings,nharms,axic  ,knset ,isopl ,strspt,ddrmm
COMMON /sdr2x7/ elesta(100)  ,bufa(100)    ,bufb(4076)
COMMON /sdr2x8/ elwork(300)
COMMON /sdr2x9/ nchk  ,isub  ,ild   ,frtmei(2)    ,two_to_p,fnchk
COMMON /zzzzzz/ z(1)
COMMON /isave / isavef,isaves
COMMON /names / rd    ,rdrew ,wrt   ,wrtrew,clsrew
COMMON /clstrs/ complx
COMMON /system /ksystm(100)
COMMON /sdr2de/ buf6  ,coef  ,deftmp,diff  ,diff1 ,device,estawd,  &
    elemid,eltype,eof   ,eofcc ,ireqx ,flag  ,fn    ,  &
    forcex,fsetno,formt ,icc   ,i     ,iedt  ,isetno,  &
    isetf ,isets ,idef  ,isymn ,sdest ,ix    ,isetnf,  &
    iseq  ,iretrn,irecx ,isave ,fdest ,ipart ,ilist ,  &
    igptta,icore ,ielem ,iesta ,buf8  ,jforc ,jstrs ,  &
    jany  ,jlist ,j     ,ktype1,khi   ,kx    ,k     ,  &
    klo   ,kn    ,ktypex,kfrq  ,kcount,lsym  ,m     ,  &
    midvec,nwdsa ,nwdstr,nlogic,nwds  ,ndef  ,n     ,  &
    n1    ,n2    ,notset,nsets ,nsetf ,nwords,nx    ,  &
    tdumm(4)     ,nwdfor,ngptt ,nesta ,nvects,nlist ,  &
    ofile ,outfl ,retx  ,setno ,stresx,SAVE  ,tload ,  &
    ugvvec,ixsets,nxsets,ixsetf,nxsetf,xsetns,xsetnf, sorc  ,tmprec,buf7  ,tgrid
COMMON /sdrett/ ieltyp,oldel ,eorflg,endid ,bufflg,itemp ,xx2(2),  &
    record,oldeid
EQUIVALENCE     (ksystm( 2),opte  ) ,(ksystm(55),iprec  ), (ksystm(56),itherm)
EQUIVALENCE     (buf(1),bufr(1)   ) ,(z(1)  ,zz(1)      ),  &
    (idstrs,ilogic(1) ) ,(idforc,ilogic(2)  ),  &
    (idlyst,ilogic(3) ) ,(idlyfr,ilogic(4)  ),  &
    (temp  ,jtemp     ) ,(nelhar,elwork(155))
DATA platit / 4HLOAD,4H fac,4HTOR , 9*0    /
DATA buf    / 50*0 /, ielold / 0 /, ielchk / 0 /

!     INITIALIZE ESTA POINTERS.

IF (stresx == 0 .AND. forcex == 0) RETURN
heat   = .false.
IF (itherm /= 0) heat = .true.
estawd = iesta
istore = 0
ix     = icc + harms
again  =.false.
oharms = z(ix)
IF (oharms < 0) oharms = nharms
isave  = ivec
isvsrc = sorc
eltype = z(estawd)
FILE   = esta
ix     = icc + istr + 2
sphase = IABS(z(ix))
ix     = icc + ielf + 2
fphase = IABS(z(ix))
two_to_p = ALOG10(2.0**mtisa)

!     POSITION TO THE PROPER THERMAL RECORD IF NECESSARY.

record = .false.
IF (tloads == 0) GO TO 18
IF (tmprec == 0) GO TO 18
CALL REWIND (gptt)
FILE = gptt
DO  i = 1,tmprec
  CALL fwdrec (*980,gptt)
END DO

!     READ AND VERIFY SET-ID  (FAILSAFE)

CALL READ (*980,*990,gptt,isetid,1,0,flag)
IF (tloads == isetid) GO TO 17
WRITE  (opte,16) sfm,tloads,isetid
16 FORMAT (a25,' 4019, SDR2E DETECTS INVALID TEMPERATURE DATA FOR ',  &
    'TEMPERATURE LOAD SET',2I10)
CALL mesage (-61,0,0)
17 record = .true.

!     INITIALIZE /SDRETT/ VARIABLES

oldeid = 0
oldel  = 0
eorflg = .false.
endid  = .true.
18 itemp  = tloads
IF (nesta /= 0) GO TO 25
CALL REWIND (esta)
20 CALL READ (*950,*990,esta,eltype,1,0,flag)

!     ELEMENT PARAMETERS FOR NEW ELEMENT TYPE

25 ielem  = (eltype-1)*incr
ieltyp = eltype
ipr    = iprec
IF (ipr /= 1) ipr = 0
jltype = 2*eltype - ipr
jcore  = icore
IF (heat .AND. eltype /= 82) GO TO 27
!                            FTUBE
nwdsa  = elem(ielem+17)
nptstr = elem(ielem+20)
nptfor = elem(ielem+21)
nwdstr = elem(ielem+18)
nwdfor = elem(ielem+19)
GO TO 28

27 nwdfor = 9
nwdstr = 0
nptfor = 0
nptstr = 0
nwdsa  = 142

!     CHOP OFF 483 WORDS FROM OPEN CORE SPACE FOR CIHEX ELEMENTS

IF (eltype < 65 .OR. eltype > 67) GO TO 28
icore  = icore - 483

28 CONTINUE

!     SETUP STRESS PRECISION CHECK.

nchk  = z(icc+146)
fnchk = nchk

!     SUBCASE ID

isub  = z(icc+1)

!     DETERMINE LOAD/MODE, EIGENVALUE/FREQ/TIME HEADER

frtmei(1) = 0.
frtmei(2) = 0.
IF (branch == 5 .OR. branch == 6)  GO TO 31
IF (branch == 2 .OR. branch == 8 .OR. branch == 9) GO TO 32

!     STATICS

i   = icc + isload
ild = z(i)
GO TO 35

!     FREQUENCY/TRANSIENT

31 i   = icc + idload
ild = z(i)
frtmei(1) = zz(jlist)
GO TO 35

!     EIGENVALUES

32 ild       =  z(jlist  )
frtmei(1) = zz(jlist+1)
frtmei(2) = zz(jlist+2)
IF (branch /= 2) GO TO 35
IF (zz(jlist+1) > 0.) frtmei(1) = SQRT(zz(jlist+1))/6.2831852
35 CONTINUE
lstres = nwdstr
lforce = nwdfor
idstrs = .false.
idforc = .false.
idlyst = .false.
idlyfr = .false.
ok2wrt = .true.
IF (ktype /= 1 .AND. nptstr == 0 .AND. nptfor == 0) GO TO 40
IF (nwdstr+nwdfor > 0) GO TO 70

!     NO STRESS OR FORCE WORDS POSSIBLE FOR THIS ELEMENT TYPE IF FALL
!     HERE

40 IF (nesta == 0) THEN
  GO TO    50
ELSE
  GO TO    60
END IF

!     FORWARD REC ON FILE TO NEXT ELEMENT TYPE

50 CALL fwdrec (*980,esta)
GO TO 20

!     FIND END OF CURRENT ELEMEMT TYPE LIST IN CORE

60 estawd = estawd + nwdsa
IF (z(estawd+1) == 0.0) THEN
  GO TO   940
ELSE
  GO TO    60
END IF

!     OK SOME STRESS AND OR FORCE REQUESTS EXIST FOR THIS ELEMENT TYPE.
!     PROCESS INDIVIDUAL ELEMENTS REQUESTED

70 IF (nesta /=     0) GO TO 90
IF (nwdsa <= icore) GO TO 80
CALL mesage (8,0,nam(1))


!     INSUFFICIENT CORE TO HOLD ESTA FOR 1 ELEMENT OF CURRENT ELEMENT
!     TYPE TRY PROCESSING THE OTHER ELEMENT TYPES IN AVAILABLE CORE.

GO TO 50

80 CALL READ (*980,*910,esta,z(iesta),nwdsa,0,flag)
estawd = iesta - 1

!     DETERMINE IF THIS PARTICULAR ELEMENT OF THE CURRENT ELEMENT TYPE
!     HAS A STRESS OR FORCE REQUEST IN THE CURRENT CASE CONTROL RECORD.

90 elemid = z(estawd+1)

!     THE FOLLOWING CODE (THRU 93) IS FOR THE COMPLEX ANALYSIS OF IHEX
!     ELEMENTS ONLY (ELEM. TYPES 65,66,67)

IF (ktype /= 2 .OR. eltype < 65 .OR. eltype > 67) GO TO 93
IF (ipart /= 2 .OR. istrpt /= (nip3+ngp1+1)) GO TO 91

!     DONE FOR THIS IHEX ELEMENT, RESET CHECKING VARIABLES

ipart  = 0
ielold = 0
ielchk = 0
GO TO 93

!     FIRST INTEGRATION POINT FOR IMAGINARY RETULS FOR THIS IHEX ELEMENT
!     SAVE ELEMENT ID AND CURRENT ESTAWD

91 IF (ipart /= 1 .OR. istrpt /= 1) GO TO 92
ielold = elemid
oldawd = estawd - nwdsa
GO TO 93

!     FIRST INTEGRATION POINT FOR REAL RESULTS FOR THIS IHEX ELEMENT,
!     SAVE ELEMENT ID TO CHECK WITH EARLIER ELEMENT ID SAVED ABOVE

92 IF (ipart == 2 .AND. istrpt == 1) ielchk = elemid

!     END OF SPECIAL TREATMENT FOR IHEX ELEMENT

93 idelem = elemid

!     DECODE ELEMID TO FIND IT IN SET

IF (.NOT. axic) GO TO 95
nelhar = elemid - (elemid/1000)*1000
elemid = elemid/1000
95 jstrs  = 0
jforc  = 0
i      = isets
IF (nwdstr == 0) GO TO 140
IF (stresx < 0.0) THEN
  GO TO   110
ELSE IF (stresx == 0.0) THEN
  GO TO   140
END IF
100 IF (i == nsets) GO TO 120
IF (z(i+1) > 0) GO TO 120
i = i + 1
IF (elemid < z(i-1) .OR. elemid > -z(i)) GO TO 130
110 jstrs = 1
GO TO 140
120 IF (elemid == z(i)) GO TO 110
130 i = i + 1
IF (i <= nsets) GO TO 100
140 i = isetf
IF (nwdfor == 0) GO TO 190
IF (forcex < 0.0) THEN
  GO TO   160
ELSE IF (forcex == 0.0) THEN
  GO TO   190
END IF
150 IF (i == nsetf) GO TO 170
IF (z(i+1) > 0) GO TO 170
i = i + 1
IF (elemid < z(i-1) .OR. elemid > -z(i)) GO TO 180
160 jforc = 1
GO TO 190
170 IF (elemid == z(i)) GO TO 160
180 i = i + 1
IF (i <= nsetf) GO TO 150
190 jany= jstrs + jforc
IF (jany == 0) IF (nesta) 890,80,890

!     OK FALL HERE AND A STRESS OR FORCE REQUEST EXISTS
!     IF THERMAL LOADING, GET THE ELEMENT THERMAL DATA.
!     IF ELEMENT DEFORMATIONS, LOOK UP THE DEFORMATION


!     ELEMENT TEMPERATURE

IF (tloads == 0) GO TO 330
n = elem(ielem+10)

!     IF NEW ELEMENTS ARE ADDED THAT HAVE SPECIAL BENDING THERMAL DATA
!     POSSIBLE THEN THE FOLLOWING TEST SHOULD BE EXPANDED TO INCLUDE
!     THEIR ELEMENT TYPE SO AS TO RECEIVE ZEROS AND ONLY THE AVERAGE
!     TEMPERATURE RATHER THAN SIMULATED GRID POINT TEMPERATURES IN THE
!     ABSENCE OF ANY USER SPECIFIED DATA.

IF (ieltyp == 34 .OR. ieltyp == 6 .OR. ieltyp == 7  .OR.  &
    ieltyp == 8 .OR. ieltyp == 15 .OR. ieltyp == 17 .OR.  &
    ieltyp == 18 .OR. ieltyp == 19) n = 0
IF (ieltyp == 74 .OR. ieltyp == 75) n = 0
IF (ieltyp == 64 .OR. ieltyp == 83) n = 0


CALL sdretd (idelem,tgrid,n)

!     SET THE AVERAGE ELEMENT TEMPERATURE CELL.

temp = tgrid(1)
GO TO 340

!     NORMALLY TGRID(1) WILL CONTAIN THE AVERAGE ELEMENT TEMPERATUE
!     AND IF GRID POINT TEMPERATURES ARE RETURNED THEY WILL BEGIN
!     IN TGRID(2).

330 jtemp = -1

!     ELEMENT DEFORMATION

340 deform = 0.0
IF (eldef == 0) GO TO 360
DO  i = idef,ndef,2
  IF (z(i) == elemid) GO TO 355
END DO
GO TO 360
355 deform = zz(i+1)

!     WRITE ID FOR STRESSES IF NOT YET WRITTEN FOR THIS ELEMENT TYPE.

360 IF (stress == 0 .OR.  nwdstr == 0 .OR. jstrs == 0) GO TO 365
IF (comps == -1 .AND. nstrop > 1) GO TO 362
IF (idstrs) GO TO 365
nlogic = 1
ofile  = oes1
device = sdest
iseq   = 4
ifltyp = dtype(iseq)
irecx  = icc + istr
nwds   = nwdstr
jcmplx = nptstr
ASSIGN 365 TO iretrn
GO TO 630

362 IF (idlyst) GO TO 365
nlogic = 3
ofile  = oes1l
device = sdest
ifltyp = 22
irecx  = icc + istr
nwds   = 10
jcmplx = 0
ok2wrt = .false.
ASSIGN 365 TO iretrn
GO TO 630

!     WRITE ID FOR FORCES IF NOT YET WRITTEN FOR THIS ELEMENT TYPE.

365 IF (force == 0 .OR.  nwdfor == 0 .OR.  jforc == 0) GO TO 375
IF (comps == -1 .AND. nstrop > 1 .AND. stress /= 0) GO TO 367
IF (idforc) GO TO 375
nlogic = 2
ofile  = oef1
device = fdest
iseq   = 5
ifltyp = dtype(iseq)
irecx  = icc + ielf
nwds   = nwdfor
jcmplx = nptfor
ASSIGN 375 TO iretrn
GO TO 630

367 IF (idlyfr) GO TO 375
nlogic = 4
ofile  = oef1l
device = fdest
ifltyp = 23
irecx  = icc + ielf
nwds   = 9
jcmplx = 0
ok2wrt = .false.
ASSIGN 375 TO iretrn
GO TO 630

!     MOVE ESTA DATA INTO /SDR2X7/

375 nsesta = estawd
IF (ielchk == 0 .OR. ipart < 2 .OR. ielchk /= ielold) GO TO 377
ipart  = 1
GO TO 380
377 ipart  = 0
380 ipart  = ipart + 1
DO  i = 1,nwdsa
  estawd = estawd + 1
  elesta(i) = z(estawd)
END DO
acstic = .false.

!     CALL APPROPRIATE ELEMENT ROUTINE FOR STRESS AND FORCE COMPUTATIONS

IF (heat) GO TO 1680
local = jltype - 100
IF (local > 0) THEN
  GO TO   395
END IF

!     PAIRED -GO TO- ENTRIES PER ELEMENT SINGLE/DOUBLE PRECISION

!             1 CROD      2 C.....    3 CTUBE     4 CSHEAR    5 CTWIST
394 GO TO( 400,  400,  610,  610,  400,  400,  420,  420,  430,  430 &
!             6 CTRIA1    7 CTRBSC    8 CTRPLT    9 CTRMEM   10 CONROD  &
,      450,  450,  460,  460,  470,  470,  480,  480,  400,  400 &
!            11 ELAS1    12 ELAS2    13 ELAS3    14 ELAS4    15 CQDPLT  &
,      490,  490,  490,  490,  490,  490,  490,  490,  500,  500 &
!            16 CQDMEM   17 CTRIA2   18 CQUAD2   19 CQUAD1   20 CDAMP1  &
,      520,  520,  450,  450,  540,  540,  540,  540,  610,  610 &
!            21 CDAMP2   22 CDAMP3   23 CDAMP4   24 CVISC    25 CMASS1  &
,      610,  610,  610,  610,  610,  610,  610,  610,  610,  610 &
!            26 CMASS2   27 CMASS3   28 CMASS4   29 CONM1    30 CONM2  &
,      610,  610,  610,  610,  610,  610,  610,  610,  610,  610 &
!            31 PLOTEL   32 C.....   33 C.....   34 CBAR     35 CCONE  &
,      610,  610,  610,  610,  610,  610,  560,  560,  570,  570 &
!            36 CTRIARG  37 CTRAPRG  38 CTORDRG  39 CTETRA   40 CWEDGE  &
,      580,  580,  590,  590,  600,  600,  601,  601,  602,  602 &
!            41 CHEXA1   42 CHEXA2   43 CFLUID2  44 CFLUID3  45 CFLUID4  &
,      603,  603,  604,  604,  610,  610,  610,  610,  610,  610 &
!            46 CFLMASS  47 CAXIF2   48 CAXIF3   49 CAXIF4   50 CSLOT3  &
,      610,  610,  605,  605,  606,  606,  607,  607,  608,  608 &
 ), jltype

!            51 CSLOT4   52 CHBDY    53 CDUM1    54 CDUM2    55 CDUM3
395 GO TO( 609,  609,  610,  610, 1614, 1614, 1615, 1615, 1616, 1616 &
!            56 CDUM4    57 CDUM5    58 CDUM6    59 CDUM7    60 CDUM8  &
,     1617, 1617, 1618, 1618, 1619, 1619, 1620, 1620, 1621, 1621 &
!            61 CDUM9    62 CQDMEM1  63 CQDMEM2  64 CQUAD4   65 CIHEX1  &
,     1622, 1622, 1623, 1623, 1624, 1624, 1625, 1625, 1626, 1626 &
!            66 CIHEX2   67 CIHEX3   68 CQUADTS  69 CTRIATS  70 CTRIAAX  &
,     1626, 1626, 1626, 1626, 1632, 1632, 1633, 1633, 1634, 1634 &
!            71 CTRAPAX  72 CAERO1   73 CTRIM6   74 CTRPLT1  75 CTRSHL  &
,     1635, 1635,  610,  610, 1640, 1640, 1645, 1645, 1650, 1650 &
!            76 CFHEX1   77 CFHEX2   78 CFTETRA  79 CFWEDGE  80 CIS2D8  &
,      610,  610,  610,  610,  610,  610,  610,  610, 1660, 1660 &
!            81 CELBOW   82 CFTUBE   83 CTRIA3  &
,      1670, 1670, 610,  610, 1630, 1630  &
 ), local

400 CALL srod2
GO TO 620
420 k = 4
GO TO 440
430 k = 5
440 CALL spanl2 (k)
GO TO 620
450 k = 3
GO TO 550
460 k = 0
GO TO 510
470 k = 3
GO TO 510
480 k = 1
GO TO 530
490 CALL selas2
GO TO 620
500 k = 4
510 CALL sbspl2 (k,tgrid(1))
GO TO 620
520 k = 2
530 CALL stqme2 (k)
GO TO 620
540 k = 4
550 CALL strqd2 (k,tgrid(1))
GO TO 620
560 CALL sbar2 (tgrid(1))
GO TO 620
570 again = .false.
CALL scone2 (sorc)
GO TO 620
580 CALL strir2 (tgrid(2))
GO TO 620
590 CALL strap2 (tgrid(2))
GO TO 620
600 CALL stord2 (tgrid(2))
GO TO 620
601 CALL ssold2 (1,tgrid(2))
GO TO 620
602 CALL ssold2 (2,tgrid(2))
GO TO 620
603 CALL ssold2 (3,tgrid(2))
GO TO 620
604 CALL ssold2 (4,tgrid(2))
GO TO 620
605 kk = 0
GO TO 611
606 kk = 1
GO TO 611
607 kk = 2
611 CALL saxif2 (kk,ipart,branch,z(jlist))
acstic = .true.
GO TO 620
608 kk = 0
GO TO 612
609 kk = 1
612 CALL sslot2 (kk,ipart,branch,z(jlist))
acstic = .true.
GO TO 620
1614 CALL sdum12
GO TO 620
1615 CALL sdum22
GO TO 620
1616 CALL sdum32
GO TO 620
1617 CALL sdum42
GO TO 620
1618 CALL sdum52
GO TO  620
1619 CALL sdum62
GO TO 620
1620 CALL sdum72
GO TO 620
1621 CALL sdum82
GO TO 620
1622 CALL sdum92
GO TO 620
1623 CALL sqdm12
GO TO 620
1624 CALL sqdm22
GO TO 620
1625 CALL squd42
GO TO 620
1626 CALL sihex2 (eltype-64,tgrid(1),nip,istrpt,istore)
ngp  = 12*(eltype-64) - 4
ngp1 = ngp + 1
IF (eltype == 67) ngp1 = 21
nip3 = nip**3
IF (istrpt < nip3+1) GO TO 905
IF (istrpt == nip3+1) GO TO 1626
IF (istrpt == nip3+1+ngp1) istore = 0
IF (ktype == 1) GO TO 620
ngpx = istrpt - (nip3+1)
nw   = 22
IF (eltype == 67) nw = 23
ist = nw*(ngpx-1)
IF (ipart >= ktype) GO TO 1628

!     STORE IMARINARY PARTS FOR THIS GRID (IHEX ELEMENTS)

ijk = ist + icore
DO  j = 1,nw
  z(j+ijk) = bufa(j)
END DO
IF (istore /= 0) GO TO 1626
ivec   = midvec
estawd = oldawd
GO TO 380

!     RETRIEVE IMAGINARY PARTS FOR THIS GRID (IHEX ELEMENTS)

1628 ijk = ist + icore
DO  j = 1,nw
  isaves(j) = z(j+ijk)
END DO
GO TO 620

1630 CALL stri32
GO TO 620
1632 CONTINUE
GO TO 620
1633 CONTINUE
GO TO 620
1634 again = .false.
CALL strax2 (sorc,tgrid(2))
GO TO 620
1635 again = .false.
CALL stpax2 (sorc,tgrid(2))
GO TO 620
1640 CALL strm62 (tgrid(1))
GO TO 620
1645 CALL strp12 (tgrid(1))
GO TO 620
1650 CALL  strsl2 (tgrid(1))
GO TO 620
1660 CALL ss2d82 (ieqex,neqex,tgrid(1))
GO TO 620
1670 CALL selbo2 (tgrid(1))
GO TO 620

!     PHASE TWO HEAT ONLY (ALL ELEMENTS)

1680 CALL sdhtf2 (ieqex,neqex)
GO TO 620
610 GO TO 900

!     CALL ELEMENT TWO TIMES FOR COMPLEX VECTOR.  IMAGINARY FIRST, REAL
!     SECOND.  CALL ELEMENT ROUTINE TWICE IF AXIC PROBLEM
!     ONCE FOR EACH OF THE 2 VECTORS IN CORE

620 IF (axic .AND. midvec /= 0 .AND. ipart == 1) GO TO 625
IF (ipart >= ktype) GO TO 615
625 ivec = midvec

!     FOR CONICAL SHELL ONLY

IF (axic .AND. ktype /= 1) GO TO 626
itemp = 1
IF (sorc == 1) itemp = 2
sorc  = itemp
626 CONTINUE
estawd = nsesta
IF (axic .AND. ktype == 1) GO TO 380

!     SAVE IMAGINARY OUTPUTS  (NOT MORE THAN 75 STRESS OR FORCE WORDS)

DO  i = 1,75
  isaves(i) = bufa(i)
  isavef(i) = bufb(i)
END DO
GO TO 380

!     SPLIT OUTPUT FROM SECOND CALL FOR ACOUSTIC ELEMENTS
!     AXIF2, AXIF3, AXIF4, SLOT3, OR SLOT4.

615 IF (.NOT. acstic) GO TO 617
IF (ipart < 2) GO TO 617
DO  i = 1,12
  isaves(i) = bufa(i   )
  bufa(i)   = bufa(i+12)
END DO


!     OUTPUT ONLY FIRST N HARMONICS REQUESTED

617 IF (.NOT. axic) GO TO 616
IF (nelhar < 0 .OR. nelhar > oharms) GO TO 880
IF (ipart == 2 .AND. ktype == 1) GO TO 880

!     OUTPUT STRESS RESULTS ON OES1 (IF REQUESTED)

616 IF (jstrs == 0 .OR. nwdstr == 0) GO TO 860
IF (ktype == 1) GO TO 850

!     COMBINE COMPLEX OUTPUT DESIRED PER FORMAT IN COMPLX ARRAY.
!          REAL PARTS ARE IN BUFA   BUFB
!          IMAG PARTS ARE IN ISAVES ISAVEF


!     COMPLEX STRESSES

iout = 0
i    = nptstr
651 npt  = complx(i)
IF (npt < 0) THEN
  GO TO   652
ELSE IF (npt == 0) THEN
  GO TO   653
ELSE
  GO TO   654
END IF
652 npt  = -npt
IF (sphase /= 3) GO TO 654

!     COMPUTE MAGNITUDE/PHASE

CALL magpha (bufa(npt),isaves(npt))
655 iout = iout + 1
elwork(iout) = bufa(npt)
i    = i + 1
GO TO 651
654 IF (npt <= lstres) GO TO 655
npt  = npt - lstres
iout = iout + 1
elwork(iout) = isaves(npt)
i    = i + 1
GO TO 651

!     TRANSFER RESULTS TO BUFA

653 DO  i = 1,iout
  bufa(i) = elwork(i)
END DO
nwdstr  = iout

!     WRITE STRESSES


!     DETERMINE DESTINATION FOR STRESS ENTRY

850 IF (stress == 0) GO TO 860
IF (.NOT.  ok2wrt) GO TO 860
id = bufa(1)
bufa(1) = 10*id + sdest
IF (xsetns < 0.0) THEN
  GO TO   858
ELSE IF (xsetns == 0.0) THEN
  GO TO   851
ELSE
  GO TO   852
END IF
851 bufa(1) = 10*id
GO TO 858
852 ix = ixsets
853 IF (ix == nxsets) GO TO 854
IF (z(ix+1) > 0) GO TO 854
IF (id >= z(ix) .AND. id <= (-z(ix+1))) GO TO 858
ix = ix + 2
GO TO 855
854 IF (id == z(ix)) GO TO 858
ix = ix + 1
855 IF (ix <= nxsets) GO TO 853
GO TO 851

!     NOW WRITE STRESS ENTRY

858 CALL WRITE (oes1,bufa(1),nwdstr,0)
bufa(1) = id

!     OUTPUT FORCE RESULTS ON OEF1 (IF REQUESTED)

860 IF (jforc == 0  .OR. nwdfor == 0) GO TO 880
IF (ktype == 1) GO TO 870

!     COMPLEX FORCES

iout = 0
i    = nptfor
951 npt  = complx(i)
IF (npt < 0) THEN
  GO TO   952
ELSE IF (npt == 0) THEN
  GO TO   953
ELSE
  GO TO   954
END IF
952 npt  = -npt
IF (fphase /= 3) GO TO 954

!     COMPUTE MAGNITUDE/PHASE FOR FORCES

CALL magpha (bufb(npt),isavef(npt))
955 iout = iout + 1
elwork(iout) = bufb(npt)
i    = i + 1
GO TO 951
954 IF (npt <= lforce) GO TO 955
npt  = npt - lforce
iout = iout + 1
elwork(iout) = isavef(npt)
i    = i + 1
GO TO 951

!     TRANSFER RESULTS TO BUFB

953 DO  i = 1,iout
  bufb(i) = elwork(i)
END DO
nwdfor  = iout

!     WRITE FORCES


!     DETERMINE DESTINATION FOR FORCE ENTRY

870 IF (force == 0) GO TO 880
IF (.NOT. ok2wrt) GO TO 880
id = bufb(1)
bufb(1) = 10*id + fdest
IF (xsetnf < 0.0) THEN
  GO TO   878
ELSE IF (xsetnf == 0.0) THEN
  GO TO   871
ELSE
  GO TO   872
END IF
871 bufb(1) = 10*id
GO TO 878
872 ix = ixsetf
873 IF (ix == nxsetf) GO TO 874
IF (z(ix+1) > 0) GO TO 874
IF (id >= z(ix) .AND. id <= (-z(ix+1))) GO TO 878
ix = ix + 2
GO TO 875
874 IF (id == z(ix)) GO TO 878
ix = ix + 1
875 IF (ix <= nxsetf) GO TO 873
GO TO 871

!     NOW WRITE FORCE ENTRY

878 CALL WRITE (oef1,bufb(1),nwdfor,0)
bufb(1) = id
880 GO TO 900
890 estawd = estawd + nwdsa
900 IF (again) GO TO 903
IF (istore == 1) GO TO 1626
IF (ktype /= 1 .OR. (axic .AND. midvec /= 0)) ivec = isave
IF (axic .AND. midvec /= 0) sorc = isvsrc
IF (.NOT. axic) GO TO 905
IF (nelhar /= nharms) GO TO 905
903 IF (eltype == 35) CALL scone3 (again)
IF (eltype == 70) CALL strax3 (again)
IF (eltype == 71) CALL stpax3 (again)
nelhar = -1
GO TO 616
905 IF (nesta == 0) GO TO 80
IF (z(estawd+1) /= 0) GO TO 90

!     END OF ESTA FOR CURRENT ELEMENT TYPE

910 IF (.NOT. idstrs) GO TO 915
CALL WRITE (oes1,0,0,1)
915 IF (.NOT. idforc) GO TO 920
CALL WRITE (oef1,0,0,1)
920 IF (.NOT. idlyst) GO TO 925
CALL WRITE (oes1l,0,0,1)
925 IF (.NOT. idlyfr) GO TO 930
CALL WRITE (oef1l,0,0,1)
930 IF (nesta == 0) GO TO 20
940 estawd = estawd + 2
IF (estawd >= nesta) GO TO 950
eltype = z(estawd)
GO TO 25

!     END OF ESTA FILE HIT

950 CONTINUE
960 CONTINUE
ivec  = isave
icore = jcore
RETURN

!     INTERNAL SUBROUTINE FOR WRITING ID RECORDS TO OUTPUT FILES

630 DO  i = 1,50
  buf(i) = 0
END DO

!     IF THE ID IS BEING WRITTEN TO A FILE WITH COMPLEX DATA,
!     CHANGE THE NUMBER OF WORDS TO REFLECT THE ACTUAL COUNT
!     OF WORDS BEING PUT TOGETHER USING THE STRING OF NUMBERS
!     IN THE 'COMPLX' ARRAY.  (SEE FORTRAN LABELS 651 THRU 654
!     AND 951 THRU 954)

IF (ktype  == 1) GO TO 645
IF (jcmplx == 0) RETURN 1
jout = 0
i    = jcmplx
638 ncmplx = complx(i)
IF (ncmplx == 0) THEN
  GO TO   642
END IF
640 jout = jout + 1
i    = i + 1
GO TO 638
642 nwds = jout

!     CHECK FOR VON MISES STRESS REQUEST.  SET WORD 11 IF
!     REQUEST IS FOUND.

645 IF (andf(nstrop,1) /= 0) buf(11) = 1

SELECT CASE ( branch )
  CASE (    1)
    GO TO 650
  CASE (    2)
    GO TO 660
  CASE (    3)
    GO TO 650
  CASE (    4)
    GO TO 650
  CASE (    5)
    GO TO 670
  CASE (    6)
    GO TO 790
  CASE (    7)
    GO TO 650
  CASE (    8)
    GO TO 660
  CASE (    9)
    GO TO 660
  CASE (   10)
    GO TO 650
END SELECT

!     NORMAL STATICS OR DIFF.STIFF. PHASE 0 OR 1 OR BUCKLING PHASE 0.

650 buf(2) = ifltyp
ix     = icc + isload
buf(5) = z(icc+1)
buf(6) = 0
buf(7) = 0
buf(8) = z(ix)
IF (branch /= 10) GO TO 840
ix     = icc + ittl + 84
z(ix)  = platit(1)
z(ix+1)= platit(2)
z(ix+2)= platit(3)
CALL int2al (ugvvec-1,z(ix+3),platit(4))
GO TO 840

!     EIGENVALUES OR BUCKLING PHASE 1.

660 buf(2) = ifltyp + ktypex
buf(5) = z(jlist)
buf(6) = z(jlist+1)
buf(7) = z(jlist+2)
buf(8) = 0
GO TO 840

!     FREQUENCY RESPONSE.

670 ix     = icc + idload
buf(8) = z(ix)
buf(6) = 0
buf(7) = 0
buf(2) = ifltyp + ktypex
671 CONTINUE

!     FIRST TIME FOR THIS LOAD VECTOR ONLY - MATCH LIST OF

IF (kfrq /= 0) GO TO 740

!     USER REQUESTED FREQS WITH ACTUAL FREQS. MARK FOR
!     OUTPUT EACH ACTUAL FREQ WHICH IS CLOSEST TO USER REQUEST.

kfrq   = 1
ix     = icc + ifrout
fsetno = z(ix)
IF (fsetno <= 0) GO TO 690
ix     = icc + ilsym
isetnf = ix+z(ix) + 1
680 isetfr = isetnf + 2
nsetfr = z(isetnf+1) + isetfr - 1
IF (z(isetnf) == fsetno) GO TO 710
isetnf = nsetfr + 1
IF (isetnf < ivec) GO TO 680
fsetno = -1
690 DO  j = ilist,nlist,2
  z(j+1) = 1
END DO
GO TO 740
710 DO  i = isetfr,nsetfr
  k      = 0
  diff   = 1.e25
  bufr(1)= zz(i)
  DO  j = ilist,nlist,2
    IF (z(j+1) /= 0) CYCLE
    diff1  = ABS(zz(j) - bufr(1))
    IF (diff1 >= diff) CYCLE
    diff = diff1
    k    = j
  END DO
  IF (k /= 0) z(k+1) = 1
END DO

!     DETERMINE IF CURRENT FREQ IS MARKED FOR OUTPUT.

740 IF (z(jlist+1) == 0) GO TO 960
buf(5) = z(jlist)
GO TO 840

!     TRANSIENT RESPONSE.

790 buf(5) = z(jlist)
buf(2) = ifltyp
ix     = icc + idload
buf(8) = z(ix)
buf(6) = 0
buf(7) = 0
GO TO 671

!     WRITE ID RECORD ON OUTPUT FILE.
!     (FOR MORE DETAIL, SEE OES1 FILE IN PROGRAMMER MANUAL P.2.3-130)

840 buf(1) = device + 10*branch
buf(3) = eltype

!     CHECK FOR TRIA1, TRIA2, TRIA3, QUAD1, QUAD2, QUAD4  ELEMENTS

IF (eltype /= 6 .AND. eltype /= 17 .AND. eltype /= 18 .AND.  &
    eltype /= 19 .AND. eltype /= 64 .AND. eltype /= 83) GO TO 845

!     CHECK FOR STRAIN OPTION

IF (buf(2) == 5 .AND. strain) buf(2) = 21
845 buf(4) = z(icc+1)
IF (ddrmm) buf(4) = 9999
buf(9) = IABS(z(irecx+2))
IF (buf(9) == 1 .AND. ktype == 2) buf(9) = 2
buf(10) = nwds
CALL WRITE (ofile,buf(1),50,0)
ix = icc + ittl
CALL WRITE (ofile,z(ix),96,1)
ilogic(nlogic) = .true.
GO TO iretrn, (365,375)

!     ERRORS

980 n = 2
GO TO 1000
990 n = 3
GO TO 1000
1000 CALL mesage (n,FILE,nam)
RETURN 1

END SUBROUTINE sdr2e
