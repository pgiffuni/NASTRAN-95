SUBROUTINE ds1a
     
!     THIS ROUTINE GENERATES THE MATRIX KGGD WHICH IS THE SECOND ORDER
!     APPROXIMATION TO THE STIFFNESS MATRIX KGG.
 
 INTEGER :: eor,clsrw,outrw,frowic,cstm,dit,ecptds,gpct,  &
     buffr1,buffr2,buffr3,FILE,bar,beam,itypi(20)
 DOUBLE PRECISION :: dz(1),dpword,dddddd
 DIMENSION        ndum(9),iz(1),inpvt(2),NAME(2),mcbkgg(7)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /BLANK /  icom
 COMMON /system/  ksystm(100)
 COMMON /zzzzzz/  z(1)
 COMMON /ds1adp/  dddddd(300)
 COMMON /ds1aet/  ecpt(112)
 COMMON /ds1aaa/  npvt,icstm,ncstm,igpct,ngpct,ipoint,npoint,  &
     i6x6k,n6x6k,cstm,mpt,dit,ecptds,gpct,kggd,  &
     inrw,outrw,eor,neor,clsrw,jmax,frowic,lrowic, nrowsc,nlinks,link(10),nogo
 COMMON /gpta1 /  nelems,last,incr,NE(1)
 COMMON /zblpkx/  dpword,dum(2),INDEX
 COMMON /condas/  pi,twopi,radeg,degra,s4pisq
 EQUIVALENCE      (ksystm( 1),   isys), (ksystm( 2),ioutpt),  &
     (ksystm(46),ndum(1)), (ksystm(55), iprec)
 EQUIVALENCE      (z(1),iz(1),dz(1))
 DATA    NAME  /  4HDS1A,4H    /, bar,beam / 4HBAR ,4HBEAM /
 
!     DEFINE VARIABLES IN COMMON /DS1AAA/
 
 cstm  = 106
 mpt   = 107
 dit   = 110
 ecptds= 301
 gpct  = 109
 kggd  = 201
 inrw  = 0
 outrw = 1
 eor   = 1
 neor  = 0
 clsrw = 1
 nlinks= 10
 nogo  = 0
 iityp = 0
 j     = 26
 NE(j) = bar
 CALL sswtch (38,l38)
 
!     DETERMINE SIZE OF VARIABLE CORE, AND SET UP BUFFERS
 
 ipr = iprec
 CALL delset
 izmax  = korsz(z)
 buffr1 = izmax  - isys
 buffr2 = buffr1 - isys
 buffr3 = buffr2 - isys
 leftt  = buffr3 - 1
 
!     READ THE CSTM INTO CORE
 
 ifile = cstm
 ncstm = 0
 icstm = 0
 CALL OPEN (*2,cstm,z(buffr1),inrw)
 CALL fwdrec (*9020,cstm)
 CALL READ (*9030,*1,cstm,z(icstm+1),leftt,eor,ncstm)
 CALL mesage (-8,0,NAME)
 1 leftt = leftt - ncstm
 
!     PRETRD SETS UP SUBSEQUENT CALLS TO TRANSD.
 
 CALL pretrs (z(icstm+1),ncstm)
 CALL pretrd (z(icstm+1),ncstm)
 CALL CLOSE  (cstm,clsrw)
 2 imat1 = ncstm
 
!     CALL PREMAT TO READ MPT AND DIT INTO CORE.
 
 CALL premat (z(imat1+1),z(imat1+1),z(buffr1),leftt,matcr,mpt,dit)
 leftt = leftt - matcr
 igpct = ncstm + matcr
 
!     OPEN KGGD, ECPTDS AND GPCT
 
 CALL gopen  (kggd,z(buffr1),outrw)
 CALL makmcb (mcbkgg,kggd,0,6,ipr)
 CALL gopen  (ecptds,z(buffr2),inrw)
 CALL gopen  (gpct,z(buffr3),inrw)
 
!     READ THE FIRST TWO WORDS OF NEXT GPCT RECORD INTO INPVT(1).
!     INPVT(1) IS THE PIVOT POINT.  INPVT(1) .GT. 0 IMPLIES THE PIVOT
!     POINT IS A GRID POINT.  INPVT(1) .LT. 0 IMPLIES THE PIVOT POINT IS
!     A SCALAR POINT.  INPVT(2) IS THE NUMBER OF WORDS IN THE REMAINDER
!     OF THIS RECORD OF THE GPCT.
 
 10 FILE = gpct
 CALL READ (*1000,*700,gpct,inpvt(1),2,neor,iflag)
 ngpct = inpvt(2)
 CALL fread (gpct,iz(igpct+1),ngpct,eor)
 IF (inpvt(1) < 0) GO TO 700
 
!     FROWIC IS THE FIRST ROW IN CORE. (1 .LE. FROWIC .LE. 6)
 
 frowic = 1
 
!     DECREMENT THE AMOUNT OF CORE REMAINING.
 
 left   = leftt - 2*ngpct
 IF (left <= 0) CALL mesage (-8,0,NAME)
 ipoint = igpct + ngpct
 npoint = ngpct
 i6x6k  = ipoint + npoint
 i6x6k  = (i6x6k - 1)/2 + 2
 
!     CONSTRUCT THE POINTER TABLE, WHICH WILL ENABLE SUBROUTINE DS1B TO
!     INSERT THE 6 X 6 MATRICES INTO KGGD.
 
 iz(ipoint+1) = 1
 i1 = 1
 i  = igpct
 j  = ipoint + 1
 30 i1 = i1 + 1
 IF (i1 > ngpct) GO TO 40
 i  = i + 1
 j  = j + 1
 inc= 6
 IF (iz(i) < 0) inc = 1
 iz(j) = iz(j-1) + inc
 GO TO 30
 
!     JMAX = NO. OF COLUMNS OF KGGD THAT WILL BE GENERATED WITH THE
!     CURRENT GRID POINT.
 
 40 inc   = 5
 ilast = igpct  + ngpct
 jlast = ipoint + npoint
 IF (iz(ilast) < 0) inc = 0
 jmax  = iz(jlast) + inc
 
!     IF 2*6*JMAX .LT. LEFT THERE ARE NO SPILL LOGIC PROBLEMS FOR
!     KGGD SINCE THE WHOLE DOUBLE PRECISION SUBMATRIX OF ORDER 6 X JMAX
!     CAN FIT IN CORE.
 
 itemp = 6*jmax
 IF (2*itemp < left) GO TO 80
 NAME(2) = inpvt(1)
 CALL mesage (30,85,NAME)
 nrowsc = 3
 70 IF (2*nrowsc*jmax < left) GO TO 90
 nrowsc = nrowsc - 1
 IF (nrowsc == 0) CALL mesage (-8,0,NAME)
 GO TO 70
 80 nrowsc = 6
 
!     LROWIC IS THE LAST ROW IN CORE. (1 .LE. LROWIC .LE. 6)
 
 90 lrowic = frowic + nrowsc - 1
 
!     ZERO OUT THE KGGD SUBMATRIX IN CORE.
 
 100 low = i6x6k + 1
 lim = i6x6k + jmax*nrowsc
 DO  i = low,lim
   dz(i) = 0.0D0
 END DO
 
!     INITIALIZE THE LINK VECTOR TO -1.
 
 DO  i = 1,nlinks
   link(i) = -1
 END DO
 
!     TURN FIRST PASS INDICATOR ON.
 
 150 ifirst = 1
 
!     READ THE 1ST WORD OF THE ECPT RECORD, THE PIVOT POINT, INTO NPVT.
!     IF NPVT .LT. 0, THE REMAINDER OF THE ECPT RECORD IS NULL SO THAT
!     1 OR 6 NULL COLUMNS MUST BE GENERATED
 
 FILE = ecptds
 CALL fread (ecptds,npvt,1,neor)
 IF (npvt < 0) GO TO 700
 
!     READ THE NEXT ELEMENT TYPE INTO THE CELL ITYPE.
 
 160 CALL READ (*9020,*500,ecptds,itype,1,neor,iflag)
 
!     READ THE ECPT ENTRY FOR THE CURRENT TYPE INTO THE ECPT ARRAY. THE
!     NUMBER OF WORDS TO BE READ WILL BE NWORDS(ITYPE).
 
 ip = iprec
 IF (ip /= 1) ip = 0
 jtyp  = 2*itype - ip
 nfree = 3
 IF (itype == 2 .OR. itype == 35 .OR. itype == 75) nfree = 6
!               BEAM             CONEAX           TRSHL
 IF (itype >= 53 .AND. itype <= 61) nfree = MOD(ndum(itype-52),10)
!                DUM1              DUM9
 idx = (itype-1)*incr
 nwords = NE(idx+12) + 2 + nfree*NE(idx+10)
 IF (itype >= 65 .AND. itype <= 67) nwords = nwords + NE(idx+10) -1
!               IHEX1             IHEX3
 IF (itype == 80) nwords = nwords + NE(idx+10)
!                 IS2D8
 IF (itype == 35) nwords = nwords + 1
!                CONEAX
 IF (NE(idx+12) <= 0) CALL mesage (-61,0,NAME)
 CALL fread (ecptds,ecpt,nwords,neor)
 itemp = NE(idx+24)
 
!     IF THIS IS THE 1ST ELEMENT READ ON THE CURRENT PASS OF THE ECPT
!     CHECK TO SEE IF THIS ELEMENT IS IN A LINK THAT HAS ALREADY BEEN
!     PROCESSED.
 
 IF (ifirst == 1) GO TO 170
 
!     THIS IS NOT THE FIRST PASS.  IF ITYPE(TH) ELEMENT ROUTINE IS IN
!     CORE, PROCESS IT.
 
 IF (itemp == lincor) GO TO 171
 
!     THE ITYPE(TH) ELEMENT ROUTINE IS NOT IN CORE.  IF THIS ELEMENT
!     ROUTINE IS IN A LINK THAT ALREADY HAS BEEN PROCESSED READ THE NEXT
!     ELEMENT.
 
 IF (link(itemp) == 1) GO TO 160
 
!     SET A TO BE PROCESSED LATER FLAG FOR THE LINK IN WHICH THE ELEMENT
!     RESIDES
 
 link(itemp) = 0
 GO TO 160
 
!     SINCE THIS IS THE FIRST ELEMENT TYPE TO BE PROCESSED ON THIS PASS
!     OF THE ECPT RECORD, A CHECK MUST BE MADE TO SEE IF THIS ELEMENT
!     IS IN A LINK THAT HAS ALREADY BEEN PROCESSED.  IF IT IS SUCH AN
!     ELEMENT, WE KEEP IFIRST = 1 AND READ THE NEXT ELEMENT.
 
 170 IF (link(itemp) == 1) GO TO 160
 
!     SET THE CURRENT LINK IN CORE = ITEMP AND IFIRST = 0
 
 lincor = itemp
 ifirst = 0
 
!     CALL THE PROPER ELEMENT ROUTINE.
 
 171 IF (itype <= 0 .OR. itype > nelems) CALL mesage (-7,0,NAME)
 
!     IF DIAG 38 IS ON, ECHO TYPE OF ELEMENT BEING PROCESSED
 
 IF (l38   == 0) GO TO 180
 IF (iityp == 0) GO TO 175
 DO  ii = 1,iityp
   IF (itype == itypi(ii)) GO TO 180
 END DO
 IF (iityp >= 20) GO TO 180
 175 iityp = iityp + 1
 itypi(iityp) = itype
 WRITE  (ioutpt,177) NE(idx+1),NE(idx+2),itype
 177 FORMAT ('0*** DS1 MODULE PROCESSING ',2A4,' ELEMENTS (ELEM.TYPE', i4,1H))
 
 180 local = jtyp - 100
 IF (local > 0) THEN
   GO TO   182
 END IF
 181 GO TO ( & 
!        1-CROD       2-CBEAM      3-CTUBE      4-CSHEAR     5-CTWIST  &
 210,  210,   220,  220,   230,  230,   240,  240,  9040, 9040, & 
!        6-CTRIA1     7-CTRBSC     8-CTRPLT     9-CTRMEM     10-CONROD  &
 260,  260,  9040, 9040,  9040, 9040,   250,  250,   210,  210, & 
!        11-CELAS1    12-CELAS2    13-CELAS3    14-CELAS4    15-CQDPLT  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,& 
!        16-CQDMEM    17-CTRIA2    18-CQUAD2    19-CQUAD1    20-CDAMP1  &
 280,  280,   270,  270,   300,  300,   290,  290,  9040, 9040, & 
!        21-CDAMP2    22-CDAMP3    23-CDAMP4    24-CVISC     25-CMASS1  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040, & 
!        26-CMASS2    27-CMASS3    28-CMASS4    29-CONM1     30-CONM2  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040, & 
!        31-PLOTEL    32-X         33-X         34-CBAR      35-CCONEAX  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,   370,  370, & 
!        36-CTRIARG   37-CTRAPRG   38-CTORDRG   39-CTETRA    40-CWEDGE  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040, & 
!        41-CHEXA1    42-CHEXA2    43-CFLUID2   44-CFLUID3   45-CFLUID4  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040, & 
!        46-CFLMASS   47-CAXIF2    48-CAXIF3    49-CAXIF4    50-CSLOT3  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040 &
  ), jtyp
 
 182 GO TO ( & 
!        51-CSLOT4    52-CHBDY     53-CDUM1     54-CDUM2     55-CDUM3  &
 9040, 9040,  9040, 9040,   321,  321,   322,  322,   323,  323, & 
!        56-CDUM4     57-CDUM5     58-CDUM6     59-CDUM7     60-CDUM8  &
 324,  324,   325,  325,   326,  326,   327,  327,   328,  328, & 
!        61-CDUM9     62-CQDMEM1   63-CQDMEM2   64-CQUAD4    65-CIHEX1  &
 329,  329,  9040, 9040,  9040, 9040,   305,  305,   310,  310, & 
!        66-CIHEX2    67-CIHEX3    68-CQUADTS   69-CTRIATS   70-CTRIAAX  &
 310,  310,   310,  310,   311,  311,   312,  312,  9040, 9040, &
!        71-CTRAPAX   72-CAERO1    73-CTRIM6    74-CTRPLT1   75-CTRSHL  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,   313,  314, & 
!        76-CFHEX1    77-CFHEX2    78-CFTETRA   79-CFWEDGE   80-CIS2D8  &
 9040, 9040,  9040, 9040,  9040, 9040,  9040, 9040,   315,  315, & 
!        81-CELBOW    82-FTUBE     83-CTRIA3    84-CPSE2     85-CPSE3  &
 9040, 9040,  9040, 9040,   275,  275,   380,  380,   385,  385, & 
!        86-CPSE4  &
 390,  390 &
  ), local
 
!     ROD
 
 210 CALL drod
 GO TO 160
 
!     BAR
 
 220 CALL dbar
 GO TO 160
 
!     TUBE
 
 230 temp = ecpt(5) - ecpt(6)
 a  = temp*ecpt(6)*pi
 fj = .25*a*(temp**2 + ecpt(6)**2)
 c  = .5*ecpt(5)
 m  = 26
 DO  i = 1,18
   m = m - 1
   ecpt(m) = ecpt(m-1)
 END DO
 GO TO 210
 
!     SHEAR
 
 240 CALL dshear
 GO TO 160
 
!     TRMEM
 
 250 CALL dtrmem (0)
 GO TO 160
 
!     TRIA1
 
 260 CALL dtria (1)
 GO TO 160
 
!     TRIA2
 
 270 CALL dtria (2)
 GO TO 160
 
!     TRIA3
 
 275 CALL dtria (3)
 GO TO 160
 
!     QDMEM
 
 280 CALL dqdmem
 GO TO 160
 
!     QUAD1
 
 290 CALL dquad (1)
 GO TO 160
 
!     QUAD2
 
 300 CALL dquad (2)
 GO TO 160
 
!     QUAD4
 
 305 CALL dquad (4)
 GO TO 160
 
!     IHEX1,IHEX2,IHEX3
 
 310 CALL dihex (itype-64)
 GO TO 160
 
!     QUADTS
 
 311 CONTINUE
 GO TO 160
 
!    TRIATS
 
 312 CONTINUE
 GO TO 160
 313 CALL dtshls
 GO TO 160
 314 CALL dtshld
 GO TO 160
 315 CALL dis2d8
 GO TO 160
 
!     DUMMY ELEMENTS
 
 321 CALL ddum1
 GO TO 160
 322 CALL ddum2
 GO TO 160
 323 CALL ddum3
 GO TO 160
 324 CALL ddum4
 GO TO 160
 325 CALL ddum5
 GO TO 160
 326 CALL ddum6
 GO TO 160
 327 CALL ddum7
 GO TO 160
 328 CALL ddum8
 GO TO 160
 329 CALL ddum9
 GO TO 160
 
!     CONE
 
 370 CALL dcone
 GO TO 160
 
!     PRESSURE STIFFNESS ELEMENTS
 
 380 CALL dpse2
 GO TO 160
 385 CALL dpse3
 GO TO 160
 390 CALL dpse4
 GO TO 160
 
!     AT STATEMENT NO. 500 WE HAVE HIT AN EOR ON THE ECPT FILE.  SEARCH
!     THE LINK VECTOR TO DETERMINE IF THERE ARE LINKS TO BE PROCESSED.
 
 500 link(lincor) = 1
 DO   i = 1,nlinks
   IF (link(i) == 0) GO TO 520
 END DO
 GO TO 525
 
!     SINCE AT LEAST ONE LINK HAS NOT BEEN PROCESSED THE ECPT FILE MUST
!     BE BACKSPACED.
 
 520 CALL bckrec (ecptds)
 GO TO 150
 
!     CHECK NOGO FLAG. IF NOGO=1, SKIP BLDPK AND PROCESS ANOTHER RECORD
!     FROM THE GPCT TABLE
 
 525 IF (nogo == 1) GO TO 10
 
!     AT THIS POINT BLDPK THE NUMBER OF ROWS IN CORE UNTO THE KGG FILE.
 
 ifile = kggd
 i1 = 0
 540 i2 = 0
 ibeg = i6x6k + i1*jmax
 CALL bldpk (2,ipr,ifile,0,0)
 550 i2 = i2 + 1
 IF (i2 > ngpct) GO TO 570
 jj = igpct + i2
 INDEX = IABS(iz(jj)) - 1
 lim = 6
 IF (iz(jj) < 0) lim = 1
 jjj = ipoint + i2
 kkk = ibeg + iz(jjj) - 1
 i3  = 0
 560 i3  = i3 + 1
 IF (i3 > lim) GO TO 550
 INDEX = INDEX + 1
 kkk = kkk + 1
 dpword = dz(kkk)
 IF (dpword /= 0.0D0) CALL zblpki
 GO TO 560
 570 CALL bldpkn (ifile,0,mcbkgg)
 i1 = i1 + 1
 IF (i1 < nrowsc) GO TO 540
 
!     TEST TO SEE IF THE LAST ROW IN CORE, LROWIC, = THE TOTAL NO. OF
!     ROWS TO BE COMPUTED = 6.  IF IT IS, WE ARE DONE.  IF NOT, THE
!     ECPTDS MUST BE BACKSPACED.
 
 IF (lrowic == 6) GO TO 10
 CALL bckrec (ecptds)
 frowic = frowic + nrowsc
 lrowic = lrowic + nrowsc
 GO TO 100
 700 IF (nogo == 1) GO TO 10
 
!     HERE WE HAVE A PIVOT POINT WITH NO ELEMENTS CONNECTED, SO THAT
!     NULL COLUMNS MUST BE OUTPUT ON THE KGGD FILE.
 
 FILE = ecptds
 lim  = 6
 IF (inpvt(1) < 0) lim = 1
 DO  i = 1,lim
   CALL bldpk  (2,ipr,kggd,0,0)
   CALL bldpkn (kggd,0,mcbkgg)
 END DO
 CALL fwdrec (*9020,ecptds)
 GO TO 10
 
!     CHECK NOGO FLAG. IF NOGO=1, TERMINATE EXECUTION
 
 1000 IF (nogo == 1) CALL mesage (-61,0,0)
 
!     WRAP UP BEFORE RETURN
 
 CALL CLOSE (ecptds,clsrw)
 CALL CLOSE (gpct,clsrw)
 CALL CLOSE (kggd,clsrw)
 mcbkgg(3) = mcbkgg(2)
 IF (mcbkgg(6) == 0) GO TO 9050
 CALL wrttrl (mcbkgg)
 j = 26
 NE(j) = beam
 RETURN
 
!     ERROR RETURNS
 
 9020 CALL mesage (-2,FILE,NAME)
 9030 CALL mesage (-3,FILE,NAME)
 9040 CALL mesage (-7,FILE,NAME)
 9050 WRITE  (ioutpt,9060) ufm
 9060 FORMAT (a23,' 2402, NULL DIFFERENTIAL STIFFNESS MATRIX ',  &
     'GENERATED IN SUBROUTINE DS1A.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE ds1a
