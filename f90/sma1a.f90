SUBROUTINE sma1a
     
!     THIS SUBROUTINE FORMERLY GENERATED THE KGG AND K4GG MATRICES FOR
!     THE SMA1 MODULE.  THESE OPERATIONS ARE NOW PERFORMED IN THE EMG
!     AND EMA MODULES AND SMA1A IS RETAINED IN SKELETAL FORM TO PROVIDE
!     A VEHICLE FOR USER-PROVIDED ELEMENTS.
 
 LOGICAL :: dodet,nogo,heat,noheat
 INTEGER :: iz(1),eor,clsrw,clsnrw,frowic,sysprt,tnrows, outrw,option
 DOUBLE PRECISION :: dz,dpword
 DIMENSION        inpvt(2),dz(1),NAME(2)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm,swm
 COMMON /BLANK /  nogenl,nok4gg,option(2)
 COMMON /system/  ksystm(65)
 COMMON /sma1ht/  heat
 COMMON /sma1io/  ifcstm,ifmpt,ifdit,idum1,ifecpt,igecpt,ifgpct,  &
     iggpct,ifgei,iggei,ifkgg,igkgg,if4gg,ig4gg,  &
     ifgpst,iggpst,inrw,outrw,clsnrw,clsrw,neor,eor, mcbkgg(7),mcb4gg(7)
 COMMON /zzzzzz/  z(1)
 COMMON /sma1bk/  icstm,ncstm,igpct,ngpct,ipoint,npoint,i6x6k,  &
     n6x6k,i6x64,n6x64
 COMMON /sma1cl/  iopt4,k4ggsw,npvt,lleft,frowic,lrowic,nrowsc,  &
     tnrows,jmax,nlinks,link(10),idetck,dodet,nogoo
 COMMON /gpta1 /  nelems,last,incr,NE(1)
 COMMON /sma1et/  ecpt(200)
 COMMON /zblpkx/  dpword,dum(2),INDEX
 EQUIVALENCE      (ksystm(2),sysprt),(ksystm(3),nogo),  &
     (ksystm(55),iprec),(z(1),iz(1),dz(1))
 DATA    NAME  /  4HSMA1, 4HA   /
 
!     FLAG FOR ERROR CHECK IF A NON-HEAT ELEMENT IS REFERENCED
!     IN A -HEAT- FORMULATION.
 
 noheat = .false.
 ipr = iprec
 
!     READ THE FIRST TWO WORDS OF NEXT GPCT RECORD INTO INPVT(1).
!     INPVT(1) IS THE PIVOT POINT.  INPVT(1) .GT. 0 IMPLIES THE PIVOT
!     POINT IS A GRID POINT.  INPVT(1) .LT. 0 IMPLIES THE PIVOT POINT IS
!     A SCALAR POINT.  INPVT(2) IS THE NUMBER OF WORDS IN THE REMAINDER
!     OF THIS RECORD OF THE GPCT.
 
 IF (nogo) WRITE (sysprt,5) swm
 5 FORMAT (a27,' 2055, NOGO FLAG IS ON AT ENTRY TO SMA1A AND IS ',  &
     'BEING TURNED OFF.')
 nogo   = .false.
 10 idetck = 0
 CALL READ (*1000,*700,ifgpct,inpvt(1),2,neor,iflag)
 ngpct = inpvt(2)
 CALL READ (*1000,*3000,ifgpct,iz(igpct+1),ngpct,eor,iflag)
 
!     FROWIC IS THE FIRST ROW IN CORE. (1 .LE. FROWIC .LE. 6)
 
 frowic = 1
 
!     DECREMENT THE AMOUNT OF CORE REMAINING.
 
 left = lleft - 2*ngpct
 IF (left <= 0) GO TO 3003
 ipoint = igpct + ngpct
 npoint = ngpct
 i6x6k  = ipoint + npoint
 i6x6k  = (i6x6k-1)/2 + 2
 
!     CONSTRUCT THE POINTER TABLE, WHICH WILL ENABLE SUBROUTINE SMA1B
!     TO ADD THE ELEMENT STRUCTURAL AND/OR DAMPING MATRICES TO KGG AND
!     K4GG.
 
 iz(ipoint+1) = 1
 i1  = 1
 i   = igpct
 j   = ipoint + 1
 30 i1  = i1 + 1
 IF (i1 > ngpct) GO TO 40
 i   = i + 1
 j   = j + 1
 inc = 6
 IF (iz(i) < 0) inc = 1
 iz(j) = iz(j-1) + inc
 GO TO 30
 
!     JMAX = THE NUMBER OF COLUMNS OF KGG THAT WILL BE GENERATED WITH
!     THE CURRENT GRID POINT.
 
 40 inc   = 5
 ilast = igpct  + ngpct
 jlast = ipoint + npoint
 IF (iz(ilast) < 0) inc = 0
 jmax  = iz(jlast) + inc
 
!     TNROWS = THE TOTAL NUMBER OF ROWS OF THE MATRIX TO BE GENERATED
!              FOR THE CURRENT PIVOT POINT.
!     TNROWS = 6 IF THE CURRENT PIVOT POINT IS A GRID POINT.
!     TNROWS = 1 IF THE CURRENT PIVOT POINT IS A SCALAR POINT.
 
 tnrows = 6
 IF (inpvt(1) < 0) tnrows = 1
 
!     IF 2*TNROWS*JMAX .LT. LEFT THERE ARE NO SPILL LOGIC PROBLEMS FOR
!     THE KGG SINCE THE WHOLE DOUBLE PRECISION SUBMATRIX OF ORDER TNROWS
!     X JMAX CAN FIT IN CORE.
 
 itemp = tnrows*jmax
 IF (2*itemp < left) GO TO 80
 NAME(2) = inpvt(1)
 CALL mesage (30,85,NAME)
 
!     THE WHOLE MATRIX CANNOT FIT IN CORE, DETERMINE HOW MANY ROWS CAN
!     FIT. IF TNROWS = 1, WE CAN DO NOTHING FURTHER.
 
 IF (tnrows == 1) GO TO 3003
 nrowsc = 3
 70 IF (2*nrowsc*jmax < left) GO TO 90
 nrowsc = nrowsc - 1
 IF (nrowsc == 0) CALL mesage (-8,0,NAME)
 GO TO 70
 80 nrowsc = tnrows
 90 frowic = 1
 
!     LROWIC IS THE LAST ROW IN CORE. (1 .LE. LROWIC .LE. 6)
 
 lrowic = frowic + nrowsc - 1
 
!     ZERO OUT THE KGG SUBMATRIX IN CORE
 
 100 low = i6x6k + 1
 lim = i6x6k + jmax*nrowsc
 DO  i = low,lim
   dz(i) = 0.0D0
 END DO
 
!     CHECK TO SEE IF THE K4GG MATRIX IS DESIRED.
 
 IF (iopt4 == 0) GO TO 137
 
!     SINCE THE K4GG MATRIX IS TO BE COMPUTED, DETERMINE IF IT TOO CAN
!     FIT INTO CORE
 
 IF (nrowsc /= tnrows) GO TO 120
 IF (4*tnrows*jmax < left)  GO TO 130
 
!     OPEN A SCRATCH FILE FOR K4GG.
 
 120 CALL mesage (-8,0,NAME)
 
!     THIS CODE TO BE FILLED IN LATER
!     ===============================
 
 130 i6x64 = i6x6k + jmax*tnrows
 low   = i6x64 + 1
 lim   = i6x64 + jmax*tnrows
 DO  i = low,lim
   dz(i) = 0.0D0
 END DO
 
!     INITIALIZE THE LINK VECTOR TO -1.
 
 137 DO  i = 1,nlinks
   link(i) = -1
 END DO
 
!     TURN FIRST PASS INDICATOR ON.
 
 150 ifirst = 1
 
!     READ THE 1ST WORD OF THE ECPT RECORD, THE PIVOT POINT, INTO NPVT.
 
 CALL fread (ifecpt,npvt,1,0)
 
!     READ THE NEXT ELEMENT TYPE INTO THE CELL ITYPE.
 
 160 CALL READ (*3025,*500,ifecpt,itype,1,neor,iflag)
 IF (itype >= 53 .OR. itype <= 61) GO TO 165
 CALL page2 (-3)
 WRITE  (sysprt,161) ufm,itype
 161 FORMAT (a23,' 2201, ELEMENT TYPE',i4,' NO LONGER SUPPORTED BY ',  &
     'SMA1 MODULE.', /5X, 'USE EMG AND EMA MODULES FOR ELEMENT MATRIX GENERATION')
 nogo = .true.
 GO TO 1000
 165 CONTINUE
 
!     READ THE ECPT ENTRY FOR THE CURRENT TYPE INTO THE ECPT ARRAY. THE
!     NUMBER OF WORDS TO BE READ WILL BE NWORDS(ITYPE).
 
 idx = (itype-1)*incr
 CALL fread (ifecpt,ecpt,NE(idx+12),0)
 itemp = NE(idx+22)
 
!     IF THIS IS THE 1ST ELEMENT READ ON THE CURRENT PASS OF THE ECPT
!     CHECK TO SEE IF THIS ELEMENT IS IN A LINK THAT HAS ALREADY BEEN
!     PROCESSED.
 
 IF (ifirst == 1) GO TO 170
 
!     THIS IS NOT THE FIRST PASS.  IF ITYPE(TH) ELEMENT ROUTINE IS IN
!     CORE, PROCESS IT.
 
 IF (itemp == lincor) GO TO 180
 
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
 itypx  = itype - 52
 
!     CALL THE PROPER ELEMENT ROUTINE.
 
 180 GO TO ( &
!                                  CDUM1   CDUM2   CDUM3   CDUM4
!                                    53      54      55      56  &
 467,    468,    469,    470, &
!          CDUM5   CDUM6   CDUM7   CDUM8   CDUM9
!            57      58      59      60      61  &
 471,    472,    473,    474,    475  ) , itypx
 
 
 467 CALL kdum1
 GO TO 160
 468 CALL kdum2
 GO TO 160
 469 CALL kdum3
 GO TO 160
 470 CALL kdum4
 GO TO 160
 471 CALL kdum5
 GO TO 160
 472 CALL kdum6
 GO TO 160
 473 CALL kdum7
 GO TO 160
 474 CALL kdum8
 GO TO 160
 475 CALL kdum9
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
 
 520 CALL bckrec (ifecpt)
 GO TO 150
 
!    CHECK NOGOO FLAG. IF 1 SKIP BKDPK AND PROCESS ANOTHER GRID POINT
!    FROM GPCT
 
 525 IF (nogoo == 1) GO TO 10
 
!     IF NO GENERAL ELEMENTS EXIST, CHECK FOR GRID POINT SINGULARITIES.
 
!WKBR IF (DODET) CALL DETCK (0)
 IF (dodet) CALL detck (0,npvt,ifgpst)
 
!     AT THIS POINT BLDPK THE NUMBER OF ROWS IN CORE UNTO THE KGG FILE.
 
 ASSIGN 580 TO iretrn
 ifile= ifkgg
 imcb = 1
 530 i1   = 0
 540 i2   = 0
 ibeg = i6x6k + i1*jmax
 CALL bldpk (2,ipr,ifile,0,0)
 550 i2  = i2 + 1
 IF (i2 > ngpct) GO TO 570
 jj  = igpct + i2
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
 570 CALL bldpkn (ifile,0,mcbkgg(imcb))
 i1 = i1 + 1
 IF (i1 < nrowsc) GO TO 540
 GO TO iretrn, (580,600)
 
!     IF THE K4GG IS CALLED FOR, BLDPK IT.
 
 580 IF (iopt4 ==  0) GO TO 600
 IF (iopt4 == -1) GO TO 590
 
!     THE K4GG MATRIX IS IN CORE.
 
 ASSIGN 600 TO iretrn
 i6x6k = i6x64
 ifile = if4gg
 imcb  = 8
 GO TO 530
 
!     HERE WE NEED LOGIC TO READ K4GG FROM A SCRATCH FILE AND INSERT.
 
 590 CONTINUE
 
!     TEST TO SEE IF THE LAST ROW IN CORE, LROWIC, = THE TOTAL NO. OF
!     ROWS TO BE COMPUTED, TNROWS.  IF IT IS, WE ARE DONE.  IF NOT, THE
!     ECPT MUST BE BACKSPACED.
 
 600 IF (lrowic == tnrows) GO TO 10
 CALL bckrec (ifecpt)
 frowic = frowic + nrowsc
 lrowic = lrowic + nrowsc
 GO TO 100
 
!     CHECK NOGOO = 1 SKIP BLDPK AND PROCESS ANOTHER RECORD
 
 700 IF (nogoo == 1) GO TO 10
 
!     HERE WE HAVE A PIVOT POINT WITH NO ELEMENTS CONNECTED, SO THAT
!     NULL COLUMNS MUST BE OUTPUT ON THE KGG AND K4GG FILES.  IF DODET
!     IS TRUE, CALL THE DETERMINANT CHECK ROUTINE TO WRITE SINGULARITY
!     INFORMATION.
 
 npvt = IABS(inpvt(1))
 IF (inpvt(1) > 0) GO TO 703
 lim  = 1
 ixx  = -1
 GO TO 706
 703 lim  = 6
 ixx  = 1
!WKBR  706 IF (DODET) CALL DETCK (IXX)
 706 IF (dodet) CALL detck (ixx,npvt,ifgpst)
 DO  i = 1,lim
   CALL bldpk (2,ipr,ifkgg,0,0)
   CALL bldpkn (ifkgg,0,mcbkgg)
   IF (iopt4 /= 1) CYCLE
   CALL bldpk (2,ipr,if4gg,0,0)
   CALL bldpkn (if4gg,0,mcb4gg)
 END DO
 CALL skprec (ifecpt,1)
 GO TO 10
 
!     RETURN SINCE AN EOF HAS BEEN HIT ON THE GPCT FILE
 
 1000 IF (.NOT.nogo .AND. nogoo == 0) RETURN
 iparm = -61
 GO TO 4010
 
!     ERROR RETURNS
 
 3000 ifile = ifgpct
 GO TO 4003
 3003 iparm = -8
 GO TO 4010
 3025 ifile = ifecpt
 iparm = -2
 GO TO 4010
 4003 iparm = -3
 4010 CALL mesage (iparm,ifile,NAME)
 CALL mesage (-30,87,itype)
 RETURN
END SUBROUTINE sma1a
