SUBROUTINE sma2a
     
!     THIS SUBROUTINE FORMERLY GENERATED THE MGG AND BGG MATRICES FOR
!     THE SMA2 MODULE.  THESE OPERATIONS ARE NOW PERFORMED IN THE EMG
!     AND EMA MODULES AND SMA2A IS RETAINED IN SKELETAL FORM TO PROVIDE
!     A VEHICLE FOR USER-SUPPLIED ELEMENTS.
 
 LOGICAL :: heat
 INTEGER :: iz(1),eor,clsrw,clsnrw,frowic,sysprt,tnrows, outrw
 DOUBLE PRECISION :: dz,dpword
 DIMENSION        inpvt(2),dz(1),NAME(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /BLANK /  wtmass,nomgg,nobgg,icmas,icmbar,icmrod,icmqd1,  &
     icmqd2,icmtr1,icmtr2,icmtub,icmqdp,icmtrp,icmtrb
 COMMON /sma2ht/  heat
 COMMON /system/  isys,isew1(53),iprec
 COMMON /sma2io/  ifcstm,ifmpt,ifdit,idum1,ifecpt,igecpt,ifgpct,  &
     iggpct,idum2,idum3,ifmgg,igmgg,ifbgg,igbgg,idum4,  &
     idum5,inrw,outrw,clsnrw,clsrw,neor,eor,mcbmgg(7), mcbbgg(7)
 COMMON /zzzzzz/  z(1)
 COMMON /sma2bk/  icstm,ncstm,igpct,ngpct,ipoint,npoint,i6x6m,  &
     n6x6m,i6x6b,n6x6b
 COMMON /sma2cl/  ioptb,bggind,npvt,lleft,frowic,lrowic,nrowsc,  &
     tnrows,jmax,nlinks,link(10),nogo
 COMMON /gpta1 /  nelems,last,incr,NE(1)
 COMMON /sma2et/  ecpt(200)
 COMMON /zblpkx/  dpword,dum(2),INDEX
 EQUIVALENCE      (z(1),iz(1),dz(1))
 DATA    NAME  /  4HSMA2,4HA   /
 
 ipr = iprec
 
!     READ THE FIRST TWO WORDS OF NEXT GPCT RECORD INTO INPVT(1).
!     INPVT(1) IS THE PIVOT POINT.  INPVT(1) .GT. 0 IMPLIES THE PIVOT
!     POINT IS A GRID POINT.  INPVT(1) .LT. 0 IMPLIES THE PIVOT POINT
!     IS A SCALAR POINT.  INPVT(2) IS THE NUMBER OF WORDS IN THE
!     REMAINDER OF THIS RECORD OF THE GPCT.
 
 10 CALL READ (*1000,*700,ifgpct,inpvt(1),2,neor,iflag)
 ngpct = inpvt(2)
 CALL READ (*1000,*3000,ifgpct,iz(igpct+1),ngpct,eor,iflag)
 
!     FROWIC IS THE FIRST ROW IN CORE. (1 .LE. FROWIC .LE. 6)
 
 frowic = 1
 
!     DECREMENT THE AMOUNT OF CORE REMAINING.
 
 left = lleft - 2*ngpct
 IF (left <= 0) GO TO 3003
 ipoint = igpct + ngpct
 npoint = ngpct
 i6x6m  = ipoint + npoint
 i6x6m  = (i6x6m-1)/2 + 2
 
!     CONSTRUCT THE POINTER TABLE, WHICH WILL ENABLE SUBROUTINE INSERT
!     TO ADD THE ELEMENT MASS AND/OR DAMPING MATRICES TO MGG AND/OR BGG.
 
 iz(ipoint+1) = 1
 i1 = 1
 i  = igpct
 j  = ipoint + 1
 30 i1 = i1 + 1
 IF (i1 > ngpct) GO TO 40
 i  = i + 1
 j  = j + 1
 inc = 6
 IF (iz(i) < 0) inc = 1
 iz(j) = iz(j-1) + inc
 GO TO 30
 
!     JMAX = THE NUMBER OF COLUMNS OF MGG THAT WILL BE GENERATED WITH
!     THE CURRENT GRID POINT.
 
 40 inc   = 5
 ilast = igpct  + ngpct
 jlast = ipoint + npoint
 IF (iz(ilast) < 0) inc = 0
 jmax = iz(jlast) + inc
 
!     TNROWS = THE TOTAL NUMBER OF ROWS OF THE MATRIX TO BE GENERATED
!              FOR THE CURRENT PIVOT POINT.
!     TNROWS = 6 IF THE CURRENT PIVOT POINT IS A GRID POINT.
!     TNROWS = 1 IF THE CURRENT PIVOT POINT IS A SCALAR POINT.
 
 tnrows = 6
 IF (inpvt(1) < 0) tnrows = 1
 
!     IF 2*TNROWS*JMAX .LT. LEFT THERE ARE NO SPILL LOGIC PROBLEMS FOR
!     THE MGG SINCE THE WHOLE DOUBLE PRECISION SUBMATRIX OF ORDER
!     TNROWS*JMAX CAN FIT IN CORE.
 
 itemp = tnrows*jmax
 IF (2*itemp < left) GO TO 80
 CALL mesage (30,86,inpvt)
 
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
 
!     ZERO OUT THE MGG SUBMATRIX IN CORE
 
 100 low = i6x6m + 1
 lim = i6x6m + jmax*nrowsc
 DO  i = low,lim
   dz(i) = 0.0D0
 END DO
 
!     CHECK TO SEE IF BGG MATRIX IS DESIRED.
 
 IF (ioptb == 0) GO TO 137
 
!     SINCE THE BGG MATRIX IS TO BE COMPUTED,DETERMINE WHETHER OR NOT IT
!     TOO CAN FIT IN CORE.
 
 IF (nrowsc /= tnrows) GO TO 120
 IF (4*tnrows*jmax < left) GO TO 130
 
!     OPEN A SCRATCH FILE FOR BGG
 
 120 CALL mesage (-8,0,NAME)
 
!     THIS CODE TO BE FILLED IN LATER
!     ===============================
 
 130 i6x6b = i6x6m + jmax*tnrows
 low = i6x6b + 1
 lim = i6x6b + jmax*tnrows
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
 IF (itype >= 53 .AND. itype <= 61) GO TO 165
 CALL page2 (-3)
 sysprt = isew1(1)
 WRITE  (sysprt,161) ufm,itype
 161 FORMAT (a23,' 2202, ELEMENT TYPE',i4,' NO LONGER SUPPORTED BY ',  &
     'SMA2 MODULE.', /5X, 'USE EMG AND EMA MODULES FOR ELEMENT MATRIX GENERATION')
 nogo = 1
 GO TO 1000
 165 CONTINUE
 
!     READ THE ECPT ENTRY FOR THE CURRENT TYPE INTO THE ECPT ARRAY. THE
!     NUMBER OF WORDS TO BE READ WILL BE NWORDS(ITYPE).
 
 idx = (itype-1)*incr
 CALL fread (ifecpt,ecpt,NE(idx+12),0)
 itemp = NE(idx+23)
 
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
 4983,   4984,   4985,   4986, &
!          CDUM5   CDUM6   CDUM7   CDUM8   CDUM9
!            57      58      59      60      61  &
 4987,   4988,   4989,   4990,   4991  ) , itypx
 
 
 4983 CALL mdum1
 GO TO 160
 4984 CALL mdum2
 GO TO 160
 4985 CALL mdum3
 GO TO 160
 4986 CALL mdum4
 GO TO 160
 4987 CALL mdum5
 GO TO 160
 4988 CALL mdum6
 GO TO 160
 4989 CALL mdum7
 GO TO 160
 4990 CALL mdum8
 GO TO 160
 4991 CALL mdum9
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
 
!     CHECK NOGO = 1 SKIP BLDPK
 
 525 IF (nogo == 1) GO TO 10
 
!     AT THIS POINT BLDPK THE NUMBER OF ROWS IN CORE UNTO THE MGG FILE.
 
 ASSIGN 580 TO iretrn
 
!     HEAT TRANSFER PROBLEM, SKIP MGG
 
 IF (heat) GO TO 580
 
 ifile = ifmgg
 imcb  = 1
 
!     MULTIPLY THE MASS MATRIX BY THE PARAMETER WTMASS IF IT IS NOT
!     UNITY
 
 IF (wtmass == 1.0) GO TO 530
 low = i6x6m + 1
 lim = i6x6m + jmax*nrowsc
 DO  i = low,lim
   dz(i) = dz(i)*wtmass
 END DO
 530 i1  = 0
 540 i2  = 0
 ibeg = i6x6m + i1*jmax
 CALL bldpk (2,ipr,ifile,0,0)
 550 i2 = i2 + 1
 IF (i2 > ngpct) GO TO 570
 jj = igpct + i2
 INDEX = IABS(iz(jj)) - 1
 lim = 6
 IF (iz(jj) < 0) lim = 1
 jjj = ipoint + i2
 kkk = ibeg + iz(jjj) - 1
 i3 = 0
 560 i3 = i3 + 1
 IF (i3 > lim) GO TO 550
 INDEX = INDEX + 1
 kkk = kkk + 1
 dpword = dz(kkk)
 IF (dpword /= 0.0D0) CALL zblpki
 GO TO 560
 570 CALL bldpkn (ifile,0,mcbmgg(imcb))
 i1 = i1 + 1
 IF (i1 < nrowsc) GO TO 540
 GO TO iretrn, (580,600)
 
!     IF THE BGG IS CALLED FOR BLDPK IT.
 
 580 IF (ioptb ==  0) GO TO 600
 IF (ioptb == -1) GO TO 590
 
!     THE BGG MATRIX IS IN CORE
 
 ASSIGN 600 TO iretrn
 i6x6m = i6x6b
 ifile = ifbgg
 imcb  = 8
 GO TO 530
 
!     HERE WE NEED LOGIC TO READ BGG FROM A SCRATCH FILE AND INSERT
 
 590 CONTINUE
 
!     TEST TO SEE IF THE LAST ROW IN CORE, LROWIC, = THE TOTAL NO. OF
!     ROWS TO BE COMPUTED, TNROWS.  IF IT IS, WE ARE DONE.  IF NOT, THE
!     ECPT MUST BE BACKSPACED.
 
 600 IF (lrowic == tnrows) GO TO 10
 CALL bckrec (ifecpt)
 frowic = frowic + nrowsc
 lrowic = lrowic + nrowsc
 GO TO 100
 
!     CHECK NOGO = 1 SKIP BLDPK
 
 700 IF (nogo == 1) GO TO 10
 
!     HERE WE HAVE A PIVOT POINT WITH NO ELEMENTS CONNECTED, SO THAT
!     NULL COLUMNS MUST BE OUTPUT ON THE MGG AND BGG FILES.
 
 lim = 6
 IF (inpvt(1) < 0) lim = 1
 DO  i = 1,lim
   IF (heat) GO TO 705
   CALL bldpk (2,ipr,ifmgg,0,0)
   CALL bldpkn (ifmgg,0,mcbmgg)
   705 IF (ioptb /= 1) CYCLE
   CALL bldpk (2,ipr,ifbgg,0,0)
   CALL bldpkn (ifbgg,0,mcbbgg)
 END DO
 CALL skprec (ifecpt,1)
 GO TO 10
 
!     RETURN SINCE AN EOF HAS BEEN HIT ON THE GPCT FILE
 
 1000 IF (nogo == 1) CALL mesage (-61,0,NAME)
 RETURN
 
!     ERROR RETURNS
 
 3000 ifile = ifgpct
 iparm = 3
 GO TO 4010
 3003 CALL mesage (-8,ifile,NAME)
 3025 ifile = ifecpt
 iparm = 2
 4010 CALL mesage (-iparm,ifile,NAME)
 CALL mesage (-30,87,itype)
 RETURN
 
END SUBROUTINE sma2a
