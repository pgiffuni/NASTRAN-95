SUBROUTINE pla42
     
!     THIS ROUTINE PROCESSES THE SCRATCH DATA BLOCK ECPTS, WHICH IS THE
!     ECPTNL DATA BLOCK APPENDED WITH THE PROPER DISPLACEMENT VECTOR
!     COMPONENTS, AND CREATES THE STIFFNESS MATRIX KGGNL AND THE UPDATED
!     ECPTNL, ECPTNL1.  ECPTNL1, NAMED ECPTO IN THIS ROUTINE, DOES NOT
!     CONTAIN DISPLACEMENT VECTOR COMPONENTS.
 
 INTEGER :: sysbuf,buffr1,buffr2,buffr3,cstm,ecpts,ecpto,gpct,  &
     dit,placnt,planos,setno,frowic,eor,clsrw,outrw, buffr4,plsetn,FILE,ecptot
 DOUBLE PRECISION :: dz,dpword,dddddd
 DIMENSION       dz(1),iz(1),inpvt(2),NAME(2),mcbkgg(7),p(4),  &
     ecptot(7),planos(2),ip(4),nwdsp2(40),tubsav(16)
 COMMON /BLANK / placnt,plsetn,plfact(2)
 COMMON /system/ sysbuf,iskpu(53),iprec
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 COMMON /zzzzzz/ z(1)
 COMMON /pla42s/ xxxxxx(325)
 COMMON /pla42d/ dddddd(300)
 COMMON /pla42e/ ecpt(100)
 COMMON /pla42c/ npvt,gamma,gammas,ipass,icstm,ncstm,igpct,ngpct,  &
     ipoint,npoint,i6x6k,n6x6k,cstm,mpt,ecpts,gpct,  &
     dit,kggnl,ecpto,inrw,outrw,eor,neor,clsrw,jmax,  &
     frowic,lrowic,nrowsc,nlinks,nwords(40),iovrly(40), link(40),nogo
 COMMON /zblpkx/ dpword,dum(2),INDEX
 COMMON /pla4es/ wordes(300)
 COMMON /pla4uv/ worduv(25)
 EQUIVALENCE     (z(1),iz(1),dz(1)) ,(p(1),ip(1))
 DATA    NAME  / 4HPLA4,4H2   /, planos / 1103,  11    /
 DATA    nwdsp2/ 20,   0,  19,   0,   0, 33,   0,   0,  27,  20,  &
     5*0, 32,  27,  32,  38,   0,  &
     13*0,     45,      6*0/
 
 
 DO  i = 1,40
   iovrly(i) = 1
 END DO
 
!     DETERMINE SIZE OF VARIABLE CORE AND SET UP BUFFERS
 
 izmax  = korsz(z)
 buffr1 = izmax  - sysbuf
 buffr2 = buffr1 - sysbuf
 buffr3 = buffr2 - sysbuf
 buffr4 = buffr3 - sysbuf
 leftt  = buffr4 - 1
 ipass  = placnt - 1
 ipr    = iprec
 
!     READ THE CSTM INTO CORE
 
 FILE  = cstm
 ncstm = 0
 icstm = 0
 CALL OPEN (*20,cstm,z(buffr1),inrw)
 CALL skprec (cstm,1)
 CALL READ (*9020,*10,cstm,z(icstm+1),leftt,eor,ncstm)
 CALL mesage (-8,0,NAME)
 10 leftt = leftt - ncstm
 
!     PRETRD SETS UP SUBSEQUENT CALLS TO TRANSD
 
 CALL pretrd (z(icstm+1),ncstm)
 CALL pretrs (z(icstm+1),ncstm)
 CALL CLOSE  (cstm,clsrw)
 20 imat = ncstm
 
!     SEARCH THE MPT FOR THE PLAFACT CARDS.
 
 FILE = mpt
 CALL preloc (*9010,z(buffr1-3),mpt)
 CALL locate (*9040,z(buffr1-3),planos,iflag)
 
!     FIND THE CORRECT PLA SET NO.
 
 30 CALL fread (mpt,setno,1,0)
 IF (setno == plsetn) GO TO 50
 40 CALL fread (mpt,nn,1,0)
 IF (nn == (-1)) GO TO 30
 GO TO 40
 
!     SKIP THE PROPER NO. OF WORDS ON THE PLFACT CARD SO THAT GAMMA AND
!     GAMMAS (GAMMA STAR) WILL BE CORRECTLY COMPUTED.
 
 50 IF (placnt <= 4) GO TO 60
 CALL fread (mpt,0,-(placnt-4),0)
 60 nwdsrd = 4
 IF (placnt < 4) nwdsrd = placnt
 CALL fread (mpt,p,nwdsrd,0)
 IF (placnt - 3 < 0.0) THEN
   GO TO    70
 ELSE IF (placnt - 3 == 0.0) THEN
   GO TO    80
 ELSE
   GO TO    90
 END IF
 70 gammas = 0.0
 plfact(1) = p(2) - p(1)
 gamma = plfact(1)/p(1)
 GO TO 100
 80 word = p(2) - p(1)
 plfact(1) = p(3) - p(2)
 gammas = word/p(1)
 gamma  = plfact(1)/word
 GO TO 100
 90 word = p(3) - p(2)
 plfact(1) = p(4) - p(3)
 gammas = word/(p(2)-p(1))
 gamma  = plfact(1)/word
 100 plfact(2) = 0.0
 CALL CLOSE (mpt,clsrw)
 
!     CALL PREMAT TO READ MPT AND DIT INTO CORE.  NOTE NEGATIVE FILE NO.
!     FOR DIT TO TRIGGER PLA FLAG IN SUBROUTINE PREMAT.
 
 CALL premat (z(imat+1),z(imat+1),z(buffr1),leftt,matcr,mpt,-dit)
 leftt = leftt - matcr
 igpct = ncstm + matcr
 
!     OPEN KGGNL, ECPTO, ECPTS, AND GPCT
 
 ifile = kggnl
 CALL gopen  (kggnl,z(buffr1),1)
 CALL makmcb (mcbkgg,kggnl,0,6,ipr)
 CALL gopen  (ecpto,z(buffr2),1)
 CALL makmcb (ecptot,ecpto,0,0,0)
 CALL gopen  (ecpts,z(buffr3),0)
 CALL gopen  (gpct,z(buffr4),0)
 
!     READ THE FIRST TWO WORDS OF NEXT GPCT RECORD INTO INPVT(1).
!     INPVT(1) IS THE PIVOT POINT.  INPVT(1) .GT. 0 IMPLIES THE PIVOT
!     POINT IS A GRID POINT.  INPVT(1) .LT. 0 IMPLIES THE PIVOT POINT
!     IS A SCALAR POINT.  INPVT(2) IS THE NUMBER OF WORDS IN THE
!     REMAINDER OF THIS RECORD OF THE GPCT.
 
 130 FILE = gpct
 CALL READ (*1000,*700,gpct,inpvt(1),2,neor,iflag)
 ngpct = inpvt(2)
 CALL fread (gpct,iz(igpct+1),ngpct,1)
 IF (inpvt(1) < 0) GO TO 700
 
!     FROWIC IS THE FIRST ROW IN CORE. (1 .LE. FROWIC .LE. 6)
 
 frowic = 1
 
!     DECREMENT THE AMOUNT OF CORE REMAINING.
 
 left = leftt - 2*ngpct
 IF (left <= 0) CALL mesage (-8,0,NAME)
 ipoint = igpct + ngpct
 npoint = ngpct
 i6x6k  = ipoint + npoint
 i6x6k  = (i6x6k - 1)/2 + 2
 
!     CONSTRUCT THE POINTER TABLE, WHICH WILL ENABLE SUBROUTINE PLA4B TO
!     INSERT THE 6 X 6 MATRICES INTO KGGNL.
 
 iz(ipoint+1) = 1
 i1 = 1
 i  = igpct
 j  = ipoint + 1
 140 i1 = i1 + 1
 IF (i1 > ngpct) GO TO 150
 i  = i + 1
 j  = j + 1
 inc= 6
 IF (iz(i) < 0) inc = 1
 iz(j) = iz(j-1) + inc
 GO TO 140
 
!     JMAX = NO. OF COLUMNS OF KGGNL THAT WILL BE GENERATED WITH THE
!     CURRENT GRID POINT.
 
 150 inc   = 5
 ilast = igpct  + ngpct
 jlast = ipoint + npoint
 IF (iz(ilast) < 0) inc = 0
 jmax  = iz(jlast) + inc
 
!     IF 2*6*JMAX .LT. LEFT, THERE ARE NO SPILL LOGIC PROBLEMS FOR KGGNL
!     SINCE THE WHOLE DOUBLE PRECISION SUBMATRIX OF ORDER 6 X JMAX CAN
!     FIT IN CORE.
 
 itemp = 6*jmax
 IF (2*itemp < left) GO TO 170
 NAME(2) = inpvt(1)
 CALL mesage (30,85,NAME)
 nrowsc = 3
 160 IF (2*nrowsc*jmax < left) GO TO 180
 nrowsc = nrowsc - 1
 IF (nrowsc == 0) CALL mesage (-8,0,NAME)
 GO TO 160
 170 nrowsc = 6
 
!     LROWIC IS THE LAST ROW IN CORE. (1 .LE. LROWIC .LE. 6)
 
 180 lrowic = frowic + nrowsc - 1
 
!     ZERO OUT THE KGGD SUBMATRIX IN CORE.
 
 185 low = i6x6k + 1
 lim = i6x6k + jmax*nrowsc
 DO  i = low,lim
   dz(i) = 0.0D0
 END DO
 
!     INITIALIZE THE LINK VECTOR TO -1.
 
 DO  i = 1,nlinks
   link(i) = -1
 END DO
 
!     TURN FIRST PASS INDICATOR ON.
 
 ifirst = 1
 
!     READ THE 1ST WORD OF THE ECPT RECORD, THE PIVOT POINT, INTO NPVT.
!     IF NPVT .LT. 0, THE REMAINDER OF THE ECPT RECORD IS NULL SO THAT
!     1 OR 6 NULL COLUMNS MUST BE GENERATED
 
 FILE = ecpts
 CALL fread (ecpts,npvt,1,0)
 IF (npvt < 0) GO TO 700
 
!     WRITE PIVOT POINT ON ECPTNL1 (ECPTO)
 
 CALL WRITE (ecpto,npvt,1,neor)
 
!     READ THE NEXT ELEMENT TYPE INTO THE CELL ITYPE.
 
 220 CALL READ (*9020,*500,ecpts,itype,1,neor,iflag)
 
!     READ THE ECPT ENTRY FOR THE CURRENT TYPE INTO THE ECPT ARRAY. THE
!     NUMBER OF WORDS TO BE READ WILL BE NWORDS(ITYPE).
 
 IF (nwords(itype) <= 0) CALL mesage (-30,61,NAME)
 CALL fread (ecpts,ecpt,nwords(itype),0)
 itemp = iovrly(itype)
 
!     IF THIS IS THE 1ST ELEMENT READ ON THE CURRENT PASS OF THE ECPT
!     CHECK TO SEE IF THIS ELEMENT IS IN A LINK THAT HAS ALREADY BEEN
!     PROCESSED.
 
 IF (ifirst == 1) GO TO 230
 
!     THIS IS NOT THE FIRST PASS.  IF ITYPE(TH) ELEMENT ROUTINE IS IN
!     CORE, PROCESS IT.
 
 IF (itemp == lincor) GO TO 235
 
!     THE ITYPE(TH) ELEMENT ROUTINE IS NOT IN CORE.  IF THIS ELEMENT
!     ROUTINE IS IN A LINK THAT ALREADY HAS BEEN PROCESSED READ THE NEXT
!     ELEMENT.
 
 IF (link(itemp) == 1) GO TO 220
 
!     SET A TO BE PROCESSED LATER FLAG FOR THE LINK IN WHICH THE ELEMENT
!     RESIDES
 
 link(itemp) = 0
 GO TO 220
 
!     SINCE THIS IS THE FIRST ELEMENT TYPE TO BE PROCESSED ON THIS PASS
!     OF THE ECPT RECORD, A CHECK MUST BE MADE TO SEE IF THIS ELEMENT
!     IS IN A LINK THAT HAS ALREADY BEEN PROCESSED.  IF IT IS SUCH AN
!     ELEMENT, WE KEEP IFIRST = 1 AND READ THE NEXT ELEMENT.
 
 230 IF (link(itemp) == 1) GO TO 220
 
!     SET THE CURRENT LINK IN CORE = ITEMP AND IFIRST = 0
 
 lincor = itemp
 ifirst = 0
 
!     CALL THE PROPER ELEMENT ROUTINE.
 
!                     ROD      BEAM      TUBE     SHEAR     TWIST
!                       1         2         3         4         5
 235 GO TO   (       240,      999,      250,      999,      999, &
!                   TRIA1     TRBSC     TRPLT     TRMEM    CONROD
!                       6         7         8         9        10  &
 260,      999,      999,      270,      240, &
!                   ELAS1     ELAS2     ELAS3     ELAS4     QDPLT
!                      11        12        13        14        15  &
 999,      999,      999,      999,      999, &
!                   QDMEM     TRIA2     QUAD2     QUAD1     DAMP1
!                      16        17        18        19        20  &
 280,      290,      300,      310,      999, &
!                   DAMP2     DAMP3     DAMP4      VISC     MASS1
!                      21        22        23        24        25  &
 999,      999,      999,      999,      999, &
!                   MASS2     MASS3     MASS4     CONM1     CONM2
!                      26        27        28        29        30  &
 999,      999,      999,      999,      999, &
!                  PLOTEL     REACT     QUAD3       BAR      CONE
!                      31        32        33        34        35  &
 999,      999,      999,      320,      999, &
!                   TRIARG    TRAPRG    CTORDRG    CORE      CAP
!                      36        37        38        39        40  &
 999,      999,      999,      999,     999),itype
 
!     ROD, CONROD
 
 240 CALL pkrod
 
!     IF THE ELEMENT IS A TUBE, RESTORE THE SAVED ECPTNL ENTRY AND STORE
!     THE UPDATED VARIABLES IN PROPER SLOTS.
 
 IF (itype /= 3) GO TO 400
 DO  i = 1,16
   ecpt(i)  = tubsav(i)
 END DO
 ecpt(17) = ecpt(18)
 ecpt(18) = ecpt(19)
 ecpt(19) = ecpt(20)
 GO TO 400
 
!     THIS IS A TUBE ELEMENT.  REARRANGE THE ECPT FOR THE TUBE SO THAT
!     IT IS IDENTICAL TO THE ONE FOR THE ROD.
 
!     SAVE THE ECPT ENTRY FOR THE TUBE EXCEPT FOR THE 3 WORDS WHICH WILL
!     BE UPDATED BY THE PKROD ROUTINE AND THE TRANSLATIONAL COMPONENTS
!     OF THE DISPLACEMENTS VECTORS.
 
 250 DO  i = 1,16
   tubsav(i) = ecpt(i)
 END DO
 
!     COMPUTE AREA, TORSIONAL INERTIA TERM AND STRESS COEFFICIENT.
 
 d = ecpt(5)
 t = ecpt(6)
 dmt = d - t
 a = dmt*t*pi
 fj= .25*a*(dmt**2 + t**2)
 c = d/2.0
 
!     MOVE THE END OF THE ECPT ARRAY DOWN ONE SLOT SO THAT ENTRIES 7
!     THROUGH  25 WILL BE MOVED TO POSITIONS 8 THROUGH 26.
 
 m = 26
 DO  i = 1,19
   ecpt(m) = ecpt(m-1)
   m = m - 1
 END DO
 ecpt(5) = a
 ecpt(6) = fj
 ecpt(7) = c
 GO TO 240
 
!     TRIA1
 
 260 CALL pktri1
 GO TO 400
 
!     TRMEM
 
 270 CALL pktrm
 GO TO 400
 
!     QDMEM
 
 280 CALL pkqdm
 GO TO 400
 
!     TRIA2
 
 290 CALL pktri2
 GO TO 400
 
!     QUAD2
 
 300 CALL pkqad2
 GO TO 400
 
!     QUAD1
 
 310 CALL pkqad1
 GO TO 400
 
!     BAR
 
 320 CALL pkbar
 
!     WRITE ELEMENT TYPE AND UPDATED ECPT ENTRY ONTO ECPTNL1 (ECPTO)
 
 400 CALL WRITE (ecpto,itype,1,neor)
 CALL WRITE (ecpto,ecpt,nwdsp2(itype),neor)
 ecptot(2) = ecptot(2) + 1
 GO TO 220
 
!     AT STATEMENT NO. 500 WE HAVE HIT AN EOR ON THE ECPT FILE.  SEARCH
!     THE LINK VECTOR TO DETERMINE IF THERE ARE LINKS TO BE PROCESSED.
 
 500 link(lincor) = 1
 DO   i = 1,nlinks
   IF (link(i) == 0) GO TO 520
 END DO
 GO TO 525
 
!     SINCE AT LEAST ONE LINK HAS NOT BEEN PROCESSED THE ECPT FILE MUST
!     BE BACKSPACED.
 
 520 CALL bckrec (ecpts)
 GO TO 150
 525 IF (nogo == 1) CALL mesage (-61,0,0)
 
!     AT THIS POINT BLDPK THE NUMBER OF ROWS IN CORE ONTO THE KGGNL FILE
 
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
 
!     WRITE AN EOR ON ECPTO
 
 CALL WRITE (ecpto,0,0,eor)
 
!     TEST TO SEE IF THE LAST ROW IN CORE, LROWIC, = THE TOTAL NO. OF
!     ROWS TO BE COMPUTED = 6.  IF IT IS, WE ARE DONE.  IF NOT, THE
!     ECPTS MUST BE BACKSPACED.
 
 IF (lrowic == 6) GO TO 130
 CALL bckrec (ecpts)
 frowic = frowic + nrowsc
 lrowic = lrowic + nrowsc
 GO TO 185
 700 IF (nogo == 1) CALL mesage (-61,0,0)
 
!     HERE WE HAVE A PIVOT POINT WITH NO ELEMENTS CONNECTED, SO THAT
!     NULL COLUMNS MUST BE OUTPUT ON THE KGGD FILE.
 
 FILE = ecpts
 lim  = 6
 IF (inpvt(1) < 0) lim = 1
 DO  i = 1,lim
   CALL bldpk  (2,ipr,ifile,0,0)
   CALL bldpkn (kggnl,0,mcbkgg)
 END DO
 CALL skprec (ecpts,1)
 
!     WRITE PIVOT POINT ON ECPTO
 
 CALL WRITE (ecpto,npvt,1,eor)
 GO TO 130
 
!     CHECK NOGO FLAG. IF NOGO = 1, TERMINATE EXECUTION
 
 1000 IF (nogo == 1) CALL mesage (-61,0,0)
 
!     WRAP UP BEFORE RETURN
 
 CALL CLOSE (ecpts,clsrw)
 CALL CLOSE (ecpto,clsrw)
 CALL CLOSE (gpct,clsrw)
 CALL CLOSE (kggnl,clsrw)
 mcbkgg(3) = mcbkgg(2)
 CALL wrttrl (mcbkgg)
 CALL wrttrl (ecptot)
 RETURN
 
!     ERROR RETURNS
 
 9010 CALL mesage (-1,FILE,NAME)
 9020 CALL mesage (-2,FILE,NAME)
 9040 CALL mesage (-4,FILE,NAME)
 999 CALL mesage (-30,92,itype)
 RETURN
END SUBROUTINE pla42
