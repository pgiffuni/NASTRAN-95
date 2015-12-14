SUBROUTINE pla1
     
!     THIS FUNCTIONAL MODULE IS THE FIRST OF FOUR FUNCTIONAL MODULES
!     UNIQUE TO THE PIECE-WISE LINEAR ANALYSIS (DISPLACEMENT METHOD)
!     RIGID FORMAT
 
 
!     PLA1   CSTM,MPT,ECPT,GPCT,DIT,CASECC,EST /KGGL,ECPTNL,ESTL,ESTNL/
!            V,N,KGGLPG/V,N,NPLALIM/V,N,ECPTNLPG/V,N,PLSETNO/
!            V,N,NONLSTR/V,N,PLFACT/ $
 
!     THE OUTPUT DATA BLOCKS AND PARAMETERS ARE DEFINED AS FOLLOWS -
 
!     KGGL   IS THE STIFFNESS MATRIX OF LINEAR (NON-STRESS DEPENDENT)
!            ELEMENTS
!     ECPTNL IS A SUBSET OF THE ECPT WHICH CONTAINS ECPT PLUS STRESS
!            INFORMATION FOR THE NON-LINEAR (STRESS DEPENDENT) ELEMENTS
!     ESTL,  A SUBSET OF THE EST, CONTAINS ALL LINEAR ELEMENTS
!     ESTNL, THE COMPLEMENT OF THE ESTL, CONTAINS INFORMATION FOR THE
!            NON-LINEAR ELEMENTS
 
!     PARAMETER NAMES BELOW ARE FORTRAN VARIABLE NAMES RATHER THAN DMAP
!     PARAMETER NAMES
 
!     KGGLPG IS THE PURGE FLAG FOR THE KGGL AND ESTL DATA BLOCKS.  IT IS
!            SET = -1 (PURGE=YES) IF ALL ELEMENTS ARE STRESS DEPENDENT
!     NPLALP IS THE NUMBER OF PASSES THAT WILL BE MADE THRU THE PLA LOOP
!     KICKOF IS SET = -1 (KICK THE USER OFF THE MACHINE = YES) IF THE
!            DIT IS PURGED OR ALL ELEMENTS ARE NON-STRESS DEPENDENT
!     PLASET IS THE SET NUMBER ON THE PLFACT CARD THAT IS OBTAINED FROM
!            THE FIRST RECORD OF CASECC
!     NONLST IS THE FLAG SUCH THAT IF IT IS A -1 THE USER DOES NOT WISH
!            TO OUTPUT HIS NON-LINEAR STRESSES.  HENCE PLA3 WILL NOT BE
!            CALLED
!     PLFACT IS THE FIRST PIECE-WISE LINEAR FACTOR TO BE USED FROM
!            PLASET
 
 LOGICAL :: phase1,all,heat
 INTEGER :: bufr1,bufr2,bufr3,bufr4,bufsz,iz(1),eor,clsrw,  &
     clsnrw,outrw,frowic,tnrows,cstm,dit,ecpt,gpct,  &
     ecptnl,estl,estnl,plaary(90),FILE,trail,plaset,  &
     planos(2),setno,outfil,iecpt(100),estltr(7), estnlt(7),yessd,est,casecc
 DOUBLE PRECISION :: dz(1),dpdum,dpword
 DIMENSION        NAME(2),inpvt(2)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm
 COMMON /BLANK /  kgglpg,nplalp,kickof,plaset,nonlst,plfact(2)
 COMMON /system/  bufsz,isyspt,isp1(37),nbpw,isp2(14),iprec
 COMMON /sma1io/  cstm,mpt,dit,bufr1,ecpt,bufr2,gpct,bufr3,bufr4,  &
     itype,kggl,estnl,ecptnl,dum1,estl,dum2,inrw,  &
     outrw,clsnrw,clsrw,neor,eor,mcbkgg(7),trail(7)
 COMMON /zzzzzz/  z(1)
 COMMON /sma1bk/  icstm,ncstm,igpct,ngpct,ipoint,npoint,i6x6k, dum4,dum5,dum6
 COMMON /sma1cl/  iopt4,dum7,npvt,left,frowic,lrowic,nrowsc,  &
     tnrows,jmax,nlinks,link(10),dum8,dum9,nogo
 COMMON /gpta1 /  nelems,last,incr,NE(1)
 COMMON /sma1et/  xecpt(100)
 COMMON /sma1dp/  dpdum(300)
 COMMON /sma1ht/  heat
 COMMON /zblpkx/  dpword,llllll(2),INDEX
 COMMON /matin /  matid,inflag,filler(4)
 COMMON /matout/  indstr,yyyyyy(19)
 EQUIVALENCE      (z(1),iz(1),dz(1)) ,(iecpt(1),xecpt(1)),  &
     (mcbkgg(1),estltr(1)), (trail(1),estnlt(1)),  &
     (trail(2),nnlel), (estltr(2),nlel), (fnn,nn), (indstr,e)
 DATA             casecc,jstset,jplset,jsym /106,23,164,200 /
 DATA             NAME  / 4HPLA1,4H    /, hmpt  / 4HMPT     /
 DATA             planos/ 1103, 11     /
 DATA             plaary/ 1,  0,  1,  0,  0,    1,  0,  0,  1,  1,  &
     0,  0,  0,  0,  0,    1,  1,  1,  1,  0,  &
     0,  0,  0,  0,  0,    0,  0,  0,  0,  0,  &
     0,  0,  0,  1,  0,    0,  0,  0,  0,  0,  &
     0,  0,  0,  0,  0,    0,  0,  0,  0,  0,  &
     0,  0,  0,  0,  0,    0,  0,  0,  0,  0, 30*0 /
 
!     IF THE DIT HAS BEEN PURGED, WE CANNOT PROCEED FURTHER
 
 ipr = iprec
 CALL delset
 trail(1) = dit
 CALL rdtrl (trail)
 IF (trail(1) < 0) CALL mesage (-1,dit,NAME)
 
!     INITIALIZE HEAT PARAMETER
 
 heat = .false.
 
!     INITIALIZE MODULE PARAMETERS AND SET IOPT4 = 0, SO THAT ELEMENT
!     ROUTINES WILL NOT CALCULATE STRUCTURAL DAMPING MATRICES.
 
 kgglpg =-1
 nplalp = 1
 kickof =-1
 plaset =-1
 nonlst = 1
 plfact(1) = 1.0
 iopt4  = 0
 estnl  = 204
 est    = 107
 phase1 = .true.
 ASSIGN 290 TO nosd
 
!     SET UP BUFFERS AND INITIALIZE FILE TRAILERS.
 
 izmax = korsz (z)
 bufr1 = izmax - bufsz
 bufr2 = bufr1 - bufsz
 bufr3 = bufr2 - bufsz
 bufr4 = bufr3 - bufsz
 left  = bufr4 - 1
 CALL makmcb (mcbkgg,kggl,0,6,ipr)
 CALL makmcb (trail,ecptnl,0,0,0)
 
!     CHECK PLAARY SIZE
 
 IF (nelems > 90) WRITE (isyspt,1) uwm
 1 FORMAT (a25,' 2151, -PLAARY- ARRAY IS SMALLER THAN MAXIMUM ',  &
     'NUMBER OF ELEMENT TYPES.')
 
!     OPEN THE KGGL FILE FOR OUTPUT
 
 CALL gopen (kggl,z(bufr1),1)
 
!     ATTEMPT TO READ THE CSTM INTO CORE
 
 icstm = 0
 ncstm = 0
 FILE  = cstm
 CALL OPEN (*20,cstm,z(bufr2),inrw)
 CALL fwdrec (*870,cstm)
 CALL READ (*870,*10,cstm,z(icstm+1),left,eor,ncstm)
 
!     INSUFFICIENT CORE - CALL MESAGE
 
 CALL mesage (-8,0,NAME)
 10 left = left - ncstm
 CALL CLOSE (cstm,clsrw)
 
!     CALL PRETRD TO SET UP FUTURE CALLS TO TRANSD.
 
 CALL pretrd (z(icstm+1),ncstm)
 CALL pretrs (z(icstm+1),ncstm)
 
!     CALL PREMAT TO READ THE MPT AND THE DIT AND TO SET UP FUTURE CALLS
!     TO SUBROUTINE MAT.  NOTE NEGATIVE ISGN FOR DIT TO TRIGGER PLA FLAG
!     IN MAT.
 
 20 imat  = ncstm
 CALL premat (iz(imat+1),z(imat+1),z(bufr2-3),left,mused,mpt,-dit)
 left  = left  - mused
 igpct = ncstm + mused
 
!     OPEN THE INPUT FILES ECPT AND GPCT AND THE OUTPUT FILE ECPTNL.
 
 CALL gopen (ecpt,  z(bufr2),0)
 CALL gopen (gpct,  z(bufr3),0)
 CALL gopen (ecptnl,z(bufr4),1)
 ileft = left
 
!     BEGIN MAIN LOOP FOR PROCESSING THE ECPT.
 
 30 CALL READ (*650,*630,gpct,inpvt(1),2,neor,iflag)
 ngpct = inpvt(2)
 left  = ileft - 2*ngpct
 IF (left <= 0) CALL mesage (-8,0,NAME)
 CALL fread (gpct,iz(igpct+1),ngpct,eor)
 
!     FROWIC IS THE FIRST ROW IN CORE (1 .LE. FROWIC .LE. 6)
 
 frowic = 1
 ipoint = igpct  + ngpct
 npoint = ngpct
 i6x6k  = ipoint + npoint
 
!     MAKE I6X6K A DOUBLE PRECISION INDEX (I6X6K POINTS TO THE 0TH
!     LOCATION OF THE 6 X 6 SUBMATRIX OF KGGL IN CORE)
 
 i6x6k = (i6x6k-1)/2 + 2
 
!     CONSTRUCT THE POINTER TABLE WHICH WILL ENABLE SUBROUTINE SMA1B TO
!     ADD THE ELEMENT STIFFNESS MATRICES TO KGGL.
 
 iz(ipoint+1) = 1
 i1 = 1
 i  = igpct
 j  = ipoint + 1
 40 i1 = i1 + 1
 IF (i1 > ngpct) GO TO 50
 i  = i + 1
 j  = j + 1
 inc= 6
 IF (iz(i) < 0) inc = 1
 iz(j) = iz(j-1) + inc
 GO TO 40
 
!     JMAX = THE NO. OF COLUMNS OF KGGL THAT WILL BE GENERATED WITH THE
!     CURRENT PIVOT POINT.
 
 50 inc   = 5
 ilast = igpct  + ngpct
 jlast = ipoint + npoint
 IF (iz(ilast) < 0) inc = 0
 jmax  = iz(jlast) + inc
 
!     TNROWS = TOTAL NO. OF ROWS OF THE MATRIX TO BE GENERATED
 
 tnrows = 6
 IF (inpvt(1) < 0) tnrows = 1
 IF (2*tnrows*jmax < left) GO TO 70
 
!     THE WHOLE SUBMATRIX CANNOT FIT IN CORE
 
 IF (tnrows == 1) CALL mesage (-8,0,NAME)
 nrowsc = 3
 plaary(39) = NAME(1)
 60 plaary(40) = npvt
 CALL mesage (30,85,plaary(39))
 IF (2*nrowsc*jmax < left) GO TO 80
 nrowsc = nrowsc - 1
 IF (nrowsc == 0) CALL mesage (-8,0,NAME)
 GO TO 60
 70 nrowsc = tnrows
 80 frowic = 1
 lrowic = frowic + nrowsc - 1
 
!     ZERO OUT THE KGGL SUBMATRIX IN CORE.
 
 low = i6x6k + 1
 lim = i6x6k + jmax*nrowsc
 DO  i = low,lim
   dz(i) = 0.0D0
 END DO
 
!     INITIALIZE THE LINK VECTOR TO -1
 
 DO  i = 1,nlinks
   link(i) = -1
 END DO
 lincor  =  1
 FILE = ecpt
 
!     TURN FIRST PASS, FIRST ELEMENT READ ON THE CURRENT PASS OF THE
!     ECPT RECORD, AND PIVOT POINT WRITTEN INDICATORS ON.
 
 ipass  = 1
 npvtwr = 0
 110 ifirst = 1
 
!     READ THE FIRST WORD OF THE ECPT RECORD, THE PIVOT POINT, INTO NPVT
 
 CALL fread (ecpt,npvt,1,neor)
 
!     READ THE NEXT ELEMENT TYPE INTO ITYPE, AND READ THE PRESCRIBED NO.
!     OF WORDS INTO THE XECPT ARRAY.
 
 120 CALL READ (*870,*520,ecpt,itype,1,neor,iflag)
 idx = (itype-1)*incr
 nn  = NE(idx+12)
 CALL fread (ecpt,xecpt,nn,neor)
 itemp = NE(idx+22)
 IF (ipass /= 1) GO TO 290
 
!     THIS IS THE FIRST PASS.  IF THE ELEMENT IS IN THE PLA SET, CALL
!     THE MAT ROUTINE TO FIND OUT IF ANY OF THE MATERIAL PROPERTIES IS
!     STRESS DEPENDENT.
 
 
!              CROD      CBEAM     CTUBE     CSHEAR    CTWIST    CTRIA1
!                 1        2         3          4         5         6
!              CTRBSC    CTRPLT    CTRMEM    CONROD    ELAS1     ELAS2
!               7           8         9        10        11        12
!              ELAS3     ELAS4     CQDPLT    CQDMEM    CTRIA2    CQUAD2
!               13        14         15        16        17        18
!              CQUAD1    CDAMP1    CDAMP2    CDAMP3    CDAMP4    CVISC
!               19        20        21         22        23        24
!              MASS1     CMASS2    CMASS3    CMASS4    CONM1     CONM2
!                25        26        27        28        29        30
!              PLOTEL    REACT     QUAD3     CBAR      CCONE
!                31        32        33       34         35
 IF (itype > 35) GO TO 120
 GO TO (130,890,150,290,290,160,290,290,170,130,  &
     290,290,290,290,290,180,190,200,210,120,  &
     120,120,120,120,120,120,120,120,120,120, 120,890,890,220,290), itype
 
!     ROD
 
 130 matid = iecpt(4)
 ASSIGN 140 TO yessd
 GO TO 230
 140 xecpt(18) = 0.0
 xecpt(19) = 0.0
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(20) = e
 IF (phase1) GO TO 145
 xecpt(21) = 0.0
 nwds = 21
 GO TO 765
 145 nwds = 20
 GO TO 260
 
!     TUBE
 
 150 matid = iecpt(4)
 ASSIGN 155 TO yessd
 GO TO 230
 155 xecpt(17) = 0.0
 xecpt(18) = 0.0
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(19) = e
 IF (phase1) GO TO 157
 xecpt(20) = 0.0
 nwds = 20
 GO TO 765
 157 nwds = 19
 GO TO 260
 
!     TRIA1
 
 160 matid = iecpt(6)
 ASSIGN 165 TO yessd
 GO TO 230
 165 DO  i = 28,33
   xecpt(i) = 0.0
 END DO
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(30) = e
 IF (phase1) GO TO 168
 DO  i = 34,38
   xecpt(i) = 0.0
 END DO
 nwds = 38
 GO TO 765
 168 nwds = 33
 GO TO 260
 
!     TRMEM
 
 170 matid = iecpt(6)
 ASSIGN 175 TO yessd
 GO TO 230
 175 DO  i = 22,27
   xecpt(i) = 0.0
 END DO
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(24) = e
 nwds = 27
 IF (phase1) GO TO 260
 GO TO 765
 
!     QDMEM
 
 180 matid = iecpt(7)
 ASSIGN 185 TO yessd
 GO TO 230
 185 DO  i = 27,32
   xecpt(i) = 0.0
 END DO
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(29) = e
 nwds = 32
 IF (phase1) GO TO 260
 GO TO 765
 
!     TRIA2
 
 190 matid = iecpt(6)
 ASSIGN 195 TO yessd
 GO TO 230
 195 DO  i = 22,27
   xecpt(i) = 0.0
 END DO
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(24) = e
 IF (phase1) GO TO 198
 DO  i = 28,32
   xecpt(i) = 0.0
 END DO
 nwds = 32
 GO TO 765
 198 nwds = 27
 GO TO 260
 
!     QUAD2
 
 200 matid = iecpt(7)
 ASSIGN 205 TO yessd
 GO TO 230
 205 DO  i = 27,32
   xecpt(i) = 0.0
 END DO
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(29) = e
 IF (phase1) GO TO 208
 DO  i = 33,37
   xecpt(i) = 0.0
 END DO
 nwds = 37
 GO TO 765
 208 nwds = 32
 GO TO 260
 
!     QUAD1
 
 210 matid = iecpt(7)
 ASSIGN 215 TO yessd
 GO TO 230
 215 DO  i = 33,38
   xecpt(i) =0.0
 END DO
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(35) = e
 IF (phase1) GO TO 218
 DO  i = 39,43
   xecpt(i) = 0.0
 END DO
 nwds = 43
 GO TO 765
 218 nwds = 38
 GO TO 260
 
!     BAR - IF COORDINATE 1 OF EITHER PT. A OR PT. B IS PINNED THE
!           ELEMENT IS TREATED AS LINEAR (NON-STRESS DEPENDENT)
 
 220 IF (iecpt(8) == 0 .AND. iecpt(9) == 0) GO TO 224
 ka = iecpt(8)
 kb = iecpt(9)
 222 IF (MOD(ka,10) == 1 .OR. MOD(kb,10) == 1) GO TO 290
 ka = ka/10
 kb = kb/10
 IF (ka <= 0 .AND. kb <= 0) GO TO 224
 GO TO 222
 224 matid = iecpt(16)
 ASSIGN 226 TO yessd
 GO TO 230
 226 xecpt(43) = 0.0
 xecpt(44) = 0.0
 inflag = 1
 CALL mat (iecpt(1))
 xecpt(45) = e
 IF (phase1) GO TO 228
 DO  i = 46,50
   xecpt(i) = 0.0
 END DO
 nwds = 50
 GO TO 765
 228 nwds = 45
 GO TO 260
 
!     TEST TO SEE IF ELEMENT IS STRESS DEPENDENT.
 
 230 inflag = 5
 CALL mat (iecpt(1))
 IF (indstr > 0) THEN
   GO TO   250
 END IF
 240 GO TO nosd,  (290,820)
 250 GO TO yessd, (140,155,165,175,185,195,205,215,226)
 
!     WRITE AN ENTRY ONTO ECPTNL
 
 260 IF (npvtwr > 0) THEN
   GO TO   280
 END IF
 270 npvtwr = 1
 CALL WRITE (ecptnl,npvt,1,neor)
 kickof = 1
 280 CALL WRITE (ecptnl,itype,1,neor)
 CALL WRITE (ecptnl,xecpt,nwds,neor)
 nnlel = nnlel + 1
 GO TO 120
 
!     IF THIS IS THE 1ST ELEMENT READ ON THE CURRENT PASS OF THE ECPT,
!     CHECK TO SEE IF THIS ELEMENT IS IN A LINK THAT HAS ALREADY BEEN
!     PROCESSED.
 
 290 kgglpg = 1
 IF (ifirst == 1) GO TO 300
 
!     THIS IS NOT THE FIRST PASS.  IF ITYPE(TH) ELEMENT ROUTINE IS IN
!     CORE, PROCESS IT.
 
 IF (itemp == lincor) GO TO 310
 
!     THE ITYPE(TH) ELEMENT ROUTINE IS NOT IN CORE.  IF THIS ELEMENT
!     ROUTINE IS IN A LINK THAT ALREADY HAS BEEN PROCESSED READ THE NEXT
!     ELEMENT.
 
 IF (link(itemp) == 1) GO TO 120
 
!     SET A TO BE PROCESSED LATER FLAG FOR THE LINK IN WHICH THE ELEMENT
!     RESIDES
 
 link(itemp) = 0
 GO TO 120
 
!     SINCE THIS IS THE FIRST ELEMENT TYPE TO BE PROCESSED ON THIS PASS
!     OF THE ECPT RECORD, A CHECK MUST BE MADE TO SEE IF THIS ELEMENT
!     IS IN A LINK THAT HAS ALREADY BEEN PROCESSED.  IF IT IS SUCH AN
!     ELEMENT, WE KEEP IFIRST = 1 AND READ THE NEXT ELEMENT.
 
 300 IF (link(itemp) == 1) GO TO 120
 
!     SET THE CURRENT LINK IN CORE = ITEMP AND IFIRST = 0
 
 lincor = itemp
 ifirst = 0
 
!     CALL THE PROPER ELEMENT ROUTINE.
 
!            CROD  CBEAM  CTUBE  CSHEAR  CTWIST  CTRIA1  CTRBSC
!              1     2      3       4       5       6       7
!          CTRPLT  CTRMEM  CONROD  ELAS1  ELAS2  ELAS3  ELAS4
!             8       9      10      11     12     13     14
!          CQDPLT  CQDMEM  CTRIA2  CQUAD2 CQUAD1  CDAMP1  CDAMP2
!            15      16      17      18      19     20      21
!          CDAMP3  CDAMP4  CVISC  CMASS1  CMASS2  CMASS3  CMASS4
!            22      23      24     25      26      27      28
!           CONM1  CONM2   PLOTEL REACT   QUAD3   CBAR    CCONE
!             29     30      31     32      33     34       35
 310 IF (itype > 35) GO TO 120
 GO TO (320,890,340,350,360,370,380,390,400,320,  &
     410,420,430,440,450,460,470,480,490,120,  &
     120,120,120,120,120,120,120,120,120,120, 120,890,890,500,510), itype
 320 CALL krod
 GO TO 120
 340 CALL ktube
 GO TO 120
 350 CALL kpanel (4)
 GO TO 120
 360 CALL kpanel (5)
 GO TO 120
 370 CALL ktriqd (1)
 GO TO 120
 380 CALL ktrbsc (0)
 GO TO 120
 390 CALL ktrplt
 GO TO 120
 400 CALL ktrmem (0)
 GO TO 120
 410 CALL kelas (1)
 GO TO 120
 420 CALL kelas (2)
 GO TO 120
 430 CALL kelas (3)
 GO TO 120
 440 CALL kelas (4)
 GO TO 120
 450 CALL kqdplt
 GO TO 120
 460 CALL kqdmem
 GO TO 120
 470 CALL ktriqd (2)
 GO TO 120
 480 CALL ktriqd (4)
 GO TO 120
 490 CALL ktriqd (3)
 GO TO 120
 500 CALL kbar
 GO TO 120
 510 IF (nbpw <= 32) CALL kconed
 IF (nbpw > 32) CALL kcones
 GO TO 120
 
!     AN END OF LOGICAL RECORD HAS BEEN HIT ON THE ECPT.  IF NPVTWR = 0,
!     THE PIVOT POINT HAS NOT BEEN WRITTEN ON ECPTNL AND NO ELEMENTS IN
!     THE CURRENT ECPT RECORD ARE PLASTIC.
 
 520 IF (ipass /= 1) GO TO 550
 IF (npvtwr > 0) THEN
   GO TO   540
 END IF
 530 CALL WRITE (ecptnl,-npvt,1,eor)
 GO TO 550
 540 CALL WRITE (ecptnl,0,0,eor)
 550 ipass = 2
 link(lincor) = 1
 DO  i = 1,nlinks
   IF (link(i) == 0) GO TO 570
 END DO
 GO TO 580
 
!     SINCE AT LEAST ONE LINK HAS NOT BEEN PROCESSED THE ECPT FILE MUST
!     BE BACKSPACED
 
 570 CALL bckrec (ecpt)
 GO TO 110
 
!     WRITE THE NO. OF ROWS IN CORE UNTO THE KGGL FILE USING ZBLPKI.
 
 580 i1  = 0
 590 i2  = 0
 ibeg = i6x6k + i1*jmax
 CALL bldpk (2,ipr,kggl,0,0)
 600 i2  = i2 + 1
 IF (i2 > ngpct) GO TO 620
 jj  = igpct + i2
 INDEX = IABS(iz(jj)) - 1
 lim = 6
 IF (iz(jj) < 0) lim = 1
 jjj = ipoint + i2
 kkk = ibeg + iz(jjj) - 1
 i3  = 0
 610 i3  = i3 + 1
 IF (i3 > lim) GO TO 600
 INDEX = INDEX + 1
 kkk = kkk + 1
 dpword = dz(kkk)
 IF (dpword /= 0.0D0) CALL zblpki
 GO TO 610
 620 CALL bldpkn (kggl,0,mcbkgg)
 i1  = i1 + 1
 IF (i1 < nrowsc) GO TO 590
 
!     IF LROWIC = TNROWS, PROCESSING OF THE CURRENT ECPT RECORD HAS BEEN
!     COMPLETED.
 
 IF (lrowic == tnrows) GO TO 30
 CALL bckrec (ecpt)
 frowic = frowic + nrowsc
 lrowic = lrowic + nrowsc
 ipass  = 2
 GO TO 110
 
!     NO ELEMENTS ARE CONNECTED TO THE PIVOT POINT.  OUTPUT ZERO
!     COLUMN(S).  ALSO, WRITE NEGATIVE PIVOT POINT ON ECPTNL.
 
 630 lim = 6
 IF (inpvt(1) < 0) lim = 1
 DO  i = 1,lim
   CALL bldpk (2,ipr,kggl,0,0)
   CALL bldpkn (kggl,0,mcbkgg)
 END DO
 CALL skprec (ecpt,1)
 CALL WRITE (ecptnl,-IABS(inpvt(1)),1,eor)
 GO TO 30
 
!     ECPT PROCESSING HAS BEEN COMPLETED SINCE AN EOF HAS BEEN READ ON
!     GPCT.
 
 650 CALL CLOSE (gpct,clsrw)
 CALL CLOSE (ecpt,clsrw)
 CALL CLOSE (kggl,clsrw)
 CALL CLOSE (ecptnl,clsrw)
 IF (kickof   == -1) GO TO 865
 IF (mcbkgg(6) /= 0) GO TO 654
 DO  i = 2,7
   mcbkgg(i) = 0
 END DO
 GO TO 656
 654 mcbkgg(3) = mcbkgg(2)
 656 CALL wrttrl (mcbkgg)
 CALL wrttrl (trail)
 
!     BEGIN EST PROCESSING
 
 left = bufr4 - 1
 icc  = ncstm +  mused
 all  = .false.
 phase1 = .false.
 
!     READ THE FIRST RECORD OF CASECC INTO CORE.
 
 FILE = casecc
 CALL gopen (casecc,z(bufr1),0)
 CALL READ (*870,*658,casecc,iz(icc+1),left,eor,ncc)
 CALL mesage (-8,0,NAME)
 658 iplset = icc + jplset
 plaset = iz(iplset)
 istset = icc + jstset
 IF (iz(istset) < 0) THEN
   GO TO   660
 ELSE IF (iz(istset) == 0) THEN
   GO TO   670
 ELSE
   GO TO   680
 END IF
 660 all = .true.
 GO TO 705
 670 nonlst = -1
 GO TO 705
 
!     THE USER HAS REQUESTED A PROPER SUBSET OF HIS SET OF ELEMENTS FOR
!     WHICH HE WANTS STRESS OUTPUT.  FIND THE SET IN OPEN CORE AND
!     DETERMINE ZERO POINTER AND LENGTH OF THE SET.
 
 680 isym = icc + jsym
 isetno = isym + iz(isym) + 1
 lset = iz(isetno+1)
 690 iset = isetno + 2
 nset = iz(isetno+1) + iset - 1
 IF (iz(isetno) == iz(istset)) GO TO 700
 isetno = nset + 1
 IF (isetno < ncc) GO TO 690
 all = .true.
 700 iz(nset+1) = 2**14 + 1
 705 CALL CLOSE (casecc,clsrw)
 IF (plaset /= -1) GO TO 706
 jj = 1
 plfact(1) = 1.0
 GO TO 731
 706 CONTINUE
 
!     SEARCH THE MPT FOR THE PLA SET
 
 FILE = mpt
 CALL preloc (*860,z(bufr1-3),mpt)
 CALL locate (*895,z(bufr1-3),planos,iflag)
 
!     READ A PLA SET NO.
 
 710 CALL READ (*895,*895,mpt,setno,1,neor,iflag)
 jj = 0
 720 CALL READ (*895,*895,mpt,nn,1,neor,iflag)
 IF (nn == -1) GO TO 730
 jj = jj + 1
 IF (jj == 1) plfact(1) = fnn
 GO TO 720
 730 IF (setno /= plaset) GO TO 710
 nplalp = jj
 plfact(2) = 0.0
 CALL CLOSE (mpt,clsrw)
 731 CONTINUE
 
!     PROCESS THE EST
 
 estltr(1) = estl
 estnlt(1) = estnl
 DO  i = 2,7
   estltr(i) = 0
   estnlt(i) = 0
 END DO
 ASSIGN 820 TO nosd
 CALL gopen (  est,z(bufr1),0)
 CALL gopen ( estl,z(bufr2),1)
 CALL gopen (estnl,z(bufr3),1)
 FILE = est
 
!     READ THE ELEMENT TYPE.  IF THE ELEMENT TYPE IS ADMISSIBLE TO
!     PIECEWISE LINEAR ANALYSIS, WRITE IT TWICE.  OTHERWISE GO TO NEXT
!     RECORD.
 
 750 CALL READ (*850,*880,est,itype,1,neor,iflag)
 IF (plaary(itype) == 1) GO TO 755
 CALL skprec (est,1)
 GO TO 750
 755 CALL WRITE (estl, itype,1,neor)
 CALL WRITE (estnl,itype,1,neor)
 
!     READ THE EST ENTRY
 
 760 idx  = (itype-1)*incr
 nwds = NE(idx+12)
 CALL READ (*870,*840,est,xecpt,nwds,neor,iflag)
 IF (plaary(itype) == 0) GO TO 820
 IF (itype > 38) GO TO 820
!              CROD      CBEAM     CTUBE     CSHEAR    CTWIST
!                1         2         3         4         5
 GO TO (   130,      820,      150,      820,      820,
!              CTRIA1    CTRBSC    CTRPLT    CTRMEM    CONROD
!                6         7         8          9       10  &
 160,      820,      820,      170,      130,
!              CELAS1    CELAS2    CELAS3    CELAS4    CQDPLT
!                11        12        13        14        15  &
 820,      820,      820,      820,      820,
!              CQDMEM    CTRIA2    CQUAD2    CQUAD1    CDAMP1
!                16        17        18        19        20  &
 180,      190,      200,      210,      820,
!              CDAMP2    CDAMP3    CDAMP4    CVISC     CMASS1
!                21        22        23        24        25  &
 820,      820,      820,      820,      820,
!              CMASS2    CMASS3    CMASS4    CONM1     CONM2
!                26        27        28        29        30  &
 820,      820,      820,      820,      820,
!              PLOTEL    REACT     QUAD3     CBAR      CCONE
!                31        32        33        34        35  &
 820,      820,      820,      220,      820,
!              CTRIARG   CTRAPRG   CTORDRG
!                36        37        38  &
 820,      820,      820), itype
 
!     THE ELEMENT IS STRESS DEPENDENT.  DETERMINE IF STRESS OUTPUT IS
!     REQUESTED.
!     AN EXAMPLE... IF WE HAVE IN CASE CONTROL
!     SET 5 = 1,2,3,98THRU100,4THRU15,81,18,82,90,92
!     THEN THE WORDS IN CASE CONTROL ARE...
!       IZ(ISETNO) = 5,12,1,2,3,4,-15,18,81,82,90,92,98,-100 = IZ(NSET)
 
 765 IF (all) GO TO 800
 IF (nonlst == -1) GO TO 760
 ielid = iecpt(1)
 i = iset
 770 IF (i    > nset) GO TO 760
 IF (iz(i+1) < 0) GO TO 780
 IF (ielid == iz(i)) GO TO 800
 IF (ielid < iz(i)) GO TO 760
 i = i + 1
 GO TO 790
 780 IF (ielid >= iz(i) .AND. ielid <= IABS(iz(i+1))) GO TO 800
 IF (ielid < iz(i)) GO TO 760
 i = i + 2
 790 IF (iz(i) > 0) GO TO 770
 all = .true.
 llllll(1) = iz(istset)
 llllll(2) = iz(i)
 CALL mesage (30,92,llllll)
 800 outfil = estnl
 nnlel  = nnlel + 1
 GO TO 830
 820 outfil = estl
 nlel   = nlel + 1
 830 CALL WRITE (outfil,xecpt,nwds,neor)
 GO TO 760
 840 CALL WRITE (estl,0,0,eor)
 CALL WRITE (estnl,0,0,eor)
 GO TO 750
 
!     WRAP UP ROUTINE
 
 850 CALL CLOSE  (est,clsrw)
 CALL CLOSE  (estl,clsrw)
 CALL CLOSE  (estnl,clsrw)
 CALL wrttrl (estltr)
 CALL wrttrl (estnlt)
 865 RETURN
 
!     FATAL ERRORS
 
 860 CALL mesage (-1,FILE,NAME)
 870 CALL mesage (-2,FILE,NAME)
 880 CALL mesage (-3,FILE,NAME)
 890 CALL mesage (-30,87,itype)
 
!     UNABLE TO FIND PLFACT CARD IN THE MPT WHICH WAS CHOSEN BY THE USER
!     IN CASECC.
 
 895 trail(1) = hmpt
 trail(2) = NAME(1)
 CALL mesage (-32,plaset,trail)
 RETURN
END SUBROUTINE pla1
