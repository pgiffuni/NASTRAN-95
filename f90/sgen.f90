SUBROUTINE sgen
     
!     THIS MODULE PREPARES THE INPUT FILES TO NASTRAN FROM A SUBSTRUCTUR
!     FORMULATION IN ORDER TO RUN THE SOLUTION PHASE OF NASTRAN.
!     3 MAJOR STEPS ARE-
 
!     1.  READ CONSTRAINT AND DYNAMICS DATA, CONVERT TO PSEUDO-STRUCTURE
!         DATA, AND OUTPUT ON GP4S AND DYNS.
 
!     2.  READ LOAD COMBO. DATA AND ASSEMBLE SCALAR LOAD SETS ON OUTPUT
!         FILE GP3S.
 
!     3.  BUILD DUMMY FILES FOR EXECUTION- CASEI, GPL, EQEXIN, GPDT,
!         BGPDT, CSTM, AND SIL.
 
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        andf,orf,complf
 LOGICAL :: nolc,nols,psuedo,stest
 INTEGER :: temp(10),temp2(10),TYPE(2),spcs1(2),spcsd(2),  &
     mpcs(2),spcs(2),loadc(2),ctypes(2,8),ctypeo(2,8),  &
     dareas(2),delays(2),dphses(2),tics(2),nlimit(3),  &
     minus(3),icode(4,9),icomp(32),ltab(4,9),mcb(7),  &
     lsload(3),lload(3),nsgen(2),ncasec(2),z(4)
 REAL :: rz,fact,rtemp(10),rtemp2(10)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / dry,NAME(2),luset,nogpdt
 COMMON /sgencm/ nono,nss,iptr,buf1,buf2,buf3,nz
 COMMON /zzzzzz/ rz(1)
 COMMON /system/ ibuf,outt
 COMMON /two   / two(2)
 COMMON /unpakx/ ity,irow,nrow,incr
 EQUIVALENCE     (rz(1),z(1)),            (temp(1),rtemp(1)),  &
     (temp2(1),rtemp2(1)), (ctypes(1,1),mpcs(1)),   (ctypes(1,2),spcs(1)),  &
     (ctypes(1,3),spcs1(1)),  (ctypes(1,4),spcsd(1)),  &
     (ctypes(1,5),dareas(1)), (ctypes(1,6),delays(1)),  &
     (ctypes(1,7),dphses(1)), (ctypes(1,8),tics(1))
 DATA    minus , nlimit /3*-1, 3*2147483647 /,  &
     eqss  / 4HEQSS /,     lods /4HLODS /
 DATA    casec , geom3 ,geom4 ,dynam        /  &
     101   , 102   ,103   ,104          /, cases , casei ,gpl   ,eqex  ,gpdt  /  &
     201   , 202   ,203   ,204   ,205   /, bgpdt , sil   ,gp3s  ,gp4s  ,dyns  /  &
     206   , 207   ,208   ,209   ,210   /, scrt  , scrt2 /  &
     201   , 202   /
 DATA    pvec  / 4HPVEC/      ,nsgen / 4HSGEN,4H    /
 
!     BULK DATA CARD CODES
 
 DATA    icode / &
!             MPCS  &
 1110  ,11    ,0     ,0      , &
!             SPCS  &
 810   ,8     ,0     ,0      , &
!             SPCS1  &
 710   ,7     ,0     ,0      , &
!             SPCSD  &
 610   ,6     ,0     ,0      , &
!             LOADC  &
 500   ,5     ,0     ,0      , &
!             DAREAS  &
 9027  ,90    ,0     ,0      , &
!             DELAYS  &
 9137  ,91    ,0     ,0      , &
!             DPHASES  &
 9277  ,92    ,0     ,0      , &
!             TICS  &
 9307  ,93    ,0     ,0      /
 
 DATA    ntypec / 4 /
 DATA    ntyped / 4 /
 DATA    ltab   / &
!             MPC  &
 4901  ,49    ,17    ,1      , &
!             SPC  &
 5501  ,55    ,16    ,2      , &
!             SPC1  &
 5481  ,58    ,12    ,3      , &
!             SPCD  &
 5110  ,51    ,256   ,4      , &
!             LOADC  &
 500   ,5     ,264   ,0      , &
!             DAREA  &
 27    ,17    ,182   ,5      , &
!             DELAY  &
 37    ,18    ,183   ,6      , &
!             DPHASE  &
 77    ,19    ,184   ,7      , &
!             TIC  &
 6607  ,66    ,137   ,8      / 
 DATA    lload / 4551  ,61    ,84    /       , lsload/ 5401  ,54    ,25    /
 DATA    mpcs  / 4HMPCS,4H    /,     spcs   / 4HSPCS,4H     /,  &
     spcs1 / 4HSPCS,4H1   /,     spcsd  / 4HSPCS,4HD    /,  &
     loadc / 4HLOAD,4HC   /,     dareas / 4HDARE,4HAS   /,  &
     delays/ 4HDELA,4HYS  /,     dphses / 4HDPHA,4HSES  /, tics  / 4HTICS,4H    /
 DATA    ncasec/ 4HCASE,4HCC  /
 DATA    ctypeo/ 4HMPC ,4H    , 4HSPC ,4H    ,  &
     4HSPC1,4H    , 4HSPCD,4H    ,  &
     4HDARE,4HA   , 4HDELA,4HY   ,  &
     4HDPHA,4HSE  , 4HTIC ,4H    /
 DATA    xxxx  / 4HXXXX       /
 
!     INITIALIZE
 
 ity  = 1
 incr = 1
 nono = 0
 large= two(2)
 nz   = korsz(z(1))
 ibs1 = nz   - ibuf + 1
 ibs2 = ibs1 - ibuf - 1
 ibs3 = ibs2 - ibuf
 buf1 = ibs3 - ibuf
 buf2 = buf1 - ibuf
 buf3 = buf2 - ibuf
 buf4 = buf3 - ibuf
 nz   = buf4 - 1
 IF (nz <= 0) GO TO 5011
 IF (NAME(1) == xxxx .AND. NAME(2) == xxxx) GO TO 3000
 
!     INITIALIZE LUSET AND NOGPDT FLAGS
 
 luset  = 0
 nogpdt = -1
 
!     FORM TABLES OF REFERENCED SID-S FOR LOAD, MPC, AND SPC
!     CASE CONTROL CARDS.
 
 
!     OPEN  SOF , GET EQSS ITEM , READ  SIL DATA INTO CORE
 
 CALL sofopn (z(ibs1),z(ibs2),z(ibs3))
 CALL sfetch (NAME,eqss,1,flag)
 item  = eqss
 IF (flag /= 1) GO TO 5001
 CALL suread (z(1),nz,nwds,flag)
 IF (flag /= 2) GO TO 5011
 nss = z(3)
 iz  = nwds + 1
 
!     READ SIL GROUP INTO CORE
 
 CALL sjump (nss)
 CALL suread (z(iz),nz-iz+1,nsil,flag)
 IF  (flag /= 2) GO TO 5011
 ipt =  iz + nsil - 2
 
!     FIND LENGTH OF VECTOR = LUSET
 
 ic = z(ipt+1)
 CALL decode (ic,icomp,nc)
 luset = z(ipt) + nc - 1
 nogpdt= luset
 
!     READ EQSS ( G ,IP, AND C AT A TIME) AND CONVERT IP TO SIL .
!     WRITE ON SCRT
 
 is   = 0
 FILE = scrt
 CALL gopen (scrt,z(buf2),1)
 CALL sfetch (NAME,eqss,1,flag)
 nj = 1
 CALL sjump (nj)
 
 50 CALL suread (temp,3,nwds,flag)
 IF (flag /= 1) GO TO 100
 ipt  = iz + 2*temp(2) - 2
 temp(2) = z(ipt)
 CALL WRITE (scrt,temp,3,0)
 GO TO 50
 100 is = is + 1
 CALL WRITE (scrt,temp,0,1)
 IF (is < nss) GO TO 50
 CALL CLOSE (scrt,1)
 
!     READ CONVERTED EQSS INTO CORE, STORE POINTERS TO THE BASIC  SUBS
!     IN  Z(IPTR) TO Z(NPTR)
!     CORE WILL CONTAIN-
!       1. 4 WORD HEADER
!       2. 2*NSS NAMES
!       3. NSS+1 POINTERS TO EACH BASIC SUBST.BLOCK
!       4. NSS BLOCKS OF G, IP, C DATA
!       5. NZB LEFT OVER
 
 
 iptr = iz
 nptr = iptr
 isub = iptr + nss + 1
 nzb  = nz - isub + 1
 FILE = scrt
 CALL gopen (scrt,z(buf2),0)
 DO  i = 1,nss
   z(nptr) = isub
   nptr = nptr + 1
   CALL READ (*9002,*110,scrt,z(isub),nzb,1,nwds)
   GO TO 5011
   110 isub = isub + nwds
   nzb  = nzb  - nwds
   IF (nzb <= 0) GO TO 5011
 END DO
 z(nptr) = isub
 CALL CLOSE (scrt,1)
 
!     ***  GEOM4 DATA CONVERSION  ***
 
!          IN  - MPCS,SPCS,SPCS1,SPCSD CARDS
!          OUT - MPC ,SPC ,SPC1 ,SCPD  ON SCRT
 
 FILE = geom4
 nog4 = 0
 CALL preloc (*400,z(buf1),geom4)
 mcb(1) = geom4
 CALL rdtrl (mcb)
 CALL gopen (scrt,z(buf2),1)
 stest = .false.
 
!     ***  MPCS CARDS  ***
 
!          IN  - NAME(2), G, C, F
!          OUT - SIL, 0, F
 
 CALL locate (*350,z(buf1),icode(1,1),idx)
 CALL WRITE (scrt,icode(1,1),3,0)
 stest = .true.
 icode(4,1) = 1
 TYPE(1) = mpcs(1)
 TYPE(2) = mpcs(2)
 ifl = 0
 lid = 0
 305 CALL READ (*9002,*346,geom4,j,1,0,nwds)
 IF (j /= lid) nsild = 0
 lid = j
 CALL WRITE (scrt,j,1,0)
 310 CALL READ (*9002,*346,geom4,temp,5,0,nwds)
 IF (temp(3) == -1) GO TO 345
 IF (temp(3) == 0) GO TO 310
 
!     FIND  REQUESTED SUBSTRUCTURE
 
 DO  i = 1,nss
   inam = 2*i + 3
   IF (z(inam) == temp(1) .AND. z(inam+1) == temp(2)) GO TO 330
 END DO
 
!     SUBSTRUCTURE NOT FOUND
 
 WRITE (outt,63290) uwm,temp(1),temp(2),TYPE,NAME
 GO TO 310
 
!     FOUND SUBSTRUCTURE NAME
 
 330 ipt  = iptr + i - 1
 igrd =  z(ipt)
 ngrd = (z(ipt+1) - z(ipt))/3
 
!     SEARCH FOR GRID POINT
 
 CALL bisloc (*334,temp(3),z(igrd),3,ngrd,igr)
 ig = igr + igrd - 1
 325 IF (z(ig-3) /= z(ig)) GO TO 331
 IF (ig <= igrd) GO TO 331
 ig = ig - 3
 GO TO 325
 331 code = z(ig+2)
 
!     FIND   THE  COMPONENT
 
 CALL decode (code,icomp,nc)
 IF (temp(4) == 0) temp(4) = 1
 DO  i = 1,nc
   IF (temp(4) /= icomp(i)+1) CYCLE
   ic = i
   GO TO 340
 END DO
 IF (z(ig+3) /= z(ig)) GO TO 334
 IF (ig+3 >= igrd+3*ngrd) GO TO 334
 ig = ig + 3
 GO TO 331
 
!     BAD COMPONENT
 
 334 nono = 1
 WRITE (outt,60220) ufm,(temp(i),i=1,4),TYPE,NAME
 GO TO 310
 
!     WRITE CONVERTED DATA ON SCRT
 
 340 temp(6) = z(ig+1) + ic - 1
 temp(7) = 0
 temp(8) = temp(5)
 CALL WRITE (scrt,temp(6),3,0)
 
!     CHECK FOR DUPLICATE DEPENDENT SIL-S
 
 IF (ifl  /=  0) GO TO 310
 IF (nsild == 0) GO TO 343
 DO  i = 1,nsild
   IF (z(isub+i-1) /= temp(6)) CYCLE
   nono = 1
   WRITE (outt,63620) ufm,j,temp(1),temp(2),temp(3),temp(4)
 END DO
 IF (nsild > nzb) GO TO 5011
 343 z(isub+nsild) = temp(6)
 nsild = nsild + 1
 ifl = 1
 GO TO 310
 
!     FINISHED ONE LOGICAL CARD, WRITE -1 FLAGS
 
 345 CALL WRITE (scrt,minus,3,0)
 ifl = 0
 GO TO 305
 
!     FINISHED ALL MPCS CARDS, WRITE EOR AND UPDATE TRAILER
 
 346 CALL WRITE (scrt,temp,0,1)
 
!     TURN OFF MPCS BIT
 
 j = (icode(2,1)-1)/16
 i =  icode(2,1)- 16*j
 mcb(j+2) = andf(complf(two(i+16)),mcb(j+2))
 
!     TURN ON MPC BIT
 
 j = (ltab(2,1)-1)/16
 i =  ltab(2,1)- 16*j
 mcb(j+2) = orf(two(i+16),mcb(j+2))
 
!     ***  SPCS CARDS  ***
 
!          IN  - SID, NAME(2), G, C, G, C, G, C, ..., -1, -1
!          OUT - SID, SIL, 0, 0 - REPEATED FOR EACH GRID
 
 350 CALL sgena (spcs,z(buf1),mcb,geom4,icode(1,2),0,scrt,ltab(1,2),1)
 
!     ***  SPCS1 CARDS  ***
 
!          IN  - SID, NAME(2), C, G, G, G, ..., -1
!          OUT - SID, 0, SIL, -1 - REPEATED FOR EACH GRID
 
 CALL sgenb (spcs1,z(buf1),mcb,geom4,icode(1,3),0,scrt,ltab(1,3),1)
 
!     ***  SPCSD CARDS  ***
 
!          IN  - SID, NAME(2), G, C, Y, ..., -1, -1, -1
!          OUT - SID, SIL, 0, Y - REPEATED FOR EACH GRID
 
 CALL sgena (spcsd,z(buf1),mcb,geom4,icode(1,4),1,scrt,ltab(1,4),1)
 
!     END OF CONSTRAINT CARD CONVERSION
 
 CALL CLOSE (geom4,1)
 CALL CLOSE (scrt,1)
 mcb(1) = gp4s
 CALL wrttrl (mcb)
 GO TO 700
 400 nog4 = 1
 
!     ***  DYNAMICS DATA CONVERSION  ***
 
!          IN  - DAREAS,DELAYS,DPHASES,TICS CARDS
!          OUT - DAREA ,DELAY ,DPHASE ,TIC  ON SCRT
 
 700 FILE  = dynam
 nodyn = 0
 CALL preloc (*750,z(buf1),dynam)
 mcb(1) = dynam
 CALL rdtrl (mcb)
 CALL gopen (scrt2,z(buf2),1)
 
!     ***  DAREAS  CARDS ***
 
!          IN  - SID, NAME(2), G, C, A, ..., -1, -1, -1
!          OUT - SID, SIL, 0, A - REPEATED FOR EACH GRID
 
 CALL sgena (dareas,z(buf1),mcb,dynam,icode(1,6),1,scrt2,ltab(1,6), 1)
 
!     ***  DELAYS CARDS  ***
 
!          IN  - SID, NAME(2), G, C, T, ..., -1, -1, -1
!          OUT - SID, SIL, 0, T - REPEATED FOR EACH GRID
 
 CALL sgena (delays,z(buf1),mcb,dynam,icode(1,7),1,scrt2,ltab(1,7), 1)
 
!    ***  DPHASES CARDS  ***
 
!         IN  - SID, NAME(2), G, C, TH, ..., -1, -1, -1
!         OUT - SID, SIL, 0, TH - REPEATED FOR EACH GRID
 
 CALL sgena (dphses,z(buf1),mcb,dynam,icode(1,8),1,scrt2,ltab(1,8), 1)
 
!     ***  TICS CARDS  ***
 
!          IN  - SID, NAME(2), G, C, U, V, ..., -1, -1, -1, -1
!          OUT - SID, SIL, 0, U, V - REPEATED FOR EACH GRID
 
 CALL sgena (tics,z(buf1),mcb,dynam,icode(1,9),2,scrt2,ltab(1,9),2)
 
!     END OF DYNAMICS CONVERSION
 
 CALL CLOSE (dynam,1)
 CALL CLOSE (scrt2,1)
 mcb(1) = dyns
 CALL wrttrl (mcb)
 GO TO 1000
 750 nodyn = 1
 
!     MERGE CONVERTED DATA WITH EXISTING DATA - GEOM4
 
 1000 IF (nog4 == 1) GO TO 1500
 CALL sgenm (ntypec,geom4,scrt,gp4s,icode(1,1),ltab(1,1),  &
     ctypes(1,1),ctypeo(1,1))
 
!     MERGE CONVERTED DATA WITH EXISTING DATA - DYNAMICS
 
 1500 IF (nodyn == 1) GO TO 2005
 CALL sgenm (ntyped,dynam,scrt2,dyns,icode(1,6),ltab(1,6),  &
     ctypes(1,1),ctypeo(1,1))
 
 
!     ***  GEOM3 PROCESSING  ***
 
!     THE LOAD VECTORS ARE COMBINED BY THE FACTORS
!     GIVEN ON THE LOADC CARDS AND MERGED WITH SLOAD CARDS
 
 2005 CONTINUE
 nolc = .true.
 nols = .true.
 CALL preloc (*2350,z(buf1),geom4)
 CALL locate (*2350,z(buf1),icode(1,5),idx)
 
!     READ FIRST GROUP OF LODS ITEM FOR SOLUTION STRUCTURE
 
 item = lods
 CALL sfetch (NAME,lods,1,flag)
 SELECT CASE ( flag )
   CASE (    1)
     GO TO 2045
   CASE (    2)
     GO TO 5001
   CASE (    3)
     GO TO 2042
   CASE (    4)
     GO TO 5001
   CASE (    5)
     GO TO 5001
 END SELECT
 
!     LODS ITEM DOES NOT EXIST
 
 2042 nolc = .true.
 GO TO 2350
 2045 CALL suread (z(1),nz,nwds,itest)
 SELECT CASE ( itest )
   CASE (    1)
     GO TO 5011
   CASE (    2)
     GO TO 2047
   CASE (    3)
     GO TO 5002
 END SELECT
 2047 nss = z(4)
 iss1= 5
 nl  = z(3)
 ipt = 2*nss + 5
 izl = 2*nss + ipt + 2
 z(ipt  ) = izl
 z(ipt+1) = 0
 IF (izl+nss+nl <= nz) GO TO 2050
 
!     INSUFFICIENT CORE
 
 CALL CLOSE (geom4,1)
 GO TO 5011
 
!     READ REMAINDER OF LODS INTO OPEN CORE AT Z(IZL)
 
 2050 DO  i = 1,nss
   CALL suread (z(izl),nz-izl+1,nwds,itest)
   SELECT CASE ( itest )
     CASE (    1)
       GO TO 5011
     CASE (    2)
       GO TO 2060
     CASE (    3)
       GO TO 5002
   END SELECT
   2060 izl= izl + nwds
   jg = ipt + 2*i
   z(jg  ) = izl
   z(jg+1) = z(jg-1) + nwds - 1
 END DO
 
!     CORE NOW CONTAINS
 
!            WORDS                 CONTENTS
!        ------------------    -----------------------------------
!        1--(IPT-1)            HEADER GROUP
!        IPT--IPT+2*(NSS+1)    LOAD DATA POINTER, NO. OF PRIOR LOAD
!                                  VECTORS (2 WORDS PER STRUCTURE)
!        IPT+2*(NSS+1)+1 --=   NO OF LOADS + LOAD SET IDS
!                                  GROUPED BY BASIC STRUCTURE
 
!     READ LOADC DATA CARDS AND CONVERT
 
!          IN  - SET ID, FACTOR, (NAME(2),SET,FACTOR) (REPEATED)
!          OUT - SET ID, FACTOR, (VECTOR NO.,FACTOR)
 
 TYPE(1) = loadc(1)
 TYPE(2) = loadc(2)
 CALL OPEN (*9001,scrt,z(buf2),1)
 2150 CALL READ (*9002,*2300,geom4,temp,2,0,nwds)
 CALL WRITE (scrt,temp,2,0)
 lid = temp(1)
 nolc=.false.
 
!     READ AN ENTRY
 
 2160 CALL fread (geom4,temp,4,0)
 IF( temp(3) == -1) GO TO 2280
 
!     FIND SUBSTRUCTURE AND SET
 
 DO  i = 1,nss
   inam = iss1 + 2*(i-1)
   IF (z(inam) == temp(1) .AND. z(inam+1) == temp(2)) GO TO 2220
 END DO
 
!     SUBSTRUCTURE NOT FOUND
 
 WRITE (outt,63290) uwm,temp(1),temp(2),TYPE,NAME
 GO TO 2160
 
!     FOUND SUBSTRUCTURE NAME
 
 2220 jpt = ipt + 2*i - 2
 
!     POINTER TO LODS DATA FOR THIS SUBSTRUCTURE
 
 ild = z(jpt)
 
!     NUMBER OF SETS IN LODS DATA FOR THIS SUBSTRUCTURE
 
 nset = z(ild)
 
!     FIND LOADC SET IN LODS DATA
 
 IF (nset == 0) GO TO 2240
 DO  i = 1,nset
   ip = ild + i
   IF (z(ip) /= temp(3)) CYCLE
   lvec = z(jpt+1) + i
   GO TO 2250
 END DO
 
!     SET NOT FOUND
 
 2240 nono = 1
 WRITE (outt,63310) ufm,NAME,lid,temp(3),temp(1),temp(2)
 GO TO 2160
 2250 temp(1) = lvec
 temp(2) = temp(4)
 CALL WRITE (scrt,temp,2,0)
 GO TO 2160
 
!     END OF LOGICAL LOADC CARD
 
 2280 CALL WRITE (scrt,temp,0,1)
 GO TO 2150
 
!     END OF LOADC RECORD
 
 2300 CALL CLOSE (scrt,1)
 2350 CALL CLOSE (geom4,1)
 
!     MERGE CONVERTED LOAD DATA WITH SLOAD DATA.
 
 
!     IF ANY ERRORS WERE DETECTED, SKIP LOAD COMPUTATION
 
 IF (nono /= 0) GO TO 3000
 CALL gopen (gp3s,z(buf4),1)
 
!     COPY LOAD CARDS TO GP3S
 
 CALL preloc (*2430,z(buf1),geom3)
 ldcd = 0
 CALL locate (*2420,z(buf1),lload,idx)
 ldcd = 1
 CALL WRITE (gp3s,lload,3,0)
 2405 CALL READ (*9002,*2410,geom3,z(1),nz,0,nwds)
 CALL WRITE (gp3s,z(1),nz,0)
 GO TO 2405
 2410 CALL WRITE (gp3s,z(1),nwds,1)
 
!     POSITION TO SLOAD CARDS
 
 2420 CALL locate (*2430,z(buf1),lsload,idx)
 nols = .false.
 2430 IF (nols) CALL CLOSE (geom3,1)
 IF (.NOT.(nols .AND. nolc)) CALL WRITE (gp3s,lsload,3,0)
 IF (nolc) GO TO 2530
 
!     COPY LOAD VECTORS TO SCRATCH FILE
 
 FILE = scrt2
 item = pvec
 IF (dry < 0) GO TO 2510
 CALL mtrxi (scrt2,NAME,pvec,z(buf3),flag)
 SELECT CASE ( flag )
   CASE (    1)
     GO TO 2520
   CASE (    2)
     GO TO 2431
   CASE (    3)
     GO TO 5001
   CASE (    4)
     GO TO 5001
   CASE (    5)
     GO TO 5001
   CASE (    6)
     GO TO 9001
 END SELECT
 2431 flag = 3
 GO TO 5001
 
!     IN DRY RUN MODE, LOADS PSEUDO-EXIST
 
 2510 psuedo =.true.
 GO TO 2530
 
!     LOADS EXIST
 
 2520 psuedo = .false.
 CALL gopen (scrt2,z(buf3),0)
 irec   = 1
 mcb(1) = scrt2
 CALL rdtrl (mcb)
 nvec  = mcb(2)
 luset = mcb(3)
 IF (2*luset < nz) GO TO 2530
 
!     INSUFFICIENT CORE
 
 CALL CLOSE (scrt2,1)
 CALL CLOSE (gp3s ,1)
 CALL CLOSE (geom3,1)
 GO TO 5011
 
!     MERGE REAL AND ARTIFICIAL SLOAD CARDS
 
 2530 sidc = 0
 irow = 1
 nrow = luset
 IF (.NOT.nolc) CALL OPEN (*9001,scrt,z(buf2),0)
 2550 IF (nols) GO TO 2560
 FILE = geom3
 CALL READ (*9002,*2560,geom3,temp2,3,0,nwds)
 GO TO 2570
 2560 IF (nolc) GO TO 2900
 temp2(1) = large
 2570 sids = temp2(1)
 IF (nolc) GO TO 2635
 IF (sidc > sids) GO TO 2600
 
!     READ THE SID AND FACTOR OF THE LOADC CARD ITSELF
 
 FILE = scrt
 CALL READ (*2580,*9003,scrt,temp,2,0,nwds)
 GO TO 2600
 2580 temp(1) = large
 nolc = .true.
 CALL CLOSE (scrt,1)
 IF (nols) GO TO 2900
 2600 CONTINUE
 DO  i = 1,luset
   rz(i) = 0.0
 END DO
 sidc = temp(1)
 fact = rtemp(2)
 IF (.NOT.nolc) GO TO 2670
 2635 IF (nols) GO TO 2900
 
!     NO MORE LOADC CARDS, WRITE ENTIRE  SLOAD  RECORD
 
 CALL WRITE (gp3s,temp2,3,0)
 FILE = geom3
 2640 CALL READ (*9002,*2650,geom3,z(1),nz,0,nwds)
 CALL WRITE (gp3s,z(1),nz,0)
 GO TO 2640
 2650 CALL WRITE (gp3s,z(1),nwds,1)
 GO TO 2900
 2670 IF (.NOT.nols) GO TO 2680
 
!     NO MORE SLOAD CARDS ARE PRESENT
 
 sids = large
 GO TO 2700
 
!     BOTH LOADC AND SLOAD CARDS ARE PRESENT
 
 2680 IF (sids < sidc) GO TO 2810
 
!     READ LOADC DATA, FIND VECTOR, UNPACK, MULT BY FACTOR, AND ADD
!     TO FIND A MATRIX COLUMN,USING FWDREC, CHANGE ON 16
 
 2700 FILE = scrt
 CALL READ (*9002,*2790,scrt,temp,2,0,nwds)
 IF (temp(1) == 0 .OR. psuedo .OR. temp(2) == 0) GO TO 2700
 n = temp(1) - irec
 IF (n < 0) THEN
   GO TO  2710
 ELSE IF (n == 0) THEN
   GO TO  2750
 ELSE
   GO TO  2720
 END IF
 2710 n = -n
 DO  i = 1,n
   CALL bckrec (scrt2)
 END DO
 GO TO 2750
 2720 DO  i = 1,n
   CALL fwdrec (*2730,scrt2)
 END DO
 GO TO 2750
 
!     CANT FIND LOAD VECTOR
 
 2730 WRITE (outt,63320) sfm,temp(1),nvec,luset,NAME
 nono = 1
 GO TO 2900
 
!     NOW SCRT2 IS POSITIONED TO THE DESIRED LOAD VECTOR.  UNPACK IT AND
!     FACTOR AND ADD IT TO VECTOR AT TOP OF OPEN CORE
 
 2750 irec = temp(1) + 1
 CALL unpack (*2700,scrt2,rz(luset+1))
 DO  i = 1,luset
   rz(i) = rtemp(2)*fact*rz(luset+i)+rz(i)
 END DO
 GO TO 2700
 
!     HERE WHEN FINISHED COMBINING VECTORS FOR ONE LOADC CARD
 
 2790 CONTINUE
 IF (sidc < sids) GO TO 2850
 2810 iz = temp2(2)
 rz(iz) = rz(iz) +rtemp2(3)
 FILE = geom3
 CALL READ (*9002,*2840,geom3,temp2,3,0,nwds)
 IF (temp2(1) == sids) GO TO 2810
 sids = temp2(1)
 GO TO 2850
 2840 nols =.true.
 
!     WRITE OUT LOAD VECTOR IN SLOAD FORMAT
 
 2850 temp(1) = MIN0(sids,sidc)
 DO  i = 1,luset
   IF (rz(i) == 0.0) CYCLE
   temp(2)  = i
   rtemp(3) = rz(i)
   CALL WRITE (gp3s,temp,3,0)
 END DO
 IF (sids /= sidc) GO TO 2570
 GO TO 2550
 
!     ALL LOADS PROCESSED
 
 2900 CALL WRITE (gp3s,0,0,1)
 CALL WRITE (gp3s,nlimit,3,1)
 CALL CLOSE (scrt,1)
 CALL CLOSE (gp3s,1)
 CALL CLOSE (scrt2,1)
 CALL CLOSE (geom3,1)
 mcb(1) = gp3s
 
!     TURN ON SLOAD BIT IN GP3S TRAILER
!     ALSO LOAD CARD BIT IF LOAD CARDS EXIST
 
 DO  i = 2,7
   mcb(i) = 0
 END DO
 j = (lsload(2)-1)/16
 i = lsload(2)-16*j
 mcb(j+2) = two(i+16)
 IF (ldcd == 0) GO TO 2920
 j = (lload(2)-1)/16
 i = lload(2)-16*j
 mcb(j+2) = orf(mcb(j+2),two(i+16))
 2920 CALL wrttrl (mcb)
 
!     SPLIT CASE CONTROL  INTO SUBSTRUCTURE AND NORMAL NASTRAN
 
 3000 CALL OPEN (*9001,casec,z(buf1),0)
 CALL OPEN (*9001,cases,z(buf2),1)
 CALL OPEN (*9001,casei,z(buf3),1)
 FILE = cases
 3250 CALL READ (*3800,*3350,casec,z(1),nz,0,nwds)
 3350 IF (z(1) == ncasec(1) .AND. z(2) == ncasec(2)) FILE = casei
 CALL WRITE (cases,z,nwds,1)
 IF (FILE == casei) CALL WRITE (casei,z(1),nwds,1)
 GO TO 3250
 3800 CONTINUE
 mcb(1) = casec
 CALL rdtrl (mcb)
 mcb(1) = cases
 CALL wrttrl (mcb)
 mcb(1) = casei
 CALL wrttrl (mcb)
 CALL CLOSE (casec,1)
 CALL CLOSE (casei,1)
 CALL CLOSE (cases,1)
 IF (NAME(1) == xxxx .AND. NAME(2) == xxxx) RETURN
 IF (nono /= 0) GO TO 4050
 
!     GENERATE  FICTITIOUS GP1 DATA BLOCKS
 
 
!     ***  GPL FILE  ***
 
!     GPL HEADER RECORD HAS 3 WORD, (SEE GP1)
!     SET THE 3RD WORD, MULTIPLIER MULT, TO 1000
 
 DO  i = 2,7
   mcb(i) = 0
 END DO
 mcb(1) = gpl
 FILE   = gpl
 n = -1
 CALL OPEN (*9200,gpl,z(buf1),1)
 CALL fname (gpl,temp(1))
 temp(3) = 1000
 CALL WRITE (gpl,temp(1),3,1)
 DO  i = 1,luset
   CALL WRITE (gpl,i,1,0)
 END DO
 CALL WRITE (gpl,i,0,1)
 DO  i = 1,luset
   temp(1) = i
   temp(2) = 1000*i
   CALL WRITE (gpl,temp,2,0)
 END DO
 CALL WRITE (gpl,i,0,1)
 CALL CLOSE (gpl,1)
 mcb(2) = luset
 CALL wrttrl (mcb)
 
!     ***  EQEXIN FILE  ***
 
 4050 mcb(1) = eqex
 CALL gopen (eqex,z(buf1),1)
 DO  i = 1,luset
   temp(1) = i
   temp(2) = i
   CALL WRITE (eqex,temp,2,0)
 END DO
 CALL WRITE (eqex,temp,0,1)
 DO  i = 1,luset
   temp(1) = i
   temp(2) = 10*i + 2
   CALL WRITE (eqex,temp,2,0)
 END DO
 CALL WRITE (eqex,temp,0,1)
 CALL CLOSE (eqex,1)
 mcb(2) = luset
 CALL wrttrl (mcb)
 
!     ***  GPDT FILE  ***
 
 mcb(1) = gpdt
 DO  i = 3,7
   temp(i) = 0
 END DO
 temp(2) = -1
 CALL gopen (gpdt,z(buf1),1)
 DO  i = 1,luset
   temp(1) = i
   CALL WRITE (gpdt,temp,7,0)
 END DO
 CALL WRITE (gpdt,temp,0,1)
 CALL CLOSE (gpdt,1)
 mcb(2) = luset
 CALL wrttrl (mcb)
 IF (nono /= 0) GO TO 4200
 
!     ***  BGPDT FILE  ***
 
 mcb(1) = bgpdt
 DO  i = 2,4
   temp(i) = 0
 END DO
 temp(1) =-1
 CALL gopen (bgpdt,z(buf1),1)
 DO  i = 1,luset
   CALL WRITE (bgpdt,temp,4,0)
 END DO
 CALL WRITE (bgpdt,temp,0,1)
 CALL CLOSE (bgpdt,1)
 mcb(2) = luset
 CALL wrttrl (mcb)
 
!     ***  SIL FILE  ***
 
 4200 mcb(1) = sil
 CALL gopen (sil,z(buf1),1)
 DO  i = 1,luset
   CALL WRITE (sil,i,1,0)
 END DO
 CALL WRITE (sil,i,0,1)
 CALL CLOSE (sil,1)
 
 
 mcb(2) = luset
 mcb(3) = luset
 CALL wrttrl (mcb)
 IF (nono /= 0) dry=-2
 CALL sofcls
 RETURN
 
!     ERRORS
 
 5001 n = 2 - flag
 GO TO 5010
 5002 n = -itest - 4
 5010 IF (dry < 0) n = IABS(n)
 dry = -2
 CALL smsg (n,item,NAME)
 RETURN
 5011 n = -8
 GO TO 9100
 9001 n = -1
 GO TO 9100
 9002 n = -2
 GO TO 9100
 9003 n = -3
 9100 CALL sofcls
 IF (dry < 0) n = IABS(n)
 dry = -2
 9200 CALL mesage (n,FILE,nsgen)
 RETURN
 
!     MESSAGE FORMATS
 
 60220 FORMAT (a23,' 6022, SUBSTRUCTURE ',2A4,', GRID POINT',i9,  &
     ', COMPONENTS',i9,1H,, /30X,'REFERENCED ON ',2A4,  &
 ' CARD, DO NOT EXIST ON SOLUTION STRUCTURE ',2A4)
   63290 FORMAT (a25,' 6329, SUBSTRUCTURE ',2A4,' REFERENCED ON ',2A4,  &
       ' CARD', /30X,'IS NOT A COMPONENT BASIC SUBSTRUCTURE OF ',  &
       'SOLUTION STRUCTURE ',2A4,/30X,'THIS CARD WILL BE IGNORED')
   63310 FORMAT (a23,' 6331, SOLUTION SUBSTRUCTURE ',2A4,' - LOADC SET',i9,  &
       ' REFERENCES UNDEFINED LOAD', /30X,'SET',i9, ' OF BASIC SUBSTRUCTURE ',2A4)
   63320 FORMAT (a25,' 6332, CANT FIND LOAD VECTOR NUMBER',i9,' IN LOAD ',  &
       'MATRIX OF',i9,' COLUMNS', /32X,'BY',i9,  &
       ' ROWS FOR SOLUTION STRUCTURE ',2A4)
   63620 FORMAT (a23,' 6362, MPCS SET',i9,' IS ILLEGAL.', //5X,  &
       'SUBSTRUCTURE ',2A4,' GRID POINT',i9,' COMPONENT',i5,  &
       ' SPECIFIES A NON-UNIQUE DEPENDENT DEGREE OF FREEDOM')
 END SUBROUTINE sgen
