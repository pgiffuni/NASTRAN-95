SUBROUTINE apdb
     
!     AERODYNAMIC POOL DISTRIBUTOR AND GEOMETRY INTERPOLATOR FOR
!     COMPRESSOR BLADES (AERODYNAMIC THEORY 6) AND SWEPT TURBOPROP
!     BLADES (AERODYNAMIC THEORY 7).
 
!     THIS IS THE DMAP DRIVER FOR APDB
 
!     DMAP CALLING SEQUENCE
 
!     APDB     EDT,USET,BGPDT,CSTM,EQEXIN,GM,GO / AEROB,ACPT,FLIST,
!              GTKA,PVECT / V,N,NK/V,N,NJ/V,Y,MINMACH/V,Y,MAXMACH/
!              V,Y,IREF/V,Y,MTYPE/V,N,NEIGV/V,Y,KINDEX $
 
!     INPUT  DATA BLOCKS CSTM, GM AND GO MAY BE PURGED
!     OUTPUT DATA BLOCK  PVECT MAY BE PURGED
!     PARAMETERS NK AND NJ ARE OUTPUT, THE OTHERS ARE INPUT
 
 
 LOGICAL :: lmkaer,first,debug
 INTEGER :: sysbuf,rd,rdrew,wrt,wrtrew,clsrew,norew,eofnrw,  &
     NAME(2),aero(3),mkaer1(3),mkaer2(3),fluttr(3),  &
     flfact(3),itrl(7),strml1(3),strml2(3),scr1,FILE,  &
     flag,name1(6,2),buf(7),edt,bgpdt,cstm,eqexin,  &
     aerob,acpt,flist,pvect,corwds,pstrm(100),typin, typout,sine,iz(6)
 REAL :: minmac,maxmac
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / nk,nj,minmac,maxmac,iref,mtype(2),neigv,kindex
 COMMON /system/ sysbuf,iout,nsys(91)
 COMMON /apdbug/ debug
 COMMON /packx / typin,typout,ii,nn,incr
 COMMON /zzzzzz/ z(1)
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,norew,eofnrw
!     NAMES  -VALUE =  2   0    3    1      1      2     3
 EQUIVALENCE     (z(1),iz(1)), (minmac,macmin), (maxmac,macmax)
 DATA    aero  / 3202,32,0/, mkaer1 /3802,38,0/, mkaer2 /3702,37,0/
 DATA    fluttr/ 3902,39,0/, flfact /4102,41,0/
 DATA    strml1/ 3292,92,0/, strml2 /3293,93,0/
 DATA    edt   , bgpdt,cstm,eqexin  / 101   , 103  ,104 ,105     /
 DATA    aerob , acpt ,flist,pvect  / 201   , 202  ,203  ,205    /
 DATA    NAME  / 4HAPDB,4H          /, scr1 /301/
 DATA    itrl  / 7*0 / , first / .true./,  sine / 4HSINE/
 DATA    name1(1,1),name1(1,2) / 4HAERO,4H      /
 DATA    name1(2,1),name1(2,2) / 4HMKAE,4HRO    /
 DATA    name1(3,1),name1(3,2) / 4HFLFA,4HCT    /
 DATA    name1(4,1),name1(4,2) / 4HFLUT,4HTER   /
 DATA    name1(5,1),name1(5,2) / 4HSTRE,4HAML1  /
 DATA    name1(6,1),name1(6,2) / 4HSTRE,4HAML2  /
 
 debug = .false.
 CALL sswtch (20,j)
 IF (j == 1) debug = .true.
 
!     SELECT AERODYNAMIC THEORY
 
!     COMPRESSOR BLADES (AERODYNAMIC THEORY 6).
!     SWEPT TURBOPROPS  (AERODYNAMIC THEORY 7).
 
!     AT PRESENT THE USER SELECTS THE THEORY VIA THE NASTRAN CARD.
!     SET SYSTEM(93)=0  FOR THEORY 6 OR SYSTEM(93)=1 FOR THEORY 7.
!     NOTE - THE DEFAULT IS THEORY 6 (SYSTEM(93)=0).
 
!     FOR EXAMPLE, TO SELECT THEORY 7, USE THE FOLLOWING CARD -
!     NASTRAN SYSTEM(93)=1
 
 IF (nsys(91) == 0) mthd = 6
 IF (nsys(91) == 1) mthd = 7
 
 IF (debug) CALL bug1 ('BLANK COMM',1,nk,9)
 nogo  = 0
 maxsl = 100
 ibuf1 = korsz(z) - sysbuf
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 last  = ibuf3 - sysbuf - 1
 IF (last <= 0) GO TO 991
 left = corwds(z(1),z(last))
 
!     CREATE AEROB DATA BLOCK
 
 CALL gopen (aerob,z(ibuf2),wrtrew)
 
!     READ AERO CARD VALUES - BREF, SYMXZ AND SYMXY
 
 FILE = edt
 CALL preloc (*992,z(ibuf1),edt)
 CALL locate (*981,z(ibuf1),aero,flag)
 CALL READ   (*993,*994,edt,z(1),6,1,flag)
 IF (debug) CALL bug1 ('AERO CARD ',2,z,6)
 iz(1) = iz(5)
 iz(2) = iz(6)
 CALL WRITE (aerob,z,3,1)
 
!     READ IN MKAERO1 CARDS
 
 lmkaer = .false.
 next = 1
 CALL locate (*60,z(ibuf1),mkaer1,flag)
 CALL READ (*993,*10,edt,z(next),left,1,nx)
 GO TO 991
 10 n1 = next
 IF (debug) CALL bug1 ('MKAERO1   ',10,z(n1),nx)
 lmkaer = .true.
 20 n2 = n1 + 7
 loop40:  DO  i = n1,n2
   IF (iz(i) == -1) EXIT loop40
   buf(1) = iz(i)
   n3 = n2 + 1
   n4 = n3 + 7
   DO  j = n3,n4
     IF (iz(j) == -1) CYCLE loop40
     buf(2) = iz(j)
     CALL WRITE (aerob,buf,2,0)
   END DO
 END DO loop40
 50 IF (n4-next+1 >= nx) GO TO 60
 n1 = n1 + 16
 GO TO 20
 
!     READ IN MKAERO2 CARDS
 
 60 CALL locate (*80,z(ibuf1),mkaer2,flag)
 CALL READ (*993,*70,edt,z(next),left,1,nx)
 GO TO 991
 70 CALL WRITE (aerob,z(next),nx,0)
 IF (debug) CALL bug1 ('MKAERO2   ',70,z(next),nx)
 lmkaer = .true.
 80 CALL WRITE (aerob,0,0,1)
 CALL CLOSE (aerob,clsrew)
 IF (.NOT.lmkaer) GO TO 982
 itrl(1) = aerob
 itrl(2) = 1
 CALL wrttrl (itrl)
 
!     CREATE FLIST TABLE
 
 CALL OPEN  (*85,flist,z(ibuf2),wrtrew)
 CALL fname (flist,iz(next))
 CALL WRITE (flist,iz(next),2,1)
 CALL locate (*981,z(ibuf1),aero,flag)
 CALL READ (*993,*90,edt,z(next),left,1,nx)
 GO TO 991
 
!     FLIST CAN BE PURGED IF THE APPROACH IS NOT AERO
 
 85 IF (IABS(nsys(19)) /= 4) GO TO 115
 FILE = flist
 GO TO 992
 90 CALL WRITE (flist,aero,3,0)
 CALL WRITE (flist,z(next),nx,1)
 IF (debug) CALL bug1 ('FLIST AERO',90,z(next),nx)
 CALL locate (*983,z(ibuf1),flfact,flag)
 CALL READ (*993,*100,edt,z(next),left,1,nx)
 GO TO 991
 100 CALL WRITE (flist,flfact,3,0)
 CALL WRITE (flist,z(next),nx,1)
 IF (debug) CALL bug1 ('FLIST FLFA',100,z(next),nx)
 CALL locate (*984,z(ibuf1),fluttr,flag)
 CALL READ (*993,*110,edt,z(next),left,1,nx)
 GO TO 991
 110 CALL WRITE (flist,fluttr,3,0)
 CALL WRITE (flist,z(next),nx,1)
 IF (debug) CALL bug1 ('FLIST FLUT',110,z(next),nx)
 CALL CLOSE (flist,clsrew)
 itrl(1) = edt
 CALL rdtrl (itrl)
 itrl(1) = flist
 CALL wrttrl (itrl)
 115 CONTINUE
 
!     CREATE ACPT TABLE
 
 CALL gopen (acpt,z(ibuf2),wrtrew)
 
!     STORE EXTERNAL NODE NUMBER, INTERNAL NODE NUMBER AND BASIC
!     COORDINATES OF ALL NODES ON BLADE ON SCR1
 
 CALL gopen (scr1,z(ibuf3),wrtrew)
 
!     READ STREAML1 AND STREAML2 CARDS. STORE IN-CORE
 
 nsl1a = next
 CALL locate (*985,z(ibuf1),strml1,flag)
 CALL READ (*993,*120,edt,z(nsl1a),left,1,nsl1l)
 GO TO 991
 120 nsl1b = nsl1a + nsl1l - 1
 IF (debug) CALL bug1 ('STREAML1  ',120,z(nsl1a),nsl1l)
 nsl2a = nsl1b + 1
 left  = corwds(z(nsl2a),z(last))
 CALL locate (*986,z(ibuf1),strml2,flag)
 CALL READ (*993,*130,edt,z(nsl2a),left,1,nsl2l)
 GO TO 991
 130 nsl2b = nsl2a + nsl2l - 1
 IF (debug) CALL bug1 ('STREAML2  ',130,z(nsl2a),nsl2l)
 CALL CLOSE (edt,clsrew)
 
!     INPUT CHECKS  (ALL ARE THEORY DEPENDENT RESTRICTIONS)
!     STREAML1 - ALL CARDS MUST HAVE THE SAME NUMBER OF NODES
!     STREAML2 - THERE MUST BE AT LEAST THREE(3) STREAML2 CARDS.
!                (THIS IS A THEORY DEPENDENT RESTRICTION,
!                SEE AMG MODULE - COMPRESSOR BLADE CODE FOR AJJL)
!              - NSTNS MUST BE THE SAME FOR ALL STREAML2 CARDS
!                AND MUST EQUAL THE NO. OF NODES ON THE STRAML1 CARD
 
!     COUNT THE NUMBER OF STREAML2 CARDS
 
 nlines = nsl2l/10
 IF (debug) CALL bug1 ('NLINES    ',131,nlines,1)
 IF (nlines >= 3) GO TO 135
 nogo = 1
 WRITE (iout,3001) ufm,nlines
 135 IF (nlines > maxsl) GO TO 988
 
!     LOCATE STREAML1 CARDS THAT CORRESPOND TO STREAML2 CARDS BY
!     MATCHING SLN VALUES
 
 nline = 0
 DO  isln = nsl2a,nsl2b,10
   nline = nline + 1
   pstrm(nline) = -iz(isln)
 END DO
 
!     LOCATE SLN AND COUNT THE NUMBER OF COMPUTING STATIONS
 
 ipos = nsl1a
 145 DO  ns = ipos,nsl1b
   IF (iz(ns) == -1) GO TO 155
 END DO
 
!     CHECK FOR VALID SLN
 
 155 DO  nline = 1,nlines
   IF (iz(ipos) == -pstrm(nline)) GO TO 165
 END DO
 GO TO 175
 165 pstrm(nline) = ipos
 nstnsx = ns - ipos - 1
 IF (.NOT.first) GO TO 170
 nstns = nstnsx
 first = .false.
 GO TO 175
 
!     ALL NSTNSX MUST BE THE SAME
 
 170 IF (nstnsx == nstns) GO TO 175
 nogo = 2
 WRITE (iout,3002) ufm,iz(ipos)
 175 ipos = ns + 1
 IF (ipos < nsl1b) GO TO 145
 
!     IS THERE A STREAML1 CARD FOR EVERY STREAML2 CARD
 
 DO  nline = 1,nlines
   IF (pstrm(nline) > 0) CYCLE
   nogo = 3
   isln = -pstrm(nline)
   WRITE (iout,3003) ufm,isln
 END DO
 IF (nogo > 0) GO TO 1000
 
!     READ BGPDT
 
 nbg1 = nsl2b + 1
 left = corwds(z(nbg1),z(last))
 FILE = bgpdt
 CALL gopen (bgpdt,z(ibuf1),rdrew)
 CALL READ (*993,*200,bgpdt,z(nbg1),left,1,nbgl)
 GO TO 991
 200 CALL CLOSE (bgpdt,clsrew)
 IF (debug) CALL bug1 ('BGPDT     ',200,z(nbg1),nbgl)
 nbg2 = nbg1 + nbgl - 1
 
!     READ EQEXIN (RECORD 1)
 
 neq1 = nbg2 + 1
 left = corwds(z(neq1),z(last))
 FILE = eqexin
 CALL gopen (eqexin,z(ibuf1),rdrew)
 CALL READ (*993,*210,eqexin,z(neq1),left,1,neql)
 GO TO 991
 210 neq2 = neq1 + neql - 1
 IF (debug) CALL bug1 ('EQEXIN R1 ',210,z(neq1),neql)
 
!     READ EQEXIN (RECORD 2)
 
 neq21 = neq2 + 1
 left = corwds(z(neq21),z(last))
 CALL READ (*993,*215,eqexin,z(neq21),left,1,neq2l)
 GO TO 991
 215 neq22 = neq2 + neq2l - 1
 IF (debug) CALL bug1 ('EQEXIN R2 ',212,z(neq21),neq2l)
 CALL CLOSE (eqexin,clsrew)
 
!     WRITE ACPT
 
!     KEY WORD = 6 FOR COMPRESSOR BLADES, I.E. METHOD ID = 6
!     KEY WORD = 7 FOR SWEPT TURBOPROPS , I.E. METHOD ID = 7
 
!     WRITE CONSTANT PARAMETERS, WORDS 1 - 6
 
 buf(1) = mthd
 buf(2) = iref
 buf(3) = macmin
 buf(4) = macmax
 buf(5) = nlines
 buf(6) = nstns
 CALL WRITE (acpt,buf,6,0)
 IF (debug) CALL bug1 ('ACPT WRT 1',216,buf,6)
 
!     WRITE STREAMLINE DATA
 
 kn = neql/2
 nline = 0
 DO  nsl = nsl2a,nsl2b,10
   
!     MAKE SURE NSTNS ON ALL STREAML2 CARDS IS THE SAME
   
   IF (iz(nsl+1) == nstns) GO TO 217
   WRITE (iout,3004) uwm,iz(nsl)
   iz(nsl+1) = nstns
   
!     WRITE STREAML2 DATA
   
   217 CALL WRITE (acpt,z(nsl),10,0)
   IF (debug) CALL bug1 ('ACPT WRT 2',217,z(nsl),10)
   
!     WRITE BASIC X, Y AND Z FOR EACH NODE ON STREAML1 CARD
   
   nline = nline + 1
   ipos  = pstrm(nline)
   ipos1 = ipos + 1
   ipos2 = ipos + nstns
   DO  igdp = ipos1,ipos2
     
!     LOCATE INTERNAL NUMBER THAT CORRESOONDS TO THIS EXTERNAL NODE
     
     CALL bisloc (*220,iz(igdp),iz(neq1),2,kn,jloc)
     GO TO 225
     
!     STREAML1 REFERNCES AN EXTERNAL ID THAT DOES NOT EXIST
     
     220 nogo = 5
     WRITE (iout,3005) ufm,iz(ipos),iz(igdp)
     CYCLE
     
!     PICK-UP BASIC GRID DATA FOR THIS NODE
     
     225 intrl  = iz(neq1+jloc)
     isilc  = iz(neq21+jloc)
     jloc   = nbg1 + (intrl-1)*4
     buf(1) = iz(igdp)
     buf(2) = intrl
     buf(3) = isilc
     buf(4) = iz(jloc  )
     buf(5) = iz(jloc+1)
     buf(6) = iz(jloc+2)
     buf(7) = iz(jloc+3)
     
!     TEST FOR SCALAR POINT (CID = -1)
     
     IF (buf(4) >= 0) GO TO 227
     nogo = 6
     WRITE (iout,3006) ufm,iz(ipos),iz(igdp)
     227 CALL WRITE (acpt,buf(5),3,0)
     CALL WRITE (scr1,buf,7,0)
     IF (debug) CALL bug1 ('ACPT WRT 3',227,buf,7)
     
!-----DETERMINE DIRECTION OF BLADE ROTATION VIA Y-COORDINATES AT TIP
!-----STREAMLINE. USE COORDINATES OF FIRST 2 NODES ON STREAMLINE.
     
     IF (nline == nlines .AND. igdp == ipos1)   ytip1 = z(jloc+2)
     IF (nline == nlines .AND. igdp == ipos1+1) ytip2 = z(jloc+2)
     
   END DO
 END DO
 
 xsign = 1.0
 IF (ytip2 < ytip1) xsign = -1.0
 IF (debug) CALL bug1 ('XSIN      ',240,xsign,1)
 CALL WRITE (acpt,0,0,1)
 CALL WRITE (scr1,0,0,1)
 CALL CLOSE (acpt,clsrew)
 CALL CLOSE (scr1,clsrew)
 itrl(1) = acpt
 itrl(2) = 1
 itrl(3) = 0
 itrl(4) = 0
 itrl(5) = 0
 itrl(6) = 0
 itrl(7) = 0
 CALL wrttrl (itrl)
 IF (nogo > 0) GO TO 1000
 
!     SET OUTPUT PARAMETERS NK AND NJ FOR APPROPRIATE THEORY.
 
!     COMPRESSOR BLADES (THEORY 6) - NK = NJ = NSTNS*NLINES.
!     SWEPT TURBOPROPS  (THEORY 7) - NK = NJ = 2*NSTNS*NLINES.
 
 IF (mthd == 6) nk = nstns*nlines
 IF (mthd == 7) nk = 2*nstns*nlines
 nj = nk
 IF (debug) CALL bug1 ('BLANK COM ',241,nk,9)
 
!     CREATE PVECT PARTITIONING VECTOR     (PVECT MAY BE PURGED)
!     PVECT IS A COLUMN PARTITIONING VECTOR TO BE USED BY MODULE PARTN
!     TO PARTITION OUT EITHER THE SINE OR COSINE COLUMNS OF MATRIX
!     PHIA WHICH IS OUTPUT BY THE CYCT2 MODULE  WHEN DOING A CYCLIC
!     NORMAL MODES ANALYSIS
!     PARAMETER MTYPE=SINE OR COSINE (DEFAULT IS COSINE)
 
!     OPEN PVECT AND WRITE HEADER
 
 CALL OPEN (*270,pvect,z(ibuf2),wrtrew)
 
!     TEST FOR VALID NEIGV AND KINDEX
 
 IF (neigv <= 0 .OR. kindex < 0) GO TO 987
 
 CALL fname (pvect,buf)
 CALL WRITE (pvect,buf,2,1)
 
!     PVECT IS TO BE GENERATED
 
 left = left - neq2
 ncol = neigv
 IF (kindex > 0) ncol = 2*ncol
 ipos1 = neq2 + 1
 ipos2 = neq2 + ncol
 DO  ipv = ipos1,ipos2
   z(ipv) = 0.0
 END DO
 IF (kindex == 0) GO TO 260
 ipos3 = ipos1
 IF (mtype(1) /= sine) ipos3 = ipos1 + 1
 DO  ipv = ipos3,ipos2,2
   z(ipv) = 1.0
 END DO
 260 typin  = 1
 typout = 1
 ii     = 1
 nn     = ncol
 incr   = 1
 CALL makmcb (itrl,pvect,ncol,2,1)
 CALL pack (z(ipos1),pvect,itrl)
 IF (debug) CALL bug1 ('PVECT     ',260,z(ipos1),ncol)
 CALL CLOSE (pvect,clsrew)
 CALL wrttrl (itrl)
 270 CONTINUE
 
!     GENERATE GTKA TRANSFORMATION MATRIX
 
!     READ CSTM INTO CORE
 
 ncstm1 = 1
 ncstml = 0
 FILE   = cstm
 itrl(1)= cstm
 CALL rdtrl (itrl)
 IF (itrl(1) /= cstm) GO TO 300
 left = corwds(z(ncstm1),z(last))
 CALL gopen (cstm,z(ibuf1),rdrew)
 CALL READ (*993,*300,cstm,z(ncstm1),left,1,ncstml)
 GO TO 991
 300 ncstm2 = ncstm1 + ncstml - 1
 IF (debug) CALL bug1 ('CSTM      ',300,z(ncstm1),ncstml)
 CALL CLOSE (cstm,clsrew)
 
!     ALLOCATE WORK STORAGE
 
 ip1  = ncstm2 + 1
 ip2  = ip1 + nstns
 ip3  = ip2 + nstns
 ip4  = ip3 + nstns
 next = ip4 + 4*nstns
 left = left - next + 1
 IF (left <= 0) GO TO 991
 
!     GENERATE GTKA TRANSFORMATION MATRIX FOR APPROPRIATE THEORY.
 
!     COMPRESSOR BLADES (AERODYNAMIC THEORY 6).
 
 IF (mthd == 6) CALL apdb1 (ibuf1,ibuf2,next,left,nstns,nlines,  &
     xsign,ncstml,z(ncstm1),z(ip1),z(ip2),z(ip3),z(ip4))
 
!     SWEPT TURBOPROPS (AERODYNAMIC THEORY 7).
 
 IF (mthd == 7) CALL apdb2 (ibuf1,ibuf2,next,left,nstns,nlines,  &
     xsign,ncstml,z(ncstm1),z(ip1),z(ip2),z(ip3),z(ip4))
 GO TO 1000
 
!     ERROR MESSAGES
 
!     NO AERO CARD FOUND
 981 kode = 1
 GO TO 989
 
!     NO MKAERO1 OR MKAERO2 CARDS FOUND
 
 982 kode = 2
 GO TO 989
 
!     NO FLFACT CARD FOUND
 
 983 kode = 3
 GO TO 989
 
!     NO FLUTTER CARD FOUND
 
 984 kode = 4
 GO TO 989
 
!     NO STREAML1 CARD FOUND
 
 985 kode = 5
 GO TO 989
 
!     NO STREAML2 CARD FOUND
 
 986 kode = 6
 GO TO 989
 
!     NEIGV OR KINDEX INVALID
 
 987 WRITE (iout,2987) ufm,neigv,kindex
 GO TO 1091
 
!     MAXIMUM NUMBER OF STREAML2 CARDS EXCEEDED FOR
!     LOCAL ARRAY PSTRM. SEE ERROR MESSAGE FOR FIX.
 
 988 WRITE (iout,3007) ufm,maxsl
 GO TO 1091
 989 WRITE (iout,2989) ufm,name1(kode,1),name1(kode,2)
 GO TO 1091
 
!     NOT ENOUGH CORE
 
 991 ip1 = -8
 GO TO 999
 
!     DATA SET NOT IN FIST
 
 992 ip1 = -1
 GO TO 999
 
!     E-O-F ENCOUNTERED
 
 993 ip1 = -2
 GO TO 999
 
!     E-O-L ENCOUNTERED
 
 994 ip1 = -3
 999 CALL mesage (ip1,FILE,NAME)
 
 1000 IF (nogo == 0) GO TO 1099
 1091 CALL mesage (-37,0,NAME)
 1099 RETURN
 
 2987 FORMAT (a23,' - APDB MODULE - INVALID PARAMETER NEIGV OR KINDEX',  &
     ' INPUT.', /40X, 'DATA BLOCK PVECT (FILE 205) CANNOT BE GENERATED.', /40X,  &
     7HNEIGV =,i8,10H, kindex =,i8)
 2989 FORMAT (a23,' - MODULE APDB - BULK DATA CARD ',2A4,  &
     ' MISSING FROM INPUT DECK.')
 3001 FORMAT (a23,' - APDB MODULE - THE NO. OF STREAML2 CARDS INPUT =',  &
     i3, /40X,'THERE MUST BE AT LEAST THREE(3) STREAML2 CARDS', ' INPUT.')
 3002 FORMAT (a23,' - APDB MODULE - ILLEGAL NO. OF NODES ON STREAML1 ',  &
     'CARD WITH SLN =',i8, /40X,  &
     'ALL STREAML1 CARDS MUST HAVE THE SAME NUMBER OF NODES.')
 3003 FORMAT (a23,' - APDB MODULE - NO STREAML1 CARD FOR THE STREAML2',  &
     ' WITH SLN =',i8)
 3004 FORMAT (a25,' - APDB MODULE - STREAML2 WITH SLN =',i8, /42X,  &
     'NSTNS INCONSISTENT WITH NO. OF NODES ON STREAML2 CARD ',  &
     'FOR BLADE ROOT.', /42X,'CORRECT VALUE OF NSTNS WILL BE ',  &
     'SUBSTITUTED ON STREAML2 CARD.')
 3005 FORMAT (a23,' - APDB MODULE - STREAML1 CARD WITH SLN =',i8,  &
     ' REFERENCES NON-EXISTENT EXTERNAL NODE =',i8)
 3006 FORMAT (a23,' - APDB MODULE - STREAML1 CARD WITH SLN =',i8,  &
     ' REFERENCES A SCALAR POINT WITH EXTERNAL ID =',i8, /40X,  &
     'SCALAR POINTS ARE ILLEGAL. USE A GRID POINT.')
 3007 FORMAT (a23,' - APDB MODULE - MAXIMUM NUMBER OF STREAML2 CARDS ',  &
     'EXCEEDED FOR LOCAL ARRAY PSTRM.', /40X,  &
     'UPDATE VARABLE MAXSL AND ARRAY PSTRM IN ROUTINE APDB.',  &
     /40X,'CURRENT VALUE OF MAXSL AND DIMENSION OF PSTRM =',i4)
END SUBROUTINE apdb
