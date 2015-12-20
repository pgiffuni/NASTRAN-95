SUBROUTINE fvrst1
!    1. ENTRY POINT - FVRST1
 
!    2. PURPOSE -  THIS MODULE IS USED FOR FORCED VIBRATION RESPONSE
!                  ANALYSIS OF ROTATING CYCLIC STRUCTURES.
!                  FVRSTR1 GENERATES DATA BLOCKS FRLX, B1GG, M1GG,
!                  M2GG, BASEXG AND PDZERO. IT ALSO COMPUTES PARAMETERS
!                  FKMAX AND NOBASEX.
 
!    3. DMAP CALLING SEQUENCE -
 
!         FVRSTR1  CASECC,BGPDT,CSTM,DIT,FRL,MGG,, / FRLX,B1GG,M1GG,
!                  M2GG,BASEXG,PDZERO,, /V,N,NOMGG/V,Y,CYCIO/V,Y,NSEGS/
!                  V,Y,KMAX/V,N,FKMAX/V,Y,BXTID=-1/V,Y,BXPTID=-1/
!                  V,Y,BYTID=-1/V,Y,BYPTID=-1/V,Y,BZTID=-1/
!                  V,Y,BZPTID=-1/V,N,NOBASEX/V,N,NOFREQ/V,N,OMEGA  $
 
!    4. INPUT DATA BLOCKS -
 
!         CASECC - CASE CONTROL
!         BGPDT  - BASIC GRID POINT DEFINITION TABLE.
!         CSTM   - COORDINATE SYSTEM TRANSFORMATION MATRICES.
!         DIT    - DIRECT INPUT TABLES.
!         FRL    - FREQUENCY RESPONSE LIST. (FREQUENCIES IN RADIANS)
!         MGG    - GLOBAL MASS MATRIX (G-SET).
 
!         NOTE   - (1) ALL INPUT DATA BLOCKS CAN BE PURGED IF ONLY
!                      PARAMETERS FKMAX AND NOBASEX ARE TO BE COMPUTED.
!                  (2) CASECC, DIT AND FRL CAN BE PURGED IF FRLX AND
!                      BASEXG ARE PURGED.
 
!    5. OUTPUT DATA BLOCKS -
 
!         FRLX    - FREQUENCY RESPONSE LIST (MODIFIED).
!         B1GG    - CORIOLIS ACCELERATION COEFFICIENT MATRIX (G-SET).
!         M1GG    - CENTRIPETAL ACCELERATION COEFFICIENT MATRIX (G-SET).
!         M2GG    - BASE ACCELERATION COEFFICIENT MATRIX (G-SET).
!         BASEXG  - BASE ACCELERATION MATRIX (G-SET).
!         PDZERO  - LOAD MODIFICATION MATRIX IN BASE ACCELERATION
!                   PROBLEMS.
 
!         NOTE    - (1) ALL OUTPUT DATA BLOCKS CAN BE PURGED IF
!                       PARAMETER NOMGG =-1.
!                   (2) B1GG AND M1GG CAN BE PURGED IF NOMGG =-1 OR
!                       IF OMEGA = 0.0.
!                   (3) FRLX AND PDZERO CAN BE PURGED IF OMEGA = 0.0.
!                   (4) FRLX, PDZERO, M2GG AND BASEXG CAN BE PURGED
!                       IF NOMGG =-1 OR NOFREQ =-1 OR CYCIO =+1 OR IF
!                       ALL PARAMETERS BXTID = BXPTID = BYTID =-1.
 
!    6. PARAMETERS -
 
!        (A) NOMGG   - INPUT-INTEGER-NO DEFAULT.  MASS MATRIX WAS NOT
!                      GENERATED IF NOMGG =-1.
!        (B) CYCIO   - INPUT-INTEGER-NO DEFAULT.  THE INTEGER VALUE
!                      OF THIS PARAMETER SPECIFIES THE FORM OF THE INPUT
!                      AND OUTPUT DATA FOR CYCLIC STRUCTURES. A VALUE
!                      OF +1 IS USED TO SPECIFY PHYSICAL SEGMENT REPRE-
!                      SENTATION AND A VALUE OF -1 FOR CYCLIC TRANSFOR-
!                      MATION REPRESENTATION.
!        (C) NSEGS   - INPUT-INTEGER-NO DEFAULT.  THE NUMBER OF
!                      IDENTICAL SEGMENTS IN THE STRUCTURAL MODEL.
!        (D) KMAX    - INPUT-INTEGER-NO DEFAULT.  THE INTEGER VALUE
!                      OF THIS PARAMETER SPECIFIES THE MAXIMUM VALUE
!                      OF THE HARMONIC INDEX.THE MAXIMUM VALUE OF
!                      KMAX IS NSEGS/2.
!        (E) FKMAX   - OUTPUT-INTEGER-NO DEFAULT.  FUNCTION OF KMAX.
!        (F) BXTID   - INPUT -INTEGER-DEFAULTS.  THE VALUES OF THESE
!        (G) BYTID     PARAMETERS DEFINE THE SET IDENTIFICATION NUMBERS
!        (H) BZTID     OF THE TABLEDI BULK DATA CARDS WHICH DEFINE THE
!        (I) BXPTID    COMPONENTS OF THE BASE ACCELERATION VECTOR. THE
!        (J) BYPTID    TABLES REFERED TO BY BXTID, BYTID AND BZTID
!        (K) BZPTID    DEFINE MAGNITUDE(LT-2) AND THE TABLES REFERED TO
!                      BY BXPTID, BYPTID AND BZPTID DEFINE PHASE(DEGREE)
!                      THE DEFAULT VALUES ARE -1 WHICH MEANS THAT THE
!                      RESPECTIVE TERMS ARE IGNORED.
!        (L) NOBASEX - OUTPUT-INTEGER-NO DEFAULT.  NOBASEX =-1 IF DATA
!                      BLOCK BASEXG IS NOT GENERATED.
!        (M) NOFREQ  - INPUT-INTEGER-NO DEFAULT. NOFREQ =-1 IF FREQUENCY
!                      WAS NOT SELECTED IN THE CASE CONTROL DECK.
!        (N) OMEGA   - INPUT-REAL-NO DEFAULT.  ROTATIONAL SPEED OF THE
!                      STRUCTURE IN RADIANS. OMEGA = 2*PI*RPS.
 
!    7. METHOD -  SEE FUNCTIONAL MODULE DESCRIPTION.
 
!    8. SUBROUTINES - FVRST1 CALLS ROUTINES FVRS1A, FVRS1B, FVRS1C,
!                     FVRS1D, FVRS1E, GMMATD, PRETRD, TRANSD, PRETAB,
!                     TAB AND OTHER STANDARD NASTRAN UTILITY ROUTINES.
!                     GINO ROUTINES.
 
!    9. DESIGN REQUIREMENTS -
 
!         (1) OPEN CORE IS DEFINED AT /ZZFVR1/.
!         (2) NO SCRATCH FILES ARE USED.
!         (3) FVRST1 RESIDES IN LINKNS07
!         (4) OPEN CORE FOR 5 BUFFERS PLUS 14*NCSTM  PLUS NTYPE*NROW OF
!             MGG IS REQUIRED.
 
!          NOTE - (1) NTYPE = 1 IF MGG IS REAL SP
!                     NTYPE = 2 IF MGG IS REAL DP
 
!   10. DIAGNOSTIC MESSAGES -
 
!         THE FOLLOWING MESSAGES MAY BE ISSUED - 3001,3002,3003,3008
!                                                AND 3031.
 
 
 LOGICAL :: modfrl
 INTEGER :: casecc,bgpdt,cstm,dit,frl,frlx,b1gg,basexg,pdzero,  &
     cycio,fkmax,bxtid,bxptid,bytid,byptid,bztid,  &
     bzptid,itlist(13),itid(6),frqset,case(14)
 DOUBLE PRECISION :: z,a(3,3),b(3,3),c(3,3),row(3),ta(3,3),avgm,  &
     dpi,dtwopi,dradeg,ddegra,d4pisq
 
 DIMENSION       mcbb1(7),mcbm1(7),mcbm2(7),mcb(7),coord(4),  &
     modnam(3),zs(1),iz(1),mcb1(7),mcb2(7),row2(3)
 
 COMMON /BLANK / nomgg,cycio,nsegs,kmax,fkmax,bxtid,bxptid,  &
     bytid,byptid,bztid,bzptid,nobasx,nofreq,omega
 COMMON /zzzzzz/ z(1)
 COMMON /system/ nbuf,nout,nerr
 COMMON /unpakx/ in1,nf1,nl1,incr
 COMMON /packx / in,iout,nf,nl,incr1
 COMMON /condad/ dpi,dtwopi,dradeg,ddegra,d4pisq
 
 EQUIVALENCE     (coord(1),ncrd),(z(1),zs(1)),(z(1),iz(1)),  &
     (mcb(1),mcb1(1)),(mcbm1(1),mcb2(1)), (itid(1),bxtid)
 
 DATA    casecc, bgpdt, cstm, dit, frl, mgg     /  &
     101,    102,   103,  104, 105, 106     /
 DATA    frlx, b1gg, m1gg, m2gg, basexg, pdzero /  &
     201,  202,  203,  204,  205   , 206    /
 DATA    modnam / 4HFRL ,  4HFVRS,4HTR1         /
 DATA    itlist / 4,  1105,11,1, 1205,12,2, 1305,13,3, 1405,14,4 /
!     LOCATE  CODES FOR -  TABLED1    TABLED2    TABLED3    TABLED4
 
!     CALCULATE PARAMETERS
 
!     TEST TO SEE IF BASEXG IS TO BE GENERATED.
 
 nobasx = -1
 IF (nomgg == -1 .OR.  cycio /= -1 .OR. nofreq == -1) GO TO 10
 IF (bxtid == -1 .AND. bytid == -1 .AND. bztid == -1) GO TO 10
 nobasx =  1
 10  CONTINUE
 
 IF (cycio /= -1) GO TO 25
 
!     DETERMINE FKMAX
 
 IF (MOD(nsegs,2) /= 0) GO TO 23
 IF (kmax == nsegs/2) GO TO 24
 23  fkmax = 2*kmax + 1
 GO TO 25
 24  fkmax = nsegs
 
!     TEST TO SEE IF ANY DATA BLOCKS ARE TO BE GENERATED.
 
 25  IF (nomgg == -1) GO TO 1000
 IF (omega == 0.0 .AND. (cycio /= -1 .OR. nofreq == -1) .AND.  &
     (bxtid == -1 .AND. bytid == -1 .AND. bztid == -1)) GO TO 1000
 
!     TEST TRAILER OF MGG TO SEE IF PURGED
 
 mcb(1) = mgg
 CALL rdtrl (mcb)
 nfile = mgg
 IF (mcb(1) <= 0) GO TO 902
 
!     COLUMN COUNT FOR MGG READ CHECK
 
 ncolc = mcb(2)
 nrowc = mcb(3)
 nform = mcb(4)
 ntype = mcb(5)
 
 nz = korsz(z)
 
!     ALLOCATE BUFFERS
 
!     MGG,CSTM (IBUF1 IS NBUF+1 LONG)
 
 ibuf1 = nz - nbuf
 
!     BGPDT
 
 ibuf2 = ibuf1 - nbuf
 
!     B1GG
 
 ibuf3 = ibuf2 - nbuf
 
!     M1GG
 
 ibuf4 = ibuf3 - nbuf
 
!     M2GG
 
 ibuf5 = ibuf4 - nbuf
 IF (omega == 0.0) ibuf5 = ibuf3
 
!     CALCULATE LENGTH OF OPEN CORE
 
 nz = ibuf5 - 1
 
!     PROCESS CSTM DATA BLOCK
 
 nfile  = cstm
 mcb(1) = cstm
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) GO TO 61
 
!     NO. OF COORDINATE SYSTEMS
 
 ncsym = mcb(3)
 lcstm = 14*ncsym
 
!     CSTM TABLE
 
 icstm = ibuf5 - lcstm
 nz    = icstm - 1
 
!     CORE FOR ENOUGH CORE FOR CSTM
 
 IF (nz < 0) GO TO 901
 
!     CORE CHECK FULL COLUMN OF MGG READ ASSUMED
 
 IF (nz-ntype*nrowc < 0) GO TO 901
 CALL gopen (cstm,zs(ibuf1),0)
 CALL READ (*903,*904,cstm,zs(icstm),lcstm,1,nwds)
 CALL pretrd (zs(icstm),lcstm)
 CALL CLOSE (cstm,1)
 GO TO 64
 
!     CORE CHECK NO CSTM
 
 61  IF (nz-ntype*nrowc < 0) GO TO 901
 64  CONTINUE
 
!     BGPDT TABLE
 
 mcb(1) = bgpdt
 CALL rdtrl (mcb)
 nfile = bgpdt
 IF (mcb(1) <= 0) GO TO 902
 
!     NO. OF GRID POINTS AND SCALAR POINTS READ CHECK FOR BGPDT
 
 ngrid = mcb(2)
 CALL gopen (bgpdt,zs(ibuf2),0)
 
!     OPEN MGG AND OUTPUT MATRICES
 
 CALL gopen (mgg,zs(ibuf1),0)
 IF (omega == 0.0) GO TO 65
 CALL gopen (b1gg,zs(ibuf3),1)
 mcbb1(1) = b1gg
 mcbb1(2) = 0
 mcbb1(3) = nrowc
 mcbb1(4) = 1
 mcbb1(5) = ntype
 mcbb1(6) = 0
 mcbb1(7) = 0
 CALL gopen (m1gg,zs(ibuf4),1)
 mcbm1(1) = m1gg
 mcbm1(2) = 0
 mcbm1(3) = nrowc
 mcbm1(4) = nform
 mcbm1(5) = ntype
 mcbm1(6) = 0
 mcbm1(7) = 0
 65  IF (nobasx == -1) GO TO 66
 CALL gopen (m2gg,zs(ibuf5),1)
 mcbm2(1) = m2gg
 mcbm2(2) = 0
 mcbm2(3) = nrowc
 mcbm2(4) = 1
 mcbm2(5) = ntype
 mcbm2(6) = 0
 mcbm2(7) = 0
 
!     SET UP PACK AND UNPACK TERMS
 
 66  in1  = ntype
 in   = 2
 iout = ntype
 incr = 1
 incr1= 1
 
!     READ INTERNAL SORT BGPDT PICK UP CID,X,Y,Z
 
 ndof = 0
 70  CALL READ (*903,*800,bgpdt,coord,4,0,m)
 ndof = ndof + 1
 IF (ncrd /= -1) GO TO 79
 
!     SCALAR POINT-UNPACK ONE COL OF MGG
!     SAVE DIAGONAL TERM
 
 nf1  = 0
 CALL unpack (*76,mgg,z)
 nrow = ndof - nf1 + 1
 nterm= nl1 - nf1 + 1
 IF (nrow < 1 .OR. nrow > nterm) GO TO 76
 IF (ntype == 1) row(1) = zs(nrow)
 IF (ntype == 2) row(1) = z(nrow)
 nf = ndof
 nl = ndof
 GO TO 77
 
!     OUT OF RANGE OF NON-ZERO BAND
 
 76  row(1) = 0.0
 nf = 1
 nl = 1
 
!     NOW PUT DIAGONAL ELEMENT INTO OUTPUT MATRICES
 
 77  IF (omega == 0.0) GO TO 78
 CALL pack (row,m1gg,mcbm1)
 CALL pack (row,b1gg,mcbb1)
 78  IF (nobasx == -1) GO TO 70
 CALL pack (row,m2gg,mcbm2)
 GO TO 70
 
!     UNPACK 3 COL OF MGG AND SAVE DIAGONAL TERMS
 
 79  DO  i = 1,3
   DO  j = 1,3
     a(i,j) = 0.0
   END DO
 END DO
 DO  i = 1,3
   nf1 = 0
   CALL unpack (*95,mgg,z)
   
!     LOCATE DIAGONAL ELEMENT IN COL-NROW
   
   nrow  = ndof - nf1 + i
   nterm = nl1  - nf1 + 1
   IF (nrow < 1 .OR. nrow > nterm) GO TO 95
   IF (ntype == 1) a(i,i) = zs(nrow)
   IF (ntype == 2) a(i,i) = z(nrow)
   CYCLE
   
!     OUT OF RANGE OF NON-ZERO ELEMENT BAND
   
   95  a(i,i) = 0.0
 END DO
 
!     NOW TRANSFORM FROM LOCAL(GLOBAL) TO BASIC
 
 IF (ncrd /= 0) GO TO 150
 
!     ALREADY IN BASIC COORDINATES
 
 avgm = (a(1,1) + a(2,2) + a(3,3))/3.0
 GO TO 161
 
!     SELECT TRANSFORMATION MATRIX-TA
 
 150  CALL transd (coord,ta)
 CALL gmmatd (ta,3,3,0,a,3,3,0,b)
 CALL gmmatd (b,3,3,0,ta,3,3,1,c)
 
!     C-IS NOW IN BASIC COORDINATES-ROW,WISE
 
 avgm = (c(1,1) + c(2,2) + c(3,3))/3.0
 
 161  IF (omega == 0.0) GO TO 307
 
!     PROCESS M1GG
 
 DO  i = 1,3
   DO  j = 1,3
     a(i,j) = 0.0
   END DO
 END DO
 
 a(2,2) = avgm
 a(3,3) = avgm
 IF (ncrd /= 0) GO TO 170
 
 DO  i = 1,3
   DO  j = 1,3
     c(i,j) = a(i,j)
   END DO
 END DO
 
 GO TO 180
 
!     TRANSFORM TO GLOBAL(LOCAL) FROM BASIC
 
 170  CALL gmmatd (ta,3,3,1,a,3,3,0,b)
 CALL gmmatd (b,3,3,0,ta,3,3,0,c)
 
!     C- IS NOW M1-11 ROW WISE
 
 180  DO  i = 1,3
   DO  k = 1,3
     row(k) = c(i,k)
   END DO
   nf = ndof
   nl = ndof + 2
   CALL pack (row,m1gg,mcbm1)
 END DO
 
!     WRITE OUT 3 NULL COLUMNS
 
 row(1) = 0.0
 DO  k = 1,3
   nf = 1
   nl = 1
   CALL pack (row,m1gg,mcbm1)
 END DO
 
!     NOW TAKE CARE OF B1GG
 
 IF (ncrd /= 0) GO TO 240
 
 DO  i = 1,3
   DO  j = 1,3
     c(i,j) = 0.0
   END DO
 END DO
 
 c(3,2) =-avgm
 c(2,3) = avgm
 GO TO 250
 
 240  DO  i = 1,3
   245  DO  j = 1,3
     a(i,j) = 0.0
   END DO
 END DO
 
 a(3,2) =-avgm
 a(2,3) = avgm
 
!     TRANSFORM TO GLOBAL(LOCAL) FROM BASIC
 
 CALL gmmatd (ta,3,3,1,a,3,3,0,b)
 CALL gmmatd (b,3,3,0,ta,3,3,0,c)
 
!     C-IS NOW B1-11 ROW WISE
 
 250  CONTINUE
 
 DO  i = 1,3
   DO  k = 1,3
     row(k) = c(i,k)
   END DO
   nf = ndof
   nl = ndof + 2
   CALL pack (row,b1gg,mcbb1)
 END DO
 
!     WRITE OUT 3 NULL COLUMNS
 
 row(1) = 0.0
 
 DO  i = 1,3
   nf = 1
   nl = 1
   CALL pack (row,b1gg,mcbb1)
 END DO
 
 307  IF (nobasx == -1) GO TO 407
 
!     NOW PROCESS M2GG
 
 IF (ncrd /= 0) GO TO 340
 
 DO  i = 1,3
   DO  j = 1,3
     c(i,j) = 0.0
   END DO
 END DO
 
 c(1,1) = avgm
 c(2,2) = avgm
 c(3,3) = avgm
 c(3,2) = avgm
 c(2,3) =-avgm
 GO TO 350
 
 340  DO  i = 1,3
   DO  j = 1,3
     a(i,j) = 0.0
   END DO
 END DO
 
 a(1,1) = avgm
 a(2,2) = avgm
 a(3,3) = avgm
 a(3,2) = avgm
 a(2,3) =-avgm
 
!     TRANSFORM TO GLOBAL(LOCAL) FROM BASIC
 
 CALL gmmatd (ta,3,3,1,a,3,3,0,c)
 
!     C-IS NOW M2-11 ROW WISE
 
 350  CONTINUE
 
 DO  i = 1,3
   DO  k = 1,3
     row(k) = c(i,k)
   END DO
   nf = ndof
   nl = ndof + 2
   CALL pack (row,m2gg,mcbm2)
 END DO
 
!     WRITE OUT 3 NULL COLUMNS
 
 row(1) = 0.0
 
 DO  i = 1,3
   nf = 1
   nl = 1
   CALL pack (row,m2gg,mcbm2)
 END DO
 
!     SPACE DOWN 3 COL IN MGG
 
 407  CONTINUE
 ndof  = ndof + 5
 nfile = mgg
 CALL fwdrec (*903,mgg)
 CALL fwdrec (*903,mgg)
 CALL fwdrec (*903,mgg)
 nfile = bgpdt
 GO TO 70
 
!     FINISH PROCESSING
 
 800  CALL CLOSE (mgg,1)
 IF (nobasx == -1) GO TO 802
 CALL CLOSE (m2gg,1)
 CALL wrttrl (mcbm2)
 802  IF (omega == 0.0) GO TO 805
 CALL CLOSE (b1gg,1)
 CALL CLOSE (m1gg,1)
 CALL wrttrl (mcbb1)
 CALL wrttrl (mcbm1)
 805  CALL CLOSE (bgpdt,1)
 
!     BEGIN PROCESSING OF FRLX, PDZERO AND BASEXG DATA BLOCKS.
 
 
!     TEST TO SEE IF BASEXG IS TO BE GENERATED.
 
 IF (nobasx == -1) GO TO 1000
 
!     RE-ESTABLISH LENGTH OF OPEN CORE FOR PHASE II PROCESSING
 
 nz = ibuf3 - 1
 
!     PROCESS FRL, FRLX AND PDZERO
 
 modfrl = .true.
 IF (omega == 0.0 .OR. bytid == -1 .AND. bztid == -1) modfrl = .false.
 
 nfile = frl
 mcb1(1) = frl
 CALL rdtrl (mcb1)
 nfsets = mcb1(2)
 ifrl = 1
 CALL OPEN (*902,frl,zs(ibuf1),0)
 
!     READ HEADER RECORD
 
 CALL READ (*903,*810,frl,iz(ifrl),nz,1,nwrds)
 GO TO 901
 
!     OPEN CASECC
 
 810  nfile = casecc
 CALL gopen (casecc,zs(ibuf2),0)
 
!     READ RECORD 1, WORD 14 (FREQUENCY SET ID)
 
 CALL READ (*903,*904,casecc,case,14,0,dummy)
 frqset = case(14)
 CALL CLOSE (casecc,1)
 
!     CHECK WHAT LOGICAL RECORD FRQSET IS IN FRL.
 
 mm = 0
 ii = ifrl + 2
 
 DO  i = ii,nwrds
   mm = mm + 1
   IF (iz(i) == frqset) GO TO 850
 END DO
 
!     FREQUENCY SET NOT FOUND.
 
 GO TO 905
 
!     MM IS LOGICAL RECORD NO. IN FRL FOR FRQSET.
 
 850  IF (.NOT.modfrl) GO TO 852
 CALL  OPEN (*902,frlx,zs(ibuf2),1)
 CALL WRITE (frlx,iz(ifrl),nwrds,1)
 CALL gopen (pdzero,zs(ibuf3),1)
 mcb2(1) = pdzero
 mcb2(2) = 0
 mcb2(3) = 0
 mcb2(4) = 1
 mcb2(5) = 1
 mcb2(6) = 0
 mcb2(7) = 0
 in      = 1
 iout    = 1
 incr1   = 1
 row2(1) = 0.0
 row2(2) = 1.0
 row2(3) = 0.0
 852  ifrl    = 1
 nfs     = 0
 nfsx    = 0
 nfile   = frl
 
 DO  i = 1,nfsets
   CALL READ (*903,*853,frl,zs(ifrl),nz,1,m)
   GO TO 901
   853  IF (i == mm) nfs = m
   IF (.NOT.modfrl .AND.i == mm) GO TO 865
   IF (.NOT.modfrl) CYCLE
   IF (i  /=  mm) GO TO 858
   
!     SET POINTERS FOR SORT INDEX ,   FRLX AND PDZERO ARRAYS.
   
   INDEX = ifrl + nfs
   ifrlx = INDEX + 3*nfs
   ipdz  = ifrlx + 3*nfs
   
!     RESET IFRL POINTER TO CONTINUE READING FRL RECORDS.
   
   ifrl  = ifrlx
   
!     CHECK CORE REQUIRED FOR EXPANDED FREQUENCY LIST AND SORT INDEX
   
   nz = nz - (ipdz + 3*nfs) + 1
   IF (nz < 0) GO TO 901
   
   ll = ifrlx - 1
   kkk= ipdz  - 1
   DO  ii = 1,nfs
     IF (zs(ii) == 0.0) GO TO 856
     DO  kk = 1,3
       kkk = kkk + 1
       zs(kkk) = row2(kk)
     END DO
     zs(ll+1) = ABS(zs(ii)-omega)
     zs(ll+2) = zs(ii)
     zs(ll+3) = ABS(zs(ii)+omega)
     ll = ll + 3
     CYCLE
     856  zs(ll+1) = 0.0
     zs(ll+2) = ABS(omega)
     kkk = kkk + 1
     zs(kkk) = row2(2)
     kkk = kkk + 1
     zs(kkk) = row2(1)
     ll  = ll + 2
   END DO
   
!     COMPUTE THE EXPANDED NUMBER OF FREQUIENCES, NFSX.
   
   nfsx = ll - ifrlx + 1
   
!     SORT EXPANDED W'S AND GET INDEX FOR SORTING BASE TABLE.
   
   CALL fvrs1e (zs(ifrlx),iz(INDEX),nfsx)
   CALL WRITE (frlx,zs(ifrlx),nfsx,1)
   CYCLE
   858  CALL WRITE (frlx,zs(ifrl),m,1)
 END DO
 IF (.NOT.modfrl) GO TO 865
 
!      FRLX IS A COPY OF FRL WITH THE SELECTED FREQUENCY SET, FRQSET,
!      EXPANDED.
 
 CALL CLOSE (frlx,1)
 mcb1(1) = frlx
 CALL wrttrl (mcb1)
 
!     SORT PDZERO BY INDEX JUST AS WAS DONE FOR FRLX
!     USE WORK   AT ZS(IFRLX)
!         INDEX  AT ZS(INDEX)   ALL NFSX LONG
!         PDZERO AT ZS(IPDZ)
 
 DO  kk = 1,nfsx
   loc = iz(INDEX+kk-1)
   zs(ifrlx+kk-1) = zs(ipdz+loc-1)
 END DO
 
!     NOW OUTPUT NFSX * FKMAX COLUMNS FOR PDZERO
 
 kkk = 0
 DO  kk = 1,fkmax
   DO  jj = 1,nfsx
     kkk = kkk + 1
     nf  = kkk
     nl  = kkk
     CALL pack (zs(ifrlx+jj-1),pdzero,mcb2)
   END DO
 END DO
 CALL CLOSE (pdzero,1)
 mcb2(3) = mcb2(2)
 CALL wrttrl (mcb2)
 865  CALL CLOSE (frl,1)
 
!     RE-ESTABLISH OPEN CORE FOR PHASE III AND
!     RESET POINTER TO ORIGINAL FREQUIENCIES.
 
 ifrl = 1
 nz   = ibuf1 - (nfs+nfsx) - 1
 
!     NFS  = THE ORIGINAL NUMBER OF FREQUIENCES
!     NFSX = THE EXPANDED NUMBER OF FREQUIENCES.
 
!     GENERATE BASE ACCELERATION MATRIX BASEXG.
 
 
!     BUILD A LIST OF UNIQUE TABLE IDS FOR PRETAB.
!     INITIALIZE THE TABLE WITH A ZERO ENTRY.
 
 itab = nfs + nfsx + 1
 ntabl=1
 k = itab + ntabl
 iz(k) = 0
 
!     WE HAVE A LIST OF TABLE ID'S TO CONSIDER
!     WE WANT ONLY A UNIQUE LIST OF TABLE ID'S GIVEN TO PRETAB
 
 loop872:  DO  i = 1,6
   iitid = itid(i)
   
!     SEARCH EXISTING LIST OF TABLE ID'S TO SEE IF IITID IS ALREADY IN
!     LIST
   
   IF (iitid <= 0 .OR. iitid > 9999999) CYCLE loop872
   DO  l = 1,ntabl
     ll = itab + l
     IF (iz(ll) == iitid) CYCLE loop872
   END DO
   
!     IITID WAS NOT AMONG EXISTING TABLE ID'S IN LIST,
!     IT'S A NEW TABLE ID,ADD IT TO LIST AND UPDATE LENGHT OF LIST
   
   ntabl = ntabl + 1
   k = itab + ntabl
   iz(k) = iitid
 END DO loop872
 
!     ALL TABLE ID'S HAVE BEEN PROCESSED,NOW PRETAB CAN BE CALLED
!     NTABL IS THE NUMBER OF TID'S IN THE LIST.
 
 iz(itab) = ntabl
 
!     ILTAB IS THE NEXT AVAILABLE LOCATION OF OPEN CORE FOR PRETAB.
 
 iltab = itab + ntabl + 1
 
!     COMPUTE LENGTH OF OPEN CORE AVAILABLE TO PRETAB.
 
 nztab = nz - ntabl - 1
 ltab  = 0
 CALL pretab (dit,zs(iltab),iz(iltab),zs(ibuf1),nztab,ltab, iz(itab),itlist)
 
!     COMPUTE LENGTH OF OPEN CORE AFTER PRETAB AND NEXT AVAILABLE LOC.
 
 nz   = nz - ltab
 next = iltab + ltab
 
!     ALLOCATE COMPLEX ARRAYS FOR BASEXG. START ON DOUBLE WORD BOUNDARY.
 
 IF (MOD(next,2) == 0) next = next + 1
 
!     DEFINE NFSX IF MODFRL IS FALSE.
 
 IF (.NOT.modfrl) nfsx = nfs
 
 n1 = next
 n2 = n1 + (3*nfsx)*2
 n3 = n2 + (3*nfsx)*2
 nt = n3 +    nrowc*2 - 1
 IF (nz < nt) GO TO 901
 CALL fvrs1a (zs(n1),zs(n2),zs(n3),zs(ifrl),zs(ibuf1),zs(INDEX),  &
     modfrl,basexg,nrowc,nfs,nfsx,fkmax,omega)
 GO TO 1000
 
!     ERROR PROCESSING
 
!     NOT ENOUGH CORE (ERROR 3008)
 
 901  ip1 = -8
 GO TO 999
 
!     DATA SET NOT IN FIST (ERROR 3001)
 
 902  ip1 = -1
 GO TO 999
 
!     EOF ENCOUNTERED (ERROR 3002)
 
 903  ip1 = -2
 GO TO 999
 
!     EOL ENCOUNTERED (ERROR 3003)
 
 904  ip1 = -3
 GO TO 999
 
!     FREQUENCY SET NOT FOUND IN FRL (ERROR 3031)
 
 905  CALL mesage (-31,frqset,modnam)
 GO TO 1000
 999  CALL mesage (ip1,nfile,modnam(2))
 
 1000 RETURN
END SUBROUTINE fvrst1
