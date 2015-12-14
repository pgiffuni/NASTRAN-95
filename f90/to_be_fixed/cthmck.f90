SUBROUTINE cthmck (nt,num,nom,io,ig,ic,ideg,idis,iw,NEW,icc,ild,  &
        ipp,jump,un,nodesl)
     
!     THIS IS THE EXECUTIVE FOR THE CUTHILL-MCKEE GRID POINT RENUMBERING
!     STRATEGY.
!     91 VERSION, WITH REVERSED NEW SEQUENCE LOGIC
 
!     IN SAN ANTONIO, TEXAS, APRIL 27, 1989, THE DOUGLAS MICHEL NASTRAN
!     ACHIEVEMENT AWARD 1989, AN ANNUAL EVENT SPONSORED BY COSMIC AND
!     NASA, WAS GIVEN TO ELIZABETH H. CUTHILL, JAMES M. McKEE AND GORDON
!     C. EVERSTINE FOR THEIR TEAMWORK THAT CREATED BANDIT, A COMPUTER
!     PROGRAM THAT MINIMIZES THE BANDWIDTHS OF NASTRAN MATRICES. THE
!     WIDOW OF DR. McKEE AND HIS FAMILY RECEIVED THE AWARD FOR HIM.
!     DRS. CUTHILL AND EVERSTINE RECEIVED THEIR AWARDS PERSONALLY.
 
!     THE PRINCIPAL INPUTS ARE THE CONNECTIVITY MATRIX IG AND THE NUMBER
!     OF GRID POINTS (NODES) NN.
 
!     INPUT   - NT,NUM,NOM,IO,IP,IG,NN,MAXGRD,ILD
!     OUTPUT  - NEW,ILD,MM,IH0,IHE,KORIG,KNEW,NCM
!     SCRATCH - IC,IDEG,IDIS,IW,ICC,IPP
 
!     SET FOLLOWING DIMENSIONS IN CALLING PROGRAM -
!     IG(II1,M),IC(L),IDEG(L),IDIS(L),IW(L),NEW(L),ICC(L),ILD(L),IP(M)
 
!     L     = HAS THE DIMENSION OF MAXGRD
!             (NEW) MAXGRD EXCEEDS NUMBER OF GRID POINTS
!     II1   = MAXGRD/(PACKING DENSITY IN INTEGERS/WORD)
!           = ROW DIMENSION OF IG
!     M     = MAX NODAL DEGREE DIVIDED BY INTEGER PACKING FACTOR
!             (NEW) EXCEEDS MAX NODAL DEGREE
!     NT    = MAX NUMBER OF STARTING NODES TO BE CONSIDERED (=80)
!     NUM AND NOM GIVE THE FRACTION OF THE RANGE FROM MIN DEGREE TO MAX
!             DEGREE TO CONSIDER FOR STARTING NODES (NUM=1, NOM=2)
!     IO    = RE-SEQUENCING CRITERION , SET BY BANDIT -
!           = 1, RMS WAVEFRONT
!           = 2, BANDWIDTH
!           = 3, PROFILE. (PROFILE IS BANDWIDTH SUM OF ALL ROWS)
!           = 4, WAVEFRONT (MAX)
!     IG(I,J) CONTAINS THE GRID POINT LABEL FOR THE JTH NODE ADJACENT
!             TO NODE I  (THE CONNECTIVITY MATRIX).
!             THE CONNECTION OF A NODE TO ITSELF IS NOT LISTED.
!     NN    = NUMBER OF GRID POINTS (NODES)
!     MM    = COLUMN DIMENSION OF IG ON INPUT,
!             MAX NODAL DEGREE ON OUTPUT
!     MAXGRD= EFFECTIVE IG ROW DIMENSION (NEGLECTING INTEGER PACKING)
!     NEW(I)= OLD LABEL FOR GRID POINT NOW LABELLED I
!     ILD(I)= NEW LABEL FOR GRID POINT ORIGINALLY LABELLED I
!             ILD AND NEW ARE INVERSES
!     ILD MUST BE INPUT TO CTHMCK TO INDICATE AN INITIAL SEQUENCE.
!             NORMALLY, ON INPUT, SET ILD(I)=I FOR ALL I.
!     JUMP  = 1 IF RESEQUENCING ATTEMPTS RESULT IN NO IMPROVEMENT
!           = 0 OTHERWISE.
!     IH0   = ORIG PROFILE
!     IHE   = NEW PROFILE
!     KORIG = ORIG BANDWIDTH
!     KNEW  = NEW BW
!     NCM   = NUMBER OF COMPONENTS
!     NODESL IS SCRATCH SPACE.
 
!     IN CALLING PROGRAM, TRY  CALL CTHMCKL (80,1,2,2,1,...)
 
!     THE FOLLOWING SUBROUTINES WERE WRITTEN BY E. CUTHILL AND J. MCKEE
!     OF NSRDC -
!     DEGREE,DIAM,IDIST,KOMPNT,MAXBND,MAXDGR,MINDEG,RELABL,CTHMCK
!     CTHMCK WAS MODIFIED BY G.C. EVERSTINE, DTRC, AND
!        PUT INTO NASTRAN BY G.C. CHAN/UNISYS
 
 
 INTEGER, INTENT(IN OUT)                  :: nt
 INTEGER, INTENT(IN OUT)                  :: num
 INTEGER, INTENT(IN)                      :: nom
 INTEGER, INTENT(IN OUT)                  :: io
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(IN)                      :: ic(1)
 INTEGER, INTENT(IN)                      :: ideg(1)
 INTEGER, INTENT(IN OUT)                  :: idis(1)
 INTEGER, INTENT(IN)                      :: iw(1)
 INTEGER, INTENT(OUT)                     :: NEW(1)
 INTEGER, INTENT(IN)                      :: icc(1)
 INTEGER, INTENT(IN OUT)                  :: ild(1)
 INTEGER, INTENT(IN OUT)                  :: ipp(1)
 INTEGER, INTENT(OUT)                     :: jump
 REAL, INTENT(IN OUT)                     :: un(1)
 INTEGER, INTENT(IN OUT)                  :: nodesl(1)
 INTEGER :: sumw
 REAL :: im1,     im2
 
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,     uwm,     uim
 COMMON /banda / ibuf1,   nompc,   nodep,    nopch,    norun, method,  icrit
 COMMON /bandb / dum3b(3),ngrid,   dumb2(2), kdim
 COMMON /bandd / korig,   knew,    ih0,      ihe,      ncm
 COMMON /bands / nn,      mm,      ih,       ib,       maxgrd
 COMMON /bandw / maxw0,   rms0,    maxw1,    rms1,     i77, brms0,   brms1
 COMMON /system/ isys,    nout,    dum6y(6), nlpp
 
!     SET UP SCRATCH SPACE NODESL.
 
 idem  = kdim
 k2    = idem + 1
 iajdim= 3*idem
 
!     DETERMINE THE DEGREE OF EACH NODE, THE NUMBER OF COMPONENTS, NCM,
!     AND THE MAXIMUM DEGREE OF ANY NODE.
 
 CALL degree (ig,ideg,un)
 ncm  = kompnt(ig,ic,ideg,iw,icc,un)
 maxd = maxdgr(0,ic,ideg)
 mmc  = maxd
 
!     INITIALIZE NEW ARRAY FROM THE ILD ARRAY.
!     ILD MUST BE INPUT TO CUTHILL.
 
 DO  i = 1,nn
   k = ild(i)
   NEW(k) = i
 END DO
 
!     COMPUTE ORIGINAL BANDWIDTH, PROFILE, WAVEFRONT AND ACTIVE COLUMN
!     IH0 = ORIGINAL PROFILE,  IS = ORIGINAL BW
 
 CALL wavey (ig,ild,NEW,0,ic,iw,is,maxw,averw,sumw,rms,brms,un)
 ih    = sumw
 maxw0 = maxw
 rms0  = rms
 brms0 = brms
 korig = is
 ih0   = ih
 CALL page1
 i = method + 2
 WRITE  (nout,20) uim,icrit,i,nompc,nodep,nopch
 20 FORMAT (a29,'S FROM RESEQUENCING PROCESSOR - BANDIT     (CRI=',i2,  &
     ',  MTH=',i2,',  MPC=',i2,',  DEP=',i2,',  PCH=',i2,')',/)
 IF (nlpp <= 50) GO TO 50
 WRITE  (nout,30)
 30 FORMAT (31X,'BEFORE RESEQUENCING - - -')
 WRITE  (nout,40) is,ih,maxw,averw,rms,brms
 40 FORMAT (40X,'BANDWIDTH',i13,     /40X,'PROFILE',i15,  &
     /40X,'MAX WAVEFRONT',i9,  /40X,'AVG WAVEFRONT',f9.3,  &
     /40X,'RMS WAVEFRONT',f9.3,/40X,'RMS BANDWIDTH',f9.3)
 
!     COMPUTE NODAL DEGREE STATISTICS.
 
 50 CALL dist (ideg,ipp,median,modd)
 IF (method == +1) RETURN
 
!     INITIALIZE ILD AND NEW ARRAYS.
 
 jump  = 0
 DO  i = 1,nn
   NEW(i) = 0
   ild(i) = 0
 END DO
 
!     GENERATE NUMBERING SCHEME FOR EACH COMPONENT, NC.
 
 DO  nc = 1,ncm
   
!     DETERMINE THE RANGE OF DEGREES (MI TO MAD) OF NODES OF INTEREST.
!     MAKE SURE MAD DOES NOT EXCEED MEDIAN
   
   mi  = mindeg(nc,ic,ideg)
   mad = mi
   IF (nom == 0) GO TO 80
   ma  = maxdgr(nc,ic,ideg)
   mad = mi + ((ma-mi)*num)/nom
   mad = MIN0(mad,median-1)
   mad = MAX0(mad,mi)
   
!     DETERMINE BANDWIDTH OR SUM CRITERION FOR EACH NODE MEETING
!     SPECIFIED CONDITION.
   
   80 CALL diam (nc,mad,nl,nodesl,idem,maxlev,ig,ic,ideg,idis,iw,icc,un)
   jmax = MIN0(nt,nl)
   jmax = MAX0(jmax,1)
   im1  = 1.e+8
   im2  = im1
   
!     CHECK SEQUENCE FOR EACH STARTING NODE SELECTED, AND
!     COMPUTE NEW BANDWIDTH,PROFILE,WAVEFRONT DATA.
!     IB = BANDWIDTH, IH = PROFILE.
   
   DO  j = 1,jmax
     CALL relabl (1,nodesl(j),ig,ic,ideg,idis,iw,NEW,icc,ild,  &
         nodesl(k2),un,iajdim)
     CALL wavey (ig,ild,NEW,nc,ic,iw,ib,maxw,averw,sumw,rms,brms,un)
     IF (ngrid == -1) RETURN
     
     ih = sumw
     SELECT CASE ( io )
       CASE (    1)
         GO TO 220
       CASE (    2)
         GO TO 230
       CASE (    3)
         GO TO 240
       CASE (    4)
         GO TO 250
     END SELECT
     220 crit1 = rms
     crit2 = ih
     GO TO 260
     230 crit1 = ib
     crit2 = ih
     GO TO 260
     240 crit1 = ih
     crit2 = ib
     GO TO 260
     250 crit1 = maxw
     crit2 = rms
     260 IF (im1-crit1 < 0) THEN
       GO TO   300
     ELSE IF (im1-crit1 == 0) THEN
       GO TO   280
     END IF
     270 im1 = crit1
     im2 = crit2
     ij  = j
     CYCLE
     280 IF (im2 <= crit2) CYCLE
     im2 = crit2
     ij  = j
     
   END DO
   
!     RECOMPUTE SEQUENCE FOR STARTING NODE WHICH IS BEST FOR CRITERION
!     SELECTED.
   
   CALL relabl (1,nodesl(ij),ig,ic,ideg,idis,iw,NEW,icc,ild,  &
       nodesl(k2),un,iajdim)
   IF (ngrid == -1) RETURN
   
 END DO
 
!     DETERMINE NODES OF ZERO DEGREE AND STACK LAST, AND
!     COMPUTE BANDWIDTH, PROFILE AND WAVEFRONT DATA.
 
 CALL stack (ideg,NEW,ild,iw)
 CALL wavey (ig,ild,NEW,0,ic,iw,ib,maxw,averw,sumw,rms,brms,un)
 ih = sumw
 
 IF (nlpp <= 50) GO TO 350
 WRITE  (nout,320)
 320 FORMAT (/31X,'AFTER RESEQUENCING BY REVERSE CUTHILL-MCKEE (CM)',  &
     ' ALGORITHM - - -')
 WRITE (nout,40) ib,ih,maxw,averw,rms,brms
 
!     CHECK CM LABELING AGAINST ORIGINAL LABELING TO SEE IF BETTER.
!     IB = BANDWIDTH,  IH = PROFILE.
 
 350 SELECT CASE ( io )
   CASE (    1)
     GO TO 400
   CASE (    2)
     GO TO 410
   CASE (    3)
     GO TO 420
   CASE (    4)
     GO TO 430
 END SELECT
 400 im1   = rms0
 im2   = ih0
 crit1 = rms
 crit2 = ih
 GO TO 440
 410 im1   = is
 im2   = ih0
 crit1 = ib
 crit2 = ih
 GO TO 440
 420 im1   = ih0
 im2   = is
 crit1 = ih
 crit2 = ib
 GO TO 440
 430 im1   = maxw0
 im2   = rms0
 crit1 = maxw
 crit2 = rms
 440 IF (crit1-im1 < 0.0) THEN
   GO TO   480
 ELSE IF (crit1-im1 == 0.0) THEN
   GO TO   450
 ELSE
   GO TO   460
 END IF
 450 IF (crit2 < im2) GO TO 480
 
!     IF NO IMPROVEMENT RETURN TO ORIGINAL SEQUENCE.
 
 460 ib   = is
 ih   = ih0
 maxw = maxw0
 rms  = rms0
 brms = brms0
 DO  i = 1,nn
   ild(i) = i
   NEW(i) = i
 END DO
 jump = 1
 
!     SET FINAL VALUES OF B, P, RMS, W.
 
 480 knew = ib
 ihe  = ih
 maxw1= maxw
 rms1 = rms
 brms1= brms
 RETURN
END SUBROUTINE cthmck
