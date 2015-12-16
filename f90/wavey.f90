SUBROUTINE wavey (ig,ild,NEW,nc,ic,kact,maxb,maxw,averw,sumw,  &
        rms,brms,jg)
     
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
!     COMPUTE WAVEFRONT AND ACTIVE COLUMN DATA -
!     MAXIMUM WAVEFRONT, AVERAGE WAVEFRONT, SUM OF ROW WAVEFRONTS,
!     SUM OF SQUARES OF ROW WAVEFRONTS, RMS WAVEFRONT, AND BANDWIDTH,
!     RMS BANDWIDTH, AND MINIMUM NODAL DEGREE.
!     DIAGONAL TERMS ARE INCLUDED.
 
!     IG     = CONNECTION TABLE
!     ILD(I) = NEW LABEL FOR NODE WITH ORIGINAL INTERNAL LABEL I
!     NEW(I) = INTERNAL LABEL CORRESPONDING TO NEW LABEL I
!              NEW AND ILD ARE INVERSES OF EACH OTHER
!     NC     = COMPONENT ID
!              IF NC.LE.0, USE ALL COMPONENTS.
!     IC(I)  = COMPONENT INDEX FOR ORIGINAL NODE I.
!     KACT(I)= LIST OF ACTIVE COLUMN FLAGS (UPDATED FOR EACH ROW)
!            = 1 IF COL I IS ACTIVE AT GIVEN ROW
!     MAXB   = BANDWIDTH
!     MAXW   = MAXIMUM WAVEFRONT
!     AVERW  = AVERAGE WAVEFRONT
!     SUMW   = SUM OF ROW WAVEFRONTS
!     SUMSQ  = SUM OF SQUARES OF ROW WAVEFRONTS
!     BSUMSQ = SUM OF SQUARES OF ROW BANDWIDTHS
!     RMS    = RMS WAVEFRONT
!     BRMS   = RMS BANDWIDTH
!     JG     = SCRATCH SPACE FOR BUNPAK
!     NN     = NUMBER OF NODES
!     MM     = MAX NODAL DEGREE
!     MINDEG = MINIMUM NODAL DEGREE
 
!     INPUT  - IG,ILD,NN,MM,NC,IC.
!     OUTPUT - NEW,KACT,MAXW,AVERW,SUMW,RMS,MAXB,BRMS,MINDEG
 
 
 INTEGER, INTENT(IN OUT)                  :: ig(1)
 INTEGER, INTENT(IN)                      :: ild(1)
 INTEGER, INTENT(OUT)                     :: NEW(1)
 INTEGER, INTENT(IN OUT)                  :: nc
 INTEGER, INTENT(IN)                      :: ic(1)
 INTEGER, INTENT(OUT)                     :: kact(1)
 INTEGER, INTENT(OUT)                     :: maxb
 INTEGER, INTENT(OUT)                     :: maxw
 REAL, INTENT(OUT)                        :: averw
 INTEGER, INTENT(OUT)                     :: sumw
 REAL, INTENT(OUT)                        :: rms
 REAL, INTENT(OUT)                        :: brms
 INTEGER, INTENT(IN)                      :: jg(1)
 
 DOUBLE PRECISION :: sumsq,    bsumsq
 
 COMMON /bands /  nn,       mm,       dum6s(6), mindeg
 
!     INITIALIZE WAVEFRONT DATA.
 
 maxb  = 0
 maxw  = 0
 sumw  = 0
 sumsq = 0.d0
 bsumsq= 0.d0
 averw = 0.
 rms   = 0.
 mindeg= MIN0(mindeg,mm)
 IF (nn*mm <= 0) RETURN
 
!     INITIALIZE NEW, THE INVERSE OF ILD
 
 IF (nc > 0) GO TO 8
 DO  i = 1,nn
   k = ild(i)
   IF (k <= 0) CYCLE
   NEW(k) = i
 END DO
 
!     INITIALIZE ACTIVE COLUMN FLAGS (1 FOR ACTIVE)
 
 8 DO  i = 1,nn
   kact(i) = 0
 END DO
 
!     COMPUTE WAVEFRONT DATA.
 
 iwave = 1
 kt = 0
 DO  i = 1,nn
   
!     COMPUTE NUMBER OF ACTIVE COLUMNS FOR ROW I
   
   k = NEW(i)
   IF (nc > 0) THEN
     GO TO    15
   ELSE
     GO TO    18
   END IF
   15 IF (k <= 0) CYCLE
   IF (nc-ic(k) == 0) THEN
     GO TO    18
   ELSE
     GO TO    40
   END IF
   18 kt = kt + 1
   CALL bunpak(ig,k,mm,jg)
   ib = 0
   DO  j = 1,mm
     l = jg(j)
     IF (l == 0) GO TO 30
     m  = ild(l)
     ib = MAX0(ib,i-m)
     IF (m <= i) CYCLE
     IF (kact(m) == 1) CYCLE
     iwave = iwave + 1
     kact(m) = 1
   END DO
   GO TO 35
   30 CONTINUE
   mindeg = MIN0(mindeg,j-1)
   35 CONTINUE
   
!     IB1 = ROW BANDWIDTH FOR ROW I (DIAGONAL INCLUDED)
   
   ib1 = ib + 1
   maxb = MAX0(maxb,ib1)
   IF (kact(i) == 1) iwave = iwave - 1
   
!     IWAVE = CURRENT NUMBER OF ACTIVE COLUMNS FOR ROW I
!             (DIAGONAL INCLUDED)
   
   maxw  = MAX0(maxw,iwave)
   sumw  = sumw + iwave
   wave  = FLOAT(iwave)
   sumsq = sumsq + wave*wave
   wave  = FLOAT(ib1)
   bsumsq= bsumsq + wave*wave
   40 CONTINUE
 END DO
 
 ann   = FLOAT(kt)
 averw = FLOAT(sumw)/ann
 rms   = SQRT(SNGL( sumsq)/ann)
 brms  = SQRT(SNGL(bsumsq)/ann)
 RETURN
END SUBROUTINE wavey
