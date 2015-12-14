SUBROUTINE shlsts (elid,pid,tlam,epsumi,epscmi)
     
!     TO PERFORM LAYER STRAIN, STRESS AND FORCE CALCULATIONS FOR THE
!     2-D SHELL ELEMENTS.
!     ONLY THE ELEMENT CENTER VALUES ARE CONSIDERED
 
!     INPUT :
!           ELID   - ELEMENT ID
!           PID    - COMPOSITE PROPERTY ID
!           TLAM   - AVERAGE ELEMENT THICKNESS
!           EPSUMI - UNCORRECTED STRAINS IN MATERIAL COORD. SYSTEM
!           EPSCMI - CORRECTED STRAINS IN MATERIAL COORD. SYSTEM
!          /CONDAS/- TRIGONOMETRIC CONSTATNTS
!          /OUTREQ/- OUTPUT REQUEST LOGICAL FLAGS
 
!     OUTPUT:
!           OUTPUT DATA ARE WRITTEN DIRECTLY TO EACH APPROPRIATE OUTPUT
!           FILE - OEF1L, OES1L/OES1AL
 
 
!     LAYER STRESS/STRAIN OUTPUT BLOCK FOR EACH CTRIA3 ELEMENT
 
!         1.    10*ELEMENT ID + DEVICE CODE
!         2.    NLAYER - NUMBER OF OUTPUT LAYERS
!         3.    TYPE OF FAILURE THEORY SELECTED
 
!         4.    LAYER ID
!        5-7.   LAYER STRESSES/STRAINS
!         8.    LAYER FAILURE INDEX, FI
!         9.    IFLAG1 = 1 IF FI.GE.0.999
!                      = 0 OTHERWISE
!       10-11.  INTERLAMINAR SHEAR STRESSES/STRAINS
!        12.    SHEAR BONDING INDEX, FB
!        13.    IFLAG2 = 1 IF FB.GE.0.999
!                      = 0 OTHERWISE
!         :
!         :     REPEAT 4-13 NLAYER TIMES FOR EACH LAYER
 
!       LAST-1. MAXIMUM FAILURE INDEX OF LAMINATE, FIMAX
!        LAST.  IFLAG3 = 1 IF FIMAX.GE.0.999
!                      = 0 OTHERWISE
 
 
!     FORCE OUTPUT BLOCK
 
!         1.    10*ELEMENT ID + DEVICE CODE
!        2-9.   FORCE RESULTANTS:
!                 MEMBRANE        BENDING     TRANSVERSE
!               -- FORCES --    - MOMENTS -  SHEAR FROCES
!               FX,  FY, FXY,   MX, MY, MXY,    VX, VY
 
 
 
 INTEGER, INTENT(IN)                      :: elid
 INTEGER, INTENT(IN OUT)                  :: pid
 REAL, INTENT(IN)                         :: tlam
 REAL, INTENT(IN)                         :: epsumi(6,1)
 REAL, INTENT(IN)                         :: epscmi(6,1)
 LOGICAL :: stress,strain,force,stsreq,stnreq,forreq,strcur,  &
     grids,vonms,layer,gridss,vonmss,layers, trnflx,nonmem,symlay,pcmp,pcmp1,pcmp2
 INTEGER :: elemid,oes1l,oes1al,oef1l,pcomp,pcomp1,  &
     pcomp2, pidloc,sym,symmem,souti,half,fthr, strinf,sdest,edest,fdest,iz(1)
 REAL :: strslr(3),trnsrr(2),epslr(3),ernsrr(2),findxr, fbondr,fimaxr,z
 REAL :: pi,twopi,raddeg, degrad,gg(9),ultstn(6),trans(9),stresl(3),  &
     epslcm(3),epslum(3),epslcf(3),epsluf(3),trnar(2),  &
     ernar(2),trnshr(2),ernshr(2),findex,fpmax,fbond,  &
     fbmax,fimax,fb(2),sb,v(2),ei(2),zbar(2), zk,zk1,zsubi,zref,theta,c,c2,s,s2,ti
 COMMON /condas/ pi,twopi,raddeg,degrad
 COMMON /outreq/ stsreq,stnreq,forreq,strcur,grids,vonms,layer  &
     ,               gridss,vonmss,layers
 COMMON /sdr2de/ ksdrde(200)
 COMMON /sdr2x2/ dum1(30),oes1l,oef1l
 COMMON /sdr2x7/ dum2(100),stres(69),dum3(31),forsul(37), dum4(163),strin(69)
 COMMON /sdr2c1/ ipcmp,npcmp,ipcmp1,npcmp1,ipcmp2,npcmp2
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1) ,iz(1)     ),(sdest,ksdrde( 26)),  &
     (fdest,ksdrde(33)),(edest,ksdrde(148)), (oes1l,oes1al    )
 DATA    symmem, mem, sym, pcomp, pcomp1, pcomp2, strinf /  &
     3     , 2  , 1  , 0    , 1     , 2     , 5      /
 
!     INITIALIZE
 
 zref  =-tlam/2.0
 findex= 0.0
 fbond = 0.0
 fpmax = 0.0
 fbmax = 0.0
 fimax = 0.0
 
 DO  ll = 1,2
   ernar(ll) = 0.0
   trnar(ll) = 0.0
   ernshr(ll)= 0.0
   trnshr(ll)= 0.0
 END DO
 
 force  = forreq .AND. layer
 stress = stsreq .AND. layer
 strain = stnreq .AND. layers
 
 itype  = -1
 lpcomp = ipcmp + npcmp + npcmp1 + npcmp2
 pcmp   = npcmp  > 0
 pcmp1  = npcmp1 > 0
 pcmp2  = npcmp2 > 0
 
 IF (.NOT.force) GO TO 20
 
!     WRITE FORCE RESULTANTS TO OEF1L IF REQUESTED
 
 elemid = 10*elid + fdest
 CALL WRITE (oef1l,elemid,1,0)
 CALL WRITE (oef1l,forsul(3),8,0)
 
!     FORCE REQUEST HAS BEEN PROCESSED. IF NO MORE REQUESTS WE ARE DONE.
!     IF NOT, PREPARE FOR OTHER REQUESTS.
!     ISSUE ERROR IF PCOMPI DATA HAS NOT BEEN READ INTO CORE.
 
 20 IF (.NOT.(stress .OR. strain)) GO TO 650
 
!     START WRITING STRESS/STRAIN OUTPUT TO OES1L/OES1AL
!     (NOTE - OES1L AND OES1AL ARE SAME FILE IN COSMIC/NASTRAN)
 
!     1.  10*ELEMENT ID + DEVICE CODE
 
 IF (lpcomp == ipcmp) GO TO 600
 elemid = 10*elid + sdest
 IF (strain) elemid = 10*elid + edest
 CALL WRITE (oes1l,elemid,1,0)
 
!     DETERMINE  IF INTERLAMINAR SHEAR STRESS CALCULATIONS ARE REQUIRED
!     BY CHECKING THE TRANSVERSE SHEAR STRESS RESULTANTS QX AND QY
 
 v(1) = forsul( 9)
 v(2) = forsul(10)
 trnflx = v(1) /= 0.0 .OR. v(2) /= 0.0
 
!     LOCATE PID BY PERFORMING A SEQUENTIAL SEARCH OF THE PCOMPI DATA
!     BLOCK WHICH IS IN CORE.
 
!     SEARCH FOR PID IN PCOMP DATA
 
 IF (.NOT.pcmp) GO TO 40
 ip = ipcmp
 IF (iz(ip) == pid) GO TO 110
 ipc11 = ipcmp1 - 1
 DO  ip = ipcmp,ipc11
   IF (iz(ip) == -1 .AND. ip < ipc11) IF (iz(ip+1)-pid) 30,100,30
 END DO
 
!     SEARCH FOR PID IN PCOMP1 DATA
 
 40 IF (.NOT.pcmp1) GO TO 60
 ip = ipcmp1
 IF (iz(ip) == pid) GO TO 140
 ipc21 = ipcmp2 - 1
 DO  ip = ipcmp1,ipc21
   IF (iz(ip) == -1 .AND. ip < ipc21) IF (iz(ip+1)-pid) 50,130,50
 END DO
 
!     SEARCH FOR PID IN PCOMP2 DATA
 
 60 IF (.NOT.pcmp2) GO TO 600
 ip = ipcmp2
 IF (iz(ip) == pid) GO TO 160
 lpc11 = lpcomp - 1
 DO  ip = ipcmp2,lpc11
   IF (iz(ip) == -1 .AND. ip < lpc11) IF (iz(ip+1)-pid) 70,150,70
 END DO
 
!     PID WAS NOT LOCATED; ISSUE ERROR
 
 GO TO 600
 
!     PID WAS LOCATED; DETERMINE TYPE
 
!     FOR PCOMP BULK DATA DETERMINE HOW MANY LAYERS HAVE THE STRESS/
!     STRAIN OUTPUT REQUEST (SOUTI).
!     FOR PCOMP1 OR PCOMP2 BULK DATA ENTRIES LAYER STRESSES/STRAINS ARE
!     OUTPUT FOR ALL LAYERS.
 
 100 ip = ip + 1
 110 itype  = pcomp
 pidloc = ip
 nlay   = iz(pidloc+1)
 nlayer = nlay
 nstrqt = 0
 DO  k = 1,nlay
   IF (iz(pidloc+8+4*k) == 1) nstrqt = nstrqt + 1
 END DO
 nlayer = nstrqt
 ipoint = pidloc + 8 + 4*nlay
 icontr = ipoint +    27*nlay
 GO TO 200
 
 130 ip = ip + 1
 140 itype  = pcomp1
 pidloc = ip
 nlay   = iz(pidloc+1)
 nlayer = nlay
 ipoint = pidloc +  8 +   nlay
 icontr = ipoint + 25 + 2*nlay
 GO TO 200
 
 150 ip = ip + 1
 160 itype  = pcomp2
 pidloc = ip
 nlay   = iz(pidloc+1)
 nlayer = nlay
 ipoint = pidloc +  8 + 2*nlay
 icontr = ipoint + 25 + 2*nlay
 
!     DETERMINE GENERAL COMPOSITE PROPERTY VALUES
 
!     LAMOPT - LAMINATION GENERATION OPTION
!            = ALL    = 0 (ALL PLYS SPECIFIED, DEFAULT)
!            = SYM    = 1 (SYMMETRIC)
!            = MEM    = 2 (MEMBRANE ONLY)
!            = SYMMEM = 3 (SYMMETRIC-MEMBRANE)
 
!     FTHR   - FAILURE THEORY
!            = 1    HILL
!            = 2    HOFFMAN
!            = 3    TSAI-WU
!            = 4    MAX-STRESS
!            = 5    MAX-STRAIN
 
!     SB     - SHEAR BONDING STRENGTH
 
 200 lamopt = iz(pidloc+8)
 fthr   = iz(pidloc+5)
 sb     =  z(pidloc+4)
 ei(1)  =  z(icontr+1)
 ei(2)  =  z(icontr+2)
 zbar(1)=  z(icontr+3)
 zbar(2)=  z(icontr+4)
 
 nonmem = lamopt /= mem .AND. lamopt /= symmem
 symlay = lamopt == sym .OR.  lamopt == symmem
 IF (symlay) nlayer = 2*nlayer
 IF (nlayer == 0) GO TO 650
 
!     CONTINUE TO WRITE LAYER-INDEPENDENT DATA TO OES1L/OES1AL
 
!     2.  NLAYER - NUMBER OF LAYERS FOR LAMINATE
!     3.  TYPE OF FAILURE THEORY SELECTED
 
 CALL WRITE (oes1l,nlayer,1,0)
 CALL WRITE (oes1l,fthr,1,0)
 
!     START THE LOOP OVER LAYERS
 
 zk = zref
 half = 1
 IF (symlay) half = 2
 
 DO  ihalf = 1,half
   DO  kk = 1,nlay
     k = kk
     IF (ihalf == 2) k = nlay + 1 - kk
     
!     OBTAIN LAYER K INFORMATION
!     - THE BOUNDARIES
!     - THE DISTANCE FROM THE REFERENCE SURFACE TO THE MIDDLE OF LAYER
!     - LAYER THICKNESS
!     - STRESS OUTPUT REQUEST (SOUTI) FOR PCOMP BULK DATA
!       (NOT SUPPORTED FOR PCOMP1 OR PCOMP2 BULK DATA)
     
     zk1 = zk
     IF (itype == pcomp ) zk = zk1 + z(pidloc + 6 + 4*k)
     IF (itype == pcomp1) zk = zk1 + z(pidloc + 7      )
     IF (itype == pcomp2) zk = zk1 + z(pidloc + 7 + 2*k)
     zsubi = (zk+zk1)/2.0
     ti = zk - zk1
     souti = 1
     IF (itype == pcomp) souti = iz(pidloc+8+4*k)
     
!     LAYER MATERIAL PROPERTIES
     
     DO  igi = 1,9
       gg(igi) = z(ipoint+igi)
     END DO
     
!     LAYER ULTIMATE STRENGTHS
     
     DO  ir = 1,6
       ultstn(ir) = z(ipoint+16+ir)
     END DO
     
!     LAYER ORIENTATION
     
     IF (itype == pcomp ) theta = z(pidloc + 7 + 4*k)
     IF (itype == pcomp1) theta = z(pidloc + 8 +   k)
     IF (itype == pcomp2) theta = z(pidloc + 8 + 2*k)
     theta = theta*degrad
     
!     BUILD THE STRAIN TENSOR TRANSFORMATION TO TRANSFORM
!     LAYER STRAINS FROM MATERIAL TO FIBER DIRECTION.
     
     c   = COS(theta)
     c2  = c*c
     s   = SIN(theta)
     s2  = s*s
     
     trans(1)  = c2
     trans(2)  = s2
     trans(3)  = c*s
     trans(4)  = s2
     trans(5)  = c2
     trans(6)  =-c*s
     trans(7)  =-2.0*c*s
     trans(8)  = 2.0*c*s
     trans(9)  = c2 - s2
     
!     CALCULATE THE CORRECTED AND UNCORRECTED STRAIN VECTORS AT ZSUBI
!     IN THE MATERIAL COORD. SYSTEM, THENTRANSFORM STRAINS FROM MATERIAL
!     TO FIBER COORD. SYSTEM AND CALCULATE THE LAYER STRESS VECTOR IN
!     THE FIBER COORD. SYSTEM
     
     DO  ir = 1,3
       epslcm(ir) = epscmi(ir,1) - zsubi*epscmi(ir+3,1)
       epslum(ir) = epsumi(ir,1) - zsubi*epsumi(ir+3,1)
     END DO
     
     CALL gmmats (trans(1),3,3,0, epslcm(1),3,1,0, epslcf(1))
     CALL gmmats (trans(1),3,3,0, epslum(1),3,1,0, epsluf(1))
     CALL gmmats (gg(1),3,3,0,    epslcf,3,1,0,    stresl(1))
     
     IF (fthr <= 0) GO TO 310
     
!     COMPUTE FAILURE INDEX FOR THIS LAYER AND THE MAXIMUM FAILURE INDEX
     
     IF (fthr == strinf) CALL failrs (fthr,ultstn,epsluf,findex)
     IF (fthr /= strinf) CALL failrs (fthr,ultstn,stresl,findex)
     IF (ABS(findex) >= ABS(fpmax)) fpmax = findex
     
     310 IF (.NOT.trnflx .OR. .NOT.nonmem) GO TO 350
     
!     CALCULATE INTERLAMINAR SHEAR STRESSES AND STRAINS
     
     IF (itype == pcomp ) icontr = ipoint + 25
     IF (itype == pcomp1) icontr = ipoint + 23 + 2*k
     IF (itype == pcomp2) icontr = ipoint + 23 + 2*k
     DO  ir = 1,2
       ernar(ir) = ernar(ir) + ti*(zbar(ir)-zsubi)
       trnar(ir) = trnar(ir) + ti*(zbar(ir)-zsubi)*z(icontr+ir)
     END DO
     
     DO  ir = 1,2
       trnshr(ir) = v(ir)*trnar(ir)/ei(ir)
       ernshr(ir) = v(ir)*ernar(ir)/ei(ir)
     END DO
     
     IF (sb <= 0.0) GO TO 350
     
!     CALCULATE SHEAR BONDING FAILURE INDEX, FB, AND THE MAX SHEAR
!     BONDING INDEX, FBMAX.
     
     DO  ir = 1,2
       fb(ir) = ABS(trnshr(ir))/sb
     END DO
     
     fbond = fb(1)
     IF (fb(2) > fb(1)) fbond = fb(2)
     IF (fbond >= fbmax) fbmax = fbond
     
     350 IF (souti == 0) GO TO 430
     
!     CONTINUE TO WRITE LAYER-DEPENDENT DATA TO OES1L AND OES1AL
     
!       4.   LAYER ID, LYRID
!     5,6,7. LAYER STRESSES/STRAINS
!       8.   LAYER FAILURE INDEX, FINDXR
!       9.   IFLAG1 (=1 IF FINDXR.GE.0.999, DEFAULT=0)
!     10,11. INTERLAMINAR SHEAR STRESSES/STRAINS
!      12.   SHEAR BONDING FAILURE INDEX, FBONDR
!      13.   IFLAG2 (=1 IF FBONDR.GE.0.999, DEFAULT=0)
!       :    REPEAT 4-13 FOR NUMBER OF LAYER WITH LAYER STRESS/STRAIN
!       :    REQUEST
     
     
     lyrid = k
     IF (ihalf == 2) lyrid = nlay + kk
     
     findxr = findex
     iflag1 = 0
     IF (ABS(findex) >= 0.999) iflag1 = 1
     
     fbondr = fbond
     iflag2 = 0
     IF (ABS(fbond ) >= 0.999) iflag2 = 1
     
     IF (.NOT.stress) GO TO 410
     DO  istr = 1,3
       strslr(istr) = stresl(istr)
     END DO
     trnsrr(1) = trnshr(1)
     trnsrr(2) = trnshr(2)
     CALL WRITE (oes1l,lyrid,1,0)
     CALL WRITE (oes1l,strslr(1),3,0)
     CALL WRITE (oes1l,findxr,1,0)
     CALL WRITE (oes1l,iflag1,1,0)
     CALL WRITE (oes1l,trnsrr(1),2,0)
     CALL WRITE (oes1l,fbondr,1,0)
     CALL WRITE (oes1l,iflag2,1,0)
     
     410 IF (.NOT.strain) GO TO 430
     DO  istr = 1,3
       epslr(istr)  = epsluf(istr)
     END DO
     ernsrr(1) = ernshr(1)
     ernsrr(2) = ernshr(2)
     CALL WRITE (oes1al,lyrid,1,0)
     CALL WRITE (oes1al,epslr(1),3,0)
     CALL WRITE (oes1al,findxr,1,0)
     CALL WRITE (oes1al,iflag1,1,0)
     CALL WRITE (oes1al,ernsrr(1),2,0)
     CALL WRITE (oes1al,fbondr,1,0)
     CALL WRITE (oes1al,iflag2,1,0)
     
!     UPDATE IPOINT FOR PCOMP BULK DATA ENTRY
     
     430 IF (itype /= pcomp) CYCLE
     IF (ihalf == 1 .AND. k /= nlay) ipoint = ipoint + 27
     IF (ihalf == 2) ipoint = ipoint - 27
   END DO
 END DO
 
!     END OF LOOP OVER LAYERS
 
 IF (fthr <= 0) GO TO 500
 
!     DETERMINE THE MAXIMUM FAILURE INDEX
 
 fimax = fpmax
 IF (fbmax > ABS(fpmax)) fimax = fbmax
 
!     CONTINUE TO OUTPUT THE MAXIMUM FAILURE INDEX TO OES1L/OES1AL
 
!     LAST-1.  MAXIMUM FAILURE INDEX OF LIMIATE, FIMAXR
!      LAST.   IFLAG3 (=1 IF FIMAXR.GE.0.999, DEFAULT=0)
 
 500 fimaxr = fimax
 iflag3 = 0
 IF (ABS(fimax) >= 0.999) iflag3 = 1
 
 CALL WRITE (oes1l,fimaxr,1,0)
 CALL WRITE (oes1l,iflag3,1,0)
 GO TO 650
 
 
!     ERROR MESSAGE
 
!     NO PCOMP, PCOMP1, PCOMP2 FOUND
 
 600 CALL mesage (-30,223,elid)
 
 650 RETURN
END SUBROUTINE shlsts
