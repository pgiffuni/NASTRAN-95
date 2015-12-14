SUBROUTINE bgrid
     
!     THIS ROUTINE COMPUTES PROBLEM SIZE, INTEGER PACKING FACTOR, AND
!     MAXGRD AND MAXDEG CONSTANTS.
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
 EXTERNAL        andf
 INTEGER :: grid(2),  seqgp,    geom1,    two,      andf,  &
     geom2,    geom4,    scr1,     rew,      sub(2), itrl(8)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,      uwm,      uim
 COMMON /machin/ machx
 COMMON /banda / ibuf1,    nompc,    nodep,    nopch,    norun,  &
     method,   icrit,    ngpts(2)
 COMMON /bandb / nbitin,   kor,      dum,      ngrid,    ipass,  &
     nw,       kdim,     nbpw,     irept
 COMMON /bandd / idum5d(5),nzero,    nel,      neq,      neqr
 COMMON /bands / nn,       mm,       dum2s(2), maxgrd,   maxdeg,  &
     kmod,     mach,     mindeg,   nedge,    mask
 COMMON /bandw / dum4w(4), i77
 COMMON /geomx / geom1,    geom2,    geom4,    scr1
 COMMON /system/ isys(100)
 COMMON /two   / two(1)
 COMMON /names / rdum(4),  rew,      norew
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (nout,isys(2))
 DATA            igeom1,   igeom2,   igeom4,   iscr1   /  &
     201,      208,      210,      301     /
 DATA            kdimx,    nelx,     neqx,     neqrx   /  &
     150,      0,        0,        0       /
 DATA            grid,     seqgp,    sub               /  &
     4501,45,  53,       4HBGRI,4HD        /
 
 IF (irept == 2) GO TO 100
 geom1 = igeom1
 geom2 = igeom2
 geom4 = igeom4
 scr1  = iscr1
 nel   = nelx
 neq   = neqx
 neqr  = neqrx
 ngrid = 0
 
!     BANDIT QUITS IF DMI CARDS ARE PRESENT. (CHK WAS DONE IN IFS2P)
!     RE-SET PROGRAM PARAMETERS IF USER REQUESTED VIA NASTRAN CARD.
 
 k = isys(i77)
 IF (k < 0) THEN
   GO TO   250
 ELSE IF (k == 0) THEN
   GO TO    30
 END IF
 10 IF (k == +9) GO TO 230
 DO  i = 1,7
   itrl(i) = MOD(k,10)
   k = k/10
 END DO
 IF (itrl(1) > 0 .AND. itrl(1) <= 4) icrit  = itrl(1)
 IF (itrl(2) > 0 .AND. itrl(2) <= 3) method = itrl(2) - 2
 nompc = itrl(3)
 IF (itrl(4) == 1) nodep = -nodep
 IF (itrl(5) == 1) nopch = -nopch
 IF (itrl(5) == 9) nopch = +9
 IF (itrl(6) == 1) norun = -norun
 IF (itrl(7) >= 2 .AND. itrl(7) <= 9) kdim = itrl(7)
 
 30 IF (norun == +1) GO TO 40
 
!     OPEN GEOM1 FILE AND CHECK THE PRESENCE OF ANY SEQGP CARD.  IF
!     ONE OR MORE IS PRESENT, ABORT BANDIT JOB.  OTHERWISE CONTINUE TO
!     COUNT HOW MANY GRID POINTS IN THE PROBLEM.
!     RESET GEOM1 TO THE BEGINNING OF GRID DATA FOR BSEQGP, AND CLOSE
!     GEOM1 WITHOUT REWINDING THE FILE
 
!     COMMENT FROM G.CHAN/SPERRY
!     IF TIME AND $ ALLOW, WE SHOULD MAKE USE OF THE SORTED GRID DATA
!     FROM GEOM1 FILE AND GET RID OF INV, INT, NORIG, ILD ARRAYS LATER.
!     THE SCATTERING TECHNEQUE (REALLY A HASHING METHOD) CAN BE REPLACED
!     BY A SIMPLE BINARY SEARCH. ROUTINES SCAT, BRIGIT, AND INTERN
!     COULD BE ELIMINATED.
 
 itrl(1) = geom1
 CALL rdtrl (itrl)
 j  = itrl(2) + itrl(3) + itrl(4) + itrl(5) + itrl(6) + itrl(7)
 IF (itrl(1) < 0 .OR. j == 0) GO TO 250
 k  = seqgp
 k1 = (k-1)/16
 k2 = k - 16*k1
 k  = andf(itrl(k1+2),two(k2+16))
 IF (k /= 0) GO TO 210
 
!     WE ASSUME THAT THE GRID POINT DATA IN GEOM1 AT THIS TIME IS NOT
!     SORTED. IF IT IS, WE CAN BLAST READ THE GRID POINT RECORD AND
!     TAKE THE LAST GRID POINT TO BE THE LARGEST GRID EXTERNAL NUMBER.
 
 40 CALL preloc (*170,z(ibuf1),geom1)
 CALL locate (*70,z(ibuf1),grid,k)
 MAX = 0
 50 CALL READ (*60,*60,geom1,itrl,8,0,k)
 ngrid = ngrid + 1
 IF (itrl(1) > MAX) MAX = itrl(1)
 GO TO 50
 60 CALL bckrec (geom1)
 70 CALL CLOSE (geom1,norew)
 
!     IF SPOINTS ARE PRESENT, ADD THEM TO THE GRID COUNT
 
 n = 0
 CALL preloc (*90,z(ibuf1),geom2)
 ngpts(1) = 5551
 ngpts(2) = 49
 CALL locate (*80,z(ibuf1),ngpts,k)
 CALL READ (*80,*80,geom2,z(1),ibuf1,1,n)
 80 CALL CLOSE (geom2,rew)
 90 ngpts(1) = ngrid
 ngpts(2) = n
 ngrid = ngrid + n
 
 IF (nopch == 9 .AND. ngrid == 1) ngrid = MAX
 100 IF (ngrid <=  0) GO TO 180
 IF (ngrid < 15) GO TO 160
 
!     SET WORD PACKING CONSTANT, NW - (NUMBER OF INTEGERS PER WORD)
!     MACHX =  1 DUMMY,   =  2 IBM 360/370, =  3 UNIVAC 1100, =  4 CDC,
!           =  5 VAX 780, =  6 DEC ULTRIX,  =  7 SUN,         =  8 AIX,
!           =  9 HP,      = 10 SILIC.GRAPH  = 11 MAC,         = 12 CRAY,
!           = 13 CONVEX,  = 14 NEC          = 15 FUJITSU,     = 16 DG,
!           = 17 AMDAHL   = 18 PRIME        = 19 486,         = 20 DUMMY
!           = 21 ALPHA    = 22 RESERVED
 
 GO TO (130,120,130,110,120,120,120,120,120,120,  &
     120,135,120,110,110,120,120,120,120,120, 120,120), machx
 110 nw = 6
 IF (ngrid >   510) nw = 5
 IF (ngrid >  2045) nw = 4
 IF (ngrid > 16380) nw = 3
 IF (ngrid > 524288) nw = 2
 GO TO 140
 120 nw = 2
 GO TO 140
 130 nw = 4
 IF (ngrid > 508) nw = 3
 IF (ngrid > 4095) nw = 2
 GO TO 140
 135 nw = 8
 IF (ngrid > 255) nw = 4
 
 140 nbitin = nbpw/nw
 mask   = 2**nbitin - 1
 
!     KDIM IS THE ARRAY DIMENSNION OF A SCRATCH ARRAY USED ONLY BY GPS
!     METHOD. IT IS 150 WORDS OR 10% OF TOTAL GRID POINT NUMBER. IF
!     USER SPECIFIED BANDTDIM = N, (WHERE N IS FROM 1 THRU 9), THE ARRAY
!     DIMENSION WILL BE N*10 PERCENT INSTEAD OF THE DEFAULT OF 10%.
 
 kdim = ngrid*kdim/10
 IF (method /= -1) kdim = MAX0(kdim,kdimx,ngrid/10)
 IF (method == -1) kdim = MIN0(kdim,kdimx,ngrid/10)
 n = ngrid
 IF (n < 10) n = 10
 
!     CALCULATE WIDTH MAXDEG AND EFFECTIVE LENGTH MAXGRD OF IG MATRIX.
 
 maxgrd = n
 kore   = kor
 150 maxdeg = ((((kore-4*kdim-8*maxgrd-5)*nw)/(maxgrd+nw))/nw)*nw
 maxdeg = MIN0(maxdeg,maxgrd-1)
 IF (maxdeg <= 0) GO TO 200
 j      = maxdeg*2.2
 kore   = kore - j
 IF (kor-j == kore) GO TO 150
 
!     INITIALIZE VARIABLES
 
 nn     = 0
 mm     = 0
 nedge  = 0
 ipass  = 0
 kmod   = 2*maxgrd - IFIX(2.3715*SQRT(FLOAT(maxgrd)))
 mindeg = 500000
 RETURN
 
!     ERROR OR QUIT
 
 160 WRITE  (nout,280) uim
 WRITE  (nout,270)
 GO TO  250
 170 CALL mesage (-1,geom1,sub)
 180 WRITE  (nout,280) uim
 WRITE  (nout,190)
 190 FORMAT (5X,25HTHE absence of grid cards)
 CALL CLOSE (geom1,rew)
 GO TO  250
 200 CALL mesage (-8,0,sub)
 210 WRITE  (nout,280) uim
 WRITE  (nout,220)
 220 FORMAT (5X,27HTHE presence of seqgp cards)
 GO TO  250
 230 WRITE  (nout,280) uim
 WRITE  (nout,240)
 240 FORMAT (5X,25HTHE presence of dmi cards)
 250 isys(i77) = 0
 IF (nopch > 0) isys(i77) = -2
 IF (isys(i77) /= -2) WRITE (nout,260)
 260 FORMAT (1H0,10X,'**NO ERRORS FOUND - EXECUTE NASTRAN PROGRAM**')
 270 FORMAT (5X,'SMALL PROBLEM SIZE')
 280 FORMAT (a29,' -  GRID-POINT RESEQUENCING PROCESSOR BANDIT IS ',  &
     'NOT USED DUE TO')
 RETURN
END SUBROUTINE bgrid
