SUBROUTINE frrd1d (pd,ull,lll,scr1,scr2,udv,ifr,nload,igood,nfreq)
     
!     ROUTINE SOLVES FOR UDV GIVEN ULL,LLL, AND PD
 
!     IF IGOOD = 1 DCOMP FAILED -- PUT ZERO SOLUTION VECTORS OUT
 
!     1. PULL LOADS FROM PD ONTO SCR1
!     2. SOLVE FOR UDV-S ON SCR2
!     3. STACK SOLVED LOADS ON UDV
 
 
 INTEGER, INTENT(IN)                      :: pd
 INTEGER, INTENT(IN)                      :: ull
 INTEGER, INTENT(IN)                      :: lll
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: udv
 INTEGER, INTENT(IN OUT)                  :: ifr
 INTEGER, INTENT(IN)                      :: nload
 INTEGER, INTENT(IN OUT)                  :: igood
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER :: sysbuf, fl,fu,fb,fx,prec,  &
     FILE,icore(1),udv1,mcb(7),NAME(2),mcore(1)
 COMMON /machin/ mach
 COMMON /unpakx/ it1,ii,jj,incr
 COMMON /packx / it2,it3,ii1,jj1,incr1
 COMMON /system/ ksystm(65)
 COMMON /gfbsx / fl(7),fu(7),fb(7),fx(7),nx,prec,ISIGN
 COMMON /fbsx  / mfl(7),mflt(7),mfb(7),mfx(7),mx,mprec,msign,iscrx
 COMMON /zzzzzz/ core(1)
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(55),iprec),  &
     (core(1),icore(1),mcore(1))
 DATA    NAME  / 4HFRRD,4H1D  /, ic / 0 /
 
 nx    = korsz(core)
 fb(1) = pd
 CALL rdtrl (fb)
 fx(1) = scr2
 IF (ifr == 1) fx(1) = udv
 fx(2) = nload
 fx(3) = fb(3)
 fx(4) = 2
 fx(5) = 2 + iprec
 it1   = fb(5)
 incr  = 1
 incr1 = 1
 it2   = it1
 it3   = 2 + iprec
 IF (igood == 1) GO TO 98
 
!     PULL LOADS FROM PD ONTO SCR1
 
 fu(1) = ull
 CALL rdtrl (fu)
 fl(1) = lll
 CALL rdtrl (fl)
 IF (nfreq == 1) GO TO 30
 nz = nx - sysbuf
 CALL gopen (pd,core(nz+1),0)
 CALL skprec (pd,ifr-1)
 nz = nz - sysbuf
 CALL gopen (scr1,core(nz+1),1)
 CALL makmcb (mcb,scr1,fb(3),2,it3)
 DO  i = 1,nload
   IF (i > 1) CALL skprec (pd,nfreq-1)
   ii = 0
   CALL unpack (*28,pd,core)
   ii1 = ii
   jj1 = jj
   22 CALL pack (core,scr1,mcb)
   CYCLE
   28 core(   1) = 0
   core(ic+2) = 0
   core(ic+3) = 0
   core(ic+4) = 0
   ii1 = 1
   jj1 = 1
   GO TO 22
 END DO
 CALL wrttrl (mcb)
 CALL CLOSE (pd,1)
 CALL CLOSE (scr1,1)
 
!     SET UP FOR GFBS
 
 fb(1) = scr1
 30 fb(2) = nload
 CALL wrttrl (fb)
 prec = 1
 IF (fb(5) == 2 .OR. fb(5) == 4) prec = 2
 ISIGN = 1
 IF (fu(1) < 0) GO TO 40
 CALL gfbs (core,core)
 CALL wrttrl (fx)
 GO TO 98
 
!     SET UP FOR FBS
 
 40 DO  i = 1,7
   mfl(i) = fl(i)
   
!     FBS DOES NOT USE THE MATRIX CONTROL BLOCK MFLT.
!     IF MFLT(1) EXISTS, SET ISCRX = MFLT(1) FILE FOR NEW FBS METHOD.
!     OTHERWISE SET ISCRX = 0, AND WE DO NOT HAVE A SCRATCH FILE FOR
!     NEW FBS. OLD FBS WILL BE USED.
   
   mfb(i) = fb(i)
   mfx(i) = fx(i)
 END DO
 mprec = prec
 msign = ISIGN
 mx    = korsz(mcore)
 iscrx = mflt(1)
 mcore(1) = mflt(1)
 CALL rdtrl (mcore(1))
 IF (mcore(1) <= 0) iscrx = 0
 CALL fbs (mcore,mcore)
 CALL wrttrl (mfx)
 98 icore(1) = 16777215
!                16777215 = '00FFFFFF'X
 iflag = 1
 
!     STACK LOADS ONTO UDV
 
 FILE = udv
 nz   = nx-sysbuf
 IF (ifr == 1) GO TO 300
 CALL OPEN (*900,udv,core(nz+1),0)
 fx(1) = udv
 CALL rdtrl (fx)
 IF (mach /= 1) GO TO 60
 50 CALL fwdrec (*51,udv)
 GO TO 50
 51 CALL bckrec (udv)
 CALL skprec (udv,1)
 GO TO 61
 60 CALL skpfil (udv,1)
 CALL skpfil (udv,-1)
 61 CALL CLOSE (udv,2)
 CALL OPEN (*900,udv,core(nz+1),3)
 
!     RESET TYPE FLAGS
 
 it1 = fx(5)
 it2 = it1
 it3 = it1
 IF (igood == 1) GO TO 101
 nz  = nz - sysbuf
 CALL gopen (scr2,core(nz+1),0)
 101 DO  i = 1,nload
   IF (igood == 1) GO TO 54
   ii  = 0
   CALL unpack (*54,scr2,core)
   ii1 = ii
   jj1 = jj
   53 CALL pack (core,udv,fx)
   CYCLE
   54 core(   1) = 0
   core(ic+2) = 0
   core(ic+3) = 0
   core(ic+4) = 0
   ii1 = 1
   jj1 = 1
   GO TO 53
 END DO
 CALL CLOSE (udv,1)
 IF (igood == 1) GO TO 56
 CALL CLOSE (scr2,1)
 56 CONTINUE
 CALL wrttrl (fx)
 350 RETURN
 
 300 IF (igood /= 1) GO TO 350
 CALL gopen (udv,core(nz+1),1)
 fx(2) = 0
 fx(6) = 0
 fx(7) = 0
 CALL wrttrl (fx)
 GO TO 101
 
!     ERROR MESAGES
 
 900 CALL mesage (-1,FILE,NAME)
 
 
 ENTRY frrd1e (udv1,udv,nload,nfreq)
!     ===================================
 
 nz = korsz(core) - sysbuf
 
!     ROUTINE REORDERS SOLUTIONS TO GET SORT BY LOADS
 
 FILE = udv1
 CALL OPEN (*900,udv1,core(nz+1),0)
 nz   = nz - sysbuf
 CALL gopen (udv,core(nz+1),1)
 FILE = udv1
 DO  i = 1,nload
   CALL skprec (udv1,i)
   DO  m = 1,nfreq
     ii  = 0
     CALL unpack (*420,udv1,core)
     ii1 = ii
     jj1 = jj
     421 CALL pack (core,udv,mcb)
     GO TO 422
     420 core(   1) = 0
     core(ic+2) = 0
     core(ic+3) = 0
     core(ic+4) = 0
     ii1 = 1
     jj1 = 1
     GO TO 421
     422 IF (m < nfreq) CALL skprec (udv1,nload-1)
   END DO
   CALL REWIND (udv1)
 END DO
 CALL CLOSE (udv1,1)
 CALL CLOSE (udv,1)
 fx(1) = udv1
 CALL rdtrl (fx)
 fx(1) = udv
 CALL wrttrl (fx)
 RETURN
END SUBROUTINE frrd1d
