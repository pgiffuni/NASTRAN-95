SUBROUTINE frrd1c (frl,frqset,mdd,bdd,kdd,ifr,ull,lll,scr1,scr2,  &
        scr3,scr4,igood)
!                                         (A)      (B)        (C)
!     THIS ROUTINE FORMS AND DECOMPOSES   KDD + I*W*BDD - W**2*MDD
!     WHERE  W = OMEGA, CYCLIC FREQ. AND I = SQUARE ROOT MINUS ONE
 
!     THE DECOMPOSITION ROUTINES ARE CALLED ACCORDING TO THE FOLLOWING
!     TABLE AS DETERMINED BY THE MATRIX RESULTING FROM THE ADDITION
 
!     IF MATRIX IS     COMPLEX SYMMETRIC    CALL SDCOMP
!                              UNSYMMETRIC  CALL CDCOMP
!                      REAL    SYMMETRIC    CALL SDCOMP
!                              UNSYMMETRIC  CALL DECOMP
 
 
 INTEGER, INTENT(IN OUT)                  :: frl
 REAL, INTENT(IN OUT)                     :: frqset
 INTEGER, INTENT(IN)                      :: mdd
 INTEGER, INTENT(IN)                      :: bdd
 INTEGER, INTENT(IN)                      :: kdd
 INTEGER, INTENT(IN OUT)                  :: ifr
 INTEGER, INTENT(IN OUT)                  :: ull
 INTEGER, INTENT(IN)                      :: lll
 INTEGER, INTENT(IN)                      :: scr1
 INTEGER, INTENT(IN)                      :: scr2
 INTEGER, INTENT(IN)                      :: scr3
 INTEGER, INTENT(IN)                      :: scr4
 INTEGER, INTENT(OUT)                     :: igood
 INTEGER :: fa,fl, fu,frq set, sr1, sr2, sysbuf, sr3,chlsky,NAME(2),  &
     mcore(1),icore(1)
 DOUBLE PRECISION :: det,minda,amcb(2),bmcb(2),cmcb(2),ddr,ddc,mindd, dett
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm
 COMMON /system/  ksystm(63)
 COMMON /cdcmpx/  fa(7),fl(7),fu(7),sr1,sr2,sr3,det(2),powr,nx, minda,ib,ibbar
 COMMON /dcompx/  ia(7),il(7),iu(7),iscr1,iscr2,iscr3,dett,ipow,  &
     ny,mindia,iib,iibb,icbr(3)
 COMMON /sfact /  mfa(7),mfl(7),mfc(7),m1fil,m2fil,mxx,ddr,ddc,  &
     power,m3fil,mindd,chlsky
 COMMON /saddx /  nomat,lcore,mcba(12),mcbb(12),mcbc(12),mcbd(12),  &
     mcbe(12),mx(7)
 COMMON /zzzzzz/  core(1)
 EQUIVALENCE      (mcore(1),core(1))
 EQUIVALENCE      (icore(1),core(1))
 EQUIVALENCE      (amcb(1),mcba(9)),(bmcb(1),mcbb(9)),  &
     (cmcb(1),mcbc(9)),(ksystm(2),nout)
 EQUIVALENCE      (ksystm(1),sysbuf),(ksystm(55),iprec)
 DATA    NAME  /  4HFRRD,4H1C  /
 
 
 nx = korsz(core)
 nz = nx - sysbuf
 
!     PICK UP CURRENT FREQUENCY
 
 CALL gopen  (frl,core(nz+1),0)
 CALL skprec (frl,frqset-1)
 CALL fread  (frl,core,ifr,1)
 w  = core(ifr)
 CALL CLOSE  (frl,1)
 
!     ADD MATRICES TOGETHER
 
 mcba(1) = kdd
 mcbb(1) = bdd
 mcbc(1) = mdd
 CALL rdtrl (mcba)
 CALL rdtrl (mcbb)
 CALL rdtrl (mcbc)
 IF (mcba(1) > 0 .AND. mcbc(1) > 0) GO TO 20
 WRITE  (nout,10) ufm
 10 FORMAT (a23,', EITHER STIFFNESS MATRIX OR MASS MATRIX IS MISSING')
 CALL mesage (-37,0,NAME)
 
 20 mcba(8) = 2
 mcbb(8) = 4
 mcbc(8) = 2
 amcb(1) = 1.0D0
 amcb(2) = 0.0D0
 bmcb(1) = 0.0D0
 bmcb(2) = w
 cmcb(1) =-w*w
 cmcb(2) = 0.0D0
 IF (mcbb(1) > 0) GO TO 30
 
!     NO BDD TO BE ADDED
 
 mcbb(1) = 0
 mcbb(8) = 0
 bmcb(2) = 0.0D0
 
 30 mx(1)   = scr3
 mx(2)   = mcba(2)
 mx(3)   = mcba(3)
 mx4a    = 6
 mx4b    = 6
 mx4c    = 6
 IF (mcba(1) > 0) mx4a = mcba(4)
 IF (mcbb(1) > 0) mx4b = mcbb(4)
 IF (mcbc(1) > 0) mx4c = mcbc(4)
 mx(4) = MIN0(mx4a,mx4b,mx4c)
 mx(5) = 2 + iprec
 IF (mcba(1) > 0 .AND. mcba(5) > 2) GO TO 40
 IF (mcbb(1) > 0) GO TO 40
 IF (mcbc(1) > 0 .AND. mcbc(5) > 2) GO TO 40
 mx(5) = iprec
 40 CONTINUE
 lcore = nx
 nomat = 3
 CALL sadd  (core,core)
 CALL wrttrl (mx)
 
!     SET UP TO DECOMPOSE MATRICES
 
 fa(1) = scr3
 CALL rdtrl (fa)
 igood = 1
 IF (fa(4) == 6) GO TO 120
 IF (fa(5) <= 2) GO TO 150
 fl(1) = lll
 fu(1) = ull
 DO  i = 2,7
   fl(i) = fa(i)
   fu(i) = fa(i)
 END DO
 fl(4) = 4
 fu(4) = 5
 sr1   = scr1
 sr2   = scr2
 sr3   = scr4
 CALL cdcomp (*100,core(1),core(1),core(1))
 igood = 0
 CALL wrttrl (fl)
 CALL wrttrl (fu)
 
!     FORCE RE-EVALUATION OF DECOMP PARAM IF W = 0.0
 
 60 IF (w /= 0.0) GO TO 70
 ib    = 0
 ibbar = 0
 70 RETURN
 
!     MATRIX SINGULR
 
 100 i = 5
 IF (w /= 0.0) i = -5
 CALL mesage (i,scr3,NAME)
 GO TO 60
 
!     USE SDCOMP TO PERFORM DECOMPOSITION
 
 120 mfa(1) = scr3
 mfl(1) = lll
 mfc(1) = ull
 DO  i = 2,7
   mfa(i) = fa(i)
   mfl(i) = fa(i)
   mfc(i) = fa(i)
 END DO
 mfl(4) = 4
 m1fil  = scr1
 m2fil  = scr2
 m3fil  = scr4
 mxx    = korsz(mcore)
 chlsky = 0
 CALL sdcomp (*100,mcore,mcore,mcore)
 igood  = 0
 
!     DIRECTION FOR FRRD1D TO USE  FBS RATHER THAN GFBS
 
 ull = -IABS(ull)
 
 CALL wrttrl (mfl)
 GO TO 60
 
!     USE DECOMP TO PERFORM DECOMPOSITION
 
 150 ia(1) = scr3
 il(1) = lll
 iu(1) = ull
 DO  i = 2,7
   ia(i) = fa(i)
   il(i) = fa(i)
   iu(i) = fa(i)
 END DO
 il(4) = 4
 iu(4) = 5
 iscr1 = scr1
 iscr2 = scr2
 iscr3 = scr4
 ny    = korsz(icore)
 CALL decomp (*100,icore,icore,icore)
 CALL wrttrl (il)
 CALL wrttrl (iu)
 igood = 0
 RETURN
END SUBROUTINE frrd1c
