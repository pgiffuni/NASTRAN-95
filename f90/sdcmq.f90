SUBROUTINE sdcmq (*,key,v1,v,dv1,dv,ic,z)
     
!     THIS SUBROUTINE CREATES A SCRATCH FILE OF QUEUED SINGULARITY
!     MESSAGES.  EACH MESSAGE IS A GINO RECORD DUE TO POSSIBLE CLOSE
!     WITHOUT REWIND.
!     THE -KEY- IS AS FOLLOWS,
!      1  - NULL COLUMN       - INPUT MATRIX
!      2  - ZERO DIAGONAL     - DECOMPOSED MATRIX.
!      3  - NEGATIVE DIAGONAL - DECOMPOSED MATRIX
!      4  - SINGULARITY TOLERANCE FAILURE - DECOMPOSED MATRIX.
!      5  - UNEXPECTED NULL COLUMN OR END OF COLUMN - ABORT IMMIDIATELY.
!      6  - NONCONSERVATIVE COLUMN  D/A.GT.1.001
!      7  - ZERO DIAGONAL     - INPUT MATRIX.
!     OTHER ARGUMENTS ARE
!      $  - NONSTANDARD RETURN IF DECOMPOSITION IS TO BE ABORTED.
!      Z  - OPEN CORE.  BUFFER LOCATIONS RELATIVE TO Z(1).
!      V  - RSP VALUE OF ENTRY IN ERROR (DV IS DOUBLE PRECISION).
!      V1 - INPUT VALUE OF DIAGONAL (DV1 IS DOUBLE PRECISION VERSION).
!      IC - COLUMN NUMBER IN ERROR.
!    THE ARGUMENTS ARE NOT CHANGED.
!    /SDCQ/ CONTAINS CONSTANT DATA.
!      FILCUR - CURRENT FILE USING BUFFER FOR SCRATCH FILE.  NEGATIVE IF
!               NONE, ZERO IF FILSCR IS TO REMAIN OPEN
!      STSCR  - GINO FILE STATUS FOR REOPENING FILCUR (1=READ,2=WRITE)
!      FILSCR - SCRATCH FILE NAME.
!      BUF    - BUFFER LOCATION RELATIVE TO Z(1).
!      NERR(2)- COUNT OF NUMBER OF ERROR CALLS ( (1)=ES, (2)=PD CHECK)
!      DIAGCK - EXIT FLAG FOR KEY=4 -- 0=NONFATAL, +N = MAX.-MESSAGES
!               WITHOUT ABORTING, -N = IMMEDIATE ABORT
!      IPREC  - 1=RSP, USE V. 2=RDP, USE DV.
!      PDEFCK - EXIT FLAG FOR KEY=3.  0 = NONFATAL IF -V, FATAL AT END
!               OF DECOMP FOR V=0.  +N = MAX-MESSAGES WITHOUT ABORTING.
!               -N = IMMEDIATE ABORT
!      NOGLEV - NOGO CODE.
!             = 0, NO FATAL ERRORS,
!             = 1, ABORT AT END OF DECOMP,
!             = 2, ABORT AT END OF PREPASS
!             = 3, ABORT,NONSTD RET.
!             = 4, INTERNAL ERRORS.  ABORT AT MAJOR CHECK-POINTS.
!-----
 
 INTEGER, INTENT(IN)                      :: key
 REAL, INTENT(IN)                         :: v1
 REAL, INTENT(IN OUT)                     :: v
 DOUBLE PRECISION, INTENT(IN)             :: dv1
 DOUBLE PRECISION, INTENT(IN OUT)         :: dv
 INTEGER, INTENT(IN)                      :: ic
 INTEGER, INTENT(IN OUT)                  :: z(1)
 LOGICAL :: opnscr,first
 INTEGER :: buf,diagck,filcur,filerr,filscr,iv(3),NAME(2), parm,pdefck,stscr
 
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm,swm
 COMMON /sdcq  /  nerr(2),noglev,buf,filscr,filcur,stscr,pdefck,  &
     diagck,diaget,iprec,parm(4),opnscr,first
 COMMON /sfact /  skpsf(32),ichly
 COMMON /names /  krd2,krr0,kwt3,kwr1, skpn,kcl2
 COMMON /system/  isb,iout
 EQUIVALENCE      (rv1,iv(2)),(rv,iv(3))
 DATA    NAME  /  4HSDCM,2HQ  /
 
 IF (filcur > 0) CALL CLOSE (filcur,kcl2)
 filerr = filscr
 IF (opnscr) GO TO 10
 IF (.NOT.first) CALL OPEN (*200,filscr,z(buf),kwt3)
 IF (first) CALL OPEN (*200,filscr,z(buf),kwr1)
 first  = .false.
 opnscr = .true.
 
 10 iv(1) = ic*10 + key
 IF (iprec == 1) GO TO 14
 rv  = dv
 rv1 = dv1
 GO TO 17
 14 rv  = v
 rv1 = v1
 17 CONTINUE
 CALL WRITE (filscr,iv,3,1)
 
!     CONVERT FILES TO ORIGINAL STATUS
 
 IF (filcur == 0) GO TO 20
 CALL CLOSE (filscr,kcl2)
 opnscr = .false.
 IF (filcur <= 0) GO TO 20
 filerr = filcur
 
!     READ MODE ON CURRENT FILE
 
 IF (stscr == 1) i = krd2
 
!     WRITE MODE ON CURRENT FILE
 
 IF (stscr == 2) i = kwt3
 CALL OPEN (*200,filcur,z(buf),i)
 
!     DETERMINE ABORT FLAG
 
 20 CONTINUE
 SELECT CASE ( key )
   CASE (    1)
     GO TO 65
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 50
   CASE (    5)
     GO TO 60
   CASE (    6)
     GO TO 50
   CASE (    7)
     GO TO 50
 END SELECT
 
!     ZERO DIAGONAL - DECOMPOSED MATRIX
 
 30 noglev = MAX0(noglev,1)
 IF (iprec == 1)  v = 1.0
 IF (iprec == 2) dv = 1.d0
 GO TO 70
 
!     NEGATIVE DIAGONAL
 
 40 CONTINUE
 IF (ichly  /= 1) GO TO 45
 IF (iprec  == 1)  v =-v
 IF (iprec  == 2) dv =-dv
 45 IF (pdefck == 0) GO TO 70
 noglev = MAX0(noglev,1)
 GO TO 70
 
!     ES SINGULARITY CHECK, DIAG-IN=0.0, NON-CONSERVATIVE MATRIX
 
 50 CONTINUE
 nerr(1) = nerr(1) + 1
 IF (diagck == 0) GO TO 100
 noglev = MAX0(noglev,1)
 IF (nerr(1) >= diagck) noglev = 3
 GO TO 100
 
!     UNEXPECTED NULL COLUMN
 
 60 noglev = 3
 GO TO 70
 
 65 noglev  = 2
 70 nerr(2) = nerr(2) + 1
 IF (nerr(2) > IABS(pdefck) .AND. pdefck /= 0)  noglev = 3
 
 100 CONTINUE
 IF (noglev == 3) RETURN 1
 RETURN
 
!     UNABLE TO USE FILES - WRITE GINO NUMBER. ABORT AT MAJOR DECOMP
!     CHCK
 
 200 CALL page2 (2)
 WRITE  (iout,210) swm,filerr,NAME,ic,key
 210 FORMAT (a27,' 2379, FILE',i8,' COULD NOT BE OPENED IN',a4,a1,  &
     '. COLUMN',i8,' SINGULAR, REASON',i3)
 parm(1) = -37
 parm(2) = filscr
 parm(3) = NAME(1)
 parm(4) = NAME(2)
 noglev  = 4
 GO TO 100
END SUBROUTINE sdcmq
