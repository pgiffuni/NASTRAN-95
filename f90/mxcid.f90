SUBROUTINE mxcid (*,z,mset,msze,nwds,uset,gpl,sil,buf1)
     
!     THIS SUBROUTINE CREATES AN ARRAY AT Z(1) OF LENGTH MSZE*NWDS
!     WHICH CONTAINS THE EXTERNAL ID*10 + COMPONENT AT Z(1,M) FOR
!     EACH DEGREE OF FREEDOM BELONGING TO SET -MSET-.
 
!     OPEN CORE IS Z(1) TO Z(BUF1-1).   TWO  BUFFERS NEEDED.
 
!     NONSTANDARD RETURN IF TASK NOT COMPLETED.
 
!     IF THIS IS A SUBSTRUCTURING PROBLEM, MXCIDS SHOULD BE CALLED
!     INSTEAD
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN OUT)                  :: z(1)
 INTEGER, INTENT(IN OUT)                  :: mset
 INTEGER, INTENT(IN)                      :: msze
 INTEGER, INTENT(OUT)                     :: nwds
 INTEGER, INTENT(IN)                      :: uset
 INTEGER, INTENT(IN)                      :: gpl
 INTEGER, INTENT(IN)                      :: sil
 INTEGER, INTENT(IN)                      :: buf1
 EXTERNAL        lshift,andf,orf
 INTEGER :: fnam(2),NAME(2),x(7)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /system/ nbufsz,koutp
 COMMON /names / krd,krdrw,kwr,kwrrw, kclrw,kcl,kweof
 COMMON /bitpos/ mask2(32),hset(32)
 COMMON /two   / itwo(32)
 DATA    nset  / 20   /
 DATA    NAME  , NONE / 4HMXCI,4HD   , 4H (no  /
 
!     ALLOCATE CORE - CHECK DATA FILE AVAILABILITY
 
 buf2 = buf1 + nbufsz
 IF (nwds <= 0) nwds = 1
 lgp  = msze*nwds + 1
 x(1) = sil
 CALL fname (sil,fnam)
 IF (fnam(1) == NONE) GO TO 220
 CALL rdtrl (x)
 ngp = x(2)
 lsil= lgp + ngp
 
!     SEVEN WORDS NEEDED IF SIL AND USET OUT OF CORE
 
 IF (lsil > buf1-7) GO TO 260
 
!     DETERMINE IF SIL (AND USET) FIT IN CORE
 
 luset = lsil + ngp
 x(1)  = uset
 CALL fname (uset,fnam)
 IF (fnam(1) == NONE) GO TO 220
 CALL rdtrl (x)
 ndof = x(3)
 l = orf(lshift(x(4),16),x(5))
 
 IF (luset+ndof > buf1) luset = 0
 IF (luset > buf1) lsil = 0
 
!     CHECK SET REQUEST
 
 DO  iset = 1,nset
   IF (hset(iset) == mset) GO TO 20
 END DO
 GO TO 240
 20 CONTINUE
 iset = mask2(iset)
 iset = itwo(iset)
 IF (andf(l,iset) == 0) GO TO 240
 
!     LOAD GPL INTO CORE
 
 x(1) = gpl
 CALL OPEN (*220,gpl,z(buf2),krdrw)
 CALL fread (gpl,0,0,1)
 CALL fread (gpl,z(lgp),ngp,0)
 CALL CLOSE (gpl,kcl)
 x(1) = sil
 CALL gopen (sil,z(buf1),krdrw)
 CALL gopen (uset,z(buf2),krdrw)
 
!     LOAD SIL AND USET IF POSSIBLE
 
 IF (lsil == 0) GO TO 30
 CALL fread (sil,z(lsil),ngp,0)
 CALL CLOSE (sil,kcl)
 sil1 = z(lsil)
 psil = lsil + 1
 i = ngp - 1
 GO TO 40
 30 CALL fread (sil,sil1,1,0)
 i = 1
 psil = lgp + ngp
 40 IF (luset == 0) GO TO 50
 CALL fread (uset,z(luset),ndof,0)
 CALL CLOSE (uset,kcl)
 puset = luset
 50 IF (luset == 0) puset = psil + i
 
!     PSIL POINTS SECOND SIL ENTRY IF SIL IN CORE, ELSE LOCATION TO USE
!     PUSET POINTS TO FIRST WORD USET, ELSE LOCATION IN Z TO USE
!     LSIL, LUSET ARE ZERO IF FILES NOT IN CORE.
!     LOOP ON NUMBER GRID POINTS - EXIT WHEN MSIZE ACHIEVED.
 
 mcount = 1
 
 DO   lll = 1,ngp
   IF (lll == ngp) GO TO 60
   IF (lsil /=   0) GO TO 70
   CALL fread (sil,z(psil),1,0)
   GO TO 70
   60 sil2 = ndof + 1
   GO TO 80
   70 sil2 = z(psil)
   IF (lsil /= 0) psil = psil + 1
   80 ndf = sil2 - sil1
   IF (ndf < 1 .OR. ndf > 6) GO TO 240
   
!     GET NDF WORDS FROM USET
   
   IF (luset == 0) CALL fread (uset,z(puset),ndf,0)
   
!     DETERMINE IF IN THE SET
   
   j = puset
   k = j + ndf - 1
   100 CONTINUE
   DO  i = j,k
     IF (andf(z(i),iset) /= 0) GO TO 120
   END DO
   GO TO 125
   
!     LOCATED A POINT IN THE SET
   
   120 CONTINUE
   ll = i - puset + 1
   l  = lgp + lll - 1
   IF (ndf == 1) ll = 0
   z(mcount) = z(l)*10 + ll
   mcount = mcount + nwds
   IF (mcount >= lgp) GO TO 310
   IF (i == k) GO TO 125
   j = i + 1
   GO TO 100
   125 IF (luset /= 0) puset = puset + ndf
   sil1 = sil2
 END DO
 
!     END OF ALL GRIDS AND MATRIX NOT FILLED - NEED IMMEDIATE MESSAGE.
 
 CALL page2 (2)
 WRITE  (koutp,210) swm,NAME
 210 FORMAT (a27,' 3016, MATRIX IS NOT IN PROPER FORM IN SUBROUTINE ', 2A4)
 GO TO 300
 
!     PURGED FILES
 
 220 CALL page2 (2)
 WRITE  (koutp,230) swm,x(1),NAME
 230 FORMAT (a27,' 3001, ATTEMPT TO OPEN DATA SET',i4,' IN SUBROUTINE',  &
     1X,2A4,' WHICH WAS NOT DEFINED IN FIST')
 GO TO 300
 
!     ILLEGAL INPUT
 
 240 CALL page2 (2)
 WRITE  (koutp,250) swm,NAME
 250 FORMAT (a27,' 3007, ILLEGAL INPUT TO SUBROUTINE ',2A4)
 GO TO 300
 
!     INSUFFICIENT CORE
 
 260 CALL page2 (2)
 WRITE  (koutp,270) swm,NAME
 270 FORMAT (a27,' 3008, INSUFFICIENT CORE AVAILABLE FOR SUBROUTINE ',  &
     2A4, 1H.)
 
 300 CONTINUE
 CALL CLOSE (sil ,kcl)
 CALL CLOSE (uset,kcl)
 RETURN 1
 310 CONTINUE
 CALL CLOSE (sil ,kcl)
 CALL CLOSE (uset,kcl)
 RETURN
END SUBROUTINE mxcid
