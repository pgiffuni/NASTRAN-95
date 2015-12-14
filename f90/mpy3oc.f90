SUBROUTINE mpy3oc (z,iz,dz)
     
!     OUT-OF-CORE PRODUCT.
 
 
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(OUT)                     :: iz(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: dz(1)
 LOGICAL :: first1,first2,first3,e
 INTEGER :: filea,filee,filec,code,prec,scr1,scr2,scr3,FILE,  &
     buf1,buf2,buf3,buf4,sysbuf,zpntrs,eol,eor,precm,  &
     typin,typout,row1,rowm,utyp,urow1,urown,uincr, buf5,signab,signc
 DOUBLE PRECISION :: da
 DIMENSION  NAME(2),nams(2)
 
!     MPYAD COMMON
 COMMON /mpyadx/  mfilea(7),mfileb(7),mfilee(7),mfilec(7),mcore,  &
     mt,signab,signc,mprec,mscr
 
!     FILES
 COMMON /mpy3tl/  filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr,lkore,code,prec,lcore,scr3(7),buf1,buf2, buf3,buf4,e
 
!     SUBROUTINE CALL PARAMETERS
 COMMON /mpy3cp/  dum1(2),n,ncb,m,nk,d,maxa,zpntrs(22),laend,  &
     first1,first2,k,k2,kcount,iflag,ka,ltbc,j,ltac
 
!     PACK
 COMMON /packx /  typin,typout,row1,rowm,incr
 
!     UNPACK
 COMMON /unpakx/  utyp,urow1,urown,uincr
 
!     TERMWISE MATRIX READ
 COMMON /zntpkx/  a(2),dum(2),irow,eol,eor
 
!     SYSTEM PARAMETERS
 COMMON /system/  sysbuf,nout
 EQUIVALENCE      (isavp,zpntrs(1)),  (nsavp,zpntrs(2)),  &
     (intbu,zpntrs(3)),  (nntbu,zpntrs(4)),  &
     (ilast,zpntrs(5)),  (nlast,zpntrs(6)),  &
     (intbu2,zpntrs(7)), (nntbu2,zpntrs(8)),  &
     (ic,zpntrs(9)),     (nc,zpntrs(10)),  &
     (ibcols,zpntrs(11)),(nbcols,zpntrs(12)),  &
     (ibcid,zpntrs(13)), (nbcid,zpntrs(14)),  &
     (ibntu,zpntrs(15)), (nbntu,zpntrs(16)),  &
     (iktbp,zpntrs(17)), (nktbp,zpntrs(18)),  &
     (iantu,zpntrs(19)), (nantu,zpntrs(20)),  &
     (iakj,zpntrs(21)),  (nakj,zpntrs(22)), (a(1),da)
 DATA    NAME  /  4HMPY3,4HOC   /
 DATA    nams  /  4HSCR3,4H     /
 
!     RECALCULATION OF NUMBER OF COLUMNS OF B ABLE TO BE PUT IN CORE.
 
 buf5  = buf4 - sysbuf
 lcore = buf5 - 1
 nk = (lcore - 4*n - prec*m - (2 + prec)*maxa)/(2 + prec*n)
 IF (nk < 1) GO TO 5008
 
!    INITIALIZATION.
 
 first1 = .true.
 first2 = .true.
 first3 = .false.
 precm  = prec*m
 
!     OPEN CORE POINTERS
 
 isavp  = 1
 nsavp  = ncb
 intbu  = nsavp + 1
 nntbu  = nsavp + ncb
 ilast  = nntbu + 1
 nlast  = nntbu + ncb
 intbu2 = nlast + 1
 nntbu2 = nlast + ncb
 ic     = nntbu2 + 1
 nc     = nntbu2 + prec*m
 ibcols = nc + 1
 nbcols = nc + prec*n*nk
 ibcid  = nbcols + 1
 nbcid  = nbcols + nk
 ibntu  = nbcid + 1
 nbntu  = nbcid + nk
 iktbp  = nbntu + 1
 nktbp  = nbntu + maxa
 iantu  = nktbp + 1
 nantu  = nktbp + maxa
 iakj   = nantu + 1
 nakj   = nantu + prec*maxa
 kf     = nsavp
 kl     = nntbu
 kn2    = nlast
 kbc    = nbcols
 kbn    = nbcid
 kt     = nbntu
 kan    = nktbp
 
!     PACK PARAMETERS
 
 typin = prec
 typout= prec
 row1  = 1
 incr  = 1
 
!     UNPACK PARAMETERS
 
 utyp  = prec
 urow1 = 1
 uincr = 1
 
!     MATRIX TRAILERS
 
 CALL makmcb (scr3,scr3(1),n,2,prec)
 IF (m == n) scr3(4) = 1
 
!     PUT B ONTO SCRATCH FILE IN UNPACKED FORM.
 
 CALL mpy3a (z,z,z)
 
!     OPEN FILES AND CHECK EXISTENCE OF MATRIX E.
 
 IF (code == 0 .OR. .NOT.e) GO TO 15
 FILE = filee(1)
 CALL OPEN (*5001,filee,z(buf5),2)
 CALL fwdrec (*5002,filee)
 15 FILE = filea(1)
 CALL OPEN (*5001,filea,z(buf1),0)
 CALL fwdrec (*5002,filea)
 FILE = scr1
 CALL OPEN (*5001,scr1,z(buf2),0)
 FILE = scr2
 CALL OPEN (*5001,scr2,z(buf3),1)
 IF (code == 0) GO TO 20
 FILE = filec(1)
 CALL gopen (filec,z(buf4),1)
 rowm = filec(3)
 GO TO 30
 20 FILE = scr3(1)
 CALL OPEN (*5001,scr3,z(buf4),1)
 CALL WRITE (scr3,nams,2,1)
 rowm = scr3(3)
 
!     PROCESS SCR2 AND SET FIRST-TIME-USED AND LAST-TIME-USED FOR EACH
!     ROW OF A.
 
 30 DO  k = 1,ncb
   iz(kf+k) = 0
   iz(kl+k) = 0
 END DO
 DO  j = 1,m
   k = 0
   CALL intpk (*80,filea,0,prec,0)
   50 CALL zntpki
   k = k + 1
   iz(kt+k) = irow
   IF (iz(kf+irow) > 0) GO TO 60
   iz(kf+irow) = j
   60 iz(kl+irow) = j
   IF (eol == 1) GO TO 70
   GO TO 50
   70 CALL WRITE (scr2,iz(iktbp),k,0)
   80 CALL WRITE (scr2,0,0,1)
 END DO
 CALL CLOSE (filea,1)
 CALL OPEN (*5001,filea,z(buf1),2)
 CALL fwdrec (*5002,filea)
 CALL CLOSE (scr2,1)
 CALL OPEN (*5001,scr2,z(buf3),0)
 
!     PROCESS COLUMNS OF A ONE AT A TIME.
 
 DO  j = 1,m
   
!     INITIALIZE SUM - ACCUMULATION MATRIX TO 0.
   
   DO  i = ic,nc
     z(i) = 0.
   END DO
   IF (code == 0 .OR. .NOT.e) GO TO 105
   urown = n
   CALL unpack (*105,filee,z(ic))
   
!     PROCESS A AND PERFORM FIRST PART OF PRODUCT BA(J).
   
   105 CALL mpy3b (z,z,z)
   
!     TEST IF PROCESSING IS COMPLETE
   
   IF (iflag == 0) GO TO 340
   
!     PROCESS REMAINING TERMS OF COLUMN J OF A.
   
!     TEST IF BCOLS IS FULL
   
   110 IF (k2 < nk) GO TO 330
   
!     CALCULATE NEW NEXT TIME USED VALUES
   
   IF (first3) GO TO 130
   first2 = .false.
   first3 = .true.
   DO  jj = 1,j
     CALL fwdrec (*5002,scr2)
   END DO
   130 FILE = scr2
   kc = 0
   kn = kf
   DO  ka = 1,ncb
     kn = kn + 1
     IF (j >= iz(kn)) GO TO 140
     kc = kc + 1
     IF (j+1 < iz(kn   )) GO TO 135
     IF (j+1 < iz(kl+ka)) GO TO 160
     iz(kn2+ka) = 99999999
     GO TO 136
     135 iz(kn2+ka) = iz(kn)
     136 kc = kc + 1
     CYCLE
     140 IF (j < iz(kl+ka)) GO TO 150
     iz(kn) = 99999999
     iz(kn2+ka) = iz(kn)
     kc = kc + 2
     CYCLE
     150 iz(kn    ) = 0
     160 iz(kn2+ka) = 0
   END DO
   IF (kc == 2*ncb) GO TO 240
   jj = j + 1
   180 CALL READ (*5002,*210,scr2,ka,1,0,kk)
   IF (iz(kn2+ka) > 0) GO TO 180
   IF (jj == j+1) GO TO 190
   iz(kn2+ka) = jj
   kc = kc + 1
   190 IF (iz(kf+ka) > 0) GO TO 200
   iz(kf+ka) = jj
   kc = kc + 1
   200 IF (kc == 2*ncb) GO TO 220
   GO TO 180
   210 jj = jj + 1
   GO TO 180
   220 mm = m - 1
   IF (j == mm) GO TO 290
   
!     POSITION SCRATCH FILE FOR NEXT PASS THROUGH
   
   jj  = jj - j
   j2  = j  + 2
   jj1 = jj - 1
   IF (j2 < jj1) GO TO 250
   IF (jj1 >  0) GO TO 270
   230 CALL fwdrec (*5002,scr2)
   GO TO 290
   240 IF (j == m) GO TO 290
   GO TO 230
   250 CALL REWIND (scr2)
   j1 = j + 1
   DO  jfwd = 1,j1
     CALL fwdrec (*5002,scr2)
   END DO
   GO TO 290
   270 DO  jbck = 1,jj1
     CALL bckrec (scr2)
   END DO
   
!     ASSIGN NEXT TIME USED TO COLUMNS OF B IN CORE
   
   290 DO  kk = 1,nk
     i = iz(kbc+kk)
     iz(kbn+kk) = iz(kf+i)
   END DO
   
!     ASSIGN NEXT TIME USED TO NON-ZERO TERMS IN COLUMN OF A
   
   DO  kk = 1,k
     IF (iz(kt+kk) == 0) GO TO 310
     i = iz(kt+kk)
     iz(kan+kk) = iz(kf+i)
     CYCLE
     310 iz(kan+kk) = 0
   END DO
   
!     PERFORM MULTIPLICATION AND SUMMATION FOR NEXT TERM OF COLUMN OF A
   
   330 CALL mpy3c (z,z,z)
   
!     TEST IF PROCESSING OF BA(J) IS COMPLETE
   
   IF (kcount == k) GO TO 340
   IF (first2) GO TO 110
   iz(kbn+ltbc) = iz(kn2+ltac)
   GO TO 330
   
!     PACK COLUMN OF C OR BA.
   
   340 IF (code == 0) GO TO 350
   CALL pack (z(ic),filec,filec)
   CYCLE
   350 CALL pack (z(ic),scr3,scr3)
 END DO
 
!     CLOSE FILES.
 
 CALL CLOSE (filea,2)
 CALL CLOSE (scr1,1)
 CALL CLOSE (scr2,1)
 IF (.NOT.e) GO TO 369
 CALL CLOSE (filee,2)
 369 IF (code == 0) GO TO 370
 CALL CLOSE (filec,1)
 GO TO 9999
 370 CALL CLOSE (scr3,1)
 CALL wrttrl (scr3)
 
!     CALL MPYAD TO FINISH PRODUCT
 
 DO  i = 1,7
   mfilea(i) = filea(i)
   mfileb(i) = scr3(i)
   mfilee(i) = filee(i)
   mfilec(i) = filec(i)
 END DO
 mt     = 1
 signab = 1
 signc  = 1
 mprec  = prec
 mscr   = scr1
 CALL mpyad (z,z,z)
 GO TO 9999
 
!     ERROR MESSAGES.
 
 5001 nerr = -1
 GO TO 6000
 5002 nerr = -2
 GO TO 6000
 5008 nerr = -8
 FILE = 0
 6000 CALL mesage (nerr,FILE,NAME)
 
 9999 RETURN
END SUBROUTINE mpy3oc
