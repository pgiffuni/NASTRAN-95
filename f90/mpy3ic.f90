SUBROUTINE mpy3ic (z,iz,dz)
     
!     IN-CORE PRODUCT.
 
 
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(IN OUT)                  :: iz(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: dz(1)
 LOGICAL :: first1,first2,e
 INTEGER :: filea,fileb,filee,filec,code,prec,scr1,scr3,FILE,  &
     buf1,buf2,buf3,buf4,d,zpntrs,typin,typout,row1,  &
     rowm,utyp,urow1,urown,uincr,eol,eor,precm
 DOUBLE PRECISION :: dd,nn,mm,pp
 DIMENSION  NAME(2)
 COMMON /mpy3tl/ filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr,lkore,code,prec,lcore,scr3(7),buf1,buf2, buf3,buf4,e
 COMMON /mpy3cp/ itrl,icore,n,ncb,m,nk,d,maxa,zpntrs(22),laend,  &
     first1,first2,k,k2,kcount,iflag,ka,ltbc,j,i,ntbu
 COMMON /packx / typin,typout,row1,rowm,incr
 COMMON /unpakx/ utyp,urow1,urown,uincr
 COMMON /zntpkx/ a(2),dum(2),irow,eol,eor
 COMMON /system/ sysbuf,nout
 EQUIVALENCE     (isavp ,zpntrs( 1)), (nsavp ,zpntrs( 2)),  &
     (ipoint,zpntrs( 3)), (npoint,zpntrs( 4)),  &
     (iacols,zpntrs( 5)), (nacols,zpntrs( 6)),  &
     (itrans,zpntrs( 7)), (ntrans,zpntrs( 8)),  &
     (ic    ,zpntrs( 9)), (nc    ,zpntrs(10)),  &
     (ibcols,zpntrs(11)), (nbcols,zpntrs(12)),  &
     (ibcid ,zpntrs(13)), (nbcid ,zpntrs(14)),  &
     (ibntu ,zpntrs(15)), (nbntu ,zpntrs(16)),  &
     (iktbp ,zpntrs(17)), (nktbp ,zpntrs(18)),  &
     (iantu ,zpntrs(19)), (nantu ,zpntrs(20)),  &
     (iakj  ,zpntrs(21)), (nakj  ,zpntrs(22))
 DATA     NAME / 4HMPY3,4HIC   /
 
 
!     INITIALIZATION.
 
 first1 = .true.
 first2 = .true.
 dd     = d
 nn     = ncb
 mm     = m
 pp     = prec
 
!     OPEN CORE POINTERS
 
 isavp  = 1
 nsavp  = ncb
 ipoint = nsavp  + 1
 npoint = nsavp  + ncb
 iacols = npoint + 1
!     NACOLS = NPOINT + D*NCB*M/10000
 nacols = npoint + (dd*nn*mm/10000.d0 + 0.5D0)
 itrans = nacols + 1
 IF (prec /= 1 .AND. MOD(itrans,2) /= 1) itrans = itrans + 1
!     NTRANS = ITRANS + PREC*D*NCB*M/10000 - 1
 ntrans = itrans + (pp*dd*nn*mm/10000.d0 + 0.5D0) - 1
 ic = ntrans + 1
 IF (prec /= 1 .AND. MOD(ic,2) /= 1) ic = ic + 1
 nc = ic + prec*m - 1
 ibcols= nc + 1
 nbcols= nc + prec*n*nk
 ibcid = nbcols + 1
 nbcid = nbcols + nk
 ibntu = nbcid  + 1
 nbntu = nbcid  + nk
 iktbp = nbntu  + 1
 nktbp = nbntu  + maxa
 iantu = nktbp  + 1
 nantu = nktbp  + maxa
 iakj  = nantu  + 1
 nakj  = nantu  + prec*maxa
 
!     PACK PARAMETERS
 
 typin = prec
 typout= prec
 row1  = 1
 incr  = 1
 
!     UNPACK PARAMETERS
 
 utyp  = prec
 urow1 = 1
 uincr = 1
 
!     PREPARE B AND A(T).
 
 CALL mpy3a (z,z,z)
 
!     OPEN FILES AND CHECK EXISTENCE OF MATRIX E.
 
 IF (.NOT.e) GO TO 20
 FILE = filee(1)
 CALL OPEN (*5001,filee,z(buf4),2)
 CALL fwdrec (*5002,filee)
 20 FILE = filea(1)
 CALL OPEN (*5001,filea,z(buf1),2)
 CALL fwdrec (*5002,filea)
 FILE = scr1
 CALL OPEN (*5001,scr1,z(buf2),0)
 FILE = filec(1)
 CALL gopen (filec,z(buf3),1)
 rowm = filec(3)
 
!     PROCESS COLUMNS OF C ONE BY ONE.
 
 DO  j = 1,m
   
!     INITIALIZE COLUMN OF C.
   
   DO  ix = ic,nc
     z(ix) = 0.
   END DO
   IF (.NOT.e) GO TO 50
   urown = m
   CALL unpack (*50,filee,z(ic))
   50 precm = prec*m
   
!     PROCESS A AND PERFORM FIRST PART OF PRODUCT.
   
   CALL mpy3b (z,z,z)
   
!     TEST IF PROCESSING IS COMPLETE
   
   IF (iflag == 0) GO TO 900
   
!     PROCESS REMAINING TERMS OF COLUMN J OF A.
   
!     TEST IF BCOLS IS FULL
   
   100 IF (k2 < nk) GO TO 150
   
!     CALCULATE NEXT TIME USED FOR COLUMNS OF B AND/OR TERMS OF A
   
   IF (.NOT.first2) GO TO 120
   first2 = .false.
   ibc = ibcid - 1
   ib  = ibntu - 1
   DO  ii = 1,nk
     ibc= ibc + 1
     i  = iz(ibc)
     CALL mpy3nu (z)
     ib = ib + 1
     iz(ib) = ntbu
   END DO
   120 ik = iktbp - 1
   ia = iantu - 1
   DO  ii = 1,k
     ik = ik + 1
     ia = ia + 1
     IF (iz(ik) == 0) GO TO 130
     i  = iz(ik)
     CALL mpy3nu (z)
     iz(ia) = ntbu
     CYCLE
     130 iz(ia) = 0
   END DO
   
!     ADD OR REPLACE COLUMN OF B INTO CORE AND PERFORM COMPUTATION
   
   150 CALL mpy3c (z,z,z)
   IF (kcount == k) GO TO 900
   IF (first2) GO TO 100
   GO TO 150
   
!     PACK COLUMN OF C.
   
   900 CALL pack (z(ic),filec,filec)
 END DO
 
!     CLOSE FILES.
 
 CALL CLOSE (filea,2)
 CALL CLOSE (scr1,1)
 CALL CLOSE (filec,1)
 IF (e) CALL CLOSE (filee,2)
 GO TO 9999
 
!     ERROR MESSAGES.
 
 5001 nerr = -1
 GO TO 6000
 5002 nerr = -2
 6000 CALL mesage (nerr,FILE,NAME)
 
 9999 RETURN
END SUBROUTINE mpy3ic
