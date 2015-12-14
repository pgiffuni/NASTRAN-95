SUBROUTINE mpy3a (z,iz,dz)
!*****
!    PREPARES B AND A(T).
!*****
 
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(OUT)                     :: iz(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dz(1)
 DOUBLE PRECISION :: da
 
 
 
 INTEGER :: filea,fileb,scr1,scr2,FILE
 INTEGER :: buf1,buf2,buf3
 INTEGER :: prec,precn
 INTEGER :: zpntrs
 INTEGER :: typin,typout,row1,rowm
 INTEGER :: utyp,urow1,urown,uincr
 INTEGER :: eol,eor
 INTEGER :: precl
 
 
 
 
 DIMENSION NAME(2)
 DIMENSION mcb(7)
 
 
! FILES
 COMMON / mpy3tl / filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr,lkore,code,prec,lcore,scr3(7),buf1,buf2, buf3,buf4,e
! SUBROUTINE CALL PARAMETERS
 COMMON / mpy3cp / itrl,icore,n,ncb,m,dumcp(3),zpntrs(22),laend
! PACK
 COMMON / packx  / typin,typout,row1,rowm,incr
! UNPACK
 COMMON / unpakx / utyp,urow1,urown,uincr
! TERMWISE MATRIX READ
 COMMON / zntpkx / a(2),dum(2),irow,eol,eor
 
 
 
 EQUIVALENCE     (ipoint,zpntrs(3)),      (npoint,zpntrs(4)),  &
     (iacols,zpntrs(5)),      (itrans,zpntrs(7)), (ibcols,zpntrs(11))
 EQUIVALENCE     (a(1),da)
 
 
 
 DATA NAME / 4HMPY3,4HA    /
!*****
!    FILE OPENING.
!*****
 FILE = scr1
 CALL OPEN(*901,scr1,z(buf2),1)
 FILE = fileb(1)
 CALL OPEN(*901,fileb,z(buf3),0)
 CALL fwdrec(*902,fileb)
!*****
!    UNPACK B AND PACK INTO SCRATCH FILE 1.
!*****
! PACK PARAMETERS
 typin = prec
 typout = prec
 row1 = 1
 rowm = n
 incr = 1
! UNPACK PARAMETERS
 utyp = prec
 urow1 = 1
 urown = n
 uincr = 1
 precn = prec*n
 mcb(1) = 301
 mcb(2) = 0
 mcb(3) = n
 mcb(4) = 1
 mcb(5) = prec
 mcb(6) = 0
 mcb(7) = 0
 DO  k=1,ncb
   CALL unpack(*20,fileb,z(ibcols))
   GO TO 40
   20 ib = ibcols - 1
   DO  l=1,precn
     ib = ib + 1
     z(ib) = 0.
   END DO
   40 CALL pack (z(ibcols),scr1,mcb)
   CALL savpos (scr1,iz(k))
 END DO
 CALL CLOSE (scr1,1)
 CALL CLOSE (fileb,1)
 IF (icore == 1) GO TO 9999
!*****
!    INITIALIZE ARRAY CONTAINING POINTERS TO ROWS OF MATRIX A TO 0.
!*****
 DO  l=ipoint,npoint
   iz(l) = 0
 END DO
!*****
!    COUNT NO. OF NON-ZERO COLUMNS IN EACH ROW OF A.
!*****
 FILE = filea(1)
 CALL OPEN(*901,filea,z(buf1),0)
 CALL fwdrec(*902,filea)
 DO  i=1,m
   CALL intpk(*120,filea,0,prec,0)
   110 CALL zntpki
   ii = ipoint + irow - 1
   iz(ii) = iz(ii) + 1
   IF (eol == 1) CYCLE
   GO TO 110
 END DO
!*****
!    CALCULATE POINTERS TO ROWS OF MATRIX A.
!*****
 jj = 1
 DO  l=ipoint,npoint
   IF (iz(l) == 0) CYCLE
   incrjj = iz(l)
   iz(l) = jj
   jj = jj + incrjj
 END DO
 laend = jj - 1
!*****
!    PROCESS A(T) MATRIX.
!*****
 FILE = filea(1)
 CALL REWIND (filea)
 CALL fwdrec(*902,filea)
 jj2 = iacols + laend - 1
 DO  jj=iacols,jj2
   iz(jj) = 0
 END DO
 DO  j=1,m
   CALL intpk(*250,filea,0,prec,0)
   210 CALL zntpki
   l = ipoint + irow - 1
   jj = iz(l)
   jjc = iacols + jj - 1
   220 IF (iz(jjc) == 0) GO TO 230
   jj = jj + 1
   jjc = jjc + 1
   GO TO 220
   230 iz(jjc) = j
   IF (prec == 2) GO TO 240
   jjt = itrans + jj - 1
   z(jjt) = a(1)
   IF (eol == 1) CYCLE
   GO TO 210
   240 jjt = (itrans - 1)/2 + jj
   dz(jjt) = da
   IF (eol == 1) CYCLE
   GO TO 210
 END DO
 precl = prec*laend
 CALL CLOSE (filea,1)
 GO TO 9999
!*****
!    ERROR MESSAGES.
!*****
 901 nerr = -1
 GO TO 1001
 902 nerr = -2
 1001 CALL mesage (nerr,FILE,NAME)
 
 9999 RETURN
END SUBROUTINE mpy3a
