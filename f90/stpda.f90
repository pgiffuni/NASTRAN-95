SUBROUTINE stpda (INPUT,ajjl,skj)
     
!     DRIVER FOR STRIP THEORY
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: ajjl
 INTEGER, INTENT(IN OUT)                  :: skj
 INTEGER :: sysbuf,iz(8), NAME(2),claf,lclaf,lcirc
 COMPLEX :: ekm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /stripc/ ns,bref,clam,fm,ncirc,nncirc,ekr(1),  &
     dum,bb(4),beta(4),ekm(4,4)
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,nout
 COMMON /condas/ pi,twopi
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk,tskj(7),isk,nsk
 COMMON /BLANK / nk,nj
 COMMON /packx / iti,it0,ii,nn,incr
 EQUIVALENCE     (iz(1),z(1))
 DATA     NAME / 4HSTPD,4HA   /
 
 icore = korsz(iz) - 4*sysbuf
 
!     BRING IN DATA AND ALLOCATE CORE
 
 CALL fread (INPUT,z,8,0)
 nnj   = iz(1)
 claf  = iz(2)
 lclaf = iz(3)
 ncirc = iz(4)
 lcirc = iz(5)
 nncirc= ncirc + 1
 nmach = iz(6)
 ns    = iz(7)
 i8    = 8
 clam  = z(i8)
 fm    = 1.0
 bref  = refc / 2.0
 ekr(1)= rfk
 idy   = 1
 ibloc = idy  + ns
 id    = ibloc+ ns
 ica   = id   + ns
 igap  = ica  + ns
 insize= igap + ns
 icla  = insize + ns
 ibm   = icla + ns
 igm   = ibm  + 16 * ns
 ipm   = igm  + 12 * ns
 ioc   = ipm  + 37 * ns
 IF (ioc > icore) CALL mesage (-8,0,NAME)
 
!     READ IN ARRAYS WHICH ARE FIXED
 
 nw = 6*ns
 CALL fread (INPUT,z,nw,0)
 
!     SET CLA ARRAY OR BB AND BETA
 
 IF (claf == 0) GO TO 40
 IF (claf < 0) GO TO 30
 
!     FIND MACH NUMBER FOR CLA
 
 DO  i = 1,nmach
   CALL fread (INPUT,rm,1,0)
   IF (rm == fmach) GO TO 20
   CALL fread (INPUT,z,-ns,0)
 END DO
 GO TO 999
 
!     MACH NUMBER NOT INPUT ON AEFACT CARD CLCAF
 
 20 CALL fread (INPUT,z(icla),ns,1)
 GO TO 90
 30 CALL fread (INPUT,rm,1,0)
 CALL fread (INPUT,z(icla),ns,1)
 DO   i = 1,ns
   z(icla+i-1) = z(icla+i-1) * SQRT((1.0-(rm*rm*clam*clam)) /  &
       (1.0-(fmach*fmach*clam*clam)))
 END DO
 GO TO 90
 40 DO  i = 1,ns
   z(icla+i-1) = twopi
 END DO
 IF (ncirc == 0) GO TO 80
 DO  i = 1,nmach
   CALL fread (INPUT,rm,1,0)
   IF (rm == fmach) GO TO 70
   CALL fread (INPUT,z,-(2*ncirc+1),0)
 END DO
 GO TO 998
 70 CALL fread (INPUT,bb(1),1,0)
 DO  i = 2,nncirc
   CALL fread (INPUT,bb(i),1,0)
   CALL fread (INPUT,beta(i),1,0)
 END DO
 80 CALL fread (INPUT,z,0,1)
 
!     OUTPUT SKJ
 
 90 iti = 1
 it0 = 3
 ii  = isk
 nsk = nsk+1
 nn  = nsk
 rm  = 1.0
 DO  i = 1,nnj
   CALL pack (rm,skj,tskj)
   ii  = ii+1
   IF (i == nnj) CYCLE
   nn  = nn+1
 END DO
 isk = ii
 nsk = nn
 iti = 3
 it0 = 3
 CALL stpbg  (z(ibm),z(igm),ns,z(ibloc),z(id),z(ica),z(insize))
 CALL stpphi (z(ica),z(ibloc),z(ipm),ns)
 CALL stpaic (z(ibloc),z(idy),z(insize),z(igap),z(ibm),z(igm),  &
     z(ipm),ns,z(icla),ajjl)
 nrow = nrow + nnj
 RETURN
 
!     ERROR MESSAGES
 
 998 n = lcirc
 GO TO 1000
 999 n = lclaf
 1000 WRITE  (nout,9999) ufm,fmach,n
 9999 FORMAT (a23,' 2426, MACH NUMBER ',f10.5,' WAS NOT FOUND ON ',  &
     'AEFACT CARD',i9)
 CALL mesage (-61,0,NAME)
 RETURN
END SUBROUTINE stpda
