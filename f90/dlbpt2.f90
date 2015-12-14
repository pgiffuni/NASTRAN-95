SUBROUTINE dlbpt2 (INPUT,w1jk,w2jk)
     
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: w1jk
 INTEGER, INTENT(IN OUT)                  :: w2jk
 INTEGER :: sysbuf,ecore,tw1jk,tw2jk,NAME(2)
 DIMENSION       a(4),iz(1)
 COMMON /packx / iti,ito,ii,nn,incr
 COMMON /amgp2 / tw1jk(7),tw2jk(7)
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,rfk
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys, z(1)
 EQUIVALENCE     (iz(1),z(1))
 DATA    NAME  / 4HDLBP,4HT2  /
 
!     GET CORE THEN SET POINTERS TO ACPT TABLE ARRAYS
 
 ecore = korsz(iz) - 4*sysbuf
 
!     READ LENGTHS OF ARRAYS
 
 CALL fread (INPUT,nj1,13,0)
 
!     COMPUTE POINTERS TO OPEN CORE
 
 IF (ntp == 0) CALL fread (INPUT,0,0,1)
 IF (ntp == 0) GO TO 50
 lns  = iz(1)
 inc  = 1
 ins  = inc
 inb  = ins + np
 inas = inb + np
 izin = inas
 iyin = izin
 inbea1 = iyin   + np
 inbea2 = inbea1 + nb
 insbea = inbea2 + nb
 izb  = insbea + nb
 iyb  = izb  + nb
 iavr = iyb  + nb
 iarb = iavr + nb
 infl = iarb + nb
 ixle = infl + nb
 ixte = ixle + nb
 int121 = ixte   + nb
 int122 = int121 + nb
 izs = int122 + nb
 n = 3*np + 12*nb
 
!     READ FIXED ARRAYS
 
 IF (n > ecore) GO TO 180
 CALL fread (INPUT,iz,n,0)
 
!     GET LENGTHS OF VARIABLE ARRAYS, PANELS THEN BODIES
 
 lnas = 0
 IF (np == 0) GO TO 20
 DO  i = 1,np
   lnas = lnas + iz(inas+i-1)
 END DO
 20 lnb  = 0
 lnsb = 0
 lnfl = 0
 lt1  = 0
 lt2  = 0
 DO  i = 1,nb
   k    = i - 1
   lnb  = lnb  + iz(inbea1+k)
   lnsb = lnsb + iz(insbea+k)
   lnfl = lnfl + iz(infl+k)
   lt1  = lt1  + iz(int121+k)
   lt2  = lt2  + iz(int122+k)
 END DO
 
!     READ VARIABLE  ARRAYS AND SET POINTERS TO CORE
 
 next = n + 1
 n = 2*nb + 5*lns + 4*ntp + 3*lnb + 4*lnsb + lnas + 2*lnfl + lt1  + lt2
 IF (next+n >= ecore) GO TO 180
 CALL READ (*190,*190,INPUT,iz(next),n,1,nw)
 next = next+ n + 1
 iys  = izs + nb + lns
 ics  = iys
 iee  = ics + nb + lns
 isg  = iee + lns
 icg  = isg + lns
 ixij = icg
 ix   = ixij+ lns
 idelx= ix  + ntp + lnb
 
!     COMPUTE TERMS AND PACK
 
 nn = ii + 1
 DO  i = 1,ntp
   a(1) = 0.0
   a(2) = 1.0
   CALL pack (a,w1jk,tw1jk)
   a(1) = -(2.0/refc)
   a(2) = z(idelx+i-1)/(2.0*refc)
   CALL pack (a,w2jk,tw2jk)
   
!     BUMP PACK INDEXES
   
   ii = ii + 2
   IF (i == ntp) CYCLE
   nn = nn + 2
 END DO
 50 ntzy = ntz + nty
 IF (ntzy == 0) GO TO 70
 nn   = ii + 1
 a(1) = 0.0
 a(2) = 0.0
 DO  i = 1,ntzy
   CALL pack (a,w1jk,tw1jk)
   CALL pack (a,w2jk,tw2jk)
 END DO
 70 ntzy = ntzs + ntys
 IF (ntzy == 0) GO TO 200
 
!     ANOTHER HARDER SHUFFLE
 
 iii = ii
 inbea2 = inbea2 - 1
 insbea = insbea - 1
 ify = ii
 IF (nbz == 0) GO TO 120
 DO  i = 1,nbz
   ibt = iz(inbea2+i)
   nbe = iz(insbea+i)
   IF (ibt == 2) GO TO 90
   a(1) = 0.0
   a(2) = 1.0
   a(3) = -2.0/refc
   a(4) = 0.0
   DO  j = 1,nbe
     nn = ii + 1
     CALL pack (a,w1jk,tw1jk)
     CALL pack (a(3),w2jk,tw2jk)
     ii  = ii + 2
     ify = ii
   END DO
   CYCLE
   90 a(1) = 0.0
   a(4) = 0.0
   DO  j = 1,nbe
     nn   = ii + 3
     a(2) = 0.0
     a(3) = 1.0
     CALL pack (a,w1jk,tw1jk)
     a(2) = -2.0/refc
     a(3) = 0.0
     CALL pack (a,w2jk,tw2jk)
     ii = ii + 4
   END DO
 END DO
 120 IF (nby == 0) GO TO 170
 ii = ify
 nbtd = nb - nby + 1
 DO  i = nbtd,nb
   ibt = iz(inbea2+i)
   nbe = iz(insbea+i)
   IF (ibt == 3) GO TO 140
   a(2) = 0.0
   a(3) = 0.0
   DO  j = 1,nbe
     nn   = ii + 3
     a(1) = 0.0
     a(4) =-1.0
     CALL pack (a,w1jk,tw1jk)
     a(1) = -2.0/refc
     a(4) = 0.0
     CALL pack (a,w2jk,tw2jk)
     ii  = ii + 4
   END DO
   CYCLE
   140 a(1) = 0.0
   a(2) =-1.0
   a(3) =-2.0/refc
   a(4) = 0.0
   DO  j = 1,nbe
     nn = ii + 1
     CALL pack (a,w1jk,tw1jk)
     CALL pack (a(3),w2jk,tw2jk)
     ii = ii + 2
   END DO
 END DO
 170 ii = iii + ntzy*2
 nn = ii  - 1
 GO TO 200
 
!     ERROR MESSAGES
 
 180 CALL mesage (-8,0,NAME)
 190 CALL mesage (-7,0,NAME)
 200 RETURN
END SUBROUTINE dlbpt2
