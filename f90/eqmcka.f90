SUBROUTINE eqmcka (ip,bgpdt,cstm,eqexin,d,iscalr)
     
!     ROUTINE FORMS D MATRIX (ACCTUALLY D TRANSPOSE)
 
 
 INTEGER, INTENT(OUT)                     :: ip
 INTEGER, INTENT(IN)                      :: bgpdt
 INTEGER, INTENT(IN)                      :: cstm
 INTEGER, INTENT(IN)                      :: eqexin
 INTEGER, INTENT(IN OUT)                  :: d
 INTEGER, INTENT(OUT)                     :: iscalr
 INTEGER :: FILE, sysbuf,iz(1),mcb(7), NAME(2)
 REAL :: tr(3,3),ti(3,3),dd(6,6),r(3)
 DIMENSION       tt(3,3)
 
 COMMON /system/ sysbuf
 COMMON /packx / it1,it2,ii,jj,incr
 COMMON /zzzzzz/ z(1)
 
 EQUIVALENCE     (iz(1),z(1))
 
 DATA    iz2   , iz3,iz4,iz5  / 2,3,4,5 /
 DATA    NAME  / 4HEQMC,4HKA  /
 
!     CONVERT  EXTERNAL IP TO INTERNAL IP
 
 ibuf = korsz(z)-sysbuf+1
 FILE = eqexin
 CALL gopen (eqexin,z(ibuf),0)
 CALL READ  (*220,*10,eqexin,iz(1),ibuf-1,0,iflag)
 GO TO 240
 10 CALL CLOSE (eqexin,1)
 DO  i = 1,iflag,2
   IF (iz(i) == ip) GO TO 40
 END DO
 CALL mesage (41,ip,NAME)
 ip = 0
 GO TO 50
 30 CALL mesage (41,ip,NAME)
 
!     SCALAR POINT
 
 GO TO 60
 40 ip = iz(i+1)
 
!     FIND RZERO FOR  IP
 
 50 FILE = bgpdt
 r(1) = 0.0
 r(2) = 0.0
 r(3) = 0.0
 CALL gopen (bgpdt,z(ibuf),0)
 IF (ip == 0) GO TO 70
 i= (ip-1)*4
 CALL fread (bgpdt,z,-i,0)
 CALL fread (bgpdt,i, 1,0)
 IF (i == -1) GO TO 30
 CALL fread (bgpdt,r,3,0)
 60 CALL REWIND (bgpdt)
 CALL skprec (bgpdt,1)
 
!     SET UP TO WRITE D
 
 70 ibuf1 = ibuf-sysbuf
 nz = ibuf1-5
 
!     BRING IN CSTM
 
 FILE = cstm
 CALL OPEN   (*90,cstm,z(ibuf1),0)
 CALL fwdrec (*220,cstm)
 CALL READ   (*220,*80,cstm,z(iz5),nz,0,ncstm)
 GO TO 240
 80 CALL CLOSE  (cstm,1)
 CALL pretrs (z(iz5),ncstm)
 90 CALL gopen  (d,z(ibuf1),1)
 CALL makmcb (mcb,d,6,2,1)
 iscalr = 0
 ii     = 1
 jj     = 6
 it1    = 1
 it2    = 1
 incr   = 1
 
!     EXAMINE BGPDT
 
 100 CALL READ (*220,*190,bgpdt,z(1),4,0,iflag)
 IF (iz(1) == -1) GO TO 170
 
!     COMPUTE  TR
 
 iscalr  = 1
 tr(1,1) = 0.0
 tr(2,2) = 0.0
 tr(3,3) = 0.0
 tr(2,1) = z(iz4) -r(3)
 tr(1,2) =-tr(2,1)
 tr(3,1) = r(2)- z(iz3)
 tr(1,3) =-tr(3,1)
 tr(3,2) = z(iz2)-r(1)
 tr(2,3) =-tr(3,2)
 DO  i = 1,3
   DO  j = 1,3
     ti(i,j) = 0.0
     IF (i == j) ti(i,j) = 1.0
   END DO
 END DO
 IF (iz(1) == 0) GO TO 130
 CALL transs (iz(1),ti)
 CALL gmmats (ti,3,3,1,tr,3,3,0,tt)
 DO  i = 1,3
   DO  j = 1,3
     tr(i,j) = tt(i,j)
   END DO
 END DO
 
!     MOVE STUFF INTO  DD
 
 130 DO  i = 1,6
   DO  j = 1,3
     IF (i > 3) GO TO 140
     dd(i  ,j  ) = ti(j,i)
     dd(i+3,j+3) = dd(i,j)
     CYCLE
     140 dd(i,j) = tr(i-3,j)
     dd(j,i) = 0.0
   END DO
 END DO
 DO  i = 1,6
   CALL pack (dd(1,i),d,mcb)
 END DO
 GO TO 100
 
!     SCALAR POINT
 
 170 DO  i = 1,6
   dd(i,1) = 0.0
 END DO
 CALL pack (dd,d,mcb)
 GO TO 100
 
!     END BGPDT
 
 190 CALL CLOSE  (bgpdt,1)
 CALL CLOSE  (d,1)
 CALL wrttrl (mcb)
 RETURN
 
!     ERROR MESAGES
 
 210 CALL mesage (ip1,FILE,NAME)
 220 ip1 = -2
 GO TO 210
 240 ip1 = -8
 GO TO 210
END SUBROUTINE eqmcka
