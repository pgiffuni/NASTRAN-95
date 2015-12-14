SUBROUTINE fa1pkv (az,amk,amb,n,e1,cz,bref,pi,vel,ibuf)
     
 
 COMPLEX, INTENT(OUT)                     :: az(1)
 REAL, INTENT(IN)                         :: amk(1)
 REAL, INTENT(IN)                         :: amb(1)
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN)                         :: e1(5)
 REAL, INTENT(IN OUT)                     :: cz(1)
 REAL, INTENT(IN OUT)                     :: bref
 REAL, INTENT(IN OUT)                     :: pi
 REAL, INTENT(IN OUT)                     :: vel
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 INTEGER :: iv(6),trl(7)
 REAL :: v(6),e(2)
 COMPLEX :: ceig,eigen,eigz
 COMMON /system/ sysbuf,nout,SPACE(6),nlpp,x(2),lines
 EQUIVALENCE     (v(1),iv(1)),(eigen,e(1))
 DATA    iscr  / 301/, ipass /0/
 
 eigz = (0.0,0.0)
 IF (n < 2) GO TO 1000
 e(1) = e1(1)
 e(2) = e1(2)
 IF (ipass /= 0) GO TO 5
 CALL OPEN (*1000,iscr,ibuf,1)
 GO TO 9
 5 CALL OPEN (*1000,iscr,ibuf,3)
 9 ipass = ipass + 1
 
!     BUILD A = IP2 + M-1B P + M-1K
 
 ceig = eigen*eigen
 k = 0
 DO  i = 1,n
   DO  j = 1,n
     k = k + 1
     az(k) = -amb(k)*eigen - amk(k)
     IF (i == j) az(k) = az(k) + ceig
   END DO
 END DO
 
!     CORE FOR EGNVCT
 
 n2 = n*2
 na = 1  + n2*n
 nb = na + n2
 nc = nb + n2
 nd = nc + n2
 CALL egnvct (az,cz(na),eigz,cz(nb),cz(nc),cz(nd),n)
 
!     BUILD ON SCR1 DATA FOR VECTOR OUTPUT
 
 iv(1) = ipass
 iv(2) = ipass
 v (3) = e1(1)
 v (4) = e1(2)
 IF (e1(2) == 0.0) GO TO 20
 v(5)  = e1(3)
 v(6)  = e1(5)
 GO TO 22
 20 v(5)  = 0.0
 v(6)  = (bref/(.34657*vel))*e1(1)
 22 CALL WRITE (iscr,iv,6,1)
 CALL WRITE (iscr,cz(nb),n2,1)
 
!     VECTOR IS IN CZ(NB)
 
 lines = nlpp
 k = 0
 DO  i = 1,n
   IF (lines < nlpp) GO TO 25
   CALL page1
   WRITE  (nout,21) eigen
   21 FORMAT (1H0,47X,30HEIGENVECTOR from the pk method, /3X,  &
       13HEIGENVALUE = ,1P,e15.5,1P,e15.5, //3X,11HEIGENVECTOR)
   lines = lines + 5
   25 lines = lines + 1
   WRITE  (nout,26) cz(nb+k),cz(nb+k+1)
   26 FORMAT (16X,1P,e15.5,1P,e15.5)
   k = k + 2
 END DO
 trl(1) = iscr
 trl(2) = 1
 CALL wrttrl (trl)
 1000 CALL CLOSE (iscr,3)
 RETURN
END SUBROUTINE fa1pkv
