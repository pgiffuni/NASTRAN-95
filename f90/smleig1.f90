SUBROUTINE smleig1(d,o,val)
     
!     COMPUTES EIGENVALUES AND VECTORS FOR 1X1 AND 2X2
 
 
 REAL, INTENT(IN OUT)                     :: d(2)
 REAL, INTENT(OUT)                        :: o(2)
 REAL, INTENT(OUT)                        :: val(2)
 REAL :: p,q
 INTEGER :: ENTRY,xentry,sysbuf,phia,mcb(7)
 DIMENSION vcom(30)
 
 COMMON/system/sysbuf
 COMMON /givn / title(150)
 COMMON /packx/it1,it2,ii,jj,incr
 COMMON /unpakx/it3,iii,jjj,incr1
 
 EQUIVALENCE  &
     (mo,title(2)),(md,title(3)),(ENTRY,title(11)),(xentry,title(20)),  &
     (vcom(1),title(101)),(n,vcom(1)),(lama,vcom(6)),(phia,vcom(12)),  &
     (nfound,vcom(10))
 
 DATA mcb/7*0/
 
!     D        ARRAY OF DIAGONALS
!     O        ARRAY OF OFF DIAGONALS
!     VAL      ARRAY OF EIGENVALUES
!     LAMA     FILE OF EIGENVALUES--HEADER,VALUES,ORDER FOUND
!     PHIA     FILE OF  VECTORS   --......,VECTORS-D.P.
!     MO       RESTART TAPE FOR MORE EIGENVALUES
!     MD       INPUT MATRIX
!     N        ORDER OF  PROBLEM
!     NFOUND   NUMBER OF EIGENVALUES/VECTOR PREVIOUSLY FOUND
 
!WKBR 2/94      IBUF1 =(KORSZ(O) - SYSBUF +1 )/2  -1
 ibuf1 =(korsz(o) - sysbuf +1 )  -1
 
!     OPEN INPUT MATRIX
 
 CALL gopen(md,o(ibuf1),0)
 
!     SETUP FOR UNPACK
 
 it3 = 1
 iii = 1
 jjj = n
 incr1= 1
 ASSIGN 101 TO itra
 CALL unpack(*1000,md,d)
 101 IF(n == 2) GO TO 110
 
!     THE MATRIX IS A 1X1
 
 o(1) = 0.0
 val(1) = d(1)
 loc = 1
 GO TO 120
 
!     THE MATRIX IS A 2X2
 
 110 o(1) = d(2)
 o(2) = 0.0
 ASSIGN 111 TO itra
 iii = 2
 CALL unpack(*1000,md,d(2))
 111 p = d(1) + d(2)
 q = SQRT( p*p -4.0*(d(1)*d(2)- o(1)**2))
 val(1) =(p+q)/2.0
 val(2) = (p-q)/2.0
 loc = 0
 
!     WRAP UP ROUTINE
 
 120 CALL CLOSE(md,1)
 
!     COPY D,O,LOC ONTO MO FOR RESTART
 
 CALL gopen(mo,o(ibuf1),1)
 
!     SETUP FOR PACK
 
 im1=1
 it1=1
 it2=1
 ii =1
 jj = n
 incr =1
 CALL pack(d,mo,mcb)
 CALL pack(o,mo,mcb)
 CALL WRITE(mo,loc,1,1)
 CALL CLOSE(mo,1)
 IF(n /= 1) GO TO 125
 
!     1X1 WRITE OUT VECTORS AND VALUES
 
 mcb(1) = phia
 mcb(2) = 0
 mcb(3) = 1
 mcb(4) = 2
 mcb(5) = 2
 mcb(6) = 0
 CALL gopen(phia,o(ibuf1),  1)
 jj = 1
 CALL pack(1.0,phia,mcb)
 CALL CLOSE(phia,1)
 CALL wrttrl(mcb(1))
 CALL gopen(lama,o(ibuf1),1 )
 IF (nfound == 0) GO TO 128
 DO  i= 1,nfound
   CALL WRITE(lama,0.0,1,0)
 END DO
 128 valx = val(1)
 CALL WRITE(lama,valx,1,1)
 IF (nfound == 0) GO TO 129
 DO  i= 1,nfound
   CALL WRITE(lama,i,1,0)
 END DO
 129 CALL WRITE(lama,nfound+1,1,1)
 CALL CLOSE(lama,1)
 mcb(1) = lama
 CALL wrttrl(mcb)
 125 xentry = -ENTRY
 RETURN
 1000 DO  i =iii,jjj
   d(i) = 0.0
 END DO
 GO TO  itra,(111,101)
END SUBROUTINE smleig1
