SUBROUTINE hess2(nrow,iden,ipv)
     
!     HESS2  WILL GENERATE  AN IDENTITY MATRIX AND A PARTIIONING VECTOR
 
 
 INTEGER, INTENT(IN)                      :: nrow
 INTEGER, INTENT(IN OUT)                  :: iden
 INTEGER, INTENT(IN OUT)                  :: ipv
 INTEGER :: mcb(7) , iz(1)
 INTEGER :: sysbuf
 
 COMMON /packx/it1,it2,ii,jj,incr
 COMMON /system/ ksystm(65)
 COMMON /zzzzzz/z(1)
 
 EQUIVALENCE ( ksystm( 1) , sysbuf )
 EQUIVALENCE ( z(1),iz(1) )
 
! ----------------------------------------------------------------------
 
 CALL makmcb( mcb, iden, nrow, 8, 1 )
 nz = korsz(z)
 ibuf1 = nz- sysbuf
 CALL gopen(iden,iz(ibuf1),1)
 it1=1
 it2=1
 incr=1
 z(1)=-1.0
 DO  i=1,nrow
   ii = i
   jj=i
   CALL pack(z,iden,mcb)
 END DO
 CALL CLOSE(iden,1)
 CALL wrttrl(mcb)
 
!     BUILD PARTITIONING VECTOR
 
 CALL makmcb( mcb, ipv, 2*nrow, 2, 1 )
 CALL gopen(ipv,iz(ibuf1),1)
 DO  i=1,nrow
   z(i)=1.0
 END DO
 ii = nrow+1
 jj= 2*nrow
 CALL pack(z,ipv,mcb)
 CALL wrttrl(mcb)
 CALL CLOSE(ipv,1)
 RETURN
END SUBROUTINE hess2
