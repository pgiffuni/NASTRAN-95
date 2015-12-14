SUBROUTINE fvrs1a (base,base1,z,w,buf,INDEX,modfrl,basexg,nrow,  &
        nf,nfx,fkmax,omega)
     
 
 COMPLEX, INTENT(IN)                      :: base(3,nfx)
 COMPLEX, INTENT(IN OUT)                  :: base1(3,nfx)
 COMPLEX, INTENT(OUT)                     :: z(nrow)
 REAL, INTENT(IN OUT)                     :: w(nf)
 REAL, INTENT(IN OUT)                     :: buf(1)
 INTEGER, INTENT(IN OUT)                  :: INDEX(1)
 LOGICAL, INTENT(IN OUT)                  :: modfrl
 INTEGER, INTENT(IN)                      :: basexg
 INTEGER, INTENT(IN)                      :: nrow
 INTEGER, INTENT(IN OUT)                  :: nf
 INTEGER, INTENT(IN)                      :: nfx
 INTEGER, INTENT(IN)                      :: fkmax
 REAL, INTENT(IN OUT)                     :: omega
 
 
 
 
 
 
 DIMENSION mcb(7)
 
 COMMON /packx/ in,iout,ns,nl,incr
 
!-----------------------------------------------------------------------
!     COMPUTE NUMBER OF GRID POINTS (SCALAR POINTS ARE NOT ALLOWED).
!-----------------------------------------------------------------------
 npts=nrow/6
!-----------------------------------------------------------------------
!     GENERATE BASE TABLE
!----------------------------------------------------------------------
 IF (modfrl) GO TO 100
 CALL fvrs1b(base,w,nf)
 GO TO 135
 100  CALL fvrs1c(base,w,omega,nf)
 135  CONTINUE
!---------------------------------------------------------------------
!     SORT BASE BY INDEX TO MAKE IT COMPATIBLE TO FRLX
 IF (.NOT. modfrl) GO TO 137
 CALL fvrs1d(base,base1,INDEX,nfx)
 137  CONTINUE
!---------------------------------------------------------------------
!     PREPARE TO OUTPUT BASEXG
!----------------------------------------------------------------------
 CALL gopen(basexg,buf,1)
!-------------------------------
!     DEFINE MCB
 mcb(1)=basexg
 mcb(2)=0
 mcb(3)=nrow
 mcb(4)=2
 mcb(5)=3
 mcb(6)=0
 mcb(7)=0
!-------------------------------
!     DEFINE PACKING CONSTANTS
 in=3
 iout=3
 ns=1
 nl=nrow
 incr=1
!-----------------------------------------------------------------------
!     GENERATE AND PACK 1ST NF COLUMNS OF BASEXG
!     BASEXG-1
!     ZERO OUT COLUMN
 DO  i=1,nrow
   z(i)=(0.0,0.0)
 END DO
 DO  i=1,nfx
   l=1
   DO  k=1,npts
     z(l)=base(1,i)
     l=l+6
   END DO
   CALL pack(z,basexg,mcb)
 END DO
 IF(fkmax < 2)GO TO 500
!----------------------------------------------------------------------
!     GENERATE AND PACK 2ND NF COLUMNS OF BASEXG
!     BASEXG-2
!     ZERO COLUMN
 DO  i=1,nrow
   z(i)=(0.0,0.0)
 END DO
 DO  i=1,nfx
   l=1
   DO  k=1,npts
     z(l+1)=base(2,i)
     z(l+2)=base(3,i)
     l=l+6
   END DO
   CALL pack(z,basexg,mcb)
 END DO
 IF(fkmax < 3) GO TO 500
!----------------------------------------------------------------------
!     GENERATE AND PACK 3RD NF COLUMNS OF BASEXG
!     BASEXG-3
 DO  i=1,nfx
   l=1
   DO  k=1,npts
     z(l+1)=base(3,i)
     z(l+2)=-base(2,i)
     l=l+6
   END DO
   CALL pack(z,basexg,mcb)
 END DO
!-----------------------------------------------------------------------
!     GENERATE 4TH THRU FKMAX NF COLUMN GROUPS-(NULL)INTO BASEXG
 IF(fkmax < 4)GO TO 500
 ns=1
 nl=1
 z(1)=(0.0,0.0)
 DO  i=4,fkmax
   DO  k=1,nfx
     CALL pack(z,basexg,mcb)
   END DO
 END DO
!----------------------------------------------------------------------
!     CLOSE OUTPUT DATA BLOCK
 500  CALL CLOSE(basexg,1)
 CALL wrttrl(mcb)
 RETURN
END SUBROUTINE fvrs1a
