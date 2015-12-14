SUBROUTINE gfsmt(mt,mmat,mrow)
     
!     ROUTINE TO ADD 1.0 TO ROW MROW AND COLUMN MROW OF MT TO PREVENT
!     SINGULARITIES IN THE MASS MATRIX FOR GIVINS
 
 
 INTEGER, INTENT(IN)                      :: mt
 INTEGER, INTENT(IN)                      :: mmat
 INTEGER, INTENT(IN)                      :: mrow
 DOUBLE PRECISION :: val
 
 INTEGER :: z        ,sysbuf   ,mcb(7)   ,NAME(2)  ,mmat  &
     , inblk(15),outblk(15)
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
!     SYSTEM COMMON
 COMMON / system /       sysbuf
 
 
!     PACK - UNPACK COMMON BLOCKS
 
 COMMON / zblpkx /       a(4)     ,irow
 
 EQUIVALENCE   ( val , a(1) )
 
 DATA NAME / 4HGFSM , 4HT    /
 DATA inblk / 15*0 /, outblk / 15*0 /
 
 
 mcb(1) = mt
 CALL rdtrl(mcb)
 nrow = mcb(2)
 
!     ALLOCATE BUFFERS
 
 nz = korsz(z(1))
 ibuf1 = nz - sysbuf
 ibuf2 = ibuf1 - sysbuf
 nz = ibuf2 - 1
 IF(nz < 100) GO TO 1008
 
!     OPEN FILES
 
 CALL makmcb(mcb,mmat,nrow,1,2)
 inblk(1) = mt
 outblk(1) = mmat
 CALL gopen(mt,z(ibuf1),0)
 CALL gopen(mmat,z(ibuf2),1)
 
!     COPY RECORDS UP TO MROW
 
 IF(mrow == 1) GO TO 310
 mr = mrow - 1
 DO  i=1,mr
   CALL cpystr(inblk,outblk,0,0)
 END DO
 
!    PACK OUT COLUMN MROW WITH A 1.0 IN ROW MROW.  THE COLUMN IS NULL
!     IN MT SO IT IS SKIPPED
 
 310 CALL bldpk(2,2,mmat,0,0)
 irow = mrow
 val = 1.0D0
 CALL zblpki
 CALL bldpkn(mmat,0,mcb)
 
 IF(mrow >= nrow) GO TO 320
 CALL fwdrec(*1002,mt)
 
!     BLAST OUT REST OF FILE
 
 CALL cpyfil(mt,mmat,z,nz,icnt)
 
!     CLOSE FILES
 
 320 CALL CLOSE(mt,1)
 CALL CLOSE(mmat,1)
 
!     COPY TRAILER OVER.  THE DENSITY WILL BE SLIGHTLY OFF BECAUSE
!     OF THE NEW TERM BUT IT:S CLOSE
 
 mcb(1) = mt
 CALL rdtrl(mcb)
 mcb(1) = mmat
 CALL wrttrl(mcb)
 RETURN
 
!     ERRORS
 
 1002 CALL mesage(-2,mt,NAME)
 1008 CALL mesage(-8,0,NAME)
 RETURN
END SUBROUTINE gfsmt
