SUBROUTINE gfsspc(nuy,pvec)
     
!     ROUTINE TO CALCULATE A PARTITIONING VECTOR TO REMOVE FIRST
!     ROW AND COLUMN OF FLUID STIFFNESS MATRIX IF NO SPC'S ARE ON
!     THE FLUID
 
 
 INTEGER, INTENT(IN)                      :: nuy
 INTEGER, INTENT(IN OUT)                  :: pvec
 INTEGER :: mcb(7)   , z        ,sysbuf   ,NAME(2)
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
!     SYSTEM COMMON
 
 COMMON / system /       sysbuf
 
!     PACK - UNPACK COMMON BLOCKS
 
 COMMON / zblpkx /       a(4)     ,irow
 
 DATA NAME / 4HGFSS , 4HPC   /
 
!     ALLOCATE CORE
 
 nz = korsz(z(1))
 ibuf = nz - sysbuf
 nz = ibuf - 1
 IF(nz < 0) GO TO 1008
 
 nuy1 = nuy - 1
 CALL makmcb(mcb,pvec,nuy,2,1)
 CALL gopen(pvec,z(ibuf),1)
 CALL bldpk(1,1,pvec,0,0)
 a(1) = 1.0
 irow = 1
 CALL zblpki
 CALL bldpkn(pvec,0,mcb)
 CALL CLOSE(pvec,1)
 CALL wrttrl(mcb)
 RETURN
 
!     ERRORS
 
 1008 CALL mesage(-8,0,NAME)
 RETURN
END SUBROUTINE gfsspc
