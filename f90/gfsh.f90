SUBROUTINE gfsh(nuy,h)
     
!     ROUTINE TO CALCULTE THE H TRANSFORMATION MATRIX USED WHEN NO
!     SPC'S ARE ON THE FLUID
 
 
 INTEGER, INTENT(IN)                      :: nuy
 INTEGER, INTENT(IN OUT)                  :: h
 REAL :: rz(2)
 
 INTEGER :: z        ,sysbuf   ,mcb(7)   , ti1 ,to1      ,NAME(2)
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
!     SYSTEM COMMON
 
 COMMON / system /       sysbuf
 
!     PACK - UNPACK COMMON BLOCKS
 
 COMMON / packx /        ti1      ,to1      ,i1       ,n1 ,incr1
 
 EQUIVALENCE   ( z(1) , rz(1) )
 
 DATA NAME / 4HGFSH , 4H     /
 
!     ALLOCATE CORE
 
 nz = korsz(z(1))
 ibuf = nz - sysbuf
 nz = ibuf - 1
 IF(nz < nuy) GO TO 1008
 nuy1 = nuy - 1
 CALL makmcb(mcb,h,nuy1,2,2)
 ti1 = 1
 to1 = 2
 i1 = 1
 n1 = nuy1
 incr1 = 1
 
 DO  i=1,nuy
   rz(i) = -1.0 / FLOAT(nuy)
 END DO
 CALL gopen(h,z(ibuf),1)
 DO  i=1,nuy
   rz(i) = FLOAT(nuy1) / FLOAT(nuy)
   CALL pack(rz(2),h,mcb)
   rz(i) = -1.0 / FLOAT(nuy)
 END DO
 CALL CLOSE(h,1)
 CALL wrttrl(mcb)
 RETURN
 
!     ERRORS
 
 1008 CALL mesage(-8,0,NAME)
 RETURN
END SUBROUTINE gfsh
