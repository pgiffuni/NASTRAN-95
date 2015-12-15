SUBROUTINE gfscom(awy,nuy,kc,ident,ac,scr)
     
!     ROUTINE TO COMPUTE THE FLUID COMPRESSIBILTY MATRIX
 
!     THIS MATRIX CONTAINS THE SPRING FACTOR WHICH COUPLES THE
!     STRUCTURE AND FREE SURFACE TO PREVENT VOLUME CHANGES
 
 
 REAL, INTENT(IN OUT)                     :: awy
 INTEGER, INTENT(IN)                      :: nuy
 INTEGER, INTENT(IN OUT)                  :: kc
 INTEGER, INTENT(IN OUT)                  :: ident
 INTEGER, INTENT(IN)                      :: ac
 INTEGER, INTENT(IN OUT)                  :: scr
 DOUBLE PRECISION :: dz(1)    ,dkcomp   ,val
 
 REAL :: kcomp    ,rz(1)
 
 INTEGER :: z        ,sysbuf   &
     ,mcb(7)   ,NAME(2)  ,ti1      ,to1      ,to2
 
 
!     MODULE PARAMETERS
 
 COMMON /BLANK/     nograv   ,nofree   ,kcomp   ,comptp ,FORM     ,lmodes
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
!     SYSTEM COMMON
 
 COMMON / system /       sysbuf
 
!     PACK - UNPACK COMMON BLOCKS
 
 COMMON / packx /        ti1      ,to1      ,i1       ,n1 ,incr1
 COMMON / unpakx /       to2      ,i2       ,n2       ,incr2
 COMMON / zblpkx /       a(4)     ,irow
 
 EQUIVALENCE   ( z(1) , rz(1) , dz(1) ) ,( val , a(1) )
 
 DATA NAME / 4HGFSC , 4HOM   /
 
!     ALLOCATE CORE
 
 nz = korsz(z(1))
 ibuf = nz - sysbuf
 nz = ibuf - 1
 IF(nz < nuy) GO TO 1008
 
!     FORM A COLUMN VECTOR OF ONES
 
 ti1 = 1
 to1 = 2
 i1 = 1
 n1 = nuy
 incr1 = 1
 DO  i=1,nuy
   rz(i) = 1.0
 END DO
 CALL makmcb(mcb,ident,nuy,2,2)
 CALL gopen(ident,z(ibuf),1)
 CALL pack(rz(1),ident,mcb)
 CALL CLOSE(ident,1)
 CALL wrttrl(mcb)
 
 CALL ssg2b(awy,ident,0,ac,0,2,1,scr)
 
!     PERFORM MULTIPLY TO GET COMPRESSIBLITY MATRIX
 
 
!     UNPACK ROW OF AC INTO CORE
 
 mcb(1) = ac
 CALL rdtrl(mcb)
 nrow = mcb(3)
 IF(nz < 2*nrow) GO TO 1008
 to2 = 2
 i2 = 1
 n2 = nrow
 incr2 = 1
 
 CALL gopen(ac,z(ibuf),0)
 CALL unpack(*40,ac,dz(1))
 GO TO 60
 
!     AC IS NULL
 
 40 DO  i=1,nrow
   dz(i) = 0.0D0
 END DO
 
!     SET UP TO CREATE KC MATRIX
 
 60 CALL CLOSE(ac,1)
 
 dkcomp = DBLE(kcomp)
 CALL gopen(kc,z(ibuf),1)
 CALL makmcb(mcb,kc,nrow,1,2)
 
!     LOOP OVER NON-ZERO TERMS OF AC TO CREATE KC
 
 DO  i=1,nrow
   CALL bldpk(2,2,kc,0,0)
   IF(dz(i) == 0.0D0) GO TO 80
   DO  j=1,nrow
     irow = j
     val = dkcomp * dz(j) * dz(i)
     CALL zblpki
   END DO
   80 CALL bldpkn(kc,0,mcb)
 END DO
 CALL CLOSE(kc,1)
 
!     WRITE TRAILER
 
 CALL wrttrl(mcb)
 RETURN
 
!     ERRORS
 
 1008 CALL mesage(-8,0,NAME)
 RETURN
END SUBROUTINE gfscom
