SUBROUTINE gfshc(awy,nuy,hc,ident,ac,mrow)
     
!     ROUTINE TO GENERATE CONSTRAINT MATRIX FOR PURELY INCOMPRESSIBLE
!     FORMULATION WHEN NO SPC'S ARE ON FLUID
 
 
 INTEGER, INTENT(IN OUT)                  :: awy
 INTEGER, INTENT(IN)                      :: nuy
 INTEGER, INTENT(IN OUT)                  :: ident
 INTEGER, INTENT(IN)                      :: ac
 INTEGER, INTENT(OUT)                     :: mrow
 DOUBLE PRECISION :: dz(1)    ,dterm    ,val
 
 REAL :: rz(1)
 
 INTEGER :: z        ,sysbuf   ,mcb(7)   ,NAME(2)  ,hc  &
     ,ti1      ,to1      ,to2      , scr
 
!     OPEN CORE
 
 COMMON / zzzzzz /        z(1)
 
!     SYSTEM COMMON
 
 COMMON / system /       sysbuf
 
!     PACK - UNPACK COMMON BLOCKS
 
 COMMON / packx /        ti1      ,to1      ,i1       ,n1 ,incr1
 COMMON / unpakx /       to2      ,i2       ,n2       ,incr2
 COMMON / zblpkx /       a(4)     ,irow
 
 EQUIVALENCE   ( z(1) , rz(1) , dz(1) ) ,( val , a(1) )
 
 DATA NAME / 4HGFSH , 4HC    /
 
 
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
 
 60 CALL CLOSE(ac,1)
 
!     LOCATE LARGEST TERM IN AC
 
 dterm = -1.0D10
 DO  i=1,nrow
   IF(dz(i) <= dterm) CYCLE
   mrow = i
   dterm = dz(i)
 END DO
 
!     GENERATE THE HC MATRIX
 
 CALL makmcb(mcb,hc,nrow,1,2)
 CALL gopen(hc,z(ibuf),1)
 
!     GENERATE COLUMNS UP TO MROW
 
 IF(mrow == 1) GO TO 230
 mr = mrow - 1
 DO  ir = 1,mr
   CALL bldpk(2,2,hc,0,0)
   irow = ir
   val = 1.0D0
   CALL zblpki
   irow = mrow
   val = -dz(ir) / dterm
   CALL zblpki
   CALL bldpkn(hc,0,mcb)
 END DO
 
!     PACK OUT NULL COLUMN FOR MROW
 
 230 CALL bldpk(2,2,hc,0,0)
 CALL bldpkn(hc,0,mcb)
 
!     GENERATE REMAINING ROWS
 
 IF(mrow >= nrow) GO TO 250
 mr = mrow + 1
 DO  ir=mr,nrow
   CALL bldpk(2,2,hc,0,0)
   irow = mrow
   val = -dz(ir) / dterm
   CALL zblpki
   irow = ir
   val = 1.0D0
   CALL zblpki
   CALL bldpkn(hc,0,mcb)
 END DO
 
 250 CALL CLOSE(hc,1)
 CALL wrttrl(mcb)
 
 RETURN
 
!     ERRORS
 
 1008 CALL mesage(-8,0,NAME)
 
 RETURN
END SUBROUTINE gfshc
