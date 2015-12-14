SUBROUTINE read7 (nr1,olama,ophia,nlama,nphia)
     
!     READ7  COPIES NR VECTORS FROM OPHIA TO NPHIA -
!     IT ALSO PLACES THE EIGENVALUES ON NLAMA
!     THIS ROUTINE HANDLES BOTH SINGLE AND DOUBLE PRECISION
 
 
 INTEGER, INTENT(IN)                      :: nr1
 INTEGER, INTENT(IN)                      :: olama
 INTEGER, INTENT(IN)                      :: ophia
 INTEGER, INTENT(IN OUT)                  :: nlama
 INTEGER, INTENT(IN OUT)                  :: nphia
 INTEGER :: sysbuf,ix(7),NAME(2),sgldbl
 REAL :: x(7)
 DOUBLE PRECISION :: dcore(2),dx
 COMMON /system/  sysbuf
 COMMON /unpakx/  itb,ii,jj,incur
 COMMON /packx /  it1,it2,iip,jjp,incrp
 COMMON /zzzzzz/  core(1)
 EQUIVALENCE      (dcore(1),core(1)) , (x(1),dx)
 DATA    NAME  /  4HREAD,4H7   /
 
!     GET ORGANIZED
 
 nr    = nr1
 lc    = korsz(core)
 ibuf1 = lc - sysbuf + 1
 ibuf2 = ibuf1 -sysbuf
 ibuf3 = ibuf2 -sysbuf
 ibuf4 = ibuf3 -sysbuf
 ix(1) = ophia
 CALL rdtrl (ix)
 nrow  = ix(3)
 ii    = 1
 jj    = nrow
 it1   = ix(5)
 it2   = it1
 itb   = it1
 dcore(1) = 0.0D0
 incrp = 1
 ASSIGN 12 TO sgldbl
 IF (itb == 2) ASSIGN 16 TO sgldbl
 incur = 1
 
!     OPEN OLD FILES
 
 CALL gopen (olama,core(ibuf1),0)
 CALL fwdrec (*3010,olama)
 CALL gopen (ophia,core(ibuf2),0)
 
!     OPEN NEW FILES TO WRITE
 
 CALL gopen (nlama,core(ibuf3),1)
 CALL gopen (nphia,core(ibuf4),1)
 
!     START COPY LOOP
 
 CALL makmcb (ix,nphia,nrow,ix(4),it2)
 DO  i = 1,nr
   CALL READ (*3010,*3020,olama,x,7,0,ifl)
   ii = 0
   CALL unpack (*150,ophia,dcore(2))
   GO TO sgldbl, (12,16)
   12 x(1) = SQRT(x(6))
   DO  j = 1,nrow
     core(j+2) = core(j+2)/x(1)
   END DO
   GO TO 20
   16 dx = SQRT(x(6))
   DO  j = 1,nrow
     dcore(j+1) = dcore(j+1)/dx
   END DO
   20 iip = ii
   jjp = jj
   CALL pack (dcore(2),nphia,ix)
   30 dx = x(3)
   CALL WRITE (nlama,dx,2,1)
   CYCLE
   
!     NULL COLUMN
   
   150 iip = 1
   jjp = 1
   CALL pack (dcore,nphia,ix)
   GO TO 30
 END DO
 CALL CLOSE (olama,1)
 CALL CLOSE (ophia,1)
 CALL CLOSE (nlama,2)
 CALL CLOSE (nphia,1)
 RETURN
 
!     ERRORS
 
 3010 nn = -2
 3012 ifile = olama
 CALL mesage (nn,ifile,NAME)
 RETURN
 3020 nn = -3
 GO TO 3012
END SUBROUTINE read7
