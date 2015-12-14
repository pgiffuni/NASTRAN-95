SUBROUTINE ctrnsp (ix,x,nx,filea,b,sr1fil)
     
!     TRANS WILL DO AN INCORE TRANSPOSE OF THE UPPER TRIANGLE OF ACTIVE
!     ELEMENTS
 
 
 INTEGER, INTENT(OUT)                     :: ix(1)
 REAL, INTENT(IN)                         :: x(1)
 INTEGER, INTENT(IN)                      :: nx
 INTEGER, INTENT(IN)                      :: filea(7)
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN)                      :: sr1fil
 EXTERNAL           lshift    ,rshift   ,orf      ,complf
 INTEGER :: typea    , eol       ,sysbuf   ,orf      ,lshift   ,  &
     NAME(2)   ,rshift   ,rdp      ,eor      , cdp       ,complf
 DOUBLE PRECISION :: di(2)
 DIMENSION  iii(6)
 COMMON   /machin/  mach      ,ihalf
 COMMON   /zntpkx/  ia(4)     ,ii       ,eol      ,eor
 COMMON   /system/  sysbuf
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      , rdp       ,csp      ,cdp
 COMMON   /TYPE  /  jprec(2)  ,nwds(4)
 EQUIVALENCE        (iii(3),di(1))
 DATA      NAME /   4HCTRN ,4HSP  /
 
 
 num   = rshift(complf(0),1)
 iobuf = nx - 4*sysbuf
 ifile = filea(1)
 
!     POSITION INPUT FILE AT START OF THE UPPER TRIANGLE
 
 CALL skprec (filea(1),b+1)
 typea = filea(5)
 ncol  = filea(2)
 no    = 0
 istor = 1
 iprec = jprec(typea)
 incr  = nwds(typea) + 1
 k     = 1
 20 CALL intpk (*70,filea(1),0,typea,0)
 30 CALL zntpki
 IF (ii > k) GO TO 50
 
!     PACK I AND J IN ONE WORD AND STORE IT AND THE NONZERO VALUE
!     IN CORE
 
 l  = orf(lshift(ii,ihalf),k+b)
 no = no + 1
 ix(istor  ) = l
 ix(istor+1) = ia(1)
 ix(istor+2) = ia(2)
 ix(istor+3) = ia(3)
 ix(istor+4) = ia(4)
 istor = istor+incr
 IF (istor+incr > iobuf) GO TO 140
 IF (eol == 0.0) THEN
   GO TO    30
 ELSE
   GO TO    70
 END IF
 50 IF (eor == 0) CALL skprec (filea(1),1)
 70 k = k + 1
 IF (k+b <= ncol) GO TO 20
 CALL REWIND (filea(1))
 
!     ALL ELEMENTS ARE IN CORE.  WRITE THEM OUT IN THE TRANSPOSED ORDER
 
 ifile = sr1fil
 CALL OPEN (*120,sr1fil,ix(iobuf),wrtrew)
 istor = istor - incr
 DO  i = 1,no
   k = num
   DO  j = 1,istor,incr
     IF (ix(j) > k) CYCLE
     kk = j
     k  = ix(j)
   END DO
   
!     UNPACK I AND J, AND WRITE OUT I,J,AND A(I,J)
   
   iii(1) = rshift(k,ihalf)
   iii(2) = k - lshift(iii(1),ihalf)
   ix(kk) = num
   IF (iprec == 2) GO TO 90
   di(1) = x(kk+1)
   di(2) = 0.d0
   IF (typea > 2) di(2) = x(kk+2)
   GO TO 100
   90 iii(3) = ix(kk+1)
   iii(4) = ix(kk+2)
   iii(5) = 0
   iii(6) = 0
   IF (typea <= 2) GO TO 100
   iii(5) = ix(kk+3)
   iii(6) = ix(kk+4)
   100 CONTINUE
   CALL WRITE (sr1fil,iii(1),6,0)
   IF (kk == istor) istor = istor - incr
 END DO
 
!     WRITE A TRAILER RECORD ON THE FILE
 
 iii(1) = -1
 CALL WRITE (sr1fil,iii(1),6,0)
 CALL CLOSE (sr1fil,rew)
 RETURN
 
 120 no = -1
 GO TO 150
 140 no = -8
 150 CALL mesage (no,ifile,NAME)
 RETURN
END SUBROUTINE ctrnsp
