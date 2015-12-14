SUBROUTINE transp (ix,x,nx,filea,b,sr1fil)
     
!     TRANSP WILL DO AN INCORE TRANSPOSE OF THE UPPER TRIANGLE OF
!     ACTIVE ELEMENTS
!     (OUT-OF-CORE TRANSPOSE IS DONE BY TRNSP)
 
 
 INTEGER, INTENT(OUT)                     :: ix(1)
 REAL, INTENT(IN)                         :: x(1)
 INTEGER, INTENT(IN)                      :: nx
 INTEGER, INTENT(IN)                      :: filea(7)
 INTEGER, INTENT(IN)                      :: b
 INTEGER, INTENT(IN)                      :: sr1fil
 EXTERNAL           lshift    ,rshift   ,orf      ,complf
 INTEGER :: typea    , eol       ,sysbuf   ,orf      ,lshift   ,  &
     NAME(2)   ,rshift   ,rdp      ,eor      , complf
 DOUBLE PRECISION :: di
 DIMENSION  iii(4)
 COMMON   /machin/  mach      ,ihalf
!     COMMON   /DESCRP/  LENGTH    ,MAJOR
 COMMON   /zntpkx/  ia(4)     ,ii       ,eol      ,eor
 COMMON   /system/  sysbuf
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      , rdp
 EQUIVALENCE        (iii(3),di)
 DATA      NAME  /  4HTRAN, 4HSP   /
 
 
 num   = rshift(complf(0),1)
 iobuf = nx - 4*sysbuf
 ifile = filea(1)
 
!     POSITION INPUT FILE AT START OF THE UPPER TRIANGLE
 
 n = b + 1
 CALL skprec (filea,n)
 typea = filea(5)
 ncol  = filea(2)
 no    = 0
 istor = 1
 k     = 1
 5 CALL intpk (*50,filea(1),0,typea,0)
 10 CALL zntpki
 IF (ii > k) GO TO 40
 
!     PACK I AND J IN ONE WORD AND STORE IT AND THE NONZERO VALUE
!     IN CORE
 
 l  = orf(lshift(ii,ihalf),k+b)
 no = no + 1
 ix(istor   ) = l
 ix(istor+ 1) = ia(1)
 istor = istor + 2
 IF (typea /= rdp) GO TO 20
 ix(istor) = ia(2)
 istor = istor + 1
 20 IF (istor+3 > iobuf) GO TO 230
 IF (eol == 0.0) THEN
   GO TO    10
 ELSE
   GO TO    50
 END IF
 40 IF (eor == 0) CALL skprec (filea,1)
 50 k = k + 1
 IF (k+b <= ncol) GO TO 5
 CALL REWIND (filea(1))
 
!     ALL ELEMENTS ARE IN CORE.  WRITE THEM OUT IN THE TRANSPOSED ORDER
 
 ifile = sr1fil
 CALL OPEN (*200,sr1fil,ix(iobuf),wrtrew)
 incr  = typea + 1
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
   IF (incr == 3) GO TO 95
   di = x(kk+1)
   GO TO 96
   95 iii(3) = ix(kk+1)
   iii(4) = ix(kk+2)
   96 CONTINUE
   CALL WRITE (sr1fil,iii(1),4,0)
   IF (kk == istor) istor = istor - incr
 END DO
 
!     WRITE A TRAILER RECORD ON THE FILE
!     NOTE - FORMAL GINO FILE TRAILER IS NOT GENERATED HERE
 
 iii(1) = -1
 CALL WRITE (sr1fil,iii(1),4,0)
 CALL CLOSE (sr1fil,rew)
 RETURN
 
 200 no = -1
 GO TO 250
 230 no = -8
 250 CALL mesage (no,ifile,NAME)
 RETURN
END SUBROUTINE transp
