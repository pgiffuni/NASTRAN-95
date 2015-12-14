SUBROUTINE cfeer3
     
!     CFEER3 IS A DRIVER ROUTINE WHICH PERFORMS THE TRIDIAGONAL
!     REDUCTION FOR THE COMPLEX FEER METHOD
 
 INTEGER :: switch   ,cdp      ,sqr      ,sysbuf     , NAME(2)
 DOUBLE PRECISION :: lambda   ,dz(1)
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew     ,  &
     rew      ,norew    ,eofnrw   ,rsp        , rdp      ,csp      ,cdp      ,sqr
 COMMON  /system/  ksystm(65)
 COMMON  /zzzzzz/  z(1)
 COMMON  /feerxc/  lambda(2),switch   ,mreduc   ,nord       ,  &
     idiag    ,epsdum(2),northo   ,nord2      ,  &
     nord4    ,nordp1   ,xcdum(12),minopn
 COMMON  /feeraa/  ik(7)    ,im(7)    ,ib(7)    ,ilam(7)    ,  &
     iphi(7)  ,dudxx    ,iscr(11) ,dumaa(84)  , mcbvec(7)
 EQUIVALENCE       (dz(1)   ,z(1)   ) ,(ksystm(55),iprec)   ,  &
     (ksystm(1),sysbuf) ,(ksystm( 2),nout )
 DATA     NAME  /  4HCFEE,4HR3  /
 
!     SCRATCH FILE AND BUFFER ALLOCATION
 
!     FILE  5  CONTAINS THE ELEMENTS OF REDUCED TRIDIAGONAL MATRIX
!     FILE  7  CONTAINS THE ORTHOGONAL VECTOR PAIRS (NUMBER OF
!              VECTOR PAIRS = NORTHO)
 
!     BUFFER Z(IBUF1) IS LOCAL SCRATCH BUFFER
!     BUFFER Z(IBUF2) IS LOCAL SCRATCH BUFFER
!     BUFFER Z(IBUF3) IS USED BY FILE 5
 
 
!     COMPUTE STORAGE ALLOCATIONS
 
 nz    = korsz(z)
 ibuf1 = nz    - sysbuf
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 itop  = ibuf3
 
!     COMPUTE LOCATIONS OF RIGHT-HANDED VECTORS
 
 iv1 = 1
 iv2 = iv1 + nord4
 iv3 = iv2 + nord4
 iv4 = iv3 + nord4
 iv5 = iv4 + nord4
 
!     TEST FOR INSUFFICIENT CORE
 
 iend = iprec*(5*nord4+1)
 IF (iend > itop) GO TO 70
 
!     COMPUTE LOCATIONS OF LEFT-HANDED VECTORS
 
 iv1l = iv1 + nord2
 iv2l = iv2 + nord2
 iv3l = iv3 + nord2
 iv4l = iv4 + nord2
 iv5l = iv5 + nord2
 iopn = itop- iend
 IF (idiag /= 0) WRITE (nout,510) iopn
 IF (iopn < minopn) minopn = iopn
 
!     INITIALIZE SCRATCH FILE TO CONTAIN TRIDIAGONAL ELEMENTS
 
 CALL gopen (iscr(5),z(ibuf3),wrtrew)
 
!     GENERATE MATRIX CONTROL BLOCK FOR SCRATCH FILE TO CONTAIN
!     ORTHOGONAL VECTORS (LEFT VECTOR PACKED IMMEDIATELY AFTER
!     RIGHT, I. E., EACH COLUMN CONTAINS RIGHT VECTOR FOLLOWED BY
!     LEFT VECTOR)
 
 jprec = iprec + 2
 CALL makmcb (mcbvec(1),iscr(7),nord2,2,jprec)
 
!     PERFORM DOUBLE PRECISION FEER
 
 IF (iprec == 2) CALL cfer3d (dz(iv1),dz(iv1l), dz(iv2),dz(iv2l),  &
     dz(iv3),dz(iv3l), dz(iv4),dz(iv4l), dz(iv5),dz(iv5l), z(ibuf1),z(ibuf2))
 
!     PERFORM SINGLE PRECISION FEER
 
 IF (iprec /= 2) CALL cfer3s (z(iv1),z(iv1l), z(iv2),z(iv2l),  &
     z(iv3),z(iv3l), z(iv4),z(iv4l), z(iv5),z(iv5l), z(ibuf1),z(ibuf2))
 
!     TERMINATE SCRATCH FILE CONTAINING TRIDIAGONAL ELEMENTS
 
 CALL CLOSE (iscr(5),norew)
 RETURN
 
 70 iend = (iend-itop)/1000 + 1
 WRITE  (nout,80) iend
 80 FORMAT (5H0NEED,i4,17HK more core words)
 CALL mesage (-8,0,NAME)
 510 FORMAT (1H ,i10,36H single PRECISION words of OPEN core,  &
     29H NOT used (SUBROUTINE cfeer3))
END SUBROUTINE cfeer3
