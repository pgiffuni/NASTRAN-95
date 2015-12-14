SUBROUTINE read4 (lama,phi,scr1,eps,mass)
     
!     READ4 WILL TEST FOR CLOSE AND EQUAL ROOTS AND MAKE SURE THE
!     CORRESPONDING VECTORS ARE ORTHOGONAL
 
 
 INTEGER, INTENT(IN)                      :: lama
 INTEGER, INTENT(IN)                      :: phi(7)
 INTEGER, INTENT(IN)                      :: scr1
 REAL, INTENT(IN)                         :: eps
 INTEGER, INTENT(IN OUT)                  :: mass
 INTEGER :: NAME(2)   , rsp      ,phi1(7)
 INTEGER :: rdrew     ,wrtrew
!WKBI ALPHA-OSF 9/94
 
 DOUBLE PRECISION :: dz(1)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm       ,uwm
 COMMON /zzzzzz/  z(1)
 COMMON /system/  ksystm(65)
 COMMON /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp
 COMMON /unpakx/  itype     ,iunpak   ,junpak   ,incr
 COMMON /packx /  itypa     ,itypb    ,ipak     ,jpak     , incrx
 EQUIVALENCE      (ksystm(1),isys)    ,(ksystm(2),iout)   , (dz(1),z(1))
 DATA    NAME  /  4HREAD,4H4   /
 
 ncol  = phi(2)
 nrow  = phi(3)
 nz    = korsz(z)
 ibuf  = nz - isys
 ibuf1 = ibuf - isys
 ibuf2 = ibuf1 - isys
 iclos = 0
 idid  = 0
 ipr   = phi(5)
 rmult = .01
 itype = rsp
 iunpak= 1
 junpak= nrow
 incr  = 1
 itypa = rsp
 itypb = rsp
 ipak  = 1
 jpak  = nrow
 incrx = 1
 epsi  = eps
 IF (eps <= 0.) epsi = .0001
 nz = nz - isys - isys - 1 - isys
 CALL makmcb (phi1,scr1,nrow,2,rsp)
 ifile = lama
 CALL gopen (lama,z(ibuf),0)
 CALL READ (*170,*10,lama,z(1),nz,1,n)
 GO TO 180
 10 CALL CLOSE (lama,rew)
 
!     REJECT ALL BUT VALUES FOR WHICH VECTORS EXIST
 
 n  = phi(2)
 nz = nz -n
 IF (nz < nrow) GO TO 180
 ifile = phi(1)
 CALL gopen (phi,z(ibuf),0)
 ipos = 1
 i    = 1
 eps1 = rmult
 20 CONTINUE
 IF (ABS(z(i))+ABS(z(i+1)) < eps1) GO TO 1111
 IF (z(i+1) == 0.0) GO TO 110
 IF (ABS(1.0-z(i)/z(i+1)) > eps1) GO TO 100
 1111 IF (iclos /= 0) GO TO 110
 iclos = i
 GO TO 110
 30 num  = i - iclos + 1
 eps1 = rmult
 
!     NUM   = NUMBER OF CLOSE ROOTS IN THIS GROUP
!     ICLOS = THE INDEX OF THE FIRST CLOSE ROOT
 
 IF (idid == 1) GO TO 40
 idid  = 1
 ifile = scr1
 CALL gopen (scr1,z(ibuf1),wrtrew)
 40 ii = n + 1
 50 IF (ipos == iclos) GO TO 70
 ifile = phi(1)
 CALL unpack (*190,phi,z(ii))
 CALL pack (z(ii),scr1,phi1)
 ipos = ipos + 1
 GO TO 50
 70 CONTINUE
 
!     CHECK FOR CORE OVERFLOW
!     EIGENVALUES + EIGENVECTORS + GEN. MASS + ACCUM.
 
 kore = ii + num*nrow + num*num + n + n + 3
 IF (kore > nz) GO TO 160
 DO  j = 1,num
   CALL unpack (*190,phi,z(ii))
   ipos = ipos + 1
   ii   = ii + nrow
   IF (ii+nrow >= nz) GO TO 180
 END DO
 ij = ii + n + n + 3
 ii = ii/2 + 1
 CALL ortck (z(n+1),mass,z(ibuf2),num,nrow,z(ij),dz(ii),epsi)
 ii = n + 1
 DO  j = 1,num
   CALL pack (z(ii),scr1,phi1)
   ii = ii + nrow
 END DO
 iclos = 0
 100 IF (iclos /= 0) GO TO 30
 110 i = i + 1
 IF (i     <   n) GO TO 20
 IF (iclos /=   0) GO TO 30
 IF (idid  ==   0) GO TO 150
 IF (ipos > ncol) GO TO 121
 DO  i = ipos,ncol
   CALL unpack (*190,phi,z)
   CALL pack (z(1),scr1,phi1)
 END DO
 121 CALL wrttrl (phi1)
 
!     COPY VECTORS FROM SCR1 TO PHI
 
 CALL CLOSE (phi,rew)
 CALL CLOSE (scr1,rew)
 CALL gopen (phi,z(ibuf),1)
 CALL gopen (scr1,z(ibuf1),rdrew )
 CALL makmcb (phi,phi,nrow,2,ipr)
 itypb = ipr
 DO  i = 1,n
   CALL unpack (*190,scr1,z)
   CALL pack (z,phi,phi)
 END DO
 CALL wrttrl (phi)
 CALL CLOSE (scr1,rew)
 150 CALL CLOSE (phi,rew)
 RETURN
 
 160 eps2 = eps1/10.
 WRITE  (iout,165) uwm,num,i,eps1,eps2
 165 FORMAT (a25,' 3142, INSUFFICIENT CORE STORAGE FOR EIGENVECTORS ',  &
     'ASSOCIATED WITH',i4,' MULTIPLE EIGENVALUES STARTING WITH',  &
     /28X,'MODE NUMBER',i4,' USING CURRENT MULTIPLE ROOT ',  &
     'CRITERIA. CRITERIA REDUCED FROM ',1P,e12.5,' TO ',e12.5)
 eps1 = eps2
 i = iclos
 GO TO 20
 170 no = -2
 GO TO 200
 180 no = -8
 GO TO 200
 190 no = -7
 200 CALL mesage (no,ifile,NAME)
 RETURN
END SUBROUTINE read4
