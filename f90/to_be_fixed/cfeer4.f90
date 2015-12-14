SUBROUTINE cfeer4
     
!     CFEER4 OBTAINS THE EIGENVALUES AND EIGENVECTORS FROM THE
!     REDUCED TRIDIAGONAL MATRIX FOR THE COMPLEX FEER METHOD
 
 LOGICAL :: no b     ,decrem   ,qpr      ,lz(1)     , dpmach
 INTEGER :: NAME(2)  ,iz(1)    ,eor
 INTEGER :: wrtrew
 DOUBLE PRECISION :: lambda   ,eps      ,dz(1)    ,d(4)      , lam1(2)
 DIMENSION         s(8)     ,dmp1(2)  ,alam(2)  ,dm(2)     ,  &
     STATUS(2),accept(2),reject(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /  ufm      ,uwm      ,uim
 COMMON  /feeraa/  ikmb(21) ,ilam(7)  ,iphi(7)  ,idmpfl    ,  &
     iscr(11) ,dumaa(84),mcbvec(7)
 COMMON  /feerxc/  lambda(2),swdum    ,mreduc   ,nord      ,  &
     idiag    ,eps      ,northo   ,nord2     ,  &
     nord4    ,nordp1   ,nswp(2)  ,no b      ,  &
     it       ,ten2mt   ,tenmht   ,nstart    ,  &
     qpr      ,jreg     ,noreg    ,nzero     , tenmtt   ,minopn
 COMMON  /zzzzzz/  z(1)
 COMMON  /unpakx/  iprc     ,ii       ,nn       ,incr
 COMMON  /system/  ksystm(65)
 COMMON  /names /  rd       ,rdrew    ,wrt      ,wrtrew    ,  &
     rew      ,norew    ,eofnrw
 EQUIVALENCE       (ksystm(2 ),nout ) ,(nrow,mreduc)       ,  &
     (ksystm(55),iprec) ,(d(1),s(1)  )       , (z(1),iz(1),lz(1),dz(1))
 DATA     NAME  /  4HCFEE,4HR4  /
 DATA     accept,  reject/4H  ac,4HCEPT,4H -re,4HJECT      /
 
!     CORE ALLOCATION FOR ALLMAT
 
!     CONTENTS                SIZE             POINTER   TYPE   NAME
!     --------                ----             -------   ----   ----
!     INPUT MATRIX--VECTORS   2*NROW*NROW      IA        COMP   A
!     EIGENVALUES             2*NROW           IL        COMP   LAM
!     H  MATRIX               2*NROW*NROW      IH        COMP   H
!     HL MATRIX               2*NROW*NROW      IHL       COMP   HL
!     VECTOR STORAGE          2*NROW           IV        COMP   VEC
!     MULTIPLIERS             2*NROW           IM        COMP   MULT
!     INTH                    NROW             INTH      INTG   INTH
!     INTQ                    NROW             INTQ      LOGL   INTQ
 
!     CORE ALLOCATION AFTER ALLMAT IS FINISHED
 
!     ALLMAT OUTPUT EIGENVECTORS               IA
!     EIGENVALUES                              IL
!     ORDER OF EXTRACTION                      IH
!     THEORETICAL ERRORS                       IHL
!     NOT USED                                 IV,IM
!     STATUS OF SOLUTIONS                      INTH
!     DISTANCES FROM CENTER                    INTQ
!     VARIABLE PRECISION PHYSICAL EIGENVECTORS IV1
!     VARIABLE PRECISION ORTHOGONAL VECTORS    IV2
 
!     DEFINITION OF INTERNAL PARAMETERS
 
!     DMP1     = D-SUB-M-PLUS-1 = EXTRANEOUS OFF-DIAGONAL ELEMENT
!                OF REDUCED TRIDIAGONAL MATRIX, USED FOR COMPUTING
!                THEORETICAL ERRORS
!     DM       = FINAL OFF-DIAGONAL ELEMENT OF REDUCED TRIDIAGONAL
!                MATRIX
!     NO B     = LOGICAL INDICATOR FOR ABSENCE OF DAMPING MATRIX B
!     DECREM   = LOGICAL INDICATOR FOR DECREMENTED SIZE OF REDUCED
!                PROBLEM
!     NROW     = SIZE OF THE REDUCED PROBLEM (EQUIVALENT TO MREDUC)
!     RMS      = ROOT-MEAN-SQUARE OF EIGENVALUES, USED IN RIGID-BODY
!                ERROR TEST
!     NOTE.....SEE LISTING OF CFCNTL FOR ADDITIONAL DEFINITIONS
 
 IF (qpr) WRITE (nout,600)
 dpmach = iprec == 2
 nord8  = 2*nord4
 decrem = .false.
 4 nrow2  = 2*nrow
 nrowsq = nrow*nrow2
 
!     ALLOCATE CORE FOR ALLMAT
 
 ia  = 1
 il  = ia  + nrowsq
 ih  = il  + nrow2
 ihl = ih  + nrowsq
 iv  = ihl + nrowsq
 im  = iv  + nrow2
 inth= im  + nrow2
 intq= inth+ nrow
 
!     ALLOCATE CORE FOR PHYSICAL EIGENVECTORS (LEFT FOLLOWS RIGHT)
 
 iv1 = intq + nrow
 iv2 = iv1  + nord8
 IF (dpmach .AND. MOD(iv2,2) == 0) iv2 = iv2 + 1
 iv1x = iv1 - 1
 
!     TEST FOR INSUFFICIENT CORE
 
 nz    = korsz(z(1))
 ibuf1 = nz    - ksystm(1)
 ibuf2 = ibuf1 - ksystm(1)
 iopn  = ibuf2 - (iv2 + nord8)
 IF (idiag /= 0) WRITE (nout,610) iopn
 IF (iopn  <= 0) CALL mesage (-8,0,NAME(1))
 IF (iopn < minopn) minopn = iopn
 IF (nswp(2)  <  0) GO TO 209
 
!     CONSTRUCT REDUCED TRIDIAGONAL MATRIX
 
 DO  i = ia,il
   z(i) = 0.
 END DO
 nrow22 = nrow2 + 2
 CALL gopen (iscr(5),z(ibuf1),rdrew)
 nw  = 4*iprec
 eor = 1
 m   = 0
 nrow1 = nrow - 1
 
!     ENTER LOOP
 
 DO  i = 1,nrow
   i1 = i - 1
   CALL READ (*420,*430,iscr(5),s(1),nw,eor,m)
   IF (qpr .AND. .NOT.dpmach) WRITE (nout,620) i,(s(j),j=1,4)
   IF (qpr .AND.      dpmach) WRITE (nout,630) i,(d(j),j=1,4)
   
!     ALLMAT ACCEPTS ONLY SINGLE PRECISION ARRAY
   
   j = ia + nrow22*i1
   IF (.NOT.dpmach) GO TO 15
   
!     LOAD MAIN DIAGONAL ELEMENT
   
   z(j  ) = d(3)
   z(j+1) = d(4)
   IF (i /= nrow1) GO TO 12
   
!     SAVE LAST OFF-DIAGONAL ELEMENT
   
   dm(1) = d(1)
   dm(2) = d(2)
   12 IF (i == nrow) CYCLE
   
!     LOAD OFF-DIAGONAL ELEMENTS
   
   z(j+2) = d(1)
   z(j+3) = d(2)
   j  = j + nrow2
   z(j  ) = d(1)
   z(j+1) = d(2)
   CYCLE
   
!     LOAD MAIN DIAGONAL ELEMENT
   
   15 z(j  ) = s(3)
   z(j+1) = s(4)
   IF (i /= nrow1) GO TO 16
   
!     SAVE LAST OFF-DIAGONAL ELEMENT
   
   dm(1) = s(1)
   dm(2) = s(2)
   16 IF (i == nrow) CYCLE
   
!     LOAD OFF-DIAGONAL ELEMENTS
   
   z(j+2) = s(1)
   z(j+3) = s(2)
   j  = j + nrow2
   z(j  ) = s(1)
   z(j+1) = s(2)
 END DO
 
!     SAVE ERROR ELEMENT FROM TRIDIAGONAL MATRIX
 
 IF (.NOT.dpmach) GO TO 25
 dmp1(1) = d(1)
 dmp1(2) = d(2)
 GO TO 26
 25 dmp1(1) = s(1)
 dmp1(2) = s(2)
 26 CONTINUE
 IF (qpr) WRITE (nout,640) (z(i),i=1,nrowsq)
 CALL CLOSE (iscr(5),rew)
 IF (decrem) GO TO 30
 
!     DECREMENT THE REDUCED PROBLEM SIZE IF THE ERROR ELEMENT IS NULL
 
 IF (dmp1(1) /= 0. .OR. dmp1(2) /= 0.) GO TO 30
 mreduc = mreduc - 1
 WRITE (nout,570) uwm,mreduc
 IF (mreduc == 0) GO TO 440
 IF (dm(1) /= 0. .OR. dm(2) /= 0.) GO TO 29
 
!     NEW ERROR ELEMENT IS ALSO NULL. RESTORE ORIGINAL REDUCED SIZE.
 
 mreduc  = mreduc + 1
 dmp1(1) = SNGL(eps)
 WRITE (nout,590) uwm,mreduc,dmp1
 GO TO 30
 29 decrem = .true.
 GO TO 4
 
 30 CALL allmat (z(ia),z(il),z(ih),z(ihl),z(iv),z(im),z(inth),z(intq),  &
     nrow,nrow,inidum)
 
!     --------------- SPECIAL PRINT -------------------------
 
 IF (.NOT.qpr) GO TO 4429
 WRITE  (nout,4408)
 4408 FORMAT (1H0,10X,15HALLMAT executed,/,1H0)
 j = ih - 1
 WRITE  (nout,4420) (z(i),i=il,j)
 4420 FORMAT (1H0,11HEIGENVALUES, //,(1H ,2E16.8))
 WRITE  (nout,4422)
 4422 FORMAT (1H0,12HEIGENVECTORS,//)
 DO  i = 1,nrow
   l = ia + nrow2*(i-1)
   k = l  + nrow2 - 1
   WRITE (nout,4424) (z(j),j=l,k)
   
!     CHECK NORMALITY
   
   sumr = 0.
   sumi = 0.
   DO  j = l,k,2
     jj   = j + 1
     sumr = sumr + z(j)**2 - z(jj)**2
     sumi = sumi + 2.*z(j)*z(jj)
   END DO
   WRITE  (nout,7770) sumr,sumi
   7770 FORMAT (//,35H self inner-product of above vector, /,1H ,6X,  &
       11HREAL part =,e16.8,8X,16HIMAGINARY part =,e16.8)
   4424 FORMAT (//,(1H ,6E16.8))
 END DO
 4429 CONTINUE
!     -------------------------------------------------------
 
!     NORMALIZE THE EIGENVECTORS OUTPUT FROM ALLMAT
 
 IF (qpr) WRITE (nout,4422)
 DO  i = 1,nrow
   l = ia + nrow2*(i-1)
   k = l  + nrow2 - 1
   sumr = 0.
   sumi = 0.
   DO  j = l,k,2
     jj   = j + 1
     sumr = sumr + z(j)**2 - z(jj)**2
     sumi = sumi + 2.*z(j)*z(jj)
   END DO
   rsqrt= SQRT(SQRT(sumr**2 + sumi**2))
   IF (rsqrt > 0.) GO TO 34
   WRITE (nout,560) uwm,NAME
   CYCLE
   34 theta2= .5*ATAN2(sumi,sumr)
   sumr  = rsqrt*COS(theta2)
   sumi  = rsqrt*SIN(theta2)
   theta2= 1./(sumr**2 + sumi**2)
   sumr  = sumr*theta2
   sumi  =-sumi*theta2
   DO  j = l,k,2
     jj    = j + 1
     theta2= z(j)
     z(j ) = sumr*z(j)   - sumi*z(jj)
     z(jj) = sumi*theta2 + sumr*z(jj)
   END DO
   
!     -------------- SPECIAL PRINT --------------------------
   
   IF (.NOT.qpr) CYCLE
   WRITE (nout,4424) (z(j),j=l,k)
   
!     CHECK NORMALITY
   
   sumr = 0.
   sumi = 0.
   DO  j = l,k,2
     jj   = j + 1
     sumr = sumr + z(j)**2 - z(jj)**2
     sumi = sumi + 2.*z(j)*z(jj)
   END DO
   WRITE (nout,7770) sumr,sumi
!     -------------------------------------------------------
   
 END DO
 
!     COMPUTE THEORETICAL EIGENVALUE ERRORS
 
 IF (qpr) WRITE (nout,650) dmp1
 ihl1 = ihl - 1
 DO  i = 1,nrow
   k = il + 2*(i-1)
   denom = SQRT(z(k)**2 + z(k+1)**2)
   IF (denom > 0.) GO TO 40
   WRITE (nout,550) uim,i
   denom = 1.e-10
   40 denom = 1./denom
   k  = ia + nrow2*i - 2
   kk = k  + 1
   j  = ihl1 + i
   z(j) = denom*SQRT((dmp1(1)*z(k) - dmp1(2)*z(kk))**2  &
       + (dmp1(1)*z(kk) + dmp1(2)*z(k))**2)
   IF (qpr) WRITE (nout,660) i,z(j),z(k),z(kk),denom
 END DO
 
!     RECOVER PHYSICAL EIGENVALUES
 
 rms = 0.
 IF (no b) GO TO 54
 alam(1) = lambda(1)
 alam(2) = lambda(2)
 GO TO 55
 54 alam(1) = lambda(1)**2 - lambda(2)**2
 alam(2) = 2.d0*lambda(1)*lambda(2)
 55 DO  i = 1,nrow
   k  = il + 2*(i-1)
   kk = k  + 1
   denom = z(k)**2 + z(kk)**2
   IF (denom == 0.) denom = 1.e-20
   denom = 1./denom
   z( k) = denom*z( k) + alam(1)
   z(kk) =-denom*z(kk) + alam(2)
   IF (no b) GO TO 60
   GO TO 70
   
!     DAMPING MATRIX ABSENT
   
   60 rsqrt  = SQRT(SQRT(z(k)**2 + z(kk)**2))
   theta2 = .5*ATAN2(z(kk),z(k))
   z( k)  = rsqrt*COS(theta2)
   z(kk)  = rsqrt*SIN(theta2)
   IF (z(kk) >= 0.) GO TO 70
   z( k)  =-z( k)
   z(kk)  =-z(kk)
   
!     COMPUTE RMS FOR RIGID-BODY ERROR TEST
   
   70 rms = rms + SQRT((z(k)**2-z(kk)**2)**2 + 4.*(z(k)*z(kk))**2)
 END DO
 rms = SQRT(rms)/FLOAT(nrow)
 IF (qpr) WRITE (nout,800) rms
 j  = ih - 1
 IF (qpr) WRITE (nout,4420) (z(i),i=il,j)
 
!     PERFORM RIGID-BODY ERROR TEST
 
 IF (rms < 1.e-20) rms = 1.e-20
 rms = 1./rms
 DO  i = 1,nrow
   k = il + 2*(i-1)
   j = ihl1 + i
   IF (rms*SQRT(z(k)**2+z(k+1)**2) <= tenmtt) z(j) = 0.
 END DO
 
!     COMPUTE DISTANCES OF EIGENVALUES TO CENTER OF NEIGHBORHOOD
 
 alam(1) = lambda(1)
 alam(2) = lambda(2)
 jj = intq - 1
 kk = ih   - 1
 ll = inth - 1
 DO  i = 1,nrow
   j  = jj + i
   k  = il + 2*(i-1)
   z(j) = SQRT((alam(1) - z(k))**2 + (alam(2)-z(k+1))**2)
   
!     LOAD ORDER OF EXTRACTION
   
   k = kk + i
   iz(k) = i
   
!     LOAD STATUS OF EACH SOLUTION
   
   k = ll + i
   lz(k) = .false.
   j = ihl1 + i
   IF (z(j) < SNGL(eps)) lz(k) = .true.
 END DO
 
!     SORT EIGENVALUES ACCORDING TO DISTANCE FROM CURRENT CENTER
 
 IF (nrow == 1) GO TO 150
 ll = nrow - 1
 DO  i = 1,ll
   k  = jj + i
   i1 = kk + i
   lll= i  + 1
   DO  j = lll,nrow
     l  = jj + j
     IF (z(k) < z(l)) CYCLE
     unidum = z(l)
     z(l)   = z(k)
     z(k)   = unidum
     i2     = kk + j
     inidum = iz(i1)
     iz(i1) = iz(i2)
     iz(i2) = inidum
   END DO
 END DO
 150 lll = il - 1
 ll  = inth - 1
 IF (idiag == 0) GO TO 170
 
!     PRINT OUT FULL SUMMARY FOR CURRENT NEIGHBORHOOD
 
 WRITE (nout,670) jreg,noreg,alam
 WRITE (nout,680)
 WRITE (nout,690)
 DO  i = 1,nrow
   k   = kk  + i
   izz = 2*iz(k) - 1
   j   = jj  + i
   l   = lll + izz
   l1  = l   + 1
   i1  = ihl1+ iz(k)
   z(i1) = 100.*z(i1)
   i2  = ll  + iz(k)
   STATUS(1) = accept(1)
   STATUS(2) = accept(2)
   IF (lz(i2)) GO TO 160
   STATUS(1) = reject(1)
   STATUS(2) = reject(2)
   160 WRITE (nout,700) i,iz(k),z(j),z(l),z(l1),z(i1),STATUS
 END DO
 
!     DECREMENT COUNTERS SO THAT ONLY ACCEPTABLE SOLUTIONS ARE RETAINED
 
 170 msave = nrow
 DO  i = 1,msave
   i2 = ll + i
   IF (lz(i2)) CYCLE
   nrow   = nrow - 1
   northo = northo - 1
   IF (nrow == 0) GO TO 450
 END DO
 nfound = nzero + nrow
 IF (nrow == msave) WRITE (nout,720) uim,msave
 IF (idiag == 0 .OR. nrow == msave) GO TO 200
 
!     PRINT OUT SUMMARY WITH REJECTED SOLUTIONS DELETED
 
 WRITE (nout,670) jreg,noreg,alam
 WRITE (nout,730)
 WRITE (nout,690)
 m  = 0
 DO  i = 1,msave
   k  = kk + i
   i2 = ll + iz(k)
   IF (.NOT.lz(i2)) CYCLE
   m  = m + 1
   izz= 2*iz(k) - 1
   j  = jj  + i
   l  = lll + izz
   l1 = l   + 1
   i1 = ihl1+ iz(k)
   WRITE (nout,700) m,iz(k),z(j),z(l),z(l1),z(i1),accept
 END DO
 200 m = msave - nrow
 IF (m > 0) WRITE (nout,740) uim,nrow,m
 
!     WRITE EIGENVALUES TO OUTPUT FILE
 
 CALL gopen (ilam(1),z(ibuf1),wrt)
 DO  i = 1,msave
   k  = kk + i
   i2 = ll + iz(k)
   IF (.NOT.lz(i2)) CYCLE
   izz = 2*iz(k) - 1
   l   = lll + izz
   lam1(1) = DBLE(z(l  ))
   lam1(2) = DBLE(z(l+1))
   CALL WRITE (ilam(1),lam1(1),4,1)
 END DO
 CALL CLOSE (ilam(1),eofnrw)
 IF (jreg < noreg .AND. nfound < nord) GO TO 214
 IF (nzero == 0) GO TO 214
 
!     IF THIS IS THE FINAL (BUT NOT THE FIRST) NEIGHBORHOOD, THEN
!     RE-WRITE THE EIGENVECTOR FILE PERTAINING TO ALL PRIOR
!     NEIGHBORHOODS (ELIMINATE LEFT-HAND VECTORS)
 
 209 IF (idiag /= 0) WRITE (nout,810) nzero,northo
 inidum = iscr(10)
 CALL OPEN  (*455,iscr(10),z(ibuf2),wrtrew)
 CALL CLOSE (iscr(10),rew)
 j = nord2
 IF (no b) j = 2*j
 inidum = iphi(1)
 CALL OPEN (*455,iphi(1),z(ibuf1),0)
 DO  i = 1,nzero
   CALL READ (*460,*211,iphi(1),z(iv2),nord8+10,0,n3)
   GO TO 470
   211 CALL gopen (iscr(10),z(ibuf2),wrt)
   CALL WRITE (iscr(10),z(iv2),j,1)
   CALL CLOSE (iscr(10),norew)
 END DO
 CALL CLOSE (iphi(1),norew)
 CALL OPEN  (*455,iphi(1),z(ibuf1),wrtrew)
 CALL CLOSE (iphi(1),rew)
 inidum = iscr(10)
 CALL OPEN  (*455,iscr(10),z(ibuf2),0)
 DO  i = 1,nzero
   CALL READ  (*460,*206,iscr(10),z(iv2),j+10,0,n3)
   GO TO 470
   206 CALL gopen (iphi(1),z(ibuf1),wrt)
   CALL WRITE (iphi(1),z(iv2),j,1)
   CALL CLOSE (iphi(1),eofnrw)
 END DO
 CALL CLOSE (iscr(10),norew)
 IF(nswp(2) < 0) GO TO 500
 
!     RECOVER PHYSICAL EIGENVECTORS, PRINT, AND WRITE TO OUTPUT FILE
 
 214 iprc = iprec + 2
 ii   = 1
 nn   = nord2
 incr = 1
 ia1  = ia - 1
 IF (qpr) WRITE (nout,750)
 ishft = nord2*iprec
 i1  = 0
 
!     ENTER LOOP
 
 DO  i = 1,msave
   k  = kk + i
   i2 = ll + iz(k)
   IF (.NOT. lz(i2)) CYCLE
   CALL gopen (iscr(7),z(ibuf2),rdrew)
   IF (nzero > 0) CALL skprec (iscr(7),nzero)
   DO  j = 1,nord8
     m = iv1x + j
     z(m) = 0.
   END DO
   
!     SET POINTER TO ALLMAT OUTPUT VECTOR
   
   ib = ia1 + 2*msave*(iz(k)-1)
   
!     CYCLE THRU ALL ORTHOGONAL VECTORS
   
   DO  j = 1,msave
     
!     NOTE.... Z(IV2) MAY BE LOADED DOUBLE-PRECISION....HIGHER DIGITS
!              ARE NOT USED
!     (HIGHER DIGITS MUST BE INCLUDED FOR THE D.P.MACHINES.  G.C/UNISYS)
     
     CALL unpack (*225,iscr(7),z(iv2))
     kr = ib + 2*j - 1
     ki = kr + 1
     DO  mm = 1,nord2,2
       mr = iv2 + (mm-1)*iprec
       mi = mr  + iprec
       jr = iv1x+ mm
       ji = jr  + 1
       IF (.NOT.dpmach) GO TO 216
       mrd = (mr+1)/2
       mid = mrd + 1
       
!     RECOVER RIGHT-HAND PHYSICAL EIGENVECTOR
       
       z(jr) = z(jr) + dz(mrd)*z(kr) - dz(mid)*z(ki)
       z(ji) = z(ji) + dz(mid)*z(kr) + dz(mrd)*z(ki)
       GO TO 217
       216 z(jr) = z(jr) + z(mr)*z(kr) - z(mi)*z(ki)
       z(ji) = z(ji) + z(mi)*z(kr) + z(mr)*z(ki)
       217 mr = mr + ishft
       mi = mr + iprec
       jr = jr + nord4
       ji = jr + 1
       IF (.NOT.dpmach) GO TO 218
       mrd = (mr+1)/2
       mid = mrd + 1
       
!     RECOVER LEFT-HAND PHYSICAL EIGENVECTOR
       
       z(jr) = z(jr) + dz(mrd)*z(kr) - dz(mid)*z(ki)
       z(ji) = z(ji) + dz(mid)*z(kr) + dz(mrd)*z(ki)
       CYCLE
       218 z(jr) = z(jr) + z(mr)*z(kr) - z(mi)*z(ki)
       z(ji) = z(ji) + z(mi)*z(kr) + z(mr)*z(ki)
     END DO
   END DO
   CALL CLOSE (iscr(7),eofnrw)
   IF (.NOT.qpr) GO TO 230
   i1  = i1 + 1
   izz = 2*iz(k) - 1
   l   = lll  + izz
   mm  = iv1x + nord8
   WRITE (nout,760) i1,iz(k),z(l),z(l+1),(z(j),j=iv1,mm)
   WRITE (nout,770)
   
!     EXPAND PHYSICAL EIGENVECTORS TO DOUBLE PRECISION FOR OUTPUT
   
   230 lim1 = iv1  + nord2
   lim2 = lim1 + nord4
   inidum = iv1x + nord4
   DO  j = 1,nord2
     ki  = lim1 - j
     mi  = 2*ki - iv1x
     mr  = mi - 1
     mrd = (mr+1)/2
     
!     EXPAND RIGHT-HAND VECTOR
     
     z(mi) = 0.
     z(mr) = z(ki)
     IF (dpmach) dz(mrd) = z(ki)
     ki  = lim2 - j
     mi  = 2*ki - inidum
     mr  = mi - 1
     mrd = (mr+1)/2
     
!     EXPAND LEFT -HAND VECTOR
     
     z(mi) = 0.
     z(mr) = z(ki)
     IF (dpmach) dz(mrd) = z(ki)
   END DO
   IF (.NOT.qpr) GO TO 250
   WRITE (nout,770)
   lim1 = iv1x + nord4
   WRITE (nout,780) (z(j),j=iv1,lim1)
   WRITE (nout,770)
   lim2 = lim1 + nord4
   lim1 = lim1 + 1
   WRITE (nout,780) (z(j),j=lim1,lim2)
   WRITE (nout,770)
   
!     PERFORM SPECIAL NORMALIZATION OF VECTORS FOR OUTPUT
   
   250 CALL cnorm1 (z(iv1),ikmb(2))
   IF (qpr) WRITE (nout,790)
   inidum = inidum + 1
   CALL cnorm1 (z(inidum),ikmb(2))
   IF (qpr) WRITE (nout,790)
   CALL gopen (iphi(1),z(ibuf1),wrt)
   IF (jreg < noreg .AND. nfound < nord) GO TO 260
   j = nord2
   IF (no b) j = 2*j
   CALL WRITE (iphi(1),z(iv1),j,1)
   CALL CLOSE (iphi(1),eofnrw)
   CYCLE
   
!     MUST USE NORD8 TO WRITE FULL RIGHT AND LEFT EIGENVECTORS
   
   260 CALL WRITE (iphi(1),z(iv1),nord8,1)
   CALL CLOSE (iphi(1),norew)
 END DO
 GO TO 500
 420 WRITE (nout,530) NAME
 GO TO 500
 430 WRITE (nout,540) m,NAME
 GO TO 500
 440 WRITE (nout,580) uwm
 IF(nzero > 0 .AND. jreg == noreg) nswp(2) = -1
 GO TO 500
 450 WRITE (nout,710) uwm,msave
 GO TO 500
 455 CALL mesage (-1,inidum,NAME)
 460 CALL mesage (-2,inidum,NAME)
 470 CALL mesage (-8,inidum,NAME)
 500 RETURN
 
 530 FORMAT (27H unexpected eof encountered,2X,2A4)
 540 FORMAT (22H unexpected word count,i5,2X,2A4)
 550 FORMAT (a29,' 3152', //5X,'SUBROUTINE ALLMAT OUTPUT EIGENVALUE',  &
     i4,' IS NULL.',//)
 560 FORMAT (a25,' 3153', //5X,'ATTEMPT TO NORMALIZE NULL VECTOR IN ',  &
     'SUBROUTINE ',a4,a2,'. NO ACTION TAKEN.',//)
 570 FORMAT (a25,' 3154', //5X,'SIZE OF REDUCED PROBLEM DECREMENTED ',  &
     'ONCE (NOW',i6,') DUE TO NULL ERROR ELEMENT.',//)
 580 FORMAT (a25,' 3155', //5X,'REDUCED PROBLEM HAS VANISHED. NO ',  &
     'ROOTS FOUND.',//)
 590 FORMAT (a25,' 3156', //5X,'SIZE OF REDUCED PROBLEM RESTORED TO',  &
     i8,' BECAUSE NEXT ERROR ELEMENT WAS ALSO NULL.', /5X,  &
     'ERROR ELEMENT SET = ',2E16.8,//)
 600 FORMAT (1H0,//7H cfeer4,//)
 610 FORMAT (1H ,i10,36H single PRECISION words of OPEN core,  &
     29H NOT used (SUBROUTINE cfeer4))
 620 FORMAT (4H row,i5,2(4X,2E16.8))
 630 FORMAT (4H row,i5,2(4X,2D16.8))
 640 FORMAT (1H0,26HREDUCED tridiagonal matrix, /(1H ,6E16.8))
 650 FORMAT (1H0,//30H theoretical eigenvalue errors,  &
     20X,18HD-sub-m-plus-one =,2E16.8,/)
 660 FORMAT (1H ,i5,e16.8,20X,2E16.8,10X,e16.8)
 670 FORMAT (1H1,27X,39H*****  f e e r  *****  (fast eigenvalue,  &
     27H extraction routine)  *****, //4X,  &
     24HSUMMARY for neighborhood,i3,3H of,i3,1H.,10X,  &
     21HNEIGHBORHOOD center =,2E16.8,/)
 680 FORMAT (4X,43HALL solutions found in current neighborhood,  &
     12H are listed.,/)
 690 FORMAT (4X,7X,8HSOLUTION,7X,8HORDER of,7X,8HDISTANCE,  &
     10X,10HEIGENVALUE,14X,11HTHEORETICAL, /4X,  &
     9X,6HNUMBER,5X,10HEXTRACTION,4X,11HFROM center,  &
     6X,4HREAL,9X,9HIMAGINARY,9X,5HERROR,12X,6HSTATUS,/)
 700 FORMAT (4X,i12,i15,1P,e18.8,1P,3E15.7,7X,2A4)
 710 FORMAT (a25,' 3163', //5X,'ALL',i6,' SOLUTIONS HAVE FAILED ',  &
     'ACCURACY TEST. NO ROOTS FOUND.',//)
 720 FORMAT (a29,' 3164',//5X,'ALL',i6,' SOLUTIONS ARE ACCEPTABLE.',//)
 730 FORMAT (4X,37HREJECTED solutions have been deleted.,/)
 740 FORMAT (a29,' 3165', //4X,i6,' SOLUTIONS HAVE BEEN ACCEPTED AND',  &
     i4,' SOLUTIONS HAVE BEEN REJECTED.',//)
 750 FORMAT (1H1,27X,39H*****  f e e r  *****  (fast eigenvalue,  &
     27H extraction routine)  *****,//  &
     42X,37HE i g e n v e c t o r   s u m m a r y,//1H , 32(4H----),2H--)
 760 FORMAT (1H ,8HSOLUTION,i4,8X,16HEXTRACTION order,i4,  &
     10X,10HEIGENVALUE,2X,1P,2E16.8, /(1H ,3(4X,1P,2E16.8)))
 770 FORMAT (3H --,32(4H----))
 780 FORMAT ((1H ,3(3X,2E16.8)))
 790 FORMAT (1H  ,12HAFTER cnorm1)
 800 FORMAT (1H  ,10X,5HRMS =,e16.8)
 810 FORMAT (1H  ,33HLEFT-hand eigenvectors eliminated,20X,2I8)
END SUBROUTINE cfeer4
