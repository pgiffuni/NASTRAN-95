SUBROUTINE feer4 (it)
     
!     FEER4 OBTAINS FROM THE REDUCED TRIDIAGONAL MATRIX THE EIGENVALUES
!     AND EIGENVECTORS
 
 
 INTEGER, INTENT(IN)                      :: it
 LOGICAL :: incore
 INTEGER :: sysbuf    ,cndflg   ,sr2fle    ,sr6fle   ,  &
     sr7fle    ,sr8fle   ,iz(1)     ,NAME(2)
!WKBNB NCL93007 11/94
 INTEGER :: sr4fle    ,sr5fle   ,rew       ,eofnrw  &
     ,                 wrtrew    ,rdrew    ,wrt
!WKBNE NCL93007 11/94
 DOUBLE PRECISION :: lambda    ,lmbda    ,dz(1)     ,b(2)     ,  &
     dsm       ,dsce
 DIMENSION         mcbc(7)   ,icr(2)   ,sb(2)
 COMMON   /machin/ mach
 COMMON   /feercx/ ifkaa(7)  ,ifmaa(7) ,iflelm(7),iflvec(7) ,  &
     sr1fle    ,sr2fle   ,sr3fle   ,sr4fle    ,  &
     sr5fle    ,sr6fle   ,sr7fle   ,sr8fle    ,  &
     dmpfle    ,nord     ,xlmbda   ,neig      ,  &
     mord      ,ibk      ,critf    ,northo    , iflrva    ,iflrvc
 COMMON   /feerxx/ lambda    ,cndflg   ,iter     ,timed     ,  &
     l16       ,ioptf    ,epx      ,errc      ,  &
     ind       ,lmbda    ,ifset    ,nzero     ,  &
     nonul     ,idiag    ,mrank    ,istart
 COMMON   /zzzzzz/ z(1)
 COMMON   /system/ sysbuf    ,io       ,systm(52),iprec
 COMMON   /opinv / mcblt(7)  ,mcbsma(7),mcbvec(7),mcbrm(7)
 COMMON   /unpakx/ iprc      ,ii       ,nn       ,incr
 COMMON   /packx / itp1      ,itp2     ,iip      ,nnp       , incrp
 COMMON   /names / rd        ,rdrew    ,wrt      ,wrtrew    ,  &
     rew       ,norew    ,eofnrw
 EQUIVALENCE       (iz(1),z(1),dz(1)), (sb(1),b(1)), (dsce,sce)
 DATA      NAME  / 4HFEER,4H4     /,    icr   / 4HPASS,4HFAIL /
 
!     SR4FLE CONTAINS THE EIGENVECTORS OF THE REDUCED PROBLEM
!     SR5FLE CONTAINS THE TRIDIAGONAL ELEMENTS AND SCRATCH IN FQRWV
!     SR6FLE CONTAINS THE G VECTORS
!     SR7FLE CONTAINS THE ORTHOGONAL VECTORS
 
 CALL sswtch (26,l26)
 mdim = mord + 1
 dsm  = 10.0D+0**(-2*it/3)
 sm   = dsm
 iprc = mcbrm(5)
 nz   = korsz(z)
 CALL makmcb (mcbc(1),sr4fle,mdim,2,iprc)
 mcbc(2) = 0
 mcbc(6) = 0
 m    = 0
 
!     INITIALIZE ALLOCATIONS
 
 ibuf1 = nz    - sysbuf
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 iv1   = 1
 iv2   = iv1 + mdim
 iv3   = iv2 + mdim
 iv4   = iv3 + mdim
 iv5   = iv4 + mdim
 iv6   = iv5 + mdim
 iv7   = iv6 + mdim
 iv8   = iv7 + mdim
 iv9   = iv8 + mdim
 ix3   = iv3 - 1
 ix4   = iv4 - 1
 iend  = iprc*(8*mdim+1) + mdim
 IF (iend > ibuf3) CALL mesage (-8,iend-ibuf3,NAME)
 CALL gopen (sr5fle,z(ibuf2),rdrew)
 IF (iprc == 2) dz(iv4+mord) = errc
 IF (iprc == 1)  z(iv4+mord) = errc
 nw = iprc*2
 DO  i = 1,mord
   CALL READ (*420,*430,sr5fle,b(1),nw,1,m)
   IF (iprc == 1) GO TO 5
   dz(ix3+i) = b(1)
   dz(ix4+i) = b(2)
   CYCLE
   5 z(ix3+i) = sb(1)
   z(ix4+i) = sb(2)
 END DO
 CALL CLOSE (sr5fle,rew)
 CALL gopen (sr4fle,z(ibuf2),wrtrew)
 IF (iprc == 1) GO TO 12
 CALL fqrwv (mord,dz(iv1),dz(iv2),dz(iv3),dz(iv4),dz(iv5),dz(iv6),  &
     dz(iv7),dz(iv8),dz(iv9),z(ibuf1),sr5fle,mcbc(1))
!                                                              SR4FLE
 GO TO 15
 12 CALL fqrw (mord,z(iv1),z(iv2),z(iv3),z(iv4),z(iv5),z(iv6),  &
     z(iv7),z(iv8),z(iv9),z(ibuf1),sr5fle,mcbc(1))
!                                                          SR4FLE
 15 CALL CLOSE (sr4fle,norew)
 
!     RECONFIGURE VECTOR INDEX TO OBTAIN PHYSICAL EIGENVECTORS
 
 ix1 = iv1 - 1
 ix2 = iv2 - 1
 ix3 = iv3 - 1
 ix4 = iv4 - 1
 ix5 = ix4 + nord
 isrv = mcbrm(1)
 iflvec(1) = iflrvc
 iflelm(1) = iflrva
 IF (nzero /= 0) GO TO 20
 
!     PREPARE FILES WHEN NO RESTART AND/OR RIGID BODY VECTORS
 
 iflvec(2) = 0
 iflvec(6) = 0
 CALL gopen (iflrvc,z(ibuf3),wrtrew)
 CALL CLOSE (iflrvc,norew)
 CALL gopen (iflrva,z(ibuf3),wrtrew)
 CALL CLOSE (iflrva,norew)
 20 itp1 = iprc
 itp2 = 1
 incrp= 1
 ii   = 1
 CALL gopen (iflrva,z(ibuf1),wrt)
 mred = 0
 mflg = 1
 DO  m = 1,mord
   IF (iprc == 1) GO TO 22
   dsce = 1.0D+0/dz(ix1+m) - lambda
   IF (l16 == 0) GO TO 24
   erf  = 0.0D+0
   IF (DABS(dsce) > dsm) erf = 100.d0*dz(ix2+m)/DABS(1.d0-dz(ix1+m)*lambda)
   dz(ix2+m) = dsce
   GO TO 23
   22 sce = 1.0/z(ix1+m) - lambda
   IF (l16 == 0) GO TO 24
   erf = 0.0D+0
   IF (ABS(sce) > sm) erf = 100.0D+0*z(ix2+m)/DABS(1.0D+0-z(ix1+m)*lambda)
   z(ix2+m) = sce
   23 IF (erf  > critf) mflg = 2
   24 IF (mflg ==     2) GO TO 25
   mred = mred + 1
   CALL WRITE (iflrva,dsce,iprec,1)
   25 IF (l16 == 0) CYCLE
   CALL page2 (1)
   IF (iprc == 2) WRITE (io,26) m,dsce,erf,icr(mflg)
   IF (iprc == 1) WRITE (io,26) m, sce,erf,icr(mflg)
   26 FORMAT (10X,'PHYSICAL EIGENVALUE',i5,1P,e16.8,  &
       '  THEOR ERROR ',e16.8,'  PERCENT',5X,a4)
 END DO
 CALL CLOSE (iflrva,eofnrw)
 IF (mord == 0) RETURN
 
 CALL gopen (isrv,z(ibuf1),rdrew)
 CALL gopen (sr4fle,z(ibuf2),rdrew)
 CALL gopen (iflrvc,z(ibuf3),wrt)
!WKBNB NCL93007 11/94
 incore = .false.
 CALL sswtch ( 43, l43 )
 IF ( l43 /= 0 ) GO TO 42
 ivw    = ix5 + nord + 1
 icreq  = nord*mord*iprc
 icavl  = ibuf3 - ivw - 1
 IF ( icavl > icreq ) incore = .true.
 IF ( .NOT. incore ) GO TO 42
 nn = nord
 DO  i = 1, mord
   ivr = ivw + (i-1)*nord
   IF ( iprc == 1 ) CALL unpack ( *41, isrv,  z(ivr+1) )
   IF ( iprc == 2 ) CALL unpack ( *41, isrv, dz(ivr+1) )
 END DO
 42 CONTINUE
!WKBNE NCL93007 11/94
 
!     IF DIAG 26 IS OFF, LIMIT EIGENSOLUTIONS TO NUMBER REQUESTED
 
 IF (mred >= neig .AND. l26 /= 0) mred = neig
 IF (iprc == 1) GO TO 200
 DO  m = 1,mred
   DO   l = 1,nord
     dz(ix5+l) = 0.0D+0
   END DO
   nn = mord
   CALL unpack (*75,sr4fle,dz(iv3))
   nn = nord
!WKBI NCL93007 11/94
   IF ( incore ) GO TO 72
   DO  i = 1,mord
     CALL unpack (*100,isrv,dz(iv4))
     DO  j = 1,nord
       dz(ix5+j) = dz(ix5+j) + dz(ix4+j)*dz(ix3+i)
     END DO
   END DO
!WKBNB NCL93007 11/94
   GO TO 73
   72 CONTINUE
   DO  i = 1, mord
     ivr = ivw + (i-1)*nord
     DO  j = 1, nord
       dz(ix5+j) = dz(ix5+j) + dz(ivr+j)*dz(ix3+i)
     END DO
   END DO
   73 CONTINUE
!WKBNE NCL93007 11/94
   75 CONTINUE
   IF (ioptf == 0) GO TO 90
   dsce = 1.0D+0/DSQRT(DABS(dz(ix1+m)))
   DO  l = 1,nord
     dz(ix5+l) = dsce*dz(ix5+l)
   END DO
   90 CONTINUE
   iip = 1
   nnp = nord
   CALL pack (dz(ix5+1),iflrvc,iflvec(1))
!WKBI NCL93007 11/94
   IF ( incore ) CYCLE
   CALL REWIND (mcbrm)
   CALL skprec (mcbrm,1)
 END DO
 GO TO 400
 200 DO  m = 1,mred
   DO  l = 1,nord
     z(ix5+l) = 0.0
   END DO
   nn = nord
   CALL unpack (*275,sr4fle,z(iv3))
   nn = nord
!WKBI NCL93007 11/94
   IF ( incore ) GO TO 272
   DO  i = 1,mord
     CALL unpack (*300,isrv,z(iv4))
     DO  j = 1,nord
       z(ix5+j) = z(ix5+j) + z(ix4+j)*z(ix3+i)
     END DO
   END DO
!WKBNB NCL93007 11/94
   GO TO 273
   272 CONTINUE
   DO  i = 1, mord
     ivr = ivw + (i-1)*nord
     DO  j = 1, nord
       z(ix5+j) = z(ix5+j) + z(ivr+j)*z(ix3+i)
     END DO
   END DO
   273 CONTINUE
!WKBNE NCL93007 11/94
   275 CONTINUE
   IF (ioptf == 0) GO TO 290
   sce = 1.0/SQRT(ABS(z(ix1+m)))
   DO  l = 1,nord
     z(ix5+l) = sce*z(ix5+l)
   END DO
   290 CONTINUE
   iip = 1
   nnp = nord
   CALL pack (z(ix5+1),iflrvc,iflvec(1))
!WKBI NCL93007 11/94
   IF ( incore ) CYCLE
   CALL REWIND (mcbrm)
   CALL skprec (mcbrm,1)
 END DO
 
 400 CALL CLOSE (iflrvc,eofnrw)
 CALL CLOSE (isrv,rew)
 CALL CLOSE (sr4fle,rew)
 mord = mred
 GO TO 500
 420 ier = 2
 GO TO 440
 430 ier = 3
 440 cndflg = 4
 CALL mesage (ier,sr5fle,NAME)
 500 iopn = ibuf3 - iend
 IF (l16 == 1) WRITE (io,510) iopn,NAME
 510 FORMAT ('  OPEN CORE NOT USED',i10,2X,2A4)
 RETURN
END SUBROUTINE feer4
