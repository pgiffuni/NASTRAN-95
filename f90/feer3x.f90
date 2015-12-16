SUBROUTINE feer3x
!                                                               T
!     FEER3 OBTAINS THE REDUCED TRIDIAGONAL MATRIX   (LI)*M*(LI)
!     WHERE M IS A SYMETRIC MATRIX AND L IS LOWER TRIANGULAR, AND (LI)
!     IS INVERSE OF L
 
!     THE TRANSFORMATION IS ALPHA = VT(L**(-1)M (L**-(1))TV
!     WHERE V IS A RECTANGULAR TRANSFORMATION.
 
!     LAST REVISED 11/91 BY G.CHAN/UNISYS, MAKE ROOM FOR NEW FBS METHOD
 
 INTEGER :: sysbuf    ,cndflg   ,mcbscl(7),sr5fle   ,  &
     sr6fle    ,sr7fle   ,sr8fle   ,sr9fle   ,  &
     sr10fl    ,srxfle   ,iz(1)    ,NAME(2)  , dashq     ,optn2
 DOUBLE PRECISION :: lambda    ,lmbda    ,dz(1)    ,dsq
 COMMON   /feercx/  ifkaa(7)  ,ifmaa(7) ,iflelm(7),iflvec(7),  &
     sr1fle    ,sr2fle   ,sr3fle   ,sr4fle   ,  &
     sr5fle    ,sr6fle   ,sr7fle   ,sr8fle   ,  &
     dmpfle    ,nord     ,xlmbda   ,neig     ,  &
     mord      ,ibk      ,critf    ,northo   , iflrva    ,iflrvc
 COMMON   /feerxx/  lambda    ,cndflg   ,iter     ,timed    ,  &
     l16       ,ioptf    ,epx      ,nochng   ,  &
     ind       ,lmbda    ,ifset    ,nzero    ,  &
     nonul     ,idiag    ,mrank    ,istart   , nzv5
 COMMON   /reigkr/  option    ,optn2
 COMMON   /TYPE  /  rc(2)     ,iwords(4)
 COMMON   /zzzzzz/  z(1)
 COMMON   /system/  sysbuf    ,io       ,systm(52),iprec    ,  &
     skip36(38),ksys94
 COMMON   /opinv /  mcblt(7)  ,mcbsma(7),mcbvec(7),mcbrm(7)
 COMMON   /unpakx/  iprc      ,ii       ,nn       ,incr
 COMMON   /packx /  itp1      ,itp2     ,iip      ,nnp      , incrp
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw
 EQUIVALENCE        (iz(1),z(1),dz(1))
 DATA      NAME  /  4HFEER,4H3   /      ,dashq    / 4H-q    /
 
!     SR5FLE CONTAINS THE TRIDIAGONAL ELEMENTS
!     SR6FLE CONTAINS THE G VECTORS
!     SR7FLE CONTAINS THE ORTHOGONAL VECTORS
!     SR8FLE CONTAINS THE CONDITIONED MAA OR KAAD MATRIX
!     SR9FLE CONTAINS MCBSMA DATA IN UNPACKED FORM = 309
!     SR10FL CONTAINS MCBLT  DATA IN UNPACKED FORM = 310
!                                              (OR = 308 IF IT IS FREE)
!     IFLVEC CONTAINS THE L OR C MATRIX FROM SDCOMP
!     IFLELM CONTAINS     KAA+ALPHA*MAA
!     IFLRVC CONTAINS THE RESTART AND/OR RIGID BODY VECTORS
 
 sr9fle = 309
 sr10fl = 308
 iprc   = mcblt(5)
 nwds   = iwords(iprc)
 nz     = korsz(z)
 CALL makmcb (mcbvec(1),sr7fle,nord,2,iprc)
 mcbvec(2) = 0
 mcbvec(6) = 0
 CALL makmcb (mcbrm(1) ,sr6fle,mord,2,iprc)
 mcbrm(2)  = 0
 mcbrm(6)  = 0
 mcbscl(1) = iflrvc
 CALL rdtrl (mcbscl(1))
 
!     INITIALIZE ALLOCATIONS
 
 ibuf1 = nz    - sysbuf
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 ibuf4 = ibuf3 - sysbuf
 iv1   = 1
 iv2   = iv1 + nord
 iv3   = iv2 + nord
 iv4   = iv3 + nord
 iv5   = iv4 + nord
 nzv5  = ibuf4 - iv5*nwds - 2
 ix2   = iv2 - 1
 iend  = nwds*(5*nord + 1) + 2
 icrq  = iend - ibuf4
 IF (icrq > 0) CALL mesage (-8,icrq,NAME)
 ifl   = mcblt(1)
 srxfle= sr8fle
 
!     CALL UNPSCR TO MOVE MCBSMA DATA INTO SR9FLE, AND MCBLT INTO SR10FL
!     (ORIGINAL MCBSMA AND MCBLT TRAILER WORDS 4,5,6,7 WILL BE CHANGED)
!     NZV5 IS THE AVAILABE SIZE OF THE WORKING SPACE FOR NEW FBS METHOD
!     USED IN FRSW/2, FRBK/2, FRMLT/D, AND FRMLTX/A ROUTINES
 
!     IF KSYS94 IS 10000 OR DIAG 41 IS ON, NEW FBS METHODS AND UNPSCR
!     ARE NOT USED
 
 IF (MOD(ksys94,100000)/10000 == 1) GO TO 10
 CALL sswtch (41,i)
 IF (i == 1) GO TO 10
 srxfle = sr9fle
 CALL unpscr (mcbsma,srxfle,z,ibuf2,ibuf1,nzv5,0,1)
 j = 2
 IF (ioptf == 1) j = 3
 CALL unpscr (mcblt,sr10fl,z,ibuf2,ibuf1,nzv5,0,j)
 nzv5 = nzv5 + 1
 ifl  = sr10fl
 
 10 CALL gopen (ifl,z(ibuf3),rdrew)
 CALL gopen (sr7fle,z(ibuf1),wrtrew)
 IF (northo == 0) GO TO 130
 
!     LOAD RESTART AND/OR RIGID BODY VECTORS
 
 CALL gopen (iflrvc,z(ibuf2),rdrew)
 incr  = 1
 incrp = 1
 itp1  = iprc
 itp2  = iprc
 
 DO  j = 1,northo
   ii  = 1
   nn  = nord
   CALL unpack (*110,iflrvc,dz(1))
   iip = ii
   nnp = nn
   IF (iprc  == 1) GO TO 60
   IF (ioptf == 0) GO TO 40
   dsq = 0.d0
   CALL frmltx (mcblt(1),dz(iv1),dz(iv2),dz(iv3))
   DO  ij = 1,nord
     dsq = dsq + dz(ix2+ij)**2
   END DO
   dsq = 1.d0/DSQRT(dsq)
   DO  ij = 1,nord
     dz(ij) = dsq*dz(ix2+ij)
   END DO
   40 IF (l16 == 0) GO TO 100
   CALL page2 (2)
   WRITE (io,50) iip,nnp,(dz(i),i=1,nord)
   50 FORMAT (10H orth vct ,2I5,  /(1X,8E16.8))
   GO TO 100
   60 IF (ioptf == 0) GO TO 90
   sq = 0.0
   CALL frmlta (mcblt(1),z(iv1),z(iv2),z(iv3))
   DO  ij = 1,nord
     sq = sq + z(ix2+ij)**2
   END DO
   sq = 1.0/SQRT(sq)
   DO  ij = 1,nord
     z(ij) = sq*z(ix2+ij)
   END DO
   90 IF (l16 == 0) GO TO 100
   CALL page2 (2)
   WRITE (io,50) iip,nnp,(z(i),i=1,nord)
   100 CALL pack (dz(1),sr7fle,mcbvec(1))
   110 CONTINUE
 END DO
 
 CALL CLOSE (iflrvc,norew)
 IF (l16 == 0) GO TO 130
 CALL page2 (1)
 WRITE  (io,120) northo,mcbvec
 120 FORMAT (5X,i5,16H orth vectors on,i5,5H FILE,5I5,i14)
 130 k = northo
 CALL CLOSE (sr7fle,norew)
 j = k
 nonul = 0
 iter  = 0
 CALL gopen (sr6fle,z(ibuf4),wrtrew)
 CALL CLOSE (sr6fle,norew)
 CALL gopen (srxfle,z(ibuf2),rdrew)
 CALL gopen (sr5fle,z(ibuf4),wrtrew)
 
!     GENERATE SEED VECTOR
 
 140 k = k + 1
 j = k
 ifn = 0
 
!     GENERATE SEED VECTOR FOR LANCZOS
 
 ss = 1.0
 IF (iprc == 1) GO TO 160
 DO  i = 1,nord
   ss =-ss
   j  = j + 1
   dsq = FLOAT(MOD(j,3)+1)/(3.0*FLOAT((MOD(j,13)+1)*(1+5*i/nord)))
   dz(ix2+i) = dsq*ss
 END DO
 IF (optn2 /= dashq) CALL fnxtvc (dz(iv1),dz(iv2),dz(iv3),  &
     dz(iv4),dz(iv5),z(ibuf1),ifn)
 GO TO 180
 
 160 DO  i = 1,nord
   ss =-ss
   j  = j + 1
   sq = FLOAT(MOD(j,3)+1)/(3.0*FLOAT((MOD(j,13)+1)*(1+5*i/nord)))
   z(ix2+i) = sq*ss
 END DO
 IF (optn2 /= dashq) CALL fnxtv  (z(iv1),z(iv2),z(iv3),z(iv4),  &
     z(iv5),z(ibuf1),ifn)
 IF (optn2 == dashq) CALL fnxtvd (z(iv1),z(iv2),z(iv3),z(iv4),  &
     z(iv5),z(ibuf1),ifn)
 
 180 IF (iter <= mord) GO TO 190
 mord = northo - nzero
 cndflg = 3
 GO TO 200
 
 190 IF (ifn < mord) GO TO 140
 200 CALL CLOSE (sr5fle,norew)
 CALL CLOSE (srxfle,rew)
 CALL CLOSE (ifl,rew)
 
!     IF NEW FBS METHOD IS USED, SR9FLE AND SR10FL FILES COULD BE VERY
!     BIG. MAKE SURE THEY ARE PHYSICALLY REDUCED TO ZERO SIZE. THIS IS
!     IMPORTANT FOR A COMPUTER SYSTEM WITH LIMITED DISC SPACE
 
 IF (ifl /= sr10fl) GO TO 210
 CALL gopen (sr9fle,z(ibuf2),wrtrew)
 CALL gopen (sr10fl,z(ibuf3),wrtrew)
 CALL CLOSE (sr9fle,rew)
 CALL CLOSE (sr10fl,rew)
 
 210 IF (l16 == 0) RETURN
 CALL page2 (1)
 i = ibuf4 - northo*nord*nwds - 2
 IF (i <  0) i = ibuf4 - iend
 WRITE  (io,220) i,NAME
 220 FORMAT (19H OPEN core NOT used,i10,2X,2A4)
 
 RETURN
END SUBROUTINE feer3x
