SUBROUTINE feer3
!                                                               T
!     FEER3 OBTAINS THE REDUCED TRIDIAGONAL MATRIX   (LI)*M*(LI)
!     WHERE M IS A SYMETRIC MATRIX AND L IS LOWER TRIANGULAR, AND (LI)
!     IS INVERSE OF L
 
!     THE TRANSFORMATION IS ALPHA = VT(L**(-1)M (L**-(1))TV
!     WHERE V IS A RECTANGULAR TRANSFORMATION.
 
!  Comments to follow refer to updates made 11/94.
!  This is a new version of FEER3.  The old version has been renamed FEER3X.
!  Diag 43 may be used to force the use of the old version.  The new version
!  uses all of available open core for storage of the orthogonal vectors,
!  the lower triangular matrix from SDCOMP, and the SMA matrix.  If
!  insufficient memory is available, only part of the lower triangular
 
 INTEGER :: sysbuf    ,cndflg   ,mcbscl(7),sr5fle   ,  &
     sr6fle    ,sr7fle   ,sr8fle   , iz(1)     ,NAME(2)  ,rew      ,wrtrew   ,  &
     optn2    ,rdrew    ,smapos
!     INTEGER            DASHQ
 DOUBLE PRECISION :: lambda    ,lmbda    ,dz(1)    ,dsq
 COMMON   /feercx/  ifkaa(7)  ,ifmaa(7) ,iflelm(7),iflvec(7),  &
     sr1fle    ,sr2fle   ,sr3fle   ,sr4fle   ,  &
     sr5fle    ,sr6fle   ,sr7fle   ,sr8fle   ,  &
     dmpfle    ,nord     ,xlmbda   ,neig     ,  &
     mord      ,ibk      ,critf    ,northo   , iflrva    ,iflrvc
 COMMON   /feerxx/  lambda    ,cndflg   ,iter     ,timed    ,  &
     l16       ,ioptf    ,epx      ,nochng   ,  &
     ind       ,lmbda    ,ifset    ,nzero    ,  &
     nonul     ,idiag    ,mrank    ,istart
 
!  NIDSMA = IN-MEMORY INDEX FOR COLUMN DATA OF SMA MATRIX
!  NIDLT  = IN-MEMORY INDEX FOR LOWER TRIANGULAR MATRIX
!  NIDORV = IN-MEMORY INDEX FOR ORTHOGONAL VECTORS
!  NLTLI  = INDEX OF LAST STRING OF LOWER TRIANGULAR MATRIX HELD IN MEMORY
!  NSMALI = INDEX OF LAST STRING OF SMA MATRIX HELD IN MEMORY
!  IBFSMA = IN-MEMORY INDEX FOR BUFFER FOR OPENING SMA MATRIX
!  IBMLT  = IN-MEMORY INDEX FOR BUFFER FOR OPENING LOWER TRIANGULAR MATRIX
!  IBFORV = IN-MEMORY INDEX FOR BUFFER FOR ORTHOGONAL VECTORS
!  SMAPOS = POSITION OF RECORD FOLLOWING LAST RECORD READ INTO MEMORY
!           AND THE LAST RECORD OF MATRIX SMA (SEE SUBROUTINE DSCPOS)
!  LTPOS  = POSITION OF RECORD FOLLOWING LAST RECORD READ INTO MEMORY
!           AND THE LAST RECORD OF THE LOWER TRIANGULAR MATRIX
 
 COMMON   /feerim/  nidsma    ,nidlt    ,nidorv   ,nltli    ,  &
     nsmali    ,ibfsma   ,ibflt    , ibforv    ,smapos(7),ltpos(7)
 COMMON   /reigkr/  option    ,optn2
 COMMON   /TYPE  /  rc(2)     ,iwords(4)
 COMMON   /zzzzzz/  z(1)
 COMMON   /system/  sysbuf    ,nout     ,systm(52),iprec    ,  &
     skip36(38),ksys94
 COMMON   /opinv /  mcblt(7)  ,mcbsma(7),mcbvec(7),mcbrm(7)
 COMMON   /unpakx/  iprc      ,ii       ,nn       ,incr
 COMMON   /packx /  itp1      ,itp2     ,iip      ,nnp      , incrp
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw
 EQUIVALENCE        (iz(1),z(1),dz(1))
 DATA      NAME  /  4HFEER,4H3   /
!     DATA      DASHQ / 4H-Q    /
 
!     SR5FLE CONTAINS THE TRIDIAGONAL ELEMENTS
!     SR6FLE CONTAINS THE G VECTORS
!     SR7FLE CONTAINS THE ORTHOGONAL VECTORS
!     SR8FLE CONTAINS THE CONDITIONED MAA OR KAAD MATRIX
!     IFLVEC CONTAINS THE L OR C MATRIX FROM SDCOMP
!     IFLELM CONTAINS     KAA+ALPHA*MAA
!     IFLRVC CONTAINS THE RESTART AND/OR RIGID BODY VECTORS
 
 CALL sswtch ( 43, l43 )
 IF ( l43 == 0 ) GO TO 1
 CALL feer3x
 GO TO 7777
 1     CONTINUE
 iprc      = mcblt(5)
 nwds      = iwords(iprc)
 nz        = korsz(z)
 CALL makmcb (mcbvec(1),sr7fle,nord,2,iprc)
 mcbvec(2) = 0
 mcbvec(6) = 0
 CALL makmcb (mcbrm(1) ,sr6fle,mord,2,iprc)
 mcbrm(2)  = 0
 mcbrm(6)  = 0
 mcbscl(1) = iflrvc
 CALL rdtrl (mcbscl(1))
 
!     INITIALIZE ALLOCATIONS
 
 ibuf1  = nz    - sysbuf
 ibuf2  = ibuf1 - sysbuf
 ibuf3  = ibuf2 - sysbuf
 ibuf4  = ibuf3 - sysbuf
 ibforv = ibuf1
 ibflt  = ibuf3
 ibfsma = ibuf2
 iv1    = 1
 iv2    = iv1 + nord
 iv2m1  = iv2 - 1
 iv3    = iv2 + nord
 iv4    = iv3 + nord
 iv5    = iv4 + nord
 iend   = nwds*(5*nord + 1) + 2
 mavail = iend - ibuf4
 IF (mavail > 0) CALL mesage (-8,mavail,NAME)
 
! COMPUTE THE MEMORY REQUIREMENT FOR ORTHOGONAL VECTORS
 
 memort = nord * ( mord+northo ) * iprc
 
! COMPUTE THE MEMORY REQUIREMENT FOR THE LOWER TRIANGULAR MATRIX
 
 CALL dssize ( mcblt, ncols, nterms, nstrgs, nwdtrm )
 memlt  = nterms*nwdtrm + nstrgs*4
 
! COMPUTE THE MEMORY REQUIREMENT FOR THE SMA MATRIX
 
 CALL dssize ( mcbsma, ncols, nterms, nstrgs, nwdtrm )
 memsma = nterms*nwdtrm + nstrgs*4
 IF ( l16 == 0 ) GO TO 2
 minnee = iend + 4*sysbuf
 memtot = memort + memlt + memsma + minnee
 WRITE ( nout, 901 ) minnee, memort, memsma, memlt, memtot, nz
 901   FORMAT(' FEER EIGENVALUE EXTRACTION NFORMATION'  &
     ,/, 5X,' THE FOLLOWING GIVES OPEN CORE REQUIREMENTS FOR KEEPING'  &
     ,/, 5X,' VARIOUS MATRICES AND VECTORS IN CORE FOR THE FEER'  &
     ,/, 5X,' EIGENVALUE EXTRACTION METHOD'  &
     ,/,10X,' MINIMUM NUMBER OF WORDS NEEDED IN OPEN CORE    =',i10  &
     ,/,10X,' NUMBER OF WORDS FOR ORTHOGONAL VECTORS         =',i10  &
     ,/,10X,' NUMBER OF WORDS FOR SMA MATRIX                 =',i10  &
     ,/,10X,' NUMBER OF WORDS FOR LOWER TRIANGULAR MATRIX    =',i10  &
     ,/,10X,' TOTAL NUMBER OF WORDS NEEDED TO ELIMINATE I/O  =',i10  &
     ,/,10X,' WORDS FOR OPEN CORE SPECIFIED IN THIS RUN      =',i10 )
 2     CONTINUE
! CHECK TO SEE IF MEMORY AVAILABLE FOR ORTHOGONAL VECTORS
 nidorv = 0
 itest  = iend + memort
 IF ( itest > ibuf4 ) GO TO 3
 nidorv = iend
 nidorv = ( nidorv/2 ) * 2 + 1
 iend   = iend + memort
 3     CONTINUE
! CHECK TO SEE IF MEMORY AVAILABLE FOR SMA MATRIX
 irmem  = ibuf4 - iend
 IF ( irmem <= 10 ) GO TO 4
 nidsma = iend
 nidsma = (nidsma/2) * 2  + 1
 memsma = memsma
 memsma = MIN0 ( memsma, irmem )
 iend   = iend + memsma
 GO TO 5
 4     CONTINUE
 nidsma = 0
 memsma = 0
 5     CONTINUE
! CHECK TO SEE IF MEMORY AVAILABLE FOR LOWER TRIANGULAR MATRIX
 irmem  = ibuf4 - iend
 IF ( irmem <= 10 ) GO TO 6
 nidlt  = iend
 nidlt  = (nidlt/2) * 2 + 1
 memlt  = memlt
 memlt  = MIN0 ( memlt, irmem )
 iend   = iend + memlt
 GO TO 7
 6     CONTINUE
 nidlt  = 0
 memlt  = 0
 7     CONTINUE
 ltpos ( 4 ) = -1
 smapos( 4 ) = -1
!      PRINT *,' FEER3, CALLING FERRDM,NIDSMA,NIDLT=',NIDSMA,NIDLT
 IF ( nidsma == 0 ) GO TO 11
 CALL ferrdm ( mcbsma,nidsma,memsma,ibfsma,nsmali,smapos)
!      PRINT *,' RETURN FROM FERRDM,MEMSMA,NSMALI=',MEMSMA,NSMALI
!      PRINT *,' SMAPOS=',SMAPOS
 11    IF ( nidlt  == 0 ) GO TO 12
 CALL ferrdm ( mcblt ,nidlt ,memlt ,ibflt ,nltli ,ltpos )
!      PRINT *,' RETURN FROM FERRDM,MEMLT,NLTLI=',MEMLT,NLTLI
!      PRINT *,' LTPOS=',LTPOS
 12    CONTINUE
 IF ( l16 == 0 ) GO TO 8
 WRITE ( nout, 902 ) 'SMA',smapos(1)
 WRITE ( nout, 902 ) 'LT ',ltpos(1)
 902   FORMAT(10X,' LAST COLUMN OF ',a3,' MATRIX IN MEMORY IS ',i4 )
!      PRINT *,' SMAPOS=',SMAPOS
!      PRINT *,' LTPOS =',LTPOS
 8     CONTINUE
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
!      PRINT *,' FERR3 CALLING FRMLTX'
   CALL frmltx (mcblt(1),dz(iv1),dz(iv2),dz(iv3))
   DO  ij = 1,nord
     dsq = dsq + dz(iv2m1+ij)**2
   END DO
   dsq = 1.d0/DSQRT(dsq)
   DO  ij = 1,nord
     dz(ij) = dsq*dz(iv2m1+ij)
   END DO
   40 CONTINUE
   GO TO 100
   60 IF (ioptf == 0) GO TO 90
   sq = 0.0
!      PRINT *,' FEER3 CALLING FRMLTA'
   CALL frmlta (mcblt(1),z(iv1),z(iv2),z(iv3))
   DO  ij = 1,nord
     sq = sq + z(iv2m1+ij)**2
   END DO
   sq = 1.0/SQRT(sq)
   DO  ij = 1,nord
     z(ij) = sq*z(iv2m1+ij)
   END DO
   90 CONTINUE
   100 CALL pack (dz(1),sr7fle,mcbvec(1))
   110 CONTINUE
 END DO
 CALL CLOSE (iflrvc,norew)
 130 k = northo
 CALL CLOSE (sr7fle,norew)
 j = k
 nonul = 0
 iter  = 0
!      PRINT *,' FEER3,SR7FLE,IFLRVC,SR6FLE=',SR7FLE,IFLRVC,SR6FLE
!      PRINT *,' FEER3,SR6FLE,SR8FLE,SR5FLE=',SR6FLE,SR8FLE,SR5FLE
!      PRINT *,' FEER3,MCBSMA,MCBLT,MCBVEC=',MCBSMA(1),MCBLT(1),MCBSMA(1)
 CALL gopen (sr6fle,z(ibuf4) ,wrtrew)
 CALL CLOSE (sr6fle,norew)
 IF ( sr8fle == mcbsma(1) ) GO TO 131
!      PRINT *,' PROBLEM IN FEER3, SR8FLE NE MCBSMA =',SR8FLE,MCBSMA(1)
 STOP
 131 CONTINUE
!      CALL GOPEN (SR8FLE,Z(IBUF2) ,RDREW )
 CALL gopen (sr5fle,z(ibuf4) ,wrtrew)
 CALL gopen (mcbsma,z(ibfsma),rdrew )
 CALL gopen (mcblt ,z(ibflt ),rdrew )
 
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
   dz(iv2m1+i) = dsq*ss
 END DO
!      PRINT *,' FEER3 CALLING FERXTD'
 CALL ferxtd (dz(iv1), dz(iv2), dz(iv3)  &
     ,            dz(iv4), dz(iv5), z(ibuf1), ifn )
 GO TO 180
 160 DO  i = 1,nord
   ss =-ss
   j  = j + 1
   sq = FLOAT(MOD(j,3)+1)/(3.0*FLOAT((MOD(j,13)+1)*(1+5*i/nord)))
   z(iv2m1+i) = sq*ss
 END DO
!      IF (OPTN2 .EQ. DASHQ) GO TO 175
 CALL ferxts ( z(iv1), z(iv2)  , z(iv3), z(iv4 )  &
     ,             z(iv5), z(ibuf1), ifn)
 GO TO 180
!  175 CALL FERXTQ ( Z(IV1), Z(IV2)  , Z(IV3), Z(IV4 )
!     1,              Z(IV5), Z(IBUF1), IFN)
 180 IF (iter <= mord) GO TO 190
 mord = northo - nzero
 cndflg = 3
 GO TO 200
 
 190 IF (ifn < mord) GO TO 140
 200 CALL CLOSE (sr5fle,norew)
 CALL CLOSE (sr8fle,rew)
 CALL CLOSE (mcblt ,rew)
 7777 CONTINUE
 RETURN
END SUBROUTINE feer3
