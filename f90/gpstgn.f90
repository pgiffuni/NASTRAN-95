SUBROUTINE gpstgn
     
!     THIS MODULE GENERATES THE GRID POINT SINGULARITY TABLE
!     BY EXAMINING THE SUBMATRICES ALONG THE LEADING DIAGONAL
!     OF THE INPUT STIFFNESS MATRIX
 
!     MODULE DMAP SEQUENCE
 
!     GPSTGEN  KGG,SIL/GPST $
 
 DIMENSION        k(3)  , mcb(7),isubnm(2)
 
 INTEGER :: sil   , gpst
 
 DOUBLE PRECISION :: b
!WKBI 8/94 SPR93026
 REAL :: bs(18)
 
 COMMON /gpstgx/  gpst  , igpst , isil  , nsing , ibuf2
 COMMON /gpstgy/  b(18)
!WKBI 8/94 SPR93026   COMMON /SYSTEM/  ISYSBF
 COMMON /system/ isysbf , nout  , dum(52),iprec
 COMMON /unpakx/  itypot, ii    , jj    , incr
 COMMON /zzzzzz/  iz(1)
!WKBI 8/94 SPR93026
 EQUIVALENCE      ( bs, b )
 
 DATA kgg, sil /101   , 102   /
 DATA isubnm   /4HGPST, 4HGN  /
 
 gpst  = 201
 igpst = 0
 nsing = 0
!WKBR 8/94 SPR93026      ITYPOT= 2
 itypot= iprec
 incr  = 1
 k(1)  = 1
 k(2)  = 1
 ibuf1 = korsz (iz) - isysbf - 2
 ibuf2 = ibuf1 - isysbf
 ifile = sil
 CALL OPEN (*120,sil,iz(ibuf1),0)
 CALL skprec (sil,1)
 mcb(1) = sil
 CALL rdtrl (mcb)
 luset = mcb(3)
 icore = luset + 1 - ibuf1
 IF (icore >= 0) GO TO 160
 CALL READ (*140,*10,sil,iz,ibuf1,0,npts)
 GO TO 160
 10 CALL CLOSE (sil,1)
 logic = 110
 IF (npts /= mcb(2)) GO TO 150
 iz(npts+1) = luset + 1
 
 ifile = kgg
 CALL OPEN (*120,kgg,iz(ibuf1),0)
 CALL skprec (kgg,1)
 mcb(1) = kgg
 CALL rdtrl (mcb)
 logic = 120
 IF (mcb(2) /= luset .OR. mcb(3) /= luset) GO TO 150
 
 DO  i = 1, npts
   ityp = 1
   isil = iz(i)
   isilnx = iz(i+1)
   IF (isilnx-isil == 1) ityp = 2
   iloop = 1
   ist   = 1
   ii = isil
   20 jj = ii + 2*(2 - ityp)
   DO  j = ii, jj
!WKBD 8/94 SPR93026      CALL UNPACK (*30,KGG,B(IST))
!WKBNB 8/94 SPR93026
     IF ( iprec == 1 ) CALL unpack (*30,kgg,bs(ist))
     IF ( iprec == 2 ) CALL unpack (*30,kgg,b(ist))
!WKBNE 8/94 SPR93026
     GO TO 50
     30 istx = ist + 2
!WKBI 8/94 SPR93026
     IF ( iprec == 1 ) GO TO 45
     DO  iii = ist, istx
       b(iii) = 0.0D0
     END DO
!WKBNB 8/94 SPR93026
     GO TO 50
     45 DO  iii = ist, istx
       bs(iii) = 0.0
     END DO
     50 ist = ist + 3
   END DO
   IF (ityp == 2) GO TO 70
   IF (iloop == 2) GO TO 90
   iloop = 2
   ii = ii + 3
   GO TO 20
!WKBD 8/94 SPR93026   70 IF (B(1).GT.0.0D0) GO TO 100
!WKBNB 8/94 SPR93026
   70 CONTINUE
   IF (iprec == 2 .AND. b(1) > 0.0D0) CYCLE
   IF (iprec == 1 .AND. bs(1) > 0.0  ) CYCLE
!WKBNE 8/94 SPR93026
   k(3) = isil
   IF (igpst == 1) GO TO 80
   igpst = 1
   CALL gopen (gpst,iz(ibuf2),1)
   80 nsing = nsing + 1
   CALL WRITE (gpst,k,3,0)
   CYCLE
!WKBD 8/94 SPR93026   90 CALL GPSTG
!WKBNB 8/94 SPR93026
   90 IF ( iprec == 1 ) CALL gpstgs
   IF ( iprec == 2 ) CALL gpstg
!WKBNE 8/94 SPR93026
 END DO
 IF (igpst == 0) GO TO 110
 CALL WRITE (gpst,0,0,1)
 CALL CLOSE (gpst,1)
 CALL makmcb (mcb,gpst,npts,luset,0)
 mcb(2) = nsing
 CALL wrttrl (mcb)
 110 CALL CLOSE (kgg,1)
 GO TO 170
 
!     ERROR MESSAGES
 
 120 n = -1
 130 CALL mesage (n,ifile,isubnm)
 140 n = -2
 GO TO 130
 150 n = -7
 GO TO 130
 160 n = -8
 ifile = icore
 GO TO 130
 
 170 RETURN
END SUBROUTINE gpstgn
