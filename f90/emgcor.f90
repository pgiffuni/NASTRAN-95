SUBROUTINE emgcor (buf)
     
!     CORE ALLOCATION AND PARAMETER INITIALIZATION FOR MAIN -EMG-
!     PROCESSOR -EMGPRO-.
 
 
 INTEGER, INTENT(OUT)                     :: buf(8)
 LOGICAL :: anycon, error, heat
 INTEGER :: z, sysbuf, outpt, subr(2), TYPE(3), bufs,  &
     buf1, buf2, rd, wrt, rdrew, wrtrew, cls, precis,  &
     clsrew, est, cstm, dit, geom2, NAME(2), eor, flags, scr4
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm, uwm, uim, sfm, swm
 COMMON /BLANK / nokmb(3), dummy(13), volume, surfac
 COMMON /system/ ksystm(65)
 COMMON /names / rd, rdrew, wrt, wrtrew, clsrew, cls
 COMMON /emgfil/ est, cstm, mpt, dit, geom2, kmbmat(3), kmbdic(3)
 COMMON /emgprm/ icore, jcore, ncore, icstm, ncstm, imat, nmat,  &
     ihmat, nhmat, idit, ndit, icong, ncong, lcong,  &
     anycon, flags(3), precis, error, heat,  &
     icmbar, lcstm, lmat, lhmat, kflags(3), l38
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm(1), sysbuf)
 EQUIVALENCE     (ksystm(2), outpt )
 DATA    TYPE  / 4HSTIF,4HMASS,4HDAMP/
 DATA    scr4  / 304   /
 DATA    subr  / 4HEMGC,4HOR  /,  eor/ 1 /
 
 IF (l38 == 0) WRITE (outpt,5) uim
 5 FORMAT (a29,' 238, TURN DIAG 38 ON FOR ADDITIONAL ELEMENT ',  &
     'PROCESSING INFORMATION',/)
 
!     DETERMINATION OF FUNCTIONS TO BE PERFORMED AND RESULTANT NUMBER
!     OF BUFFERS NEEDED.
 
 bufs = 1
 DO  i = 1,3
   flags(i)  = 0
   kflags(i) = 0
   IF (nokmb(i) == -1) CYCLE
   flags(i) = -1
   bufs = bufs + 2
 END DO
 IF (volume > 0.0 .OR. surfac > 0.0) bufs = bufs + 1
 
!     ALLOCATE BUFFERS
 
 n = ncore
 DO  i = 1,bufs
   buf(i) = n - sysbuf - 2
   n = buf(i)
 END DO
 ncore = n - 1
 IF (ncore < jcore) CALL mesage (-8,jcore-ncore,subr)
 
!  OPEN REQUIRED DATA BLOCKS.
 
 buf1 = buf(1)
 CALL OPEN (*60,est,z(buf1),rdrew)
 CALL skprec (est,1)
 ibuf = 1
 
!     K, M, OR B MATRIX DATA BLOCKS
 
 DO  i = 1,3
   IF (flags(i) == 0) CYCLE
   buf1 = buf(ibuf+1)
   buf2 = buf(ibuf+2)
   CALL OPEN (*30,kmbmat(i),z(buf1),wrtrew)
   CALL OPEN (*30,kmbdic(i),z(buf2),wrtrew)
   CALL fname (kmbmat(i),NAME)
   CALL WRITE (kmbmat(i),NAME,2,eor)
   CALL fname (kmbdic(i),NAME)
   CALL WRITE (kmbdic(i),NAME,2,eor)
   ibuf = ibuf + 2
   kflags(i) = 1
   CYCLE
   
!     FILE REQUIRED IS MISSING
   
   30 flags(i) = 0
   CALL page2 (2)
   WRITE  (outpt,40) uwm,kmbmat(i),kmbdic(i),TYPE(i)
   40 FORMAT (a25,' 3103, EMGCOR OF EMG MODULE FINDS EITHER OF DATA ',  &
       'BLOCKS ',i4,4H OR ,i4,' ABSENT AND THUS,', /5X,a4,  &
       ' MATRIX WILL NOT BE FORMED.')
 END DO
 
!     IF VOLUME OR SURFACE COMPUTATION IS REQUESTED BY USER FOR THE 2-D
!     AND 3-D ELEMENTS, OPEN SCR4 FILE. (ONLY TO BE CLOSED BY EMGFIN)
 
 IF (volume <= 0.0 .AND. surfac <= 0.0) GO TO 55
 ibuf = ibuf + 1
 buf1 = buf(ibuf)
 CALL OPEN (*80,scr4,z(buf1),wrtrew)
 
!     ALL FILES READY TO GO.
 
 55 ncore = buf(ibuf) - 1
 RETURN
 
!     EST MISSING
 
 60 CALL page2 (2)
 WRITE  (outpt,70) swm,est
 70 FORMAT (a27,' 3104, EMGCOR FINDS EST (ASSUMED DATA BLOCK',i5,  &
     ') MISSING.  EMG MODULE COMPUTATIONS LIMITED.')
 flags(1) = 0
 flags(2) = 0
 flags(3) = 0
 RETURN
 
 80 CALL mesage (-1,scr4,subr)
 RETURN
END SUBROUTINE emgcor
