SUBROUTINE mtrxo (FILE,NAME,item,dumbuf,itest)
     
!     COPIES MATRIX ITEM OF SUBSTRUCTURE NAME FROM THE NASTRAN FILE
!     TO THE SOF
!     ITEST VALUES RETURNED ARE
!        1 - ITEM ALREADY EXISTS ON THE SOF
!        2 - THE ITEM WAS PESUDO WRITTEN
!        3 - NORMAL RETURN
!        4 - NAME DOES NOT EXIST
!        5 - ITEM IS NOT ONE OF THE ALLOWABLE MATIX ITEMS
!        6 - NASTRAN FILE HAS BEEN PURGED
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN OUT)                  :: item
 REAL, INTENT(IN OUT)                     :: dumbuf
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        lshift,orf,andf
 LOGICAL :: mdiup
 INTEGER :: nmsbr(2),buf(1), trail(7),oldbuf, blksiz,first,orf,andf
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /sof   / ditdum(6),io,iopbn,iolbn,iomode,ioptr,iosind,  &
     ioitcd,ioblk,sofdum(20),mdiup
 COMMON /system/ nbuf,nout
 COMMON /sys   / blksiz
 COMMON /zzzzzz/ nstrn
 EQUIVALENCE     (buf(1),nstrn)
 DATA    nmsbr / 4HMTRX, 4HO   /
 DATA    iwrt  / 2 /
 DATA    ifetch/ -2/
 DATA    idle  / 0 /
 
!     CHECK IF ITEM IS ONE OF THE FOLLOWING ALLOWABLE NAMES.
!     KMTX, MMTX, PVEC, POVE, UPRT, HORG, UVEC, QVEC, PAPP, POAP, LMTX
 
 CALL chkopn (nmsbr(1))
 itm = ittype(item)
 IF (itm  /= 1) GO TO 1030
 IF (FILE > 0) GO TO 20
 
!     THE MATRIX ITEM IS TO BE PESUDO WRITTEN
 
 itest = 2
 CALL sfetch (NAME(1),item,ifetch,itest)
 GO TO 100
 
!     CHECK IF NASTRAN FILE HAS BEEN PURGED
 
 20 trail(1) = FILE
 CALL rdtrl (trail)
 IF (trail(1) <= 0) GO TO 1020
 
!     OPEN ITEM TO WRITE AND FETCH FIRST BLOCK FOR SOF
 
 itest = 3
 CALL sfetch (NAME(1),item,ifetch,itest)
 IF (itest /= 3) GO TO 100
 
!     OPEN NASTRAN FILE
!     MAKE SURE BUFFER IS  DOUBLE WORD ALIGNED.
 
 idisp = locfx(buf(io-2)) - locfx(nstrn)
 IF (andf(idisp,1) /= 0) io = io + 1
 iopt = 0
 CALL OPEN (*1020,FILE,buf(io-2),iopt)
 
!     ADJUST SOF BUFFER TO COINCIDE WITH GINO BUFFER
 
 oldbuf = io
 
!IBMD 6/93 IF (MACH .GT. 2) GO TO 50
!IBMD 6/93 IF (BUF(IO-2) .EQ. FILE) GO TO 50
!IBMD 6/93 IO = IO + 1
!IBMD 6/93 IF (BUF(IO-2) .NE. FILE) GO TO 1000
 
!     BEGIN COPYING DATA TO SOF
 
!     FIRST CHECK IF CALL TO OPEN OBTAINED ONLY BLOCK IN FILE
 
 50 first = 1
 CALL rdblk (*70,FILE,first,left)
 first = 0
 
!     WRITE OUT BLOCK IN BUFFER TO SOF AND OBTAIN A FREE SOF BLOCK
 
 60 CALL sofio (iwrt,iopbn,buf(io-2))
 CALL getblk (iopbn,j)
 IF (j == -1) GO TO 120
 iopbn = j
 iolbn = iolbn + 1
 
!     OBTAIN A NEW BLOCK FROM THE GINO FILE
 
 CALL rdblk (*70,FILE,first,left)
 GO TO 60
 
!     THE LAST BUFFER OF THE GINO FILE HAS BEEN FOUND - DETERMINE
!     IF SUFFICIENT SPACE IN BUFFER REMAINS FOR TRAILER
 
 70 CONTINUE
 IF (left >= 6) GO TO 80
 
!     INSUFFICIENT SPACE - OBTAIN NEW SOF BLOCK
!     SET BLOCK NUMBER OF CURRENT BLOCK TO ZERO TO INDICATE TRAILER
!     IS STORED IN NEXT BLOCK
 
 buf(io+1) = 0
 CALL sofio (iwrt,iopbn,buf(io-2))
 CALL getblk (iopbn,j)
 IF (j == -1) GO TO 120
 iopbn = j
 iolbn = iolbn + 1
 
!     STORE TRAILER IN LAST SIX WORDS OF BLOCK
!     SET BLOCK NUMBER NEGATIVE TO INDICATE LAST BLOCK AND
!     WRITE OUT FINAL BLOCK TO SOF
 
 80 DO  i = 2,7
   buf(io+blksiz-7+i) = trail(i)
 END DO
 buf(io+1) = -iolbn
 CALL sofio (iwrt,iopbn,buf(io-2))
 
!     CLOSE FILE AND UPDATE MDI
 
 CALL CLOSE (FILE,1)
 CALL fmdi (iosind,imdi)
 buf(imdi+ioitcd) = orf(andf(ioblk,jhalf),lshift(iolbn,ihalf))
 mdiup = .true.
 
!     RETURN
 
 itest  = 3
 io     = oldbuf
 iomode = idle
 100 RETURN
 
!     THERE ARE NO MORE FREE BLOCKS ON THE SOF.
!     RETURN THE BLOCKS THAT HAVE BEEN USED SO FAR BY THE ITEM BEING
!     WRITTEN, CLOSE THE SOF AND ISSUE A FATAL MESSAGE
 
 120 CALL retblk (ioblk)
 CALL sofcls
 GO TO 1010
 
!     ERROR RETURNS
 
 
!     BUFFER ALIGNMENT ERROR
 
!IBMD 6/93 1000 CALL SOFCLS
!IBMD 6/93      CALL MESAGE (-8,0,NMSBR)
!IBMD 6/93      GO TO 100
 
!     NO MORE FREE BLOCKS ON THE SOF
 
 1010 WRITE  (nout,1011) ufm
 1011 FORMAT (a23,' 6223, THERE ARE NO MORE FREE BLOCKS AVAILABLE ON ',  &
     'THE SOF.')
 CALL mesage (-37,0,nmsbr)
 
!     GINO FILE PURGED
 
 1020 itest = 6
 GO TO 100
 
!     INVALID ITEM NAME
 
 1030 itest = 5
 GO TO 100
 
END SUBROUTINE mtrxo
