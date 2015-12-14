SUBROUTINE mtrxi (FILE,NAME,item,dumbuf,itest)
     
!     COPIES MATRIX ITEM OF SUBSTRUCTURE NAME FROM THE SOF TO THE
!     NASTRAN FILE
!     ITEST VALUES RETURNED ARE
!        1 - NORMAL RETURN
!        2 - ITEM PSEUDO-EXISTS ONLY ON THE SOF
!        3 - ITEM DOES NOT EXIST ON THE SOF
!        4 - NAME DOES NOT EXIT
!        5 - ITEM IS NOT ONE OF THE ALLOWABLE MATRIX ITEMS
!        6 - THE NASTRAN FILE HAS BEEN PURGED
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN OUT)                  :: item
 REAL, INTENT(IN OUT)                     :: dumbuf
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        rshift, andf
 INTEGER :: rshift, andf, trail(7), buf(1), eof, nmsbr(2), oldbuf, blksiz
 COMMON /machin/ mach, ihalf, jhalf
 COMMON /sof   / ditdum(6), io, iopbn, iolbn, iomode
 COMMON /sys   / blksiz
 COMMON /zzzzzz/ nstrn
 EQUIVALENCE     (buf(1), nstrn)
 DATA    nmsbr / 4HMTRX,4HI   /
 DATA    ird   / 1 /,   idle  / 0 /,  ifetch / -1 /
 
!     CHECK IF ITEM IS ONE OF THE FOLLOWING ALLOWABLE NAMES.
!     KMTX, MMTX, PVEC, POVE, UPRT, HORG, UVEC, QVEC, PAPP, POAP, LMTX
 
 CALL chkopn (nmsbr(1))
 itm = ittype(item)
 IF (itm /= 1) GO TO 1030
 
!     MAKE SURE BUFFER IS DOUBLE WORD ALIGNED, OPEN NASTRAN FILE, AND
!     ADUST SOF BUFFER TO COINCIDE WITH GINO
!     ALSO DETERMINE PLACEMENT OF MATRIX NAME IN FIRST BUFFER
 
 idisp = locfx(buf(io-2)) - locfx(nstrn)
 IF (andf(idisp,1) /= 0) io = io + 1
 iopt = 1
 CALL OPEN (*1000,FILE,buf(io-2),iopt)
 oldbuf = io
 
 in = 4
 IF (mach > 2) GO TO 40
 in = 7
!IBMD 6/93 IF (BUF(IO-2) .EQ. FILE) GO TO 40
!IBMD 6/93 IO = IO + 1
!IBMD 6/93 IF (BUF(IO-2) .NE. FILE) GO TO 1010
 
!     OPEN ITEM TO READ AND FETCH FIRST BLOCK FROM SOF
 
 40 CONTINUE
 CALL sfetch (NAME(1),item,ifetch,itest)
 IF (itest /= 1) GO TO 1050
 
!     INSERT CORRECT MATRIX NAME INTO BUFFER
 
 CALL fname (FILE,buf(io+in))
 
!     WRITE BLOCK ON NASTRAN FILE
 
 ASSIGN 50 TO ijump
 eof = 0
 50 IF (buf(io+1) <= 0) GO TO 90
 CALL wrtblk (FILE,eof)
 
!     READ NEXT SOF BLOCK
 
 60 CALL fnxt (iopbn,inxt)
 IF (MOD(iopbn,2) == 1) GO TO 70
 next = andf(rshift(buf(inxt),ihalf),jhalf)
 GO TO 80
 70 next = andf(buf(inxt),jhalf)
 80 IF (next == 0) GO TO 1020
 iopbn = next
 iolbn = iolbn + 1
 CALL sofio (ird,iopbn,buf(io-2))
 GO TO ijump, (50,100)
 
!     LAST DATA BLOCK HAS BEEN READ FROM SOF
 
 90 itrail = buf(io+1)
 buf(io+1) = iolbn
 IF (itrail >= 0) GO TO 97
 trail(1) = FILE
 DO  i = 2,7
   trail(i) = buf(io+blksiz-7+i)
 END DO
 CALL wrttrl (trail)
 97 eof = 1
 CALL wrtblk (FILE,eof)
 CALL CLOSE (FILE,1)
 IF (itrail /= 0) GO TO 120
 
!     TRAILER IS STORED ON NEXT SOF BLOCK - READ IT
 
 ASSIGN 100 TO ijump
 GO TO 60
 
!     WRITE TRAILER OF NASTRAN DATA BLOCK
 
 100 trail(1) = FILE
 DO  i = 2,7
   trail(i) = buf(io+blksiz-7+i)
 END DO
 CALL wrttrl (trail)
 120 itest = 1
 io = oldbuf
 RETURN
 
!     ERROR RETURNS
 
 
!     NASTRAN FILE PURGED
 
 1000 itest = 6
 iomode = idle
 RETURN
 
!     BUFFER ALIGNMENT ERROR
 
!IBMD 6/93 1010 CALL SOFCLS
!IBMD 6/93 CALL MESAGE (-8,0,NMSBR)
 
!     SOF CHAINING ERROR
 
 1020 CALL errmkn (19,9)
 RETURN
 
!     INVALID ITEM NAME
 
 1030 itest = 5
 RETURN
 
!     ERROR IN SFETCH CALL
 
 1050 CALL CLOSE (FILE,1)
 io = oldbuf
 RETURN
 
END SUBROUTINE mtrxi
