SUBROUTINE sjump (n)
     
!     JUMP OVER N GROUPS WITHIN AN ITEM WHEN IN READ MODE.  N WILL BE
!     RETURNED AS -1 IF THE END OF ITEM IS REACHED BEFORE JUMPING OVER
!     N GROUPS.
 
 
 INTEGER, INTENT(OUT)                     :: n
 EXTERNAL        andf,rshift
 INTEGER :: andf,rshift,buf,eog,eoi,blksiz,dirsiz,nmsbr(2)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / ditdum(6),io,iopbn,iolbn,iomode,ioptr,iosind, ioitcd,ioblk
 COMMON /sys   / blksiz,dirsiz
 DATA    ird   / 1   /
 DATA    eog   , eoi /  4H$eog ,4H$eoi       /
 DATA    indsbr/ 17  /, nmsbr  /4HSJUM,4HP   /
 
 CALL chkopn (nmsbr(1))
 IF (n <= 0) RETURN
 icount = 0
 IF (iomode == ird) GO TO 20
 n = -2
 RETURN
 
!     SEARCH THROUGH SOF FOR END OF ITEM AND END OF GROUP.
 
 10 ioptr = ioptr + 1
 20 IF (ioptr > blksiz+io) GO TO 50
 30 IF (buf(ioptr) /=  eoi) GO TO 40
 n = -1
 RETURN
 
 40 IF (buf(ioptr) /= eog) GO TO 10
 icount = icount + 1
 IF (icount /= n) GO TO 10
 ioptr = ioptr + 1
 RETURN
 
!     REACHED END OF BLOCK.  REPLACE THE BLOCK CURRENTLY IN CORE BY ITS
!     LINK BLOCK.
 
 50 CALL fnxt (iopbn,inxt)
 IF (MOD(iopbn,2) == 1) GO TO 60
 next = andf(rshift(buf(inxt),ihalf),jhalf)
 GO TO 70
 60 next = andf(buf(inxt),jhalf)
 70 IF (next == 0) GO TO 510
 iopbn = next
 iolbn = iolbn + 1
 CALL sofio (ird,iopbn,buf(io-2))
 ioptr = io + 1
 GO TO 30
 510 CALL errmkn (indsbr,9)
 RETURN
END SUBROUTINE sjump
