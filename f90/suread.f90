SUBROUTINE suread (ia,nd,nout,itest)
     
!     READS DATA FROM THE SOF INTO THE ARRAY IA.  ND IS AN INPUT
!     PARAMETER INDICATING THE NUMBER OF WORDS DESIRED.  ND=-1 MEANS
!     READ UNTIL END OF GROUP. NOUT IS AN OUTPUT PARAMETER INDICATING
!     THE NUMBER OF WORD THAT HAVE BEEN READ.  ITEST IS AN OUTPUT
!     PARAMETER WHERE ITEST=3 MEANS END OF ITEM ENCOUNTERED, ITEST=2
!     MEANS END OF GROUP ENCOUNTERED, AND ITEST=1 OTHERWISE.
 
 
 INTEGER, INTENT(OUT)                     :: ia(1)
 INTEGER, INTENT(IN)                      :: nd
 INTEGER, INTENT(OUT)                     :: nout
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        andf,rshift
 INTEGER :: andf,rshift,buf,blksiz,dirsiz, nmsbr(2)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / ditdum(6),io,iopbn,iolbn,iomode,ioptr,iosind, ioitcd,ioblk
 COMMON /sys   / blksiz,dirsiz
 DATA    idle  , ird / 0,1    /
 DATA    ieog  , ieoi/ 4H$eog ,4H$eoi       /
 DATA    indsbr/ 19  /, nmsbr /4HSURE,4HAD  /
 
 CALL chkopn (nmsbr(1))
 icount = 0
 IF (iomode == ird) GO TO 20
 itest = 4
 nout  = 0
 RETURN
 
 10 icount = icount + 1
 ia(icount) = buf(ioptr)
 ioptr = ioptr + 1
 IF (icount == nd) GO TO 35
 20 IF (ioptr > blksiz+io) GO TO 80
 
!     READ SOF INTO ARRAY IA, BUT WATCH FOR END OF GROUP AND END OF ITEM
 
 30 IF (buf(ioptr) == ieoi) GO TO 50
 IF (buf(ioptr) == ieog .AND. nd /= -2) GO TO 40
 GO TO 10
 
!     READ THE REQUIRED NUMBER OF WORDS.
 
 35 itest = 1
 GO TO 70
 
!    REACHED END OF GROUP.
 
 40 itest = 2
 GO TO 60
 
!    REACHED END OF ITEM.
 
 50 itest  = 3
 iomode = idle
 60 ioptr  = ioptr + 1
 70 nout   = icount
 RETURN
 
!     REACHED END OF BLOCK.  REPLACE THE BLOCK CURRENTLY IN CORE BY ITS
!     LINK BLOCK.
 
 80 CALL fnxt (iopbn,inxt)
 IF (MOD(iopbn,2) == 1) GO TO 90
 next = andf(rshift(buf(inxt),ihalf),jhalf)
 GO TO 100
 90 next = andf(buf(inxt),jhalf)
 100 IF (next == 0) GO TO 510
 iopbn = next
 iolbn = iolbn + 1
 CALL sofio (ird,iopbn,buf(io-2))
 ioptr = io + 1
 GO TO 30
 
!     ERROR MESSAGES.
 
 510 CALL errmkn (indsbr,9)
 RETURN
END SUBROUTINE suread
