SUBROUTINE suwrt (ia,nwords,itest)
     
!     COPIES DATA FROM THE ARRAY IA ON THE SOF.  NWORD IS AN INPUT
!     PARAMETER INDICATING THE NUMBER OF WORDS TO BE COPIED.  ITEST IS
!     AN INPUT PARAMETER WHERE ITEST=1 MEANS MORE TO COME, ITEST=2 MEANS
!     WRITE END OF GROUP, AND ITEST=3 MEANS WRITE END OF ITEM.
 
 
 INTEGER, INTENT(IN)                      :: ia(1)
 INTEGER, INTENT(IN OUT)                  :: nwords
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        lshift,andf,orf
 LOGICAL :: mdiup
 INTEGER :: buf,mdi,mdipbn,mdilbn,mdibl,blksiz,dirsiz,andf,orf
 DIMENSION  nmsbr(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / ditdum(6),  &
     io,iopbn,iolbn,iomode,ioptr,iosind,ioitcd,ioblk,  &
     mdi,mdipbn,mdilbn,mdibl,nxtdum(15),ditup,mdiup
 COMMON /sys   / blksiz,dirsiz
 COMMON /system/ nbuff,nout
 DATA    idle  , iwrt / 0,2    /
 DATA    ieog  , ieoi / 4H$eog ,4H$eoi /, nmsbr /4HSUWR,4HT   /
 
 CALL chkopn (nmsbr(1))
 icount = 0
 IF (iomode == iwrt) GO TO 10
 itest = 4
 RETURN
 
!     KEEP COPYING DATA FROM THE ARRAY IA INTO THE INPUT/OUTPUT BUFFER
!     UNTIL THE BUFFER IS FULL, OR UNTIL THE REQUESTED NUMBER OF WORDS
!     HAS BEEN COPIED.
 
 10 IF (ioptr > blksiz+io) GO TO 30
 20 IF (icount == nwords) SELECT CASE ( itest )
   CASE (    1)
     GO TO 80
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 50
 END SELECT
 icount = icount + 1
 buf(ioptr) = ia(icount)
 ioptr = ioptr + 1
 GO TO 10
 
!     THE BUFFER IS FULL.  OUTPUT IT ON THE SOF.
 
 30 CALL sofio (iwrt,iopbn,buf(io-2))
 CALL getblk (iopbn,j)
 IF (j == -1) GO TO 40
 iopbn = j
 iolbn = iolbn + 1
 ioptr = io + 1
 GO TO 20
 
!     THERE ARE NO MORE FREE BLOCKS ON THE SOF.  RETURN THE BLOCKS THAT
!     HAVE BEEN USED SO FAR BY THE ITEM BEING WRITTEN, AND CLOSE THE SOF
!     THEN ISSUE A FATAL ERROR MESSAGE.
 
 40 CALL retblk (ioblk)
 CALL sofcls
 GO TO 90
 
!     WRITE END OF ITEM, OUTPUT THE INPUT/OUTPUT BUFFER ON THE SOF, AND
!     UPDATE THE MDI.
 
 50 buf(ioptr) = ieoi
 CALL sofio (iwrt,iopbn,buf(io-2))
 CALL fmdi (iosind,imdi)
 buf(imdi+ioitcd) = ioblk
 buf(imdi+ioitcd) = orf(andf(buf(imdi+ioitcd),jhalf), lshift(iolbn,ihalf))
 mdiup  = .true.
 iomode = idle
 GO TO 70
 
!     WRITE END OF GROUP.
 
 60 buf(ioptr) = ieog
 70 ioptr = ioptr + 1
 80 RETURN
 
!     ERROR MESSAGES.
 
 90 WRITE  (nout,100) ufm
 100 FORMAT (a23,' 6223, THERE ARE NO MORE FREE BLOCKS AVAILABLE ON',  &
     ' THE SOF FILE.')
 CALL sofcls
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE suwrt
