SUBROUTINE fdsub (NAME,i)
!                                                           ** PRETTIED
!     SEARCHES IF THE SUBSTRUCTURE NAME HAS AN ENTRY IN THE DIT. IF IT
!     DOES, THE OUTPUT VALUE OF I WILL INDICATE THAT NAME IS THE ITH
!     SUBSTRUCTURE IN THE DIT.  I WILL BE SET TO -1 IF NAME DOES NOT
!     HAVE AN ENTRY IN THEDIT.
 
 
 INTEGER, INTENT(IN)                      :: NAME(2)
 INTEGER, INTENT(OUT)                     :: i
 LOGICAL :: ditup
 INTEGER :: buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl, blksiz,dirsiz
 DIMENSION  nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,iodum(8),  &
     mdidum(4),nxtdum(15),ditup
 COMMON /sys   / blksiz,dirsiz
 DATA    nmsbr / 4HFDUB,4HB   /
 
!     NNMS IS THE NUMBER OF NAMES ON ONE BLOCK OF THE DIT, AND NBLKS IS
!     THE SIZE OF THE DIT IN NUMBER OF BLOCKS.
 
 CALL chkopn (nmsbr(1))
 IF (ditnsb == 0) GO TO 70
 nnms  = blksiz/2
 nblks = ditsiz/blksiz
 IF (ditsiz == nblks*blksiz) GO TO 30
 nblks = nblks + 1
 
!     START LOOKING FOR THE SUBSTRUCTURE NAME.
 
 30 MAX = blksiz
 DO  j = 1,nblks
   i = 1 + (j-1)*nnms
   CALL fdit (i,dummy)
   IF (j /= nblks) GO TO 40
   MAX = ditsiz - (nblks-1)*blksiz
   
!     SEARCH THE BLOCK OF THE DIT WHICH IS PRESENTLY IN CORE.
   
   40 DO  k = 1,MAX,2
     IF (buf(dit+k) /= NAME(1) .OR. buf(dit+k+1) /= NAME(2)) CYCLE
     kk = k
     GO TO 80
   END DO
 END DO
 
!     DID NOT FIND NAME IN THE DIT.
 
 70 i = -1
 RETURN
 
!     DID FIND NAME IN THE DIT.  RETURN NAME INDEX NUMBER
 
 80 i = (ditlbn-1)*nnms + (kk+1)/2
 RETURN
END SUBROUTINE fdsub
