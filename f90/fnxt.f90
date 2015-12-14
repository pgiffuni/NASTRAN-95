SUBROUTINE fnxt (ii,j)
     
!     FETCHES FROM THE RANDOM ACCESS STORAGE DEVICE THE BLOCK OF THE
!     ARRAY NXT CONTAINING THE ENTRY FOR BLOCK I.  IT STORES THE FETCHED
!     BLOCK IN THE ARRAY BUF, STARTING AT LOCATION NXT.  THE OUTPUT J
!     INDICATES THAT BLOCK I HAS THE JTH ENTRY IN THE ARRAY BUF.
 
 
 INTEGER, INTENT(IN)                      :: ii
 INTEGER, INTENT(OUT)                     :: j
 LOGICAL :: ditup,nxtup
 INTEGER :: buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     blksiz,dirsiz,supsiz,filsiz,filnum,filsup
 DIMENSION       nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,iodum(8),  &
     mdidum(4),nxt,nxtpbn,nxtlbn,nxttsz,nxtfsz(10),  &
     nxtcur,ditup,mdiup,nxtup,nxtrst
 COMMON /sys   / blksiz,dirsiz,supsiz
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10)
 DATA    ird   , iwrt,indsbr  / 1,2, 9 /
 DATA    nmsbr / 4HFNXT,4H    /
 
!     FILNUM IS THE NUMBER OF THE DEVICE TO WHICH BLOCK I BELONGS.
 
 CALL chkopn (nmsbr(1))
 INDEX = ii
 DO  l = 1,nfiles
   IF (INDEX > filsiz(l)) GO TO 2
   filnum = l
   GO TO 10
   2 INDEX = INDEX - filsiz(l)
 END DO
 GO TO 500
 
!     INDEX IS THE INDEX OF BLOCK I WITHIN FILE FILNUM.
!     FILSUP IS THE NUMBER OF THE SUPERBLOCK WITHIN FILE FILNUM TO WHICH
!     BLOCK I BELONGS, AND SUPSIZ IS THE SIZE OF A SUPERBLOCK.
 
 10 filsup = (INDEX-1)/supsiz
 IF (INDEX-1 == filsup*supsiz) GO TO 20
 filsup = filsup + 1
 
!     COMPUTE THE LOGICAL BLOCK NUMBER, WITHIN THE ARRAY NXT, IN WHICH
!     THE ITH BLOCK HAS AN ENTRY, ALSO COMPUTE THE INDEX OF THIS ENTRY
!     RELATIVE TO THE ARRAY BUF.  STORE THE BLOCK NUMBER IN IBLOCK, AND
!     THE INDEX IN J.
 
 20 iblock = 0
 MAX = filnum - 1
 IF (MAX < 1) GO TO 26
 DO  i = 1,MAX
   iblock = iblock + nxtfsz(i)
 END DO
 26 iblock = iblock + filsup
 j = (INDEX-(filsup-1)*supsiz)/2 + 1 + nxt
 IF (iblock == nxtlbn) RETURN
 IF (iblock > nxttsz) GO TO 500
 
!     THE DESIRED NXT BLOCK IS NOT PRESENTLY IN CORE, MUST THEREFORE
!     FETCH IT.
 
 IF (ditpbn == 0) GO TO 40
 
!     THE IN CORE BLOCK SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY ONE BLOCK OF THE DIT.
 
 IF (.NOT.ditup) GO TO 30
 
!     THE DIT BLOCK NOW IN CORE HAS BEEN UPDATED.  MUST THEREFORE WRITE
!     IT OUT BEFORE READING IN THE DESIRED NXT BLOCK.
 
 CALL sofio (iwrt,ditpbn,buf(dit-2))
 ditup  = .false.
 30 ditpbn = 0
 ditlbn = 0
 GO TO 50
 40 IF (nxtpbn == 0) GO TO 50
 
!     THE IN CORE BLOCK SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY ONE BLOCK OF NXT.
 
 IF (.NOT.nxtup) GO TO 50
 
!     THE NEXT BLOCK CURRENTLY IN CORE HAS BEEN UPDATED.  MUST THEREFORE
!     WRITE IT OUT BEFORE READING IN A NEW BLOCK.
 
 CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 nxtup = .false.
 
!     READ THE DESIRED NXT BLOCK INTO CORE.
 
 50 nxtlbn = iblock
 nxtpbn = 0
 IF (MAX < 1) GO TO 70
 DO  i = 1,MAX
   nxtpbn = nxtpbn+filsiz(i)
 END DO
 70 nxtpbn = nxtpbn + (filsup-1)*supsiz + 2
 CALL sofio (ird,nxtpbn,buf(nxt-2))
 RETURN
 
!     ERROR MESSAGES.
 
 500 CALL errmkn (indsbr,1)
 RETURN
END SUBROUTINE fnxt
