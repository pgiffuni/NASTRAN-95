SUBROUTINE retblk (ibl)
     
!     RETURNS BLOCK I AND ALL BLOCKS LINKED TO IT TO THE LIST OF FREE
!     BLOCKS IN THE SUPERBLOCK TO WHICH BLOCK I BELONGS   IF SOME OF
!     THE BLOCKS THAT ARE LINKED TO BLOCK I DO NOT BELONG TO THE SAME
!     SUPERBLOCK, THEY ARE RETURNED TO THE FREE LIST OF THEIR OWN
!     RESPECTIVE SUPERBLOCKS.
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN)                      :: ibl
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: ditup,nxtup,repeat
 DIMENSION       nmsbr(2)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl, iodum(8),mdidum(4),  &
     nxt,nxtpbn,nxtlbn,nxttsz,nxtfsz(10),nxtcur, ditup,mdiup,nxtup
 COMMON /sys   / blksiz,dirsiz,supsiz,avblks
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10)
 DATA    ird   , iwrt/ 1,2   /
 DATA    indsbr/ 13 /, nmsbr /4HRETB,4HLK  /
 
 CALL chkopn (nmsbr(1))
 i = ibl
 IF (i <= 0) GO TO 500
 lmask = lshift(jhalf,ihalf)
 
!     COMPUTE THE NUMBER OF THE FILE TO WHICH BLOCK I BELONGS,
!     THE INDEX OF BLOCK I WITHIN THAT FILE, THE NUMBER WITHIN THE
!     FILE OF THE SUPERBLOCK TO WHICH BLOCK I BELONGS, AND THE LOGICAL
!     BLOCK NUMBER OVER THE SYSTEM OF THAT SUPERBLOCK.
 
 5 left = i
 DO  l = 1,nfiles
   IF (left > filsiz(l)) GO TO 6
   filnum = l
   GO TO 10
   6 left = left - filsiz(l)
 END DO
 GO TO 500
 10 filind = left
 filsup = (filind-1)/supsiz
 IF (filind-1 == filsup*supsiz) GO TO 20
 filsup = filsup + 1
 20 ilbn = 0
 MAX  = filnum - 1
 IF (MAX < 1) GO TO 28
 DO  l = 1,MAX
   ilbn = ilbn + nxtfsz(l)
 END DO
 28 ilbn = ilbn + filsup
 IF (ilbn == nxtlbn) GO TO 60
 
!     THE DESIRED BLOCK OF THE ARRAY NXT IS NOT IN CORE.
 
 IF (nxtlbn == 0) GO TO 30
 
!     THE IN CORE BUFFER SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY A BLOCK OF NXT.  IF THAT BLOCK HAS BEEN UPDATED,
!     MUST WRITE IT OUT BEFORE READING IN THE NEW BLOCK.
 
 IF (.NOT.nxtup) GO TO 50
 CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 nxtup = .false.
 GO TO 50
 
!     THE IN CORE BUFFER SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY A BLOCK OF THE DIT.  IF THAT BLOCK HAS BEEN UPDATED,
!     MUST WRITE IT OUT BEFORE READING IN THE NEW BLOCK.
 
 30 IF (.NOT.ditup) GO TO 40
 CALL sofio (iwrt,ditpbn,buf(dit-2))
 ditup  = .false.
 40 ditpbn = 0
 ditlbn = 0
 
!     READ IN THE DESIRED BLOCK OF NXT.
 
 50 nxtlbn = ilbn
 nxtpbn = 0
 MAX = filnum - 1
 IF (MAX < 1) GO TO 58
 DO  l = 1,MAX
   nxtpbn = nxtpbn + filsiz(l)
 END DO
 58 nxtpbn = nxtpbn + (filsup-1)*supsiz + 2
 CALL sofio (ird,nxtpbn,buf(nxt-2))
 
!     THE DESIRED BLOCK OF NXT IS IN CORE.
 
 60 btfree = andf(buf(nxt+1),jhalf)
 tpfree = rshift(buf(nxt+1),ihalf)
 IF (btfree == 0) GO TO 90
 
!     CHECK IF BLOCK I IS ALREADY IN THE LIST OF FREE BLOCKS.
 
 j = tpfree
 70 IF (j == i) GO TO 220
 IF (j == 0) GO TO 90
 ind = (j-nxtpbn+2)/2 + 1
 IF (MOD(j,2) == 1) GO TO 80
 j = rshift(buf(nxt+ind),ihalf)
 GO TO 70
 80 j = andf(buf(nxt+ind),jhalf)
 GO TO 70
 
!     BLOCK I IS NOT IN THE LIST OF FREE BLOCKS.
!     SET TPFREE TO I
 
 90 buf(nxt+1) = lshift(i,ihalf)
 
!     EXAMINE THE BLOCKS THAT ARE LINKED TO BLOCK I.
 
 repeat = .false.
 IF (filsup /= nxtfsz(filnum)) GO TO 123
 lstblk = nxtpbn + filsiz(filnum) - (filsup-1)*supsiz - 2
 GO TO 126
 123 lstblk = nxtpbn + supsiz - 1
 126 avblks = avblks + 1
 ind = (i-nxtpbn+2)/2 + 1
 IF (MOD(i,2) == 1) GO TO 130
 isv = rshift(buf(nxt+ind),ihalf)
 GO TO 140
 130 isv = andf(buf(nxt+ind),jhalf)
 140 IF (isv == 0) GO TO 160
 IF (isv < nxtpbn .OR. isv > lstblk) GO TO 150
 i = isv
 GO TO 126
 150 repeat = .true.
 
!     ALL THE BLOCKS IN THIS SUPERBLOCK HAVE BEEN FOUND.
!     SET POINTER OF I TO VALUE OF OLD TPFREE.
 
 160 IF (MOD(i,2) == 1) GO TO 170
 buf(nxt+ind) = orf(andf(buf(nxt+ind),jhalf),lshift(tpfree,ihalf))
 GO TO 180
 170 buf(nxt+ind) = orf(andf(buf(nxt+ind),lmask),tpfree)
 180 IF(btfree == 0) btfree = i
 
!     SET BTFREE TO LAST BLOCK IN CHAIN.
 
 buf(nxt+1) = orf(andf(buf(nxt+1),lmask),btfree)
 nxtup = .true.
 IF (.NOT. repeat) GO TO 220
 
!     ISV BELONGS TO A DIFFERENT SUPERBLOCK, REPEAT
!     SUBROUTINE FOR BLOCK ISV.
 
 i = isv
 GO TO 5
 
!     NO MORE BLOCKS LINKED TO BLOCK I, RETURN.
 
 220 CONTINUE
 RETURN
 
!     ERROR MESSAGE.
 
 500 CALL errmkn (indsbr,2)
 RETURN
END SUBROUTINE retblk
