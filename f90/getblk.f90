SUBROUTINE getblk (iold,inew)
     
!     FINDS A FREE BLOCK INEW.  IF IOLD IS NOT ZERO IOLD POINTER WILL
!     BE SET TO INEW.
 
 
 INTEGER, INTENT(IN OUT)                  :: iold
 INTEGER, INTENT(OUT)                     :: inew
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: ditup,nxtup,nxtrst
 INTEGER :: buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,orf,  &
     blksiz,dirsiz,supsiz,filsiz,avblks,andf,rshift,  &
     filnum,tpfree,btfree,filind,filsup,nmsbr(2)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl, iodum(8),mdidum(4),  &
     nxt,nxtpbn,nxtlbn,nxttsz,nxtfsz(10),nxtcur, ditup,mdiup,nxtup,nxtrst
 COMMON /sys   / blksiz,dirsiz,supsiz,avblks
 COMMON /sofcom/ nfiles,filnam(10),filsiz(10)
 COMMON /system/ nbuff
 DATA   ird    , iwrt / 1, 2    /
 DATA   indsbr / 11   /,  nmsbr / 4HGETB,4HLK  /
 
!     CHECK IF THE SUPERBLOCK NXTCUR HAS A FREE BLOCK.
 
 CALL chkopn (nmsbr(1))
 lmask = lshift(jhalf,ihalf)
 5 IF (nxtcur == nxtlbn) GO TO 40
 
!     THE SUPERBLOCK NXTCUR IS NOT IN CORE.
 
 IF (nxtlbn == 0) GO TO 10
 
!     THE IN CORE BUFFER SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY A BLOCK OF NXT.
 
 IF (.NOT.nxtup) GO TO 20
 
!     THE BLOCK OF THE ARRAY NXT WHICH IS NOW IN CORE HAS BEEN UPDATED,
!     MUST THEREFORE WRITE IT OUT BEFORE READING IN A NEW BLOCK.
 
 CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 nxtup = .false.
 GO TO 20
 10 IF (ditpbn == 0) GO TO 20
 
!     THE IN CORE BUFFER SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY A BLOCK OF DIT.
 
 IF (.NOT.ditup) GO TO 15
 
!     THE DIT BLOCK WHICH IS NOW IN CORE HAS BEEN UPDATED, MUST
!     THEREFORE WRITE IT OUT BEFORE READING IN A NEW BLOCK.
 
 CALL sofio (iwrt,ditpbn,buf(dit-2))
 ditup  = .false.
 15 ditpbn = 0
 ditlbn = 0
 
!     READ INTO CORE THE DESIRED BLOCK OF THE ARRAY NXT.
 
 20 nxtlbn = nxtcur
 nxtpbn = 0
 left   = nxtlbn
 DO  i = 1,nfiles
   IF (left > nxtfsz(i)) GO TO 23
   filnum = i
   GO TO 30
   23 nxtpbn = nxtpbn + filsiz(i)
   left   = left - nxtfsz(i)
 END DO
 GO TO 510
 30 nxtpbn = nxtpbn + (left-1)*supsiz + 2
 CALL sofio (ird,nxtpbn,buf(nxt-2))
 
!     CHECK THE FREE LIST OF SUPERBLOCK NXTCUR.
 
 40 tpfree = rshift(buf(nxt+1),ihalf)
 IF (tpfree > 0) GO TO 180
 
!     THE SUPERBLOCK NXTCUR DOES NOT HAVE ANY FREE BLOCKS.
 
 IF (nxtcur == nxttsz) GO TO 50
 nxtcur = nxtcur + 1
 GO TO 5
 
!     NXTCUR IS THE LAST SUPERBLOCK.
 
 50 IF (nxtrst) GO TO 60
 nxtcur = 1
 nxtrst = .true.
 GO TO 5
 
!     MUST START A BRAND NEW SUPERBLOCK.
 
 60 nxtrst = .false.
 IF (nxtup) CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 nxtup  = .false.
 70 nxtcur = nxtcur + 1
 left   = nxtcur
 DO  i = 1,nfiles
   IF (left > nxtfsz(i)) GO TO 75
   filnum = i
   GO TO 85
   75 left   = left-nxtfsz(i)
 END DO
 nxtcur = nxtcur - 1
 GO TO 500
 85 last   = nbuff - 4
 DO  i = 1,last
   buf(nxt+i) = 0
 END DO
 IF (left == 1) GO TO 110
 
!     NXTCUR IS NOT THE FIRST SUPERBLOCK ON FILE FILNUM.
 
 nxtpbn = nxtpbn + supsiz
 IF (left /= nxtfsz(filnum)) GO TO 120
 
!     NXTCUR IS THE LAST BLOCK ON FILE FILNUM.
 
 lstsiz = MOD(filsiz(filnum)-2,supsiz) + 1
 IF (lstsiz > 1) GO TO 90
 
!     THE SIZE OF THE LAST BLOCK ON FILE FILNUM IS EQUAL TO 1.
!     THERE ARE THEREFORE NO FREE BLOCKS AVAILABLE ON SUPERBLOCK NXTCUR.
!     SET TPFREE AND BTFREE OF NXTCUR EQUAL TO ZERO.
 
 buf(nxt+1) = 0
 avblks = avblks - 1
 CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 GO TO 70
 
!     THE SIZE OF SUPERBLOCK NXTCUR IS LARGER THAN 1.
 
 90 IF (lstsiz > 2) GO TO 100
 
!     THE SIZE OF SUPERBLOCK NXTCUR IS EQUAL TO 2.  THERE IS THEREFORE
!     ONLY ONE FREE BLOCK IN NXTCUR.  SET TPFREE AND BTFREE TO ZERO.
 
 buf(nxt+1) = 0
 GO TO 170
 
!     THE SIZE OF SUPERBLOCK NXTCUR IS LARGER THAN 2.
 
 100 btfree = nxtpbn + lstsiz - 1
 GO TO 130
 
!     NXTCUR IS THE FIRST SUPERBLOCK ON FILE FILNUM.
 
 110 lstsiz = MOD(filsiz(filnum-1)-2,supsiz) + 1
 nxtpbn = nxtpbn + lstsiz + 1
 avblks = avblks - 1
 IF (filsiz(filnum) >= supsiz+1) GO TO 120
 btfree = nxtpbn + filsiz(filnum) - 2
 GO TO 130
 120 btfree = nxtpbn + supsiz - 1
 
!     INITIALIZE THE NEW SUPERBLOCK.
 
 130 tpfree = nxtpbn + 2
 
!     PUT THE VALUES OF BTFREE AND TPFREE IN THE FIRST WORD OF THE ARRAY
!     NXT BELONGING TO SUPERBLOCK NXTCUR.
 
 buf(nxt+1) = btfree
 buf(nxt+1) = orf(buf(nxt+1),lshift(tpfree,ihalf))
 IF (MOD(btfree,2) == 1) GO TO 140
 
!     BTFREE IS AN EVEN INTEGER.
 
 MAX = (btfree-nxtpbn+2)/2
 buf(nxt+MAX+1) = 0
 GO TO 150
 
!     BTFREE IS AN ODD INTEGER.
 
 140 MAX = (btfree-nxtpbn+1)/2
 buf(nxt+MAX+1) = lshift(btfree,ihalf)
 
!     SET UP THE THREAD THROUGH THE BLOCKS OF SUPERBLOCK NXTCUR.
 
 150 IF (MAX < 3) GO TO 170
 DO  i = 3,MAX
   buf(nxt+i) = 2*i + nxtpbn - 2
   buf(nxt+i) = orf(buf(nxt+i),lshift(2*i+nxtpbn-3,ihalf))
 END DO
 
!     SETUP VARIABLES RELATED TO THE SUPERBLOCK NXTCUR.
 
 170 buf(nxt+2) = 0
 inew   = nxtpbn + 1
 avblks = avblks - 2
 nxtlbn = nxtcur
 nxttsz = nxtcur
 GO TO 230
 
!     SUPERBLOCK NXTCUR DOES HAVE A FREE BLOCK.
 
 180 inew   = tpfree
 avblks = avblks - 1
 
!     COMPUTE THE INDEX OF TPFREE ENTRY IN THE BLOCK OF ARRAY NXT
!     BELONGING TO SUPERBLOCK NXTCUR.
 
 filind = tpfree
 DO  i = 1,nfiles
   IF (filind <= filsiz(i)) EXIT
   filind = filind - filsiz(i)
 END DO
 187 filsup = (filind-1)/supsiz
 IF (filind-1 == filsup*supsiz) GO TO 190
 filsup = filsup + 1
 190 INDEX  = (filind-(filsup-1)*supsiz)/2 + 1
 IF (MOD(tpfree,2) == 1) GO TO 200
 
!     TPFREE IS AN EVEN INTEGER.  THE ENTRY FOR TPFREE IS THEREFORE
!     IN BITS (IHALF+1) TO (2*IHALF-1) OF THE WORD.  SAVE TPFREE ENTRY
!     IN NXTBLK AND THEN SET IT TO ZERO.
 
 nxtblk = rshift(buf(nxt+INDEX),ihalf)
 buf(nxt+INDEX) = andf(buf(nxt+INDEX),jhalf)
 GO TO 210
 
!     TPFREE IS AN ODD INTEGER.  THE ENTRY FOR TPFREE IS THEREFORE
!     IN BITS 0 TO IHALF OF THE WORD.  SAVE TPFREE ENTRY IN NXTBLK
!     AND THEN SET IT TO ZERO.
 
 200 nxtblk = andf(buf(nxt+INDEX),jhalf)
 buf(nxt+INDEX) = andf(buf(nxt+INDEX),lmask)
 210 btfree = andf(buf(nxt+1),jhalf)
 IF (tpfree == btfree) GO TO 220
 
!     SET TPFREE TO NXTBLK.
 
 buf(nxt+1) = orf(andf(buf(nxt+1),jhalf),lshift(nxtblk,ihalf))
 GO TO 230
 
!     SET TPFREE AND BTFREE TO ZERO.
 
 220 buf(nxt+1) = 0
 230 IF (iold == 0) GO TO 250
 
!     WANT TO SET IOLD POINTER TO INEW.
 
 nxtup =.true.
 CALL fnxt (iold,ind)
 IF (MOD(iold,2) == 1) GO TO 240
 
!     IOLD IS AN EVEN INTEGER
 
 buf(ind) = orf(andf(buf(ind),jhalf),lshift(inew,ihalf))
 GO TO 250
 
!     IOLD IS AN ODD INTEGER
 
 240 buf(ind) = orf(andf(buf(ind),lmask),inew)
 250 nxtup = .true.
 RETURN
 
!     ERROR MESSAGES.
 
 500 inew = -1
 RETURN
 510 CALL errmkn (indsbr,4)
 RETURN
END SUBROUTINE getblk
