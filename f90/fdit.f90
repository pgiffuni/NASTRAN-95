SUBROUTINE fdit (i,k)
     
!     FETCHES FROM THE RANDOM ACCESS STORAGE DEVICE THE BLOCK OF THE
!     DIT CONTAINING THE ITH SUBSTRUCTURE NAME, AND STORES IT IN THE
!     ARRAY BUF STARTING AT LOCATION (DIT+1) AND EXTENDING TO LOCATION
!     (DIT+BLKSIZ).  THE OUTPUT K INDICATES THAT THE SUBSTRUCTURE HAS
!     THE KTH ENTRY IN BUF.
 
 
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(OUT)                     :: k
 EXTERNAL        rshift,andf
 LOGICAL :: ditup,nxtup,newblk
 INTEGER :: buf,dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl,  &
     blksiz,dirsiz,andf,rshift,nmsbr(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / dit,ditpbn,ditlbn,ditsiz,ditnsb,ditbl, iodum(8),mdidum(4),  &
     nxt,nxtpbn,nxtlbn,nxttsz,nxtfsz(10),nxtcur, ditup,mdiup,nxtup,nxtrst
 COMMON /sys   / blksiz,dirsiz
 COMMON /system/ nbuff,nout
 COMMON /machin/ mach,ihalf,jhalf
 DATA    ird   , iwrt  / 1,2   /
 DATA    iempty/ 4H    /
 DATA    indsbr/ 5     /, nmsbr /4HFDIT,4H    /
 
 CALL chkopn (nmsbr(1))
 
!     NDIR IS THE NUMBER OF SUBSTRUCTURE NAMES IN ONE BLOCK OF THE DIT
 
 ndir = blksiz/2
 
!     COMPUTE THE LOGICAL BLOCK NUMBER, AND THE WORD NUMBER WITHIN
!     BUF IN WHICH THE ITH SUBSTRUCTURE NAME IS STORED.  STORE THE BLOCK
!     NUMBER IN IBLOCK, AND THE WORD NUMBER IN K.
 
 iblock = i/ndir
 IF (i == iblock*ndir) GO TO 10
 iblock = iblock + 1
 10 k = 2*(i-(iblock-1)*ndir) - 1 + dit
 IF (ditlbn == iblock) RETURN
 
!     THE DESIRED DIT BLOCK IS NOT PRESENTLY IN CORE, MUST THEREFORE
!     FETCH IT.
 
 newblk = .false.
 
!     FIND THE PHYSICAL BLOCK NUMBER OF THE BLOCK ON WHICH THE LOGICAL
!     BLOCK IBLOCK IS STORED.
 
 j = ditbl
 icount = 1
 30 IF (icount == iblock) GO TO 40
 icount = icount + 1
 CALL fnxt (j,inxt)
 IF (MOD(j,2) == 1) GO TO 33
 ibl = rshift(buf(inxt),ihalf)
 GO TO 36
 33 ibl = andf(buf(inxt),jhalf)
 36 IF (ibl == 0) GO TO 70
 j = ibl
 GO TO 30
 40 IF (ditpbn == 0) GO TO 43
 
!     THE IN CORE BLOCK SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY THE DIT.  WRITE IT OUT IF IT HAS BEEN UPDATED.
 
 IF (.NOT.ditup) GO TO 50
 CALL sofio (iwrt,ditpbn,buf(dit-2))
 GO TO 50
 43 IF (nxtpbn == 0) GO TO 50
 
!     THE IN CORE BLOCK SHARED BY THE DIT AND THE ARRAY NXT IS NOW
!     OCCUPIED BY NXT.  WRITE OUT NXT IF IT HAS BEEN UPDATED.
 
 IF (.NOT.nxtup) GO TO 46
 CALL sofio (iwrt,nxtpbn,buf(nxt-2))
 nxtup  = .false.
 46 nxtpbn = 0
 nxtlbn = 0
 
!     READ THE DESIRED DIT BLOCK INTO CORE.
 
 50 ditpbn = j
 ditlbn = iblock
 IF (newblk) GO TO 60
 CALL sofio (ird,j,buf(dit-2))
 RETURN
 
 60 istart = dit + 1
 iend   = dit + blksiz
 DO  ll = istart,iend
   buf(ll) = iempty
 END DO
 RETURN
 
!     WE NEED A FREE BLOCK FOR THE DIT.
 
 70 CALL getblk (j,ibl)
 IF (ibl == -1) GO TO 80
 newblk = .true.
 j = ibl
 IF (icount == iblock) GO TO 40
 
!     ERROR MESSAGES.
 
 CALL errmkn (indsbr,7)
 80 WRITE  (nout,85) ufm
 85 FORMAT (a23,' 6223, SUBROUTINE FDIT - THERE ARE NO MORE FREE ',  &
     'BLOCKS AVAILABLE ON THE SOF')
 CALL sofcls
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE fdit
