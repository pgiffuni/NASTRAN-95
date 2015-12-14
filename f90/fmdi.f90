SUBROUTINE fmdi (i,j)
     
!     THE SUBROUTINE FETCHES FROM THE RANDOM ACCESS STORAGE DEVICE THE
!     BLOCK OF MDI CONTAINING THE I-TH DIRECTORY, AND STORES THAT BLOCK
!     IN THE ARRAY BUF STARTING AT LOCATION (MDI+1) AND EXTENDING TO
!     LOCATION (MDI+BLKSIZ).  IT ALSO RETURNS IN J THE (INDEX-1) OF THE
!     DIRECTORY IN BUF.
 
 
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(OUT)                     :: j
 EXTERNAL        rshift,andf
 LOGICAL :: mdiup,newblk
 INTEGER :: buf,mdi,mdipbn,mdilbn,mdibl,blksiz,dirsiz, andf,rshift,nmsbr(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / ditdum(6),iodum(8),mdi,mdipbn,mdilbn,mdibl,  &
     nxtdum(15),ditup,mdiup
 COMMON /sys   / blksiz,dirsiz
 COMMON /system/ nbuff,nout
 DATA    ird   , iwrt / 1, 2    /
 DATA    indsbr/ 7    /,  nmsbr /4HFMDI,4H    /
 
!     NDIR IS THE NUMBER OF DIRECTORIES ON ONE BLOCK OF THE MDI.
 
 CALL chkopn (nmsbr(1))
 ndir = blksiz/dirsiz
 
!     COMPUTE THE LOGICAL BLOCK NUMBER, AND THE WORD NUMBER WITHIN
!     BUF IN WHICH THE ITH SUBSTRUCTURE DIRECTORY IS STORED.  STORE THE
!     BLOCK NUMBER IN IBLOCK, AND THE WORD NUMBER IN J.
 
 iblock = i/ndir
 IF (i == iblock*ndir) GO TO 10
 iblock = iblock + 1
 10 j = dirsiz*(i-(iblock-1)*ndir-1) + mdi
 IF (mdilbn == iblock) RETURN
 IF (mdipbn == 0) GO TO 20
 IF (.NOT.mdiup) GO TO 20
 
!     THE MDI BLOCK CURRENTLY IN CORE HAS BEEN UPDATED.  MUST THEREFORE
!     WRITE IT OUT BEFORE READING IN A NEW BLOCK.
 
 CALL sofio (iwrt,mdipbn,buf(mdi-2))
 mdiup = .false.
 
!     THE DESIRED MDI BLOCK IS NOT PRESENTLY IN CORE, MUST THEREFORE
!     FETCH IT.
 
 20 newblk = .false.
 
!     FIND THE PHYSICAL BLOCK NUMBER OF THE BLOCK ON WHICH THE LOGICAL
!     BLOCK IBLOCK IS STORED.
 
 k = mdibl
 icount = 1
 30 IF (icount == iblock) GO TO 35
 icount = icount + 1
 CALL fnxt (k,nxtk)
 IF (MOD(k,2) == 1) GO TO 32
 ibl = rshift(buf(nxtk),ihalf)
 GO TO 34
 32 ibl = andf(buf(nxtk),jhalf)
 34 IF (ibl == 0) GO TO 60
 k = ibl
 GO TO 30
 35 IF (mdipbn == k) GO TO 500
 
!     READ THE DESIRED MDI BLOCK INTO CORE.
 
 mdipbn = k
 mdilbn = iblock
 IF (newblk) RETURN
 CALL sofio (ird,mdipbn,buf(mdi-2))
 RETURN
 
!     WE NEED A FREE BLOCK FOR THE MDI.
 
 60 CALL getblk (k,ibl)
 IF (ibl == -1) GO TO 1000
 newblk = .true.
 k   = ibl
 MIN = mdi + 1
 MAX = mdi + blksiz
 DO  ll = MIN,MAX
   buf(ll) = 0
 END DO
 CALL sofio (iwrt,k,buf(mdi-2))
 GO TO 30
 
!     ERROR IN UPDATING EITHER MDIPBN OR MDILBN.
 
 500 CALL errmkn (indsbr,6)
 
!     ERROR MESSAGES.
 
 1000 WRITE  (nout,1001) ufm
 1001 FORMAT (a23,' 6223, SUBROUTINE FMDI - THERE ARE NO MORE FREE ',  &
     'BLOCKS AVAILABLE ON THE SOF.')
 CALL sofcls
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE fmdi
