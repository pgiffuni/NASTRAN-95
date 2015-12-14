SUBROUTINE DELETE (NAME,itemx,itest)
     
!     DELETES ITEM WHICH BELONGS TO THE SUBSTRUCTURE NAME.  THE MDI IS
!     UPDATED ACCORDINGLY AND THE BLOCKS ON WHICH ITEM WAS WRITTEN ARE
!     RETURNED TO THE LIST OF FREE BLOCKS.  ITEST IS AN OUTPUT PARAMETER
!     WHICH TAKES ON ONE OF THE FOLLOWING VALUES
 
!              1  IF ITEM DOES EXIST
!              2  IF ITEM PSEUDO-EXISTS
!              3  IF ITEM DOES NOT EXIST
!              4  IF NAME DOES NOT EXIST
!              5  IF ITEM IS AN ILLEGAL ITEM NAME
 
!     THE BLOCKS OCCUPIED BY THE ITEM ARE RETURNED TO THE LIST OF FREE
!     BLOCKS IF THEY BELONG TO THE SPECIFIED SUBSTRUCTURE
 
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN OUT)                  :: itemx
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        rshift,andf
 LOGICAL :: mdiup
 INTEGER :: buf,mdi,mdipbn,mdilbn,mdibl,blksiz,dirsiz,ps,ss, andf,rshift
 DIMENSION  nmsbr(2)
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / ditdum(6),iodum(8),mdi,mdipbn,mdilbn,mdibl,  &
     nxtdum(15),ditup,mdiup
 COMMON /sys   / blksiz,dirsiz,sys(3),ifrst
 COMMON /itemdt/ nitem,item(7,1)
 DATA    is,ps , ss/ 1,1,1    /
 DATA    nmsbr / 4HDELE,4HTE  /
 
 CALL chkopn (nmsbr(1))
 CALL fdsub  (NAME(1),k)
 IF (k == -1) GO TO 500
 CALL fmdi (k,imdi)
 ii  = itcode(itemx)
 IF (ii == -1) GO TO 510
 itm = ii - ifrst + 1
 ibl = andf(buf(imdi+ii),65535)
!                             55535 = 2**16 - 1
 IF (ibl /= 0) GO TO 10
 
!     ITEM DOES NOT EXIST.
 
 itest = 3
 RETURN
 
 10 buf(imdi+ii) = 0
 mdiup = .true.
 IF (ibl /= 65535) GO TO 20
 
!     ITEM PSEUDO-EXISTS.
 
 itest = 2
 GO TO 30
 
!     ITEM DOES EXIST.
 
 20 itest = 1
 30 IF (andf(buf(imdi+is),1073741824) == 0) GO TO 35
!                           1073741824 = 2**30
 
!     IMAGE SUBSTRUCTURE
 
 IF (itest /= 1) RETURN
 IF (item(4,itm) == 0) GO TO 32
 CALL retblk (ibl)
 32 RETURN
 
!     NAME IS A SECONDARY OR A PRIMARY SUBSTRUCTURE
 
 35 isvps = andf(buf(imdi+ps),1023)
!                               1023 = 2**10 - 1
 IF (isvps == 0) GO TO 39
 
!     SECONDARY SUBSTRUCTURE
 
 IF (itest /= 1) RETURN
 IF (item(5,itm) == 0) GO TO 37
 CALL retblk (ibl)
 37 RETURN
 
!     PRIMARY SUBSTRUCTURE
 
 39 IF (itest == 1) CALL retblk (ibl)
 40 isvss = rshift(andf(buf(imdi+ss),1048575),10)
!                                      1048575 = 2*20 - 1
 IF (isvss == 0) RETURN
 CALL fmdi (isvss,imdi)
 IF (andf(buf(imdi+ii),65535) /= ibl) GO TO 40
 buf(imdi+ii) = 0
 mdiup = .true.
 GO TO 40
 
!     NAME DOES NOT EXIST.
 
 500 itest = 4
 RETURN
 
!     ITEM IS AN ILLEGAL ITEM NAME.
 
 510 itest = 5
 RETURN
END SUBROUTINE DELETE
