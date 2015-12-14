SUBROUTINE sfetch (NAME,item,irw,itest)
     
!     POSITIONS THE SOF TO READ OR WRITE DATA ASSOCIATED WITH ITEM OF
!     SUBSTRUCTURE NAME.
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN OUT)                  :: item
 INTEGER, INTENT(IN OUT)                  :: irw
 INTEGER, INTENT(OUT)                     :: itest
 EXTERNAL        andf
 LOGICAL :: mdiup
 INTEGER :: andf,buf,mdi,mdipbn,mdilbn,mdibl,blksiz,dirsiz
 DIMENSION  nmsbr(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /zzzzzz/ buf(1)
 COMMON /sof   / ditdum(6),io,iopbn,iolbn,iomode,ioptr,iosind,  &
     ioitcd,ioblk,mdi,mdipbn,mdilbn,mdibl,nxtdum(15), ditup,mdiup
 COMMON /sys   / blksiz,dirsiz
 COMMON /system/ nbuff,nout
 DATA    idle  , ird,iwrt /0,1,2/, nmsbr /4HSFET,4HCH  /
 
 CALL chkopn (nmsbr(1))
 CALL fdsub  (NAME(1),iosind)
 IF (iosind == -1) GO TO 500
 ioitcd = itcode(item)
 IF (ioitcd == -1) GO TO 510
 
!     CHECK IF ITEM IS A TABLE ITEM UNLESS SPECIAL CALL FROM MTRXO OR
!     MTRXI
 
 IF (irw < 0) GO TO 10
 itm = ittype(item)
 IF (itm /= 0) GO TO 530
 10 CALL fmdi (iosind,imdi)
 iolbn = 1
 ioptr = io + 1
 ibl   = andf(buf(imdi+ioitcd),65535)
 irdwrt= IABS(irw)
 SELECT CASE ( irdwrt )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 80
   CASE (    3)
     GO TO 30
 END SELECT
 
!     READ OPERATION.
 
 30 IF (ibl ==     0) GO TO 50
 IF (ibl /= 65535) GO TO 60
 
!     ITEM WAS PSEUDO-WRITTEN.
 
 itest = 2
 GO TO 520
 
!     ITEM HAS NOT BEEN WRITTEN.
 
 50 itest = 3
 GO TO 520
 
!     UPDATE THE COMMON BLOCK SOF, AND BRING INTO CORE THE DESIRED BLOCK
 
 60 itest = 1
 IF (irdwrt == 3) GO TO 520
 iopbn  = ibl
 iomode = ird
 CALL sofio (ird,iopbn,buf(io-2))
 RETURN
 
!     WRITE OPERATION.
 
 80 IF (ibl == 0 .OR. ibl == 65535) GO TO 90
 
!     ITEM HAS ALREADY BEEN WRITTEN.
 
 itest = 1
 GO TO 520
 90 itest1 = itest - 1
 SELECT CASE ( itest1 )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 110
 END SELECT
 
!     ITEM IS TO BE PSEUDO-WRITTEN.
 
 100 buf(imdi+ioitcd) = 65535
 mdiup = .true.
 RETURN
 
!     ITEM IS TO BE WRITTEN.  GET A FREE BLOCK AND UPDATE THE COMMON
!     BLOCK SOF.
 
 110 CALL getblk (0,ioblk)
 IF (ioblk == -1) GO TO 1000
 iopbn  = ioblk
 iomode = iwrt
 RETURN
 
!     NAME DOES NOT EXIST.
 
 500 itest = 4
 GO TO 520
 
!     ITEM IS AN ILLEGAL ITEM NAME.
 
 510 itest  = 5
 520 iomode = idle
 RETURN
 
!     ATTEMPT TO OPERATE ON A MATRIX ITEM
 
 530 WRITE  (nout,540) sfm,item,NAME
 540 FORMAT (a25,' 6227, AN ATTEMPT HAS BEEN MADE TO OPERATE ON THE ',  &
     'MATRIX ITEM ',a4,' OF SUBSTRUCTURE ',2A4,' USING SFETCH.')
 GO TO 1010
 
!     NO MORE BLOCKS ON SOF
 
 1000 WRITE  (nout,1001) ufm
 1001 FORMAT (a23,' 6223, SUBROUTINE SFETCH - THERE ARE NO MORE FREE ',  &
     'BLOCKS AVAILABLE ON THE SOF.')
 1010 CALL sofcls
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE sfetch
