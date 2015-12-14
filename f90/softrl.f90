SUBROUTINE softrl (NAME,item,mcb)
     
!     UTILITY SUBROUTINE TO OBTAIN THE MATRIX TRAILER FOR A MATRIX
!     STORED ON THE SOF
!     STATUS OF THE SOF ITEM IS RETURNED IN WORD ONE OF THE MATRIX
!     CONTROL BLOCK
 
!         1 NORMAL RETURN - THE TRAILER IS STORED IN WORDS 2 THRU 7
!         2 ITEM WAS PESUDO WRITTEN
!         3 ITEM DOES NOT EXIST
!         4 SUBSTRUCTURE NAME DOES NOT EXIST
!         5 ILLEGAL ITEM NAME
 
 
 INTEGER, INTENT(IN OUT)                  :: NAME(2)
 INTEGER, INTENT(IN OUT)                  :: item
 INTEGER, INTENT(OUT)                     :: mcb(7)
 EXTERNAL        andf,rshift
 INTEGER :: andf,rshift,buf,blksiz,nmsbr(2)
 COMMON /machin/ mach,ihalf,jhalf
 COMMON /sof   / ditdum(6),io,iopbn,iolbn,iomode,ioptr,iosind, ioitcd,ioblk
 COMMON /sys   / blksiz
 COMMON /zzzzzz/ buf(1)
 DATA    ird   / 1/
 DATA    nmsbr / 4HSOFT,4HRL  /
 
 
!     CHECK IF ITEM IS ONE OF THE FOLLOWING ALLOWABLE NAMES.
!     KMTX,MMTX,PVEC,POVE,UPRT,HORG,UVEC,QVEC,PAPP,POAP,LMTX
 
 CALL chkopn (nmsbr(1))
 ioitcd = itcode(item)
 itm = ittype(item)
 IF (itm /= 1) GO TO 1000
 
!     FIND SUBSTRUCTURE NAME AND MDI BLOCK
 
 CALL fdsub (NAME,iosind)
 IF (iosind < 0) GO TO 1010
 CALL fmdi (iosind,imdi)
 
!     GET BLOCK NUMBER OF FIRST BLOCK
 
 iopbn = andf(buf(imdi+ioitcd),jhalf)
 IF (iopbn ==     0) GO TO 1020
 IF (iopbn == jhalf) GO TO 1030
 iolbn = 1
 
!     GET NEXT BLOCK IN CHAIN
 
 20 CALL fnxt (iopbn,inxt)
 IF (MOD(iopbn,2) == 1) GO TO 30
 next = andf(rshift(buf(inxt),ihalf),jhalf)
 GO TO 40
 30 next = andf(buf(inxt),jhalf)
 40 CONTINUE
 IF (next == 0) GO TO 50
 iopbn = next
 iolbn = iolbn + 1
 GO TO 20
 
!     WE HAVE HIT END OF CHAIN - READ THE LAST BLOCK
 
 50 CALL sofio (ird,iopbn,buf(io-2))
 i1 = io - 2
 i2 = i1 + blksiz + 4
 
!     EXTRACT TRAILER FROM BLOCK
 
 DO  i = 1,6
   mcb(i+1) = buf(io+blksiz-6+i)
 END DO
 mcb(1  ) = 1
 RETURN
 
 
!     ERRORS
 
!     ILLEGAL ITEM
 
 1000 mcb(1) = 5
 RETURN
 
!     SUBSTRUCTURE DOES NOT EXIST
 
 1010 mcb(1) = 4
 RETURN
 
!     ITEM DOES NOT EXIST
 
 1020 mcb(1) = 3
 RETURN
 
!     ITEM IS PESUDO WRITTEN
 
 1030 mcb(1) = 2
 RETURN
END SUBROUTINE softrl
