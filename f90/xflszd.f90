SUBROUTINE xflszd (FILE,iblock,filnam)
     
!     XFLSZD (EXECUTIVE FILE SIZE DETERMINATOR) ACCUMULATES THE
!     NUMBER OF BLOCKS USED FOR A FILE (FILE LT 0) IN THE FIAT OR
!     FOR A FILE (FILE GT 0) IN THE DATA POOL FILE.
!     IF FILE GT 0 IT IS THE INDEX OF THE FILE ON THE DATA POOL FILE
!     IF FILE = 0 THE NUMBER OF WORDS PER BLOCK IS RETURNED IN IBLOCK
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: iblock
 INTEGER, INTENT(IN OUT)                  :: filnam
 EXTERNAL         rshift,andf
 COMMON / machin/ mach
 COMMON / xfiat / fiat(1)
 COMMON / xfist / nfist,lfist,ifist(1)
 COMMON / xdpl  / pool(1)
 COMMON / system/ kystem
 
 DATA     mask  / 32767 /
 
 IF (FILE) 10,150,100
 
!     FILE IS IN THE FIAT
 
!     COMMENTS FROM G.CHAN/UNIVAC 8/90
!     VAX AND VAX-DERIVED MACHINES DO NOT SAVE ANY INFORMATION OF BLOCKS
!     USED IN FIAT 7TH AND 8TH WORDS. THEREFORE, IBLOCK IS ALWAYS ZERO.
 
 10 CONTINUE
 
 lim = 2*lfist
 DO  i = 1,lim,2
   IF (filnam /= ifist(i)) CYCLE
   IF (ifist(i+1)  <=   0) EXIT
   indx   = ifist(i+1)
   iblock = rshift(fiat(indx+7),16) + andf(mask,fiat(indx+8)) +  &
       rshift(fiat(indx+8),16)
!            = BLOCK COUNT ON PRIMARY, SECONDARY AND TERTIARY FILES ??
   GO TO 200
 END DO
 50 iblock = 0
 GO TO 200
 
!     FILE IS ON THE DATA POOL FILE
 
 100 indx   = FILE*3 + 3
 iblock = rshift(pool(indx),16)
 GO TO 200
 
!     USER WANTS THE NUMBER OF WORDS PER BLOCK
 
 150 CONTINUE
 IF (mach == 2 .OR. mach >= 5) iblock = kystem - 4
 
 200 RETURN
END SUBROUTINE xflszd
