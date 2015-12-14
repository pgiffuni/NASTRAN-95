SUBROUTINE xychar (row,col,CHAR)
     
 
 INTEGER, INTENT(IN OUT)                  :: row
 INTEGER, INTENT(IN)                      :: col
 INTEGER, INTENT(IN)                      :: CHAR
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf,complf
 LOGICAL :: pass,exceed
 DIMENSION       mask(4)
 COMMON /machin/ mach
 COMMON /system/ dum(38),bperch,bperwd
 COMMON /xypppp/ iframe,titlec(32),titlel(14),titler(14),  &
     xtitle(32),id(300),maxplt,xmin,xinc,exceed, i123,maxrow
 COMMON /zzzzzz/ z(1)
 DATA    pass  / .false. /
 
 IF (row <= maxrow) GO TO 1
 exceed = .true.
 RETURN
 
 1 IF (col > 119 .OR. col < 1 .OR. row < 1) RETURN
 
!     CHAR COMING IN IS ASSUMED LEFT ADJUSTED
 
 IF (pass) GO TO 20
 pass = .true.
 
!     SET UP MASKS FIRST TIME THROUGH AFTER LOADING
 
 n = 2**bperch  -  1
 ishift = bperwd - bperch
 n = lshift(n,ishift)
 nmask = n
 DO  i = 1,4
   mask(i) = complf(n)
   n = rshift(n,bperch)
 END DO
 
!     COMPUTE WORD AND CHARACTER OF WORD
 
 20 iword = (col-1)/4 + 1
 ICHAR = col - (iword-1)*4
 iword = (row-1)*30 + iword
 
!     PACK THE CHARACTER
 
 IF (mach == 5 .OR. mach == 6 .OR. mach == 21) GO TO 30
 let = rshift(andf(CHAR,nmask),bperch*(ICHAR-1))
 z(iword) = orf(andf(z(iword),mask(ICHAR)),let)
 RETURN
 
!     VAX, ULTRIX, AND ALPHA
 
 30 z(iword) = khrfn1(z(iword),ICHAR,CHAR,1)
 RETURN
END SUBROUTINE xychar
