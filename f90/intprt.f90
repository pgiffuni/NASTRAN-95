SUBROUTINE intprt (a,cr,o,NAME)
     
 
 REAL, INTENT(IN OUT)                     :: a(1)
 INTEGER, INTENT(IN OUT)                  :: cr
 INTEGER, INTENT(IN OUT)                  :: o
 REAL, INTENT(IN OUT)                     :: NAME(2)
 INTEGER :: colnum,crfmt(3),cropt(2,2)
 
 COMMON /system/ skip,mo
 DATA    crfmt / 4H(60X , 4H,2A4 , 4H,i5) /
 DATA    cropt / 4HCOLU , 4HMN   , 4HROW  , 4H     /
 
!     CR   = 0  IF MATRIX BY COLUMNS.
!          = 1  IF MATRIX BY ROWS.
!     IF O = 0, THE MATRIX WILL NOT BE PRINTED.
!     NAME = 8  CHARACTER BCD NAME OF THE MATRIX.
 
 IF (cr /= 0) GO TO 100
 icropt = 1
 GO TO 110
 100 icropt = 2
 
 110 CALL matprt (*120,*130,a,-1,colnum)
 GO TO 150
 120 WRITE  (mo,125) NAME(1),NAME(2)
 125 FORMAT (50X,24HINTERMEDIATE matrix ... ,2A4//)
 130 WRITE  (mo,crfmt) (cropt(i,icropt),i=1,2),colnum
 CALL prtmat (*120,*130,colnum)
 150 RETURN
 
END SUBROUTINE intprt
