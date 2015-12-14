SUBROUTINE makmcb (mcb,ifile,irow,IF,it)
     
 
 
 INTEGER, INTENT(OUT)                     :: mcb(7)
 INTEGER, INTENT(IN)                      :: ifile
 INTEGER, INTENT(IN)                      :: irow
 INTEGER, INTENT(IN)                      :: IF
 INTEGER, INTENT(IN)                      :: it
 
 
 mcb(1) = ifile
 mcb(2) = 0
 mcb(3) = irow
 mcb(4) = IF
 mcb(5) = it
 mcb(6) = 0
 mcb(7) = 0
 
 RETURN
 
END SUBROUTINE makmcb
