SUBROUTINE ofpcc1 (ix, l1, l2, l3, l4, l5, ipoint)
!*****
!     SETS HEADER LINE FORMATS FOR COMPLEX ELEMENT STRESSES IN
!     MATERIAL COORDINATE SYSTEM  --  SORT 1 OUTPUT
!*****
 
 INTEGER, INTENT(OUT)                     :: ix
 INTEGER, INTENT(OUT)                     :: l1
 INTEGER, INTENT(OUT)                     :: l2
 INTEGER, INTENT(OUT)                     :: l3
 INTEGER, INTENT(OUT)                     :: l4
 INTEGER, INTENT(OUT)                     :: l5
 INTEGER, INTENT(IN OUT)                  :: ipoint
 DIMENSION idata(48)
 
 DATA idata/3951,104, 139, 125, 0, 432, 3977,104, 139, 126, 0, 432,  &
     3951,104, 140, 125, 0, 432, 3977,104, 140, 126, 0, 432,  &
     3951,104, 135, 125, 0, 432, 3977,104, 135, 126, 0, 432,  &
     3951,104, 134, 125, 0, 432, 3977,104, 134, 126, 0, 432/
 
 ix = idata(ipoint  )
 l1 = idata(ipoint+1)
 l2 = idata(ipoint+2)
 l3 = idata(ipoint+3)
 l4 = idata(ipoint+4)
 l5 = idata(ipoint+5)
 
 RETURN
END SUBROUTINE ofpcc1
