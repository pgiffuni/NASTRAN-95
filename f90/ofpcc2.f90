SUBROUTINE ofpcc2 (ix, l1, l2, l3, l4, l5, ipoint)
!*****
!     SETS HEADER LINE FORMATS FOR COMPLEX ELEMENT STRESSES IN
!     MATERIAL COORDINATE SYSTEM  --  SORT 2 OUTPUT
!*****
 
 INTEGER, INTENT(OUT)                     :: ix
 INTEGER, INTENT(OUT)                     :: l1
 INTEGER, INTENT(OUT)                     :: l2
 INTEGER, INTENT(OUT)                     :: l3
 INTEGER, INTENT(OUT)                     :: l4
 INTEGER, INTENT(OUT)                     :: l5
 INTEGER, INTENT(IN OUT)                  :: ipoint
 DIMENSION idata(48)
 
 DATA idata/4003,108, 139, 125, 0, 433, 4031,108, 139, 126, 0, 433,  &
     4003,108, 140, 125, 0, 433, 4031,108, 140, 126, 0, 433,  &
     4003,108, 135, 125, 0, 433, 4031,108, 135, 126, 0, 433,  &
     4003,108, 134, 125, 0, 433, 4031,108, 134, 126, 0, 433/
 
 ix = idata(ipoint  )
 l1 = idata(ipoint+1)
 l2 = idata(ipoint+2)
 l3 = idata(ipoint+3)
 l4 = idata(ipoint+4)
 l5 = idata(ipoint+5)
 
 RETURN
END SUBROUTINE ofpcc2
