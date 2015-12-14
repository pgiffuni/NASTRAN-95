SUBROUTINE ofrf2s(ix,l1,l2,l3,l4,l5,point)
!*****
!  SETS HEADER LINE FORMATS FOR REAL FORCES SORT2 - STATICS
!*****
 
 INTEGER, INTENT(OUT)                     :: ix
 INTEGER, INTENT(OUT)                     :: l1
 INTEGER, INTENT(OUT)                     :: l2
 INTEGER, INTENT(OUT)                     :: l3
 INTEGER, INTENT(OUT)                     :: l4
 INTEGER, INTENT(OUT)                     :: l5
 INTEGER, INTENT(IN OUT)                  :: point
 INTEGER :: c
 COMMON/ofpb7s/ c(10)
!*****
 ix = c(point)
 l1 = c(point+1)
 l2 = c(point+2)
 l3 = c(point+3)
 l4 = c(point+4)
 l5 = c(point+5)
 RETURN
END SUBROUTINE ofrf2s
