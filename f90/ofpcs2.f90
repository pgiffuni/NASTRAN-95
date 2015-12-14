SUBROUTINE ofpcs2(ix,l1,l2,l3,l4,l5,point)
!*****
!  SETS HEADER LINE FORMATS FOR COMPLEX STRESSES SORT2
!*****
 
 INTEGER, INTENT(OUT)                     :: ix
 INTEGER, INTENT(OUT)                     :: l1
 INTEGER, INTENT(OUT)                     :: l2
 INTEGER, INTENT(OUT)                     :: l3
 INTEGER, INTENT(OUT)                     :: l4
 INTEGER, INTENT(OUT)                     :: l5
 INTEGER, INTENT(IN OUT)                  :: point
 INTEGER :: c
 COMMON/ofpb4/ c(10)
 
 ix = c(point)
 l1 = c(point+1)
 l2 = c(point+2)
 l3 = c(point+3)
 l4 = c(point+4)
 l5 = c(point+5)
 RETURN
END SUBROUTINE ofpcs2
