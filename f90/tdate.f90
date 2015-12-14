SUBROUTINE tdate (date)
     
!     VAX VERSION
!     ===========
!     (ALSO SiliconGraphics, DEC/ultrix, and SUN.
!      CRAY AND HP DO NOT HAVE IDATE)
 
!     THIS ROUTINE OBTAINS THE MONTH, DAY AND YEAR, IN INTEGER FORMAT
 
 
 
 INTEGER, INTENT(OUT)                     :: date(3)
 INTEGER :: date1(3)
 
 CALL idate (date1)
!                 DAY   MONTH     YEAR
!     THESE DATES HAD TO BE INTERCHANGED FOR THE SUN
 date(1)=date1(2)
 date(2)=date1(1)
!WUT  BAD FIX
 date(3)=date1(3)-2000
 RETURN
END SUBROUTINE tdate
