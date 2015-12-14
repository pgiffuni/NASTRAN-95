SUBROUTINE encode( ii )
     
!     THIS SUBROUTINE CONVERTS THE DEGREE OF FREEDOM CODES AS GIVEN
!     IN BULK DATA FORM ( INTEGERS FROM 1-6 ) TO THE BIT PATTERN
!     USED IN SUBSTRUCTURE ANALYSIS.
 
 
 INTEGER, INTENT(IN OUT)                  :: ii
 DIMENSION idiv(6)
 DATA idiv/ 100000 , 10000 , 1000 , 100 , 10 , 1 /
 
 isum = 0
 DO  i=1,6
   j = ii/idiv(i)
   IF( j == 0 ) CYCLE
   isum = isum + 2 ** (j-1)
   ii = ii - j*idiv(i)
 END DO
 ii = isum
 RETURN
END SUBROUTINE encode
