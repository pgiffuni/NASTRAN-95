BLOCK DATA ofssbd
!OFSSBD
!     C ARRAY FOR REAL STRESSES SORT1 (IN MATERIAL COORDINATES)
 
 INTEGER :: c1
 COMMON /ofss1/ c1(30)
 DATA c1  / 2542,  0, 70,390,376,377   ,2542,  0, 71,390,376,377  &
     ,2542,  0, 69,390,376,377   ,2542,  0, 68,390,376,377  &
     ,2542,  0,378,390,387,377   /
END SUBROUTINE ofsplt