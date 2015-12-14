SUBROUTINE setval
     
 EXTERNAL        andf,rshift
 INTEGER :: andf,rshift,p,oscar,vps,subnam(2)
 COMMON /BLANK / p(2,5)
 COMMON /system/ ksystm(65)
 COMMON /xvps  / vps(1)
 COMMON /oscent/ oscar(1)
 EQUIVALENCE     (ksystm(40),nbpw)
 DATA    subnam/ 4HSETV,4HAL   /
 
 j = 12
 DO  i = 1,5
   
!     CHECK ODD PARAMETERS TO FIND VARIABLE ONES
   
   IF (andf(rshift(oscar(j+1),nbpw-1),1) == 0) GO TO 200
   
!     PARAMETER IS VARIABLE
   
   k = andf(oscar(j+1),65535)
   p(1,i) = p(2,i)
   vps(k) = p(1,i)
   j = j + 2
   IF (andf(rshift(oscar(j),nbpw-1),1) == 0) j = j + 1
 END DO
 GO TO 500
 
 200 CONTINUE
 IF (i > 1) GO TO 500
 CALL mesage (-7,0,subnam)
 
 500 CONTINUE
 RETURN
END SUBROUTINE setval
