FUNCTION zj ( arg )
     
!     ZERO ORDER BESSEL FUNCTION OF FIRST KIND
 
 dbslj = 1.0E-10
 a  =  - ( arg / 2.0 ) ** 2
 zj  =  1.0
 pf  =  1.0
 an  =  1.0
 DO      i = 1 , 20
   an  =  an * a / pf ** 2
   pf  =  pf + 1.0
   IF ( ABS ( an ) <= dbslj )   RETURN
   zj  =  zj + an
 END DO
 RETURN
END FUNCTION zj
