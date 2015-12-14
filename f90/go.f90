FUNCTION GO ( r , etar , etal , ekm )
     
 
 REAL, INTENT(IN)                         :: r
 REAL, INTENT(IN)                         :: etar
 REAL, INTENT(IN)                         :: etal
 REAL, INTENT(IN)                         :: ekm
 DIMENSION  as(2) , c(2) , s(2) , s0(2)
 DIMENSION bsl(23)
 
 dbslj = 1.0E-10
 s(1)  =  etar
 s(2)  =  etal
 DO      i = 1 , 2
   IF ( ABS ( s(i) ) >= r )   GO TO  200
   s(i)  =  s(i) / r
   c(i)  =  SQRT ( 1.0 - s(i) ** 2 )
   as(i)  =  2.0 * ATAN ( s(i) / ( 1.0 + c(i) ) )
   s(i)  =  2.0 * s(i) * c(i)
   c(i)  =  2.0 * c(i) ** 2 - 1.0
   GO TO  300
   
   200  as(i)  =  SIGN ( 1.570796   , s(i) )
   s(i)  =  0.0
   
   300  s0(i)  =  0.0
 END DO
 
 GO  =  as(1) - as(2)
 IF ( ABS ( GO ) <= dbslj )   GO TO  700
 
 arg  =  ekm * r
 IF ( arg == 0.0 )   RETURN
 CALL mbbslj(arg,n,bsl)
 
 GO  =  bsl(1) * GO
 f  =  1.0
 fi  =  1.0
 DO      j = 2 , n
   GO  =  bsl(j) * ( s(1) - s(2) ) / fi - GO
   
   DO      i = 1 , 2
     s4  =  2.0 * s(i) * c(i) - s0(i)
     s0(i)  =  s(i)
     s(i)  =  s4
   END DO
   
   f  =  -f
   fi  =  fi + 1.0
 END DO
 
 IF ( f < 0.0 )   GO  =  -GO
 RETURN
 
 700  GO  =  0.0
 RETURN
END FUNCTION GO
