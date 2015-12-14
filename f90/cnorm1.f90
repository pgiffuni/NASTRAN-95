SUBROUTINE cnorm1 (x,n)
     
!     CNORM1 WILL SEARCH A VECTOR FOR THE LARGEST VALUE AND NORMALIZE
!     THE VECTOR TO LARGEST ELEMENT EQUAL TO ONE
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: x(1)
 INTEGER, INTENT(IN)                      :: n
 INTEGER :: NAME(2)
 DOUBLE PRECISION :: dum,MAX,div(2)
 COMMON /system/  ibuf,nout
 DATA    NAME  /  4HCNOR,4HM1   /
 
 nn   = n + n
 MAX  = 0.d0
 INDEX= 0
 DO  i = 1,nn,2
   dum = x(i)*x(i) + x(i+1)*x(i+1)
   IF (dum <= MAX) CYCLE
   MAX   = dum
   INDEX = i
 END DO
 IF (INDEX == 0) GO TO 30
 div(1) = x(INDEX  )
 div(2) = x(INDEX+1)
 MAX    = div(1)*div(1) + div(2)*div(2)
 DO  i = 1,nn,2
   dum    = (x(i)*div(1) + x(i+1)*div(2))/MAX
   x(i+1) = (x(i+1)*div(1) - x(i)*div(2))/MAX
   x(i)   = dum
 END DO
 RETURN
 
 30 WRITE  (nout,40)
 40 FORMAT (//5X,37HCNORM1 received a vector of all zeros)
 CALL mesage (-37,0,NAME)
 RETURN
END SUBROUTINE cnorm1
