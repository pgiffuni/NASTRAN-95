SUBROUTINE stppt2(INPUT,w1jk,w2jk)
     
 INTEGER, INTENT(IN OUT)                  :: INPUT
 INTEGER, INTENT(IN OUT)                  :: w1jk
 INTEGER, INTENT(IN OUT)                  :: w2jk
 
 COMPLEX :: one,zero
 COMMON /packx/ iti,it0,ii,nn,incr
 COMMON /amgp2/ tw1jk(7),tw2jk(7)
 
 one = (1.0,0.0)
 zero = (0.0,0.0)
 CALL fread(INPUT,nj,1,1)
 DO  i=1,nj
   nn = ii
   CALL pack(one,w1jk,tw1jk)
   CALL pack(zero,w2jk,tw2jk)
   ii = ii+1
 END DO
 RETURN
END SUBROUTINE stppt2
