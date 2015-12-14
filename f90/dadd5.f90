SUBROUTINE dadd5
     
!     DMAP DRIVER FOR SADD (MATRIX ADD) ROUTINE
!     THE DMAP CALL FOR THIS MODULE IS
!     ADD5 A,B,C,D,E / X / V,N,P1 / V,N,P2 / V,N,P3 / V,N,P4 / V,N,P5 $
!     THE PARAMETERS ARE ALL COMPLEX SINGLE-PRECISION.
 
 DIMENSION       inx(5),amcbs(1),mc(5)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ibuf,nout
 COMMON /saddx / nomat,lcore,mcbs(67)
 COMMON /zzzzzz/ core(1)
 COMMON /BLANK / alpha(10)
 EQUIVALENCE     (mcbs(1),amcbs(1)),(mcbs(61),mc(1))
 DATA    inx   / 101,102,103,104,105 /, iout /201/
 
 lcore = korsz(core)
 
 DO  i = 1,67
   mcbs(i) = 0
 END DO
 
!     SETUP MATRIX CONTROL BLOCKS OF THE INPUT MATRICES
 
 i = 1
 k = 0
 
 mc(5) = 1
 DO  j = 1,5
   mcbs(i) = inx(j)
   CALL rdtrl (mcbs(i))
   
!     EXCLUDE NULL MATRICES FROM MCBS ARRAY
   
   IF (mcbs(i) <= 0) CYCLE
   
!     MOVE MULTIPLIERS TO MCBS ARRAY
   
   mcbs (i+7) = 1
   amcbs(i+8) = alpha(2*j-1)
   amcbs(i+9) = alpha(2*j)
   IF (amcbs(i+9) /= 0.0) mcbs(i+7) = 3
   
!     DETERMINE THE PRECISION AND TYPE OF THE OUTPUT MATRIX
   
   mc(5) = MAX0(mc(5),mcbs(i+4),mcbs(i+7))
   IF (mcbs(i+4) == 2) k = 1
   i = i + 12
 END DO
 
 mc(1) = iout
 nomat = i/12
 IF (nomat == 0) RETURN
 IF (nomat == 1) GO TO 60
 
!     CHECK TO ENSURE THAT THE MATRICES BEING ADDED ARE OF THE SAME
!     ORDER
 
 i = 14
 DO  j = 2, nomat
   IF (mcbs(2) == mcbs(i) .AND. mcbs(3) == mcbs(i+1)) GO TO 40
   WRITE  (nout,30) ufm
   30 FORMAT (a23,' 4149, ATTEMPT TO ADD MATRICES OF UNEQUAL ORDER IN ',  &
       'MODULE ADD5.')
   CALL mesage (-61,0,0)
   40 i = i + 12
 END DO
 60 mc(2) = mcbs(2)
 mc(3) = mcbs(3)
 mc(4) = mcbs(4)
 IF (mc(5) == 3 .AND. k /= 0) mc(5) = 4
 mc(5) = MIN0(4,mc(5))
 
!     ADD MATRICES
 
 CALL sadd   (core,core)
 CALL wrttrl (mc(1))
 RETURN
 
END SUBROUTINE dadd5
