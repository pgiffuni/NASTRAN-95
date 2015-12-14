SUBROUTINE kpltst (g1,g2,g3,g4)
     
!     THIS ROUTINE WILL VERIFY THAT THE 4 GRID POINTS IN 3 SPACE LIE IN
!     AN APPROXIMATE PLANE. IF NOT THE NOGO FLAG IS SET TRUE AND A
!     MESSAGE IS WRITEN.
 
 
 REAL, INTENT(IN)                         :: g1(3)
 REAL, INTENT(IN)                         :: g2(3)
 REAL, INTENT(IN)                         :: g3(3)
 REAL, INTENT(IN)                         :: g4(3)
 LOGICAL :: nogo
 INTEGER :: out
 
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ sysbuf,out,nogo
 COMMON /sma1et/ id
 COMMON /sma1dp/ r13(3),r24(3),rxr(3),r(3)
 
 r13(1) = g3(1) - g1(1)
 r13(2) = g3(2) - g1(2)
 r13(3) = g3(3) - g1(3)
 r24(1) = g4(1) - g2(1)
 r24(2) = g4(2) - g2(2)
 r24(3) = g4(3) - g2(3)
 CALL saxb (r13,r24,rxr)
 
!     NORMALIZE
 
 dl     = SQRT(rxr(1)**2 + rxr(2)**2 + rxr(3)**2)
 IF (dl > 0.0) THEN
   GO TO    10
 ELSE
   GO TO    20
 END IF
 10 rxr(1) = rxr(1)/dl
 rxr(2) = rxr(2)/dl
 rxr(3) = rxr(3)/dl
 r1l    = SQRT(r13(1)**2 + r13(2)**2 + r13(3)**2)
 r2l    = SQRT(r24(1)**2 + r24(2)**2 + r24(3)**2)
 dl     = AMIN1(r1l,r2l)
 r(1)   = g2(1) - g1(1)
 r(2)   = g2(2) - g1(2)
 r(3)   = g2(3) - g1(3)
 dh     = sadotb(r,rxr)
 IF (dl > 0.0) THEN
   GO TO    15
 ELSE
   GO TO    20
 END IF
 15 IF (ABS(dh/dl) <= 0.10) RETURN
 
!     NOT PLANER
 
 20 CALL page2 (-2)
 WRITE  (out,30) uwm,id
 30 FORMAT (a25,' 4000, ONE SIDE OF ELEMENT',i10,  &
     ' CONNECTING FOUR POINTS IS NOT APPROXIMATELY PLANER.')
 RETURN
END SUBROUTINE kpltst
