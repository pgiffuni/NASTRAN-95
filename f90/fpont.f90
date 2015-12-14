SUBROUTINE fpont
     
!     DOES DIRECT,TPONT,FPONT,AND SCALAR LOADS
 
 INTEGER :: gpid,slt,pont(5),swload(2)
 DIMENSION       igpco(4,5),gpco1(3),gpco2(3),gpco3(3),gpco4(3),  &
     vect1(3),vect2(3),iord(5),vect(3),gridp(7)
 COMMON /loadx / lc,slt,bg,OLD,n(12),ifm
 COMMON /zzzzzz/ core(1)
 EQUIVALENCE     (gpid,gridp(2)), (gridp(4),ip1), (gridp(5),ip2),  &
     (gridp(6),ip3) , (gridp(7),ip4),  &
     (igpco(2,1),gpco1(1)), (igpco(2,2),gpco2(1)),  &
     (igpco(2,3),gpco3(1)), (igpco(2,4),gpco4(1)), (icosyt,gridp(3))
 DATA    swload/ 4HFPON,4HT   /
 
 nr = 6
 np = 5
 minus = 5
 10 CALL READ (*120,*130,slt,gridp(2),nr,0,flag)
 scale   = gridp(3)
 pont(1) = ip1
 pont(2) = ip2
 IF (np == 3) GO TO 20
 pont(3) = ip3
 pont(4) = ip4
 20 pont(np)= gpid
 CALL permut (pont(1),iord(1),np,OLD)
 DO  i = 1,np
   l = iord(i)
   CALL fndpnt (igpco(1,l),pont(l))
 END DO
 IF (np == 3) GO TO 50
 DO  i = 1,3
   vect1(i) = gpco2(i) - gpco1(i)
   vect2(i) = gpco4(i) - gpco3(i)
 END DO
 CALL cross (vect1(1),vect2(1),vect(1))
 GO TO 70
 50 DO  i = 1,3
   vect(i) = gpco2(i) - gpco1(i)
 END DO
 70 CALL norm (vect(1),xl)
 80 IF (igpco(1,np) == 0) THEN
   GO TO   100
 END IF
 90 CALL basglb (vect(1),vect(1),igpco(2,np),igpco(1,np))
 100 CALL fndsil (gpid)
 gpid = gpid + (ifm-minus)*3 - 1
 DO  i = 1,3
   in = gpid + i
   core(in) = core(in) + vect(i)*scale
 END DO
 GO TO 150
 120 n1 = -2
 GO TO 140
 130 n1 = -3
 140 iparm = slt
 CALL mesage (n1,iparm,swload)
 150 RETURN
 
 
 ENTRY tpont
!     ===========
 
!     TPONT PROCESSES FORCE1 AND MOMENT1 CARDS
 
 nr = 4
 np = 3
 minus = 3
 GO TO 10
 
 
 ENTRY DIRECT
!     ============
 
!     DIRECT PROCESSES FORCE+ MOMENT CARDS
 
 np = 1
 minus = 1
 CALL READ (*120,*130,slt,gridp(2),6,0,flag)
 DO  i = 1,3
   vect(i) = gridp(i+4)
 END DO
 CALL fndpnt (igpco(1,1),gpid)
 scale = gridp(4)
 IF (icosyt == igpco(1,np)) GO TO 100
 IF (icosyt == 0) THEN
   GO TO    80
 END IF
 180 CALL glbbas (vect(1),vect(1),igpco(2,1),icosyt)
 GO TO 80
 
 
 ENTRY sload
!     ===========
 
!     SLOAD PROCESSES SLOAD CARDS
 
 CALL READ (*120,*130,slt,gridp(2),2,0,flag)
 CALL fndsil (gpid)
 core(gpid) = core(gpid) + gridp(3)
 GO TO 150
END SUBROUTINE fpont
