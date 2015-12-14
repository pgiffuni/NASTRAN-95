SUBROUTINE qhbdy
!*****
 
!  THIS ROUTINE APPLIES THE LOADS DUE TO A SELECTED HEAT FLUX
!  LOADING CONDITION.
 
!  DATA CARD IS...
 
!  QHBDY  SETID  FLAG  Q0  AF  G1  G2  G3  G4
!                ============================
!                ABOVE FIELDS AVAILABLE TO THIS ROUTINE ONLY.
!                GRIDS ARE IN INTERNAL NOTATION AT THIS POINT.
!*****
 INTEGER :: map(15)  ,card(7)  ,igrids(5),slt      ,OLD      ,bg
 INTEGER :: sils(4)  ,grids(4) ,order(4) ,subr(2)
 
 REAL :: bgpdt(4,4)         ,x(4)     ,y(4)     ,z(4)
 REAL :: data4(4) ,p(4)     ,r12(3)   ,r13(3)   ,length
 
 COMMON /condas/ consts(5)
 COMMON/loadx / lc, slt, bg, OLD, n(12), ifm
 COMMON/zzzzzz/ core(1)
 
 EQUIVALENCE ( consts(1) , pi     )
 EQUIVALENCE(x(1),bgpdt(1,2)), (y(1),bgpdt(1,3)), (z(1),bgpdt(1,4))
 EQUIVALENCE(iflag,card(1)), (q0,card(2)), (af,card(3))
 EQUIVALENCE(grids(1),card(4),sils(1))
 
 DATA map/ 1,2,3,  1,2,4,  2,3,1,  3,4,2,  4,1,3 /
 DATA igrids/ 1,2,2,3,4 /
 DATA subr/ 4HQHBD,4HY    /
!*****
!  READ AND PROCESS ONE QHBDY IMAGE PER CALL TO THIS ROUTINE.
!*****
 CALL READ(*902,*903,slt,card(1),7,0,flag)
 ngrids = igrids(iflag)
!*****
!  OBTAIN A GRID (INTERNAL) POINT SORT VECTOR SO AS TO CALL FOR BGPDT
!  DATA EFFICIENTLY.
!*****
 IF( ngrids <= 1 ) GO TO 35
 CALL permut( grids(1), order(1), ngrids, OLD )
 GO TO 38
 35 order(1) = 1
!*****
!  PICK UP BGPDT FOR THE 1 TO 4 POINTS AND OBTAIN THE SILS.
!*****
 38 DO  i = 1,ngrids
   l = order(i)
   CALL fndpnt( data4(1), grids(l) )
   bgpdt(l,1) = data4(1)
   bgpdt(l,2) = data4(2)
   bgpdt(l,3) = data4(3)
   bgpdt(l,4) = data4(4)
   CALL fndsil( grids(l) )
 END DO
!*****
!  ALL DATA IS AT HAND FOR LOAD CALCULATIONS
!*****
 af = af * q0
 SELECT CASE ( iflag )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 300
   CASE (    4)
     GO TO 400
   CASE (    5)
     GO TO 500
 END SELECT
!*****
!  IFLAG=1   A POINT...
!*****
 100 p(1) = af
 GO TO 700
!*****
!  IFLAG=2  A LINE...
!*****
 200 length = SQRT( (x(2)-x(1))**2 + (y(2)-y(1))**2 +  (z(2)-z(1))**2 )
 p(1) = af * length * 0.50E0
 p(2) = p(1)
 GO TO 700
!*****
!  IFLAG=3  A LINE OF REVOLUTION...
!*****
 300 fact = pi*q0*SQRT( (x(2)-x(1))**2 + (z(2)-z(1))**2 ) / 3.0E0
 p(1) = fact * (2.0E0*x(1) + x(2))
 p(2) = fact * (x(1) + 2.0E0*x(2))
 GO TO 700
!*****
!  IFLAG=4  A TRIANGLE...
!*****
 400 fact = q0 / 6.0E0
 imap = 1
 nmap = 3
 GO TO 600
!*****
!  IFLAG=5  A QUADRILATERAL...
!*****
 500 fact = q0 / 12.0E0
 imap = 4
 nmap = 15
!*****
!  MAP 1 OR 4 TRIANGLES INTO 3 OR 4 POINTS.
!*****
 600 p(1) = 0.0E0
 p(2) = 0.0E0
 p(3) = 0.0E0
 p(4) = 0.0E0
 DO  i = imap,nmap,3
   i1 = map(i)
   i2 = map(i+1)
   i3 = map(i+2)
   r12(1) = x(i2) - x(i1)
   r12(2) = y(i2) - y(i1)
   r12(3) = z(i2) - z(i1)
   r13(1) = x(i3) - x(i1)
   r13(2) = y(i3) - y(i1)
   r13(3) = z(i3) - z(i1)
   CALL saxb( r12(1), r13(1), r12(1) )
   factx= fact * SQRT( r12(1)**2  + r12(2)**2  + r12(3)**2 )
   p(i1) = p(i1) + factx
   p(i2) = p(i2) + factx
   p(i3) = p(i3) + factx
 END DO
!*****
!  LOAD VALUES COMPLETE.
!*****
 700 DO  i = 1,ngrids
   isil = sils(i)
   core(isil  ) = core(isil  ) + p(i)
 END DO
 RETURN
!*****
!  END OF FILE OR END OF RECORD HIT ERROR.
!*****
 902 CALL mesage(-2,slt,subr)
 903 CALL mesage(-3,slt,subr)
 GO TO 902
END SUBROUTINE qhbdy
