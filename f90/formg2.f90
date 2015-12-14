SUBROUTINE formg2(ig,jr,jd,ir,id)
     
!     FORMGG FORMS THE GG MATRIX FOR EACH RIGID ELEMENT DEGREE OF
!     FREEDOM.  IG IS THE START OF THE ROW STORED GG MATRIX - 1
!     JR IS THE START OF THE TA MATRIX - 1.
!     JD IS THE START OF THE TB MATRIX - 1.
!     IR IS THE START OF THE BGPDT INFORMATION FOR REFERENCE POINT
!     ID IS THE START OF THE BGPDT INFORMATION FOR DEPENDENT POINT
 
 
 INTEGER, INTENT(IN OUT)                  :: ig
 INTEGER, INTENT(IN OUT)                  :: jr
 INTEGER, INTENT(IN OUT)                  :: jd
 INTEGER, INTENT(IN OUT)                  :: ir
 INTEGER, INTENT(IN OUT)                  :: id
 DOUBLE PRECISION :: xd,yd,zd,zz(1)
 DIMENSION zr(1)
 INTEGER :: z
 COMMON/zzzzzz/z(1)
 EQUIVALENCE (zz(1),zr(1))
 EQUIVALENCE (zz(1),z(1))
 
!     CALCULATE THE X,Y,AND Z DIRECTED DISTANCES WITH RESPECT TO THE
!     REFERENCE GRID POINT
 
 xd = zr(id+1) - zr(ir+1)
 yd = zr(id+2) - zr(ir+2)
 zd = zr(id+3) - zr(ir+3)
 
!     IF NO TRANSFORMATION IS NECESSARY, GO TO 30
 
 IF (z(ir) == 0.AND.z(id) == 0) GO TO 30
 
!     IF ONLY DEPENDENT GRID POINT HAS A TRANSFORMATION, GO TO 20
 
 IF (z(ir) == 0) GO TO 20
 
!     IF BOTH HAVE TRANSFORMATIONS, GO TO 10
 
 
 IF (z(id) /= 0) GO TO 10
 
!     ONLY REFERENCE GRID POINT HAS A TRANSFORMATION
 
 zz(ig+ 1) = zz(jr+1)
 zz(ig+ 2) = zz(jr+2)
 zz(ig+ 3) = zz(jr+3)
 zz(ig+ 4) =zd * zz(jr+4) - yd * zz(jr+7)
 zz(ig+ 5) =zd * zz(jr+5) - yd * zz(jr+8)
 zz(ig+ 6) =zd * zz(jr+6) - yd * zz(jr+9)
 zz(ig+ 7) = zz(jr+4)
 zz(ig+ 8) = zz(jr+5)
 zz(ig+ 9) = zz(jr+6)
 zz(ig+10) =xd * zz(jr+7) - zd * zz(jr+1)
 zz(ig+11) =xd * zz(jr+8) - zd * zz(jr+2)
 zz(ig+12) =xd * zz(jr+9) - zd * zz(jr+3)
 zz(ig+13) = zz(jr+7)
 zz(ig+14) = zz(jr+8)
 zz(ig+15) = zz(jr+9)
 zz(ig+16) =yd * zz(jr+1) - xd * zz(jr+4)
 zz(ig+17) =yd * zz(jr+2) - xd * zz(jr+5)
 zz(ig+18) =yd * zz(jr+3) - xd * zz(jr+6)
 zz(ig+19) = 0.0
 zz(ig+20) = 0.0
 zz(ig+21) = 0.0
 zz(ig+22) = zz(ig+ 1)
 zz(ig+23) = zz(ig+ 2)
 zz(ig+24) = zz(ig+ 3)
 zz(ig+25) = 0.0
 zz(ig+26) = 0.0
 zz(ig+27) = 0.0
 zz(ig+28) = zz(ig+ 7)
 zz(ig+29) = zz(ig+ 8)
 zz(ig+30) = zz(ig+ 9)
 zz(ig+31) = 0.0
 zz(ig+32) = 0.0
 zz(ig+33) = 0.0
 zz(ig+34) = zz(ig+13)
 zz(ig+35) = zz(ig+14)
 zz(ig+36) = zz(ig+15)
 RETURN
 10 CONTINUE
 
!     BOTH HAVE TRANSFORMATIONS
 
 zz(ig+ 1) = zz(jd+1)*zz(jr+1) + zz(jd+4)*zz(jr+4) + zz(jd+7)*zz(jr+7)
 zz(ig+ 2) = zz(jd+1)*zz(jr+2) + zz(jd+4)*zz(jr+5) + zz(jd+7)*zz(jr+8)
 zz(ig+ 3) = zz(jd+1)*zz(jr+3) + zz(jd+4)*zz(jr+6) + zz(jd+7)*zz(jr+9)
 zz(ig+ 4) = zz(jd+1)*zd*zz(jr+4)-zz(jd+1)*yd*zz(jr+7) +  &
     zz(jd+4)*xd*zz(jr+7) - zz(jd+4)*zd*zz(jr+1) +  &
     zz(jd+7)*yd*zz(jr+1) - zz(jd+7)*xd*zz(jr+4)
 zz(ig+ 5) = zz(jd+1)*zd*zz(jr+5)-zz(jd+1)*yd*zz(jr+8) +  &
     zz(jd+4)*xd*zz(jr+8) - zz(jd+4)*zd*zz(jr+2) +  &
     zz(jd+7)*yd*zz(jr+2) - zz(jd+7)*xd*zz(jr+5)
 zz(ig+ 6) = zz(jd+1)*zd*zz(jr+6)-zz(jd+1)*yd*zz(jr+9) +  &
     zz(jd+4)*xd*zz(jr+9) - zz(jd+4)*zd*zz(jr+3) +  &
     zz(jd+7)*yd*zz(jr+3) - zz(jd+7)*xd*zz(jr+6)
 zz(ig+ 7) = zz(jd+2)*zz(jr+1) + zz(jd+5)*zz(jr+4) + zz(jd+8)*zz(jr+7)
 zz(ig+ 8) = zz(jd+2)*zz(jr+2) + zz(jd+5)*zz(jr+5) + zz(jd+8)*zz(jr+8)
 zz(ig+ 9) = zz(jd+2)*zz(jr+3) + zz(jd+5)*zz(jr+6) + zz(jd+8)*zz(jr+9)
 zz(ig+10) = zz(jd+2)*zd*zz(jr+4)-zz(jd+2)*yd*zz(jr+7) +  &
     zz(jd+5)*xd*zz(jr+7) - zz(jd+5)*zd*zz(jr+1) +  &
     zz(jd+8)*yd*zz(jr+1) - zz(jd+8)*xd*zz(jr+4)
 zz(ig+11) = zz(jd+2)*zd*zz(jr+5)-zz(jd+2)*yd*zz(jr+8) +  &
     zz(jd+5)*xd*zz(jr+8) - zz(jd+5)*zd*zz(jr+2) +  &
     zz(jd+8)*yd*zz(jr+2) - zz(jd+8)*xd*zz(jr+5)
 zz(ig+12) = zz(jd+2)*zd*zz(jr+6)-zz(jd+2)*yd*zz(jr+9) +  &
     zz(jd+5)*xd*zz(jr+9) - zz(jd+5)*zd*zz(jr+3) +  &
     zz(jd+8)*yd*zz(jr+3) - zz(jd+8)*xd*zz(jr+6)
 zz(ig+13) = zz(jd+3)*zz(jr+1) + zz(jd+6)*zz(jr+4) + zz(jd+9)*zz(jr+7)
 zz(ig+14) = zz(jd+3)*zz(jr+2) + zz(jd+6)*zz(jr+5) + zz(jd+9)*zz(jr+8)
 zz(ig+15) = zz(jd+3)*zz(jr+3) + zz(jd+6)*zz(jr+6) + zz(jd+9)*zz(jr+9)
 zz(ig+16) = zz(jd+3)*zd*zz(jr+4)-zz(jd+3)*yd*zz(jr+7) +  &
     zz(jd+6)*xd*zz(jr+7) - zz(jd+6)*zd*zz(jr+1) +  &
     zz(jd+9)*yd*zz(jr+1) - zz(jd+9)*xd*zz(jr+4)
 zz(ig+17) = zz(jd+3)*zd*zz(jr+5)-zz(jd+3)*yd*zz(jr+8) +  &
     zz(jd+6)*xd*zz(jr+8) - zz(jd+6)*zd*zz(jr+2) +  &
     zz(jd+9)*yd*zz(jr+2) - zz(jd+9)*xd*zz(jr+5)
 zz(ig+18) = zz(jd+3)*zd*zz(jr+6)-zz(jd+3)*yd*zz(jr+9) +  &
     zz(jd+6)*xd*zz(jr+9) - zz(jd+6)*zd*zz(jr+3) +  &
     zz(jd+9)*yd*zz(jr+3) - zz(jd+9)*xd*zz(jr+6)
 zz(ig+19) = 0.0
 zz(ig+20) = 0.0
 zz(ig+21) = 0.0
 zz(ig+22) = zz(ig+ 1)
 zz(ig+23) = zz(ig+ 2)
 zz(ig+24) = zz(ig+ 3)
 zz(ig+25) = 0.0
 zz(ig+26) = 0.0
 zz(ig+27) = 0.0
 zz(ig+28) = zz(ig+ 7)
 zz(ig+29) = zz(ig+ 8)
 zz(ig+30) = zz(ig+ 9)
 zz(ig+31) = 0.0
 zz(ig+32) = 0.0
 zz(ig+33) = 0.0
 zz(ig+34) = zz(ig+13)
 zz(ig+35) = zz(ig+14)
 zz(ig+36) = zz(ig+15)
 RETURN
 20 CONTINUE
 
!     DEPENDENT GRID POINT HAS TRANSFORMATION
 
 zz(ig+ 1) = zz(jd+1)
 zz(ig+ 2) = zz(jd+4)
 zz(ig+ 3) = zz(jd+7)
 zz(ig+ 4) = zz(jd+7)*yd - zz(jd+4)*zd
 zz(ig+ 5) = zz(jd+1)*zd - zz(jd+7)*xd
 zz(ig+ 6) = zz(jd+4)*xd - zz(jd+1)*yd
 zz(ig+ 7) = zz(jd+2)
 zz(ig+ 8) = zz(jd+5)
 zz(ig+ 9) = zz(jd+8)
 zz(ig+10) = zz(jd+8)*yd - zz(jd+5)*zd
 zz(ig+11) = zz(jd+2)*zd - zz(jd+8)*xd
 zz(ig+12) = zz(jd+5)*xd - zz(jd+2)*yd
 zz(ig+13) = zz(jd+3)
 zz(ig+14) = zz(jd+6)
 zz(ig+15) = zz(jd+9)
 zz(ig+16) = zz(jd+9)*yd - zz(jd+6)*zd
 zz(ig+17) = zz(jd+3)*zd - zz(jd+9)*xd
 zz(ig+18) = zz(jd+6)*xd - zz(jd+3)*yd
 zz(ig+19) = 0.0
 zz(ig+20) = 0.0
 zz(ig+21) = 0.0
 zz(ig+22) = zz(ig+ 1)
 zz(ig+23) = zz(ig+ 2)
 zz(ig+24) = zz(ig+ 3)
 zz(ig+25) = 0.0
 zz(ig+26) = 0.0
 zz(ig+27) = 0.0
 zz(ig+28) = zz(ig+ 7)
 zz(ig+29) = zz(ig+ 8)
 zz(ig+30) = zz(ig+ 9)
 zz(ig+31) = 0.0
 zz(ig+32) = 0.0
 zz(ig+33) = 0.0
 zz(ig+34) = zz(ig+13)
 zz(ig+35) = zz(ig+14)
 zz(ig+36) = zz(ig+15)
 RETURN
 30 CONTINUE
 
!     NO TRANSFORMATIONS
 
 DO  i = 1,36
   zz(ig+i) = 0.0
 END DO
 zz(ig+ 1) = 1.0
 zz(ig+ 8) = 1.0
 zz(ig+15) = 1.0
 zz(ig+22) = 1.0
 zz(ig+29) = 1.0
 zz(ig+36) = 1.0
 zz(ig+ 5) =  zd
 zz(ig+ 6) = -yd
 zz(ig+10) = -zd
 zz(ig+12) =  xd
 zz(ig+16) =  yd
 zz(ig+17) = -xd
 RETURN
END SUBROUTINE formg2
