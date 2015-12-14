SUBROUTINE grav (ngrav,gvect,nlist,ilist,nloop)
     
 
 INTEGER, INTENT(OUT)                     :: ngrav
 REAL, INTENT(OUT)                        :: gvect(1)
 INTEGER, INTENT(IN)                      :: nlist
 INTEGER, INTENT(IN OUT)                  :: ilist(1)
 INTEGER, INTENT(IN)                      :: nloop
 INTEGER :: NAME(2)
 DIMENSION  gl(5),x(3)
 COMMON /tranx/ nsys,tysys,ro(3),TO(3,3)
 COMMON /loadx/ lcore,slt,n(14),nobld
 EQUIVALENCE    (igl,gl(2))
 DATA    NAME / 4HGRAV,4H    /
 
!     CONVERTS GRAV CARD TO BASIC AND STORES
!     GB = G*TON*V
 
 CALL READ (*30,*40,slt,gl(1),5,0,flag)
 GO TO 50
 20 RETURN
 
 30 CONTINUE
 40 CALL mesage (-7,NAME,NAME)
 50 ngrav = ngrav + 1
 IF (gl(1) == 0.0) THEN
   GO TO    70
 END IF
 60 CALL fdcstm (gl(1))
 CALL mpyl (TO,gl(3),3,3,1,x(1))
 DO  i = 1,3
   gl(i+2) = x(i)
 END DO
 70 DO  i = 1,3
   j = (ngrav-1)*3 + i
   gvect(j) = gl(i+2)*gl(2)
 END DO
 nl1 = nloop - ngrav + 1
 IF (nl1 == nlist) GO TO 20
 nsave  = ilist(nl1)
 nlist1 = nlist - 1
 DO  i = nl1,nlist1
   ilist(i) = ilist(i+1)
 END DO
 ilist(nlist) = nsave
 GO TO 20
END SUBROUTINE grav
