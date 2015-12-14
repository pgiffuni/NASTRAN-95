SUBROUTINE jacobs (elid,shp,dshp,gpth,bgpdt,gpnorm,jacob)
     
!     THIS SUBROUTINE CALCULATES JACOBIAN AT EACH GIVEN INTEGRATION
!     POINT FOR QUAD4 POTVIN TYPE ELEMENTS.
 
!     SINGLE PRECISION VERSION
 
 
 INTEGER, INTENT(IN OUT)                  :: elid
 REAL, INTENT(IN)                         :: shp(1)
 REAL, INTENT(IN)                         :: dshp(1)
 REAL, INTENT(IN)                         :: gpth(1)
 REAL, INTENT(IN)                         :: bgpdt(4,1)
 REAL, INTENT(IN)                         :: gpnorm(4,1)
 REAL, INTENT(OUT)                        :: jacob(3,3)
 LOGICAL :: badj
 INTEGER :: INDEX(3,3), nogo,nout
 
 REAL :: psitrn(9), tgrid(3,8),sk(3),tk(3),enk(3),v1(3),v2(3),v3(3),  &
     val,hzta,thick ,detj,dum(3),eps
 COMMON /q4dt  / detj,hzta,psitrn,nnode,badj,n1
 COMMON /system/ ibuf,nout,nogo
 
 EQUIVALENCE     (psitrn(1),v1(1))
 EQUIVALENCE     (psitrn(4),v2(1))
 EQUIVALENCE     (psitrn(7),v3(1))
 
 DATA   eps    / 1.0E-15 /
 
!     INITIALIZE BADJ LOGICAL
 
 badj=.false.
 
!     COMPUTE THE JACOBIAN AT THIS GAUSS POINT,
!     ITS INVERSE AND ITS DETERMINANT.
 
 DO  i=1,nnode
   thick=gpth(i)
   tgrid(1,i)=bgpdt(2,i)+hzta*thick*gpnorm(2,i)
   tgrid(2,i)=bgpdt(3,i)+hzta*thick*gpnorm(3,i)
   tgrid(3,i)=bgpdt(4,i)+hzta*thick*gpnorm(4,i)
 END DO
 DO  i=1,2
   ipoint=n1*(i-1)
   DO  j=1,3
     jacob(i,j)=0.0
     DO  k=1,nnode
       jacob(i,j)=jacob(i,j)+dshp(k+ipoint)*tgrid(j,k)
     END DO
   END DO
 END DO
 DO  j=1,3
   jacob(3,j)=0.0
   DO  k=1,nnode
     jtemp=j+1
     jacob(3,j)=jacob(3,j)+0.5*gpth(k)*gpnorm(jtemp,k)*shp(k)
   END DO
 END DO
 
!     SAVE THE S,T, AND N VECTORS FOR CALCULATING PSI LATER.
 
 DO  i=1,3
   IF (ABS(jacob(1,i)) <= eps) jacob(1,i)=0.0
   sk(i)=jacob(1,i)
   IF (ABS(jacob(2,i)) <= eps) jacob(2,i)=0.0
   tk(i)=jacob(2,i)
   IF (ABS(jacob(3,i)) <= eps) jacob(3,i)=0.0
   enk(i)=jacob(3,i)
 END DO
 
!     THE INVERSE OF THE JACOBIAN WILL BE STORED IN
!     JACOB AFTER THE SUBROUTINE INVERS HAS EXECUTED.
 
 CALL invers (3,jacob,3,dum,0,detj,ising,INDEX)
 IF (ising == 1 .AND. detj > 0.0) GO TO 350
 WRITE (nout,550) elid
 nogo=1
 badj=.true.
 GO TO 500
 350 CALL saxb (sk,tk,v3)
 val=SQRT(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
 v3(1)=v3(1)/val
 v3(2)=v3(2)/val
 v3(3)=v3(3)/val
 
!     CROSS ELEMENT Y DIRECTION WITH UNIT VECTOR V3 IN ORDER
!     TO BE CONSISTENT WITH THE ELEMENT COORDINATE SYSTEM.
 
!     NOTE - THIS IS IMPORTANT FOR THE DIRECTIONAL REDUCED
!            INTEGRATION CASES.
 
 
 
 v2(1)=0.0
 v2(2)=1.0
 v2(3)=0.0
 
 CALL saxb (v2,v3,v1)
 val=SQRT(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
 v1(1)=v1(1)/val
 v1(2)=v1(2)/val
 v1(3)=v1(3)/val
 CALL saxb (v3,v1,v2)
 
!     REMEMBER THAT V1(1) IS EQUIVALENCED TO PSITRN(1), AND SO ON.
 
!     ELIMINATE SMALL NUMBERS
 
 DO  i = 1,3
   IF (ABS(v1(i)) <= eps) v1(i)=0.0
   IF (ABS(v2(i)) <= eps) v2(i)=0.0
   IF (ABS(v3(i)) <= eps) v3(i)=0.0
 END DO
 
 500 CONTINUE
 RETURN
 
 550 FORMAT ('0*** USER FATAL ERROR, ELEMENT ID =',i10,  &
     '  HAS BAD OR REVERSE GEOMETRY')
END SUBROUTINE jacobs
