SUBROUTINE dipole(buf,ibuf,xx,yy,zz,hc1,hc2,hc3)
     
! DIPOLE COMPUTES THE MAGNETIC INTENSITY AT THE POINT (XX,YY,ZZ) DUE
! TO A MAGNETIC DIPOLE MOMENT DEFINED ON AN MDIPOLE CARD STORED IN BUF.
! THE FORMULATION COMES FROM DARRELL NIXONS REPORT 27-23 MARCH 1972
 
 
 REAL, INTENT(IN)                         :: buf(50)
 INTEGER, INTENT(IN)                      :: ibuf(50)
 REAL, INTENT(IN OUT)                     :: xx
 REAL, INTENT(IN OUT)                     :: yy
 REAL, INTENT(IN OUT)                     :: zz
 REAL, INTENT(OUT)                        :: hc1
 REAL, INTENT(OUT)                        :: hc2
 REAL, INTENT(OUT)                        :: hc3
 REAL :: mx,my,mz,MIN,MAX,mxa,myb,mzc
 
 DATA fpi/12.566371/
 
 hc1=0.
 hc2=0.
 hc3=0.
 
! ICID IS 0 FOR NOW AND WILL NOT BE USED. COORDS. AND MOMENT
! ARE ASSUMED TO BE IN BASIC COORDS
 
 icid=ibuf(1)
 
! COORDS OF POINT DIPOLE
 
 cx=buf(2)
 cy=buf(3)
 cz=buf(4)
 mx=buf(5)
 my=buf(6)
 mz=buf(7)
 MIN=buf(8)
 MAX=buf(9)
 
! H WILL BE COMPUTED ONLY IF DISTANCE FROM (CX,CY,CZ) TO (XX,YY,ZZ) IS
! BETWEEN MIN AND MAX- IF MAX IS 0, COMPUTE FOR ALL POINTS GREATER THAN
! MIN
 
 rmr1=SQRT((xx-cx)**2+(yy-cy)**2+(zz-cz)**2)
 IF(MIN <= 1.e-6)GO TO 5
 IF(rmr1 < MIN)GO TO 20
 5 IF(MAX <= 1.e-6)GO TO 10
 IF(rmr1 > MAX)GO TO 20
 
 10 mxa=3.*mx*(xx-cx)
 myb=3.*my*(yy-cy)
 mzc=3.*mz*(zz-cz)
 
 r3=rmr1**3
 r5=r3*rmr1**2
 xnum=(mxa+myb+mzc)/r5
 
 hc1=-mx/r3+xnum*(xx-cx)
 hc1=hc1/fpi
 
 hc2=-my/r3+xnum*(yy-cy)
 hc2=hc2/fpi
 
 hc3=-mz/r3+xnum*(zz-cz)
 hc3=hc3/fpi
 
 20 RETURN
END SUBROUTINE dipole
