SUBROUTINE geloop(rbuf,buf,xx,yy,zz,hc1,hc2,hc3)
     
! GELOOP COMPUTES MAGNETIC FIELD COMPONENTS HC1,HC2,HC3(IN BASIC
! COORDS. AT XX,YY,ZZ DUE TO GEMLOOP CARD. DATA FIELDS(EXCEPT SET ID)
! OF GEMLOOP ARE IN RBUF=REAL AND BUF=INTEGER
 
 
 REAL, INTENT(IN)                         :: rbuf(50)
 INTEGER, INTENT(IN)                      :: buf(50)
 REAL, INTENT(IN)                         :: xx
 REAL, INTENT(IN)                         :: yy
 REAL, INTENT(IN)                         :: zz
 REAL, INTENT(OUT)                        :: hc1
 REAL, INTENT(OUT)                        :: hc2
 REAL, INTENT(OUT)                        :: hc3
 INTEGER :: ti1,ti2
 DIMENSION  zi(3),zj(3),zk(3),zjxi(3)
 DATA fpi/12.566371/
 
 hc1=0.
 hc2=0.
 hc3=0.
 
 xi=rbuf(1)
 
! ICID IS 0 FOR NOW AND UNUSED
 
 icid=buf(2)
 npts=buf(3)
 nptsm1=npts-1
 DO  i=1,nptsm1
   
! 2 CONSECUTIVE POINTS DEFINE A SEGMENT OF A COIL. LET ZI BE THE VECTOR
! FROM 1ST POINT OF SEGMENT TO 2ND. LET ZJ BE VECTOR FROM FILED POINT
! XX,YY,ZZ TO 1ST POINT OF SEGMENT. ZK IS VECTOR FROM FILED POINT
! TO 2ND POINT. IF THE FILED POINT LIES ON A SEGMENT, IGNORE THE
! COMPUTATION FOR THAT SEGMENT FOR THAT POINT
   
   ti1=3*i+3
   ti2=3*(i+1)+3
   zi(1)=rbuf(ti2-2)-rbuf(ti1-2)
   zi(2)=rbuf(ti2-1)-rbuf(ti1-1)
   zi(3)=rbuf(ti2)-rbuf(ti1)
   zj(1)=rbuf(ti1-2)-xx
   zj(2)=rbuf(ti1-1)-yy
   zj(3)=rbuf(ti1)-zz
   zk(1)=rbuf(ti2-2)-xx
   zk(2)=rbuf(ti2-1)-yy
   zk(3)=rbuf(ti2)-zz
   
   zkl=SQRT(zk(1)**2+zk(2)**2+zk(3)**2)
   IF(zkl < 1.e-8)CYCLE
   zjl=SQRT(zj(1)**2+zj(2)**2+zj(3)**2)
   IF(zjl < 1.e-8)CYCLE
   zdot=0.
   DO  ii=1,3
     zdot=zdot+zi(ii)*(zk(ii)/zkl-zj(ii)/zjl)
   END DO
   zjxi(1)=zj(2)*zi(3)-zj(3)*zi(2)
   zjxi(2)=zj(3)*zi(1)-zj(1)*zi(3)
   zjxi(3)=zj(1)*zi(2)-zj(2)*zi(1)
   zlen2=zjxi(1)**2+zjxi(2)**2+zjxi(3)**2
   IF(zlen2 < 1.e-8)CYCLE
   factor=xi*zdot/fpi/zlen2
   hc1=hc1+zjxi(1)*factor
   hc2=hc2+zjxi(2)*factor
   hc3=hc3+zjxi(3)*factor
 END DO
 RETURN
END SUBROUTINE geloop
