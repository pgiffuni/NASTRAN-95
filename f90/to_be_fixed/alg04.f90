SUBROUTINE alg04(h,s,vw,r1,r2,x1,x2,vm,eps,sclfac,g,ej,hmin,vmin,  &
        psmid,nstrms,log2,lnct,ifail)
     
 
 REAL, INTENT(IN OUT)                     :: h(1)
 REAL, INTENT(IN OUT)                     :: s(1)
 REAL, INTENT(IN OUT)                     :: vw(1)
 REAL, INTENT(IN OUT)                     :: r1(1)
 REAL, INTENT(IN)                         :: r2(1)
 REAL, INTENT(IN)                         :: x1(1)
 REAL, INTENT(IN)                         :: x2(1)
 REAL, INTENT(IN)                         :: vm(1)
 REAL, INTENT(IN)                         :: eps
 REAL, INTENT(IN)                         :: sclfac
 REAL, INTENT(IN)                         :: g
 REAL, INTENT(IN)                         :: ej
 REAL, INTENT(IN OUT)                     :: hmin
 REAL, INTENT(IN)                         :: vmin
 REAL, INTENT(IN)                         :: psmid
 INTEGER, INTENT(IN)                      :: nstrms
 INTEGER, INTENT(IN OUT)                  :: log2
 INTEGER, INTENT(IN OUT)                  :: lnct
 INTEGER, INTENT(OUT)                     :: ifail
 
 DIMENSION vzup(21),vzdn(21),psdn(21),hup(21),vwup(21),sup(21),rmid  &
     (20),delr(20),psup(21),hsup(21),vwfun(21),vzfun(21),sdn(21),hsdn(2 1),vwdn(21),hdn(21),xx1(21),xx2(21),r(21)
 
 DO  j=1,nstrms
   r(j)=(r1(j)+r2(j))*0.5
 END DO
 q1=r(nstrms)-r(1)
 q2=vm(1)
 DO  j=2,nstrms
   IF(r(j)-r(j-1) < q1)q1=r(j)-r(j-1)
   IF(vm(j) < q2)q2=vm(j)
 END DO
 delz=q2*q1**2/(eps*sclfac)*0.25
 q1=(x2(1)+x2(nstrms)-x1(1)-x1(nstrms))*0.5
 istep=q1/delz+1.0
 delz=q1/FLOAT(istep)
 vm2=vmin**2
 itub=nstrms-1
 imid=nstrms/2+1
 DO  j=1,nstrms
   psup(j)=psmid
   hup(j)=h(j)
   vwup(j)=vw(j)
   sup(j)=s(j)
 END DO
 DO  j=1,itub
   rmid(j)=(r(j)+r(j+1))*0.5
   delr(j)=r(j+1)-r(j)
 END DO
 ifail=0
 kstep=1
 130   CALL alg29(vwup,r,xx2,nstrms)
 DO  j=1,nstrms
   vwfun(j)=eps/r(j)*(xx2(j)-vwup(j)/r(j))*sclfac
 END DO
 IF(kstep > 1)GO TO 280
 jstep=1
 j1=imid
 190   j2=j1+jstep
 jj=j1
 IF(jstep == -1)jj=j2
 q1=((vwup(j1)+vwup(j2))*0.5)**2/rmid(jj)
 q1=delr(jj)*q1*FLOAT(jstep)
 x3=(sup(j1)+sup(j2))*0.5
 k=1
 200   q2=alg2(x3,(psup(j1)+psup(j2))*0.5)
 IF(q2 >= hmin)GO TO 210
 ifail=1
 GO TO 600
 210   q2=alg5(q2,x3)/g
 x4=psup(j2)
 psup(j2)=psup(j1)+q1*q2
 IF(ABS(x4/psup(j2)-1.0) <= 1.0E-5)GO TO 220
 k=k+1
 IF(k <= 10)GO TO 200
 ifail=2
 GO TO 600
 220   IF(j2 == 1)GO TO 240
 IF(j2 == nstrms)GO TO 230
 j1=j2
 GO TO 190
 230   jstep=-1
 j1=imid
 GO TO 190
 240   DO  j=1,nstrms
   hsup(j)=alg2(sup(j),psup(j))
   IF(hsup(j) >= hmin)GO TO 250
   ifail=3
   GO TO 600
   250   q1=2.0*g*ej*(hup(j)-hsup(j))-vwup(j)**2
   IF(q1 >= vm2)GO TO 260
   ifail=4
   GO TO 600
   260   vzup(j)=SQRT(q1)
 END DO
 flow=0.0
 DO  j=1,itub
   flow=flow+(r(j+1)**2-r(j)**2)*(vzup(j)+vzup(j+1))*alg5((hsup(j)+hs  &
       up(j+1))*0.5,(sup(j)+sup(j+1))*0.5)
 END DO
 280   CALL alg29(hsup,r,xx2,nstrms)
 DO  j=1,nstrms
   hdn(j)=hup(j)+delz/vzup(j)*eps/r(j)*xx2(j)*sclfac
 END DO
 DO  j=1,nstrms
   vwdn(j)=vwup(j)+delz/vzup(j)*vwfun(j)
 END DO
 CALL alg29(vzup,r,vzfun,nstrms)
 DO  j=1,nstrms
   vzfun(j)=delz*eps*sclfac*vzfun(j)/r(j)
   sdn(j)=sup(j)
   psdn(j)=psup(j)
 END DO
 kk=1
 340   j1=imid
 jstep=1
 350   j2=j1+jstep
 jj=j1
 IF(jstep == -1)jj=j2
 q1=((vwdn(j1)+vwdn(j2))*0.5)**2/rmid(jj)
 q1=delr(jj)*q1*FLOAT(jstep)
 x3=(sdn(j1)+sdn(j2))*0.5
 k=1
 360   q2=alg2(x3,(psdn(j1)+psdn(j2))*0.5)
 IF(q2 >= hmin)GO TO 370
 ifail=5
 GO TO 600
 370   q2=alg5(q2,x3)/g
 x4=psdn(j2)
 psdn(j2)=psdn(j1)+q1*q2
 IF(ABS(x4/psdn(j2)-1.0) <= 1.0E-5)GO TO 380
 k=k+1
 IF(k <= 10)GO TO 360
 ifail=6
 GO TO 600
 380   IF(j2 == 1)GO TO 400
 IF(j2 == nstrms)GO TO 390
 j1=j2
 GO TO 350
 390   j1=imid
 jstep=-1
 GO TO 350
 400   DO  j=1,nstrms
   vzdn(j)=vzup(j)+(vzfun(j)-(psdn(j)-psup(j))/alg5(hsup(j),sup(j))*g  &
       )/vzup(j)
   hsdn(j)=hdn(j)-(vzdn(j)**2+vwdn(j)**2)/(2.0*g*ej)
   IF(hsdn(j) >= hmin)GO TO 410
   ifail=7
   GO TO 600
   410   sdn(j)=alg3(psdn(j),hsdn(j))
 END DO
 xx1(1)=0.0
 DO  j=1,itub
   xx1(j+1)=xx1(j)+(r(j+1)**2-r(j)**2)*(vzdn(j+1)+vzdn(j))*alg5((hsdn  &
       (j)+hsdn(j+1))*0.5,(sdn(j)+sdn(j+1))*0.5)
 END DO
 q1=xx1(nstrms)
 IF(ABS(q1/flow-1.0) <= 1.0E-5.AND.kk > 1)GO TO 450
 IF(kk <= 15)GO TO 430
 ifail=8
 GO TO 600
 430   q2=alg9(hsdn(imid),sdn(imid),vzdn(imid)**2)
 q1=(q1-flow)*psdn(imid)*q2/(flow*(1.0-q2))
 DO  j=1,nstrms
   psdn(j)=psdn(j)+q1
 END DO
 kk=kk+1
 GO TO 340
 450   IF(kstep == istep)GO TO 510
 DO  j=1,nstrms
   psup(j)=psdn(j)
   hsup(j)=hsdn(j)
   vzup(j)=vzdn(j)
   vwup(j)=vwdn(j)
   hup(j)=hdn(j)
   sup(j)=sdn(j)
 END DO
 kstep=kstep+1
 GO TO 130
 510   DO  j=1,nstrms
   h(j)=hdn(j)
   s(j)=sdn(j)
   vw(j)=vwdn(j)
 END DO
 RETURN
 600   CALL alg03(lnct,1)
 WRITE(log2,610)ifail
 610   FORMAT(5X,30HMIXING calculation failure no.,i2)
 RETURN
END SUBROUTINE alg04
