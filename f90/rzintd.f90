DOUBLE PRECISION FUNCTION rzintd(ip,iq,r,z,ngrids)
     
 INTEGER, INTENT(IN)                      :: ip
 INTEGER, INTENT(IN)                      :: iq
 DOUBLE PRECISION, INTENT(IN)             :: r(4)
 DOUBLE PRECISION, INTENT(IN)             :: z(4)
 INTEGER, INTENT(IN OUT)                  :: ngrids
 DOUBLE PRECISION :: pt(3),h(3)
 DOUBLE PRECISION :: xint,rrp,zzq,drdxi,dzdxi,drdeta,dzdeta,detj
 DOUBLE PRECISION :: rr,zz
 
 IF(ngrids == 3)GO TO 200
 npt=3
 pt(1)=-.7745966692D0
 pt(2)=0.d0
 pt(3)=-pt(1)
 h(1)=5.d0/9.d0
 h(2)=8.d0/9.d0
 h(3)=h(1)
 xint=0.d0
 DO  iii=1,npt
   DO  jjj=1,npt
     rr=.25D0*((1.d0-pt(iii))*(1.d0-pt(jjj))*r(1)  &
         +(1.d0+pt(iii))*(1.d0-pt(jjj))*r(2)  &
         +(1.d0+pt(iii))*(1.d0+pt(jjj))*r(3) +(1.d0-pt(iii))*(1.d0+pt(jjj))*r(4))
     zz=.25D0*((1.d0-pt(iii))*(1.d0-pt(jjj))*z(1)  &
         +(1.d0+pt(iii))*(1.d0-pt(jjj))*z(2)  &
         +(1.d0+pt(iii))*(1.d0+pt(jjj))*z(3) +(1.d0-pt(iii))*(1.d0+pt(jjj))*z(4))
     rrp=rr**ip
     zzq=zz**iq
     drdxi=.25D0*((1.d0-pt(jjj))*(r(2)-r(1))+(1.d0+pt(jjj))*(r(3)-r(4)) )
     dzdxi=.25D0*((1.d0-pt(jjj))*(z(2)-z(1))+(1.d0+pt(jjj))*(z(3)-z(4)) )
     drdeta=.25D0*((1.d0-pt(iii))*(r(4)-r(1))+(1.d0+pt(iii))*(r(3)-r(2) ))
     dzdeta=.25D0*((1.d0-pt(iii))*(z(4)-z(1))+(1.d0+pt(iii))*(z(3)-z(2) ))
     detj=drdxi*dzdeta-dzdxi*drdeta
     detj=DABS(detj)
     xint=xint+rrp*zzq*h(iii)*h(jjj)*detj
   END DO
 END DO
 rzintd=xint
 200 RETURN
END FUNCTION rzintd
