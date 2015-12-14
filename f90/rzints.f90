FUNCTION rzints(ip,iq,r,z,ngrids)
     
 
 INTEGER, INTENT(IN)                      :: ip
 INTEGER, INTENT(IN)                      :: iq
 REAL, INTENT(IN)                         :: r(4)
 REAL, INTENT(IN)                         :: z(4)
 INTEGER, INTENT(IN OUT)                  :: ngrids
 DIMENSION  pt(3),h(3)
 IF(ngrids == 3)GO TO 200
 npt=3
 pt(1)=-.7745966692
 pt(2)=0.
 pt(3)=-pt(1)
 h(1)=5./9.
 h(2)=8./9.
 h(3)=h(1)
 xint=0.
 DO  iii=1,npt
   DO  jjj=1,npt
     rr=.25*((1.-pt(iii))*(1.-pt(jjj))*r(1) +(1.+pt(iii))*(1.-pt(jjj))*r(2)  &
         +(1.+pt(iii))*(1.+pt(jjj))*r(3) +(1.-pt(iii))*(1.+pt(jjj))*r(4))
     zz=.25*((1.-pt(iii))*(1.-pt(jjj))*z(1) +(1.+pt(iii))*(1.-pt(jjj))*z(2)  &
         +(1.+pt(iii))*(1.+pt(jjj))*z(3) +(1.-pt(iii))*(1.+pt(jjj))*z(4))
     rrp=rr**ip
     zzq=zz**iq
     drdxi=.25*((1.-pt(jjj))*(r(2)-r(1))+(1.+pt(jjj))*(r(3)-r(4)))
     dzdxi=.25*((1.-pt(jjj))*(z(2)-z(1))+(1.+pt(jjj))*(z(3)-z(4)))
     drdeta=.25*((1.-pt(iii))*(r(4)-r(1))+(1.+pt(iii))*(r(3)-r(2)))
     dzdeta=.25*((1.-pt(iii))*(z(4)-z(1))+(1.+pt(iii))*(z(3)-z(2)))
     detj=drdxi*dzdeta-dzdxi*drdeta
     detj=ABS(detj)
     xint=xint+rrp*zzq*h(iii)*h(jjj)*detj
   END DO
 END DO
 rzints=xint
 200 RETURN
END FUNCTION rzints
