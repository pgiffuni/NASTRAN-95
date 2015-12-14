SUBROUTINE linein(x1,y1,z1,x2,y2,z2,hcdl)
     
! PERFORMS LINE INTEGRAL FROM (X1,Y1,Z1) TO (X2,Y2,Z2) OF BIOT-SAVART
! FILED DOTTED INTO THE LINE, IE INT(HC.DL)
 
 
 REAL, INTENT(IN)                         :: x1
 REAL, INTENT(IN)                         :: y1
 REAL, INTENT(IN)                         :: z1
 REAL, INTENT(IN)                         :: x2
 REAL, INTENT(IN)                         :: y2
 REAL, INTENT(IN)                         :: z2
 REAL, INTENT(OUT)                        :: hcdl
 DIMENSION xi(4),w(4)
 DATA xi/.06943184,.33000948,.66999052,.93056816/
 DATA w/.17392742,2*.32607258,.173927423/
 
! COMPONENTS OF LINE SEGMENT
 
 hcdl=0.
 segx=x2-x1
 segy=y2-y1
 segz=z2-z1
 segl=SQRT(segx**2+segy**2+segz**2)
 IF(segl == 0.)RETURN
 
! 4 POINT INTEGRATION OVER LINE SEGMENT(XI= / TO +1)
 
 DO  i=1,4
   xx=x1+segx*xi(i)
   yy=y1+segy*xi(i)
   zz=z1+segz*xi(i)
   CALL biotsv(xx,yy,zz,hcx,hcy,hcz)
   hcdl=hcdl+(hcx*segx+hcy*segy+hcz*segz)*w(i)
 END DO
 RETURN
END SUBROUTINE linein
