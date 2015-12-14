SUBROUTINE q4shps (xi,eta, shp, sshp)
!     &    ENTRY Q4SHPD (DI,ETD,DSHP,DSSHP)
 
!*****
!     COMPUTES SHAPE FUNCTIONS AND THEIR DERIVATIVES
!     FOR THE QUAD4 ELEMENT
!*****
 
 
 REAL, INTENT(IN)                         :: xi
 REAL, INTENT(IN)                         :: eta
 REAL, INTENT(OUT)                        :: shp(4)
 REAL, INTENT(OUT)                        :: sshp(8)
 REAL :: clc(2,4)
 DOUBLE PRECISION :: di,etd,dshp(4),dsshp(8),dld(2,4)
 DATA   clc /-1.0  ,-1.0  ,1.0  ,-1.0  ,1.0  ,1.0  ,-1.0  ,1.0  /
 DATA   dld /-1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0,1.0D0/
 
!     SINGLE PRECISION -
 
 DO  i=1,4
   shp (i  ) = 0.25 * (1.0+xi *clc(1,i)) * (1.0+eta*clc(2,i))
   sshp(i  ) = 0.25 * (1.0+eta*clc(2,i)) * clc(1,i)
   sshp(i+4) = 0.25 * (1.0+xi *clc(1,i)) * clc(2,i)
 END DO
 RETURN
 
 ENTRY q4shpd (di,etd,dshp,dsshp)
!     ================================
 
!     DOUBLE PRECISION -
 
 DO  i=1,4
   dshp (i  ) = .25D0 * (1.d0+di *dld(1,i)) * (1.d0+etd*dld(2,i))
   dsshp(i  ) = .25D0 * (1.d0+etd*dld(2,i)) * dld(1,i)
   dsshp(i+4) = .25D0 * (1.d0+di *dld(1,i)) * dld(2,i)
 END DO
 RETURN
END SUBROUTINE q4shps
