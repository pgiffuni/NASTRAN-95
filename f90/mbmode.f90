SUBROUTINE mbmode(INPUT,out,icor,ncor,z,ni,nd,xd,yd,is,cr)
     
!     MBMODE BUILDS THE MODE LIKE DATA ON OUT FROM SURFACE SPLINE INTER
 
 
 INTEGER, INTENT(IN OUT)                  :: INPUT
 REAL, INTENT(IN OUT)                     :: out
 INTEGER, INTENT(IN)                      :: icor
 INTEGER, INTENT(IN)                      :: ncor
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(IN)                      :: ni
 INTEGER, INTENT(IN)                      :: nd
 REAL, INTENT(IN)                         :: xd(1)
 REAL, INTENT(IN)                         :: yd(1)
 INTEGER, INTENT(IN OUT)                  :: is
 REAL, INTENT(IN)                         :: cr
 DIMENSION  NAME(2)
 DATA NAME /4HMBMO,4HDE  /
 
 nni = ni*2
 nnd = nd*2
 IF(icor+nni+nnd > ncor) CALL mesage(-8,0,NAME)
 CALL fread(INPUT,z(icor),nni,0)
 idp = icor + nni
 l = 0
 DO  i=1,nd
   z(idp+l) = xd(i)
   z(idp+l+1) = yd(i)
   l = l+2
 END DO
 icc = idp+l
 ncore = ncor-icc
 
!     CALL SSPLIN TO INTERPOLATE
 
 CALL ssplin(ni,z(icor),nd,z(idp),0,0,1,0,0.0,z(icc),ncore,is)
 IF(is == 2) GO TO 1000
 
!     REORDER INTO MACH BOX ORDER
 
 m = idp+nd
 icc = icc-1
 DO  i=1,ni
   l = 0
   DO  j=1,nnd,2
     z(idp+l) = z(icc+j)
     z(m+l) = z(icc+j+1) * cr
     l = l+1
   END DO
   CALL WRITE(out,z(idp),nnd,0)
   icc = icc + nnd
 END DO
 1000 RETURN
END SUBROUTINE mbmode
