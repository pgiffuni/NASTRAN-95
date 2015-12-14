SUBROUTINE dzymat (d,nfb,nlb,ntzys,idzdy,ntape,xp,beta,iprnt,ns,  &
        nc,yp,zp,sg,cg,yb,zb,nbea)
     
!     CALCULATION OF DZ AND DY MATRICES SLENDER BODY CALCULATIONS
 
!     D         WORKING ARRAY USED TO STORE A ROW OF DZ OR DY
!     NFB       NUMBER OF THE FIRST BODY WITH THE ORIENTATION REQUESTED
!     NLB       NUMBER OF THE LAST BODY WITH THE ORIENTATION
!               REQUESTED
!     NTZYS     NUMBER OF Z OR Y ORIENTED SLENDER BODY ELE.
!     NTAPE     I/O UNIT NUMBER WHICH THE OUTPUT MATRIX IS TO
!               BE WRITTEN ON
!     XP        X-CONTROL POINT COORDINATE OF LIFTING SURFACE
!               BOXES
!     BETA      SQRT(1.0 - M**2)
 
 
 REAL, INTENT(IN OUT)                     :: d(2,ntzys)
 INTEGER, INTENT(IN OUT)                  :: nfb
 INTEGER, INTENT(IN OUT)                  :: nlb
 INTEGER, INTENT(IN OUT)                  :: ntzys
 INTEGER, INTENT(IN OUT)                  :: idzdy
 INTEGER, INTENT(IN OUT)                  :: ntape
 REAL, INTENT(IN)                         :: xp(1)
 REAL, INTENT(IN OUT)                     :: beta
 INTEGER, INTENT(IN OUT)                  :: iprnt
 INTEGER, INTENT(IN)                      :: ns(1)
 INTEGER, INTENT(IN)                      :: nc(1)
 REAL, INTENT(IN)                         :: yp(1)
 REAL, INTENT(IN)                         :: zp(1)
 REAL, INTENT(IN)                         :: sg(1)
 REAL, INTENT(IN)                         :: cg(1)
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: zb(1)
 INTEGER, INTENT(IN)                      :: nbea(1)
 INTEGER :: by,bz,c,c1,p,s,s1,yt,zt
 
 COMMON /dlbdy / nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,  &
     izb,iyb,iavr,iarb,infl,ixle,ixte,int121,int122,  &
     izs,iys,ics,iee,isg,icg,ixij,ix,idelx,ixic,ixlam,  &
     ia0,ixis1,ixis2,ia0p,iria,inasb,ifla1,ifla2,ith1a,  &
     ith2a,ecore,next,scr1,scr2,scr3,scr4,scr5
 COMMON /zzzzzz/ z(1)
 COMMON /amgmn / mcb(7),nrow,nd,NE,refc,fmach,kr
 COMMON /system/ sysbuf,npot
 
 c1 = 0
 s1 = 0
 nfyb = nb - nby  + 1
 IF (np == 0) GO TO 410
 
!     THIS LOOP IS FOR EACH LIFTING SURF. PANEL
 
 isn = 0
 DO  p = 1,np
   nsp = ns(p)
   ncp = nc(p)
   nsp = (nsp-isn)/ncp
   isn = ns(p)
   
!     LOOP FOR EACH STRIP IN PANEL -P-
   
   DO  s = 1,nsp
     s1  = s1 + 1
     
!     Y AND Z COORDINATE OF STRIP
     
     dy  = yp(s1)
     dz  = zp(s1)
     sgr = sg(s1)
     cgr = cg(s1)
     
!     LOOP FOR EACH CHORDWISE ELEMENT IN STRIP
     
     DO  c = 1,ncp
       c1  = c1 + 1
       dx  = xp(c1)
       
!     - ROWDYC -  CALCULATES ROW -C1- OF DZ OR DY
       
       CALL rowdyz (nfb,nlb,c1,ntzys,d,dx,dy,dz,beta,idzdy,ntape,sgr,  &
           cgr,iprnt,yb,zb,z(iarb),z(insbea),z(ixis1),z(ixis2), z(ia0))
       
     END DO
   END DO
 END DO
 
!     WE HAVE NOW CALCULATED -C1- ROWS WHICH ARE THE LIFTING SURFACES.
!     NOW, LOOP FOR THE -Z- ORIENTED BODIES
 
 410 CONTINUE
 IF (nbz <= 0 .OR. ntz <= 0) GO TO 510
 sgr = 0.0
 cgr = 1.0
 DO  bz = 1,nbz
   dy  = yb(bz)
   dz  = zb(bz)
   nbez = nbea(bz)
   
!     LOOP FOR EACH ELEMENT OF BODY -BZ-
   
   DO  zt = 1,nbez
     c1  = c1 + 1
     dx  = xp(c1)
     
     CALL rowdyz (nfb,nlb,c1,ntzys,d,dx,dy,dz,beta,idzdy,ntape,sgr,  &
         cgr,iprnt,yb,zb,z(iarb),z(insbea),z(ixis1),z(ixis2), z(ia0))
   END DO
 END DO
 
!     NOW, LOOP FOR THE -Y- ORIENTED BODIES
 
 510 IF (nb < nfyb .OR. nty <= 0) GO TO 650
 ixp = ntp
 IF (nfyb <= 1) GO TO 530
 nfybm1 = nfyb - 1
 DO  i = 1,nfybm1
   ixp = ixp + nbea(i)
 END DO
 530 CONTINUE
 sgr =-1.0
 cgr = 0.0
 DO  by = nfyb,nb
   dy  = yb(by)
   dz  = zb(by)
   nbey = nbea(by)
   
!     LOOP FOR EACH ELEMENT OF BODY -BY-
   
   DO  yt = 1,nbey
     c1  = c1  + 1
     ixp = ixp + 1
     dx  = xp(ixp)
     
     CALL rowdyz (nfb,nlb,c1,ntzys,d,dx,dy,dz,beta,idzdy,ntape,sgr,  &
         cgr,iprnt,yb,zb,z(iarb),z(insbea),z(ixis1),z(ixis2), z(ia0))
     
   END DO
 END DO
 650 CONTINUE
 RETURN
END SUBROUTINE dzymat
