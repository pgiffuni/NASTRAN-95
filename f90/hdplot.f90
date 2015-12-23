SUBROUTINE hdplot (gplst,nmax,maxsf,iopcor,ib)
     
 IMPLICIT INTEGER (a-z) 
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 INTEGER, INTENT(IN)                      :: nmax
 INTEGER, INTENT(IN OUT)                  :: maxsf
 INTEGER, INTENT(IN)                      :: iopcor
 INTEGER, INTENT(IN OUT)                  :: ib
 LOGICAL :: debug
 INTEGER :: NAME(2),isys(100),ptrs(29)
 REAL :: dv,psi,phi,theta,scf,p,x(20),y(20),z(20)
 COMMON /BLANK / ngp,nsil,nsets,skp1(7), skp2(2),elset,skp22(7),  &
     merr,idum(3),nscr1,nscr2,nscr3
 COMMON /system/ skps,iout
 COMMON /pltscr/ nnn,g(3)
 COMMON /hdrec / nofsur,ns,elid,lid,npers,p(3,13)
 COMMON /zzzzzz/ rz(1)
 COMMON /hdptrs/ xdum,xcc,xasolv,yasolv,zasolv,x1skt,y1skt,z1skt,  &
     zcoef1,zcoef,icount,irct,x21,y21,z21,iia,xe,ye,  &
     xu,yu,xi,yi,zi,di,ibeg,iend,ict,icct,work
 COMMON /hdsc  / scf,psi,phi,theta,mne,dv,mnp,icore
 COMMON /plothd/ used
 EQUIVALENCE     (isys(1),skps), (ptrs(1),xdum)
 DATA    NAME  / 4HHDPL,4HOT  /
 DATA    debug / .false.      /
 
!     CALL SSWTCH (47,J)
!     IF (J .EQ. 1) DEBUG = .TRUE.
 
!     SET MNE EQUAL TO THE MAXIMUM NUMBER OF EDGES IN ANY ONE POLYGON.
 
 mne = nmax
 
!     MNP=MNE+2+2*NHOLES   WHERE NHOLES IS THE NUMBER OF HOLES,IF ANY
 
 nholes = 0
 mnp = mne + 2 + 2*nholes
 
!     SET DISTANCE FROM VIEWER, AND SET SCALING FACTOR = 1 UNITS/INCH
 
 dv  = 99999.
 scf = 1.00
 
!     SET MAX. LINES OF INTERSECTION ALLOWED IN HDSOLV (DIMEN. OF XCC)
 
 lintc = 800
 IF (isys(85) /= 0) lintc = isys(85)
 
!     DEFINE EULERIAN ANGLES IN DEGREES.
 
 psi = 0.
 phi = 0.
 theta = 0.
 
!     INITIALIZE ARRAY POINTERS IN OPEN CORE SPACE (USED, SET BY PLOT,
!     IS NO. OF WORDS ALREADY IN USE)
 
 xdum  = 1
 xcc   = xdum  + used
 xasolv= xcc   + lintc
 yasolv= xasolv+ 50
 zasolv= yasolv+ 50
 x1skt = zasolv+ 50
 y1skt = x1skt + 160
 z1skt = y1skt + 160
 zcoef1= z1skt + 160
 zcoef = zcoef1+ 150
 icount= zcoef + 150
 irct  = icount+ 150
 x21   = irct  + 100
 y21   = x21   + 200
 z21   = y21   + 200
 iia   = z21   + 200
 xe    = iia   + 200
 ye    = xe    + 150
 xu    = ye    + 150
 yu    = xu    + 150
 ibeg  = yu    + 150
 iend  = ibeg  + 100
 ict   = iend  + 100
 icct  = ict   + 100
 xi    = icct  + 100
 icore = (25+5*mne+4*mnp)*(maxsf+1)
 j     = (iopcor-icore-xi)/5
 yi    = xi    + j
 zi    = yi    + j
 di    = zi    + j
 work  = di    + j
 IF (debug .OR. j < 300) WRITE (iout,55) nmax,maxsf,icore,used,  &
     lintc,iopcor,ib,nsets,j,ptrs
 IF (j >= 300) GO TO 5
 j = 300*5 + xi + icore - iopcor
 CALL mesage (-8,j,NAME)
 
 5 CALL gopen (nscr2,gplst(ib),0)
 CALL line (0.,0.,0.,0.,1,-1)
 10 CONTINUE
 CALL READ (*25,*25,nscr2,nofsur,44,0,m)
 nps = npers
 DO  i = 1,nps
   x(i) = p(1,i)
   y(i) = p(2,i)
   z(i) = p(3,i)
 END DO
 IF (debug) WRITE (iout,65) nofsur,ns,elid,lid,nps,(x(n),y(n),z(n),n=1,nps)
 nc = 0
 CALL hdsket (x,y,z,nps,nc)
 GO TO 10
 25 CALL CLOSE (nscr2,1)
 nc = 1
 CALL hdsket (x,y,z,nps,nc)
 IF (nc == 0) GO TO 40
 WRITE  (iout,30) nc,icore,dv
 30 FORMAT (22H code for hidden error,i3,6H icore,i9,3H dv,f15.5)
 40 CALL line (0.,0.,0.,0.,1,+1)
 IF (debug) WRITE (iout,60)
 RETURN
 
 55 FORMAT (1X,10HIN hdplot ,9I8, /,(5X,15I8))
 60 FORMAT (1X,10HOUT hdplot)
 65 FORMAT (1X,5I10/(1X,3G20.4))
END SUBROUTINE hdplot
