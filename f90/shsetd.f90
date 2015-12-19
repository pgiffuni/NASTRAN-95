SUBROUTINE shsetd (*,mm,sil,bgpdt,igpdt,gpth,elth,gptemp,bgpdm,  &
        egpdt,dgpth,gpnorm,epnorm,nnode,mmn,nsil,  &
        iorder,iordrn,teb,tub,cente,avgthk,tce,elid)
     
!     TO SET UP FOR ISOPARAMETRIC SHELL ELEMENTS, CALLED ONLY BY SHHMGD
 
!     DOUBLE PRECISION VERSION
 
!     INPUT :
!           MM       - MAXIMUM NO. OF NODES PER THIS TYPE ELEMENT
!           SIL      - ARRAY OF SIL NUMBERS
!           BGPDT    - BGPDT DATA FROM EST (REAL ARRAY)
!           IGPDT    - BGPDT DATA FROM EST (INTEGER ARRAY)
!           GPTH     - GRID POINT THICKNESS DATA
!           ELTH     - ELEMENT THICKNESS FROM EPT
!           GPTEMP   - GRID POINT TEMPERATURE DATA
!           ELID     - ELEMENT ID
!     OUTPUT:
!           SIL      - ARRAY OF SIL NUMBERS        (REARRANGED)
!           BGPDT    - BGPDT DATA (REAL ARRAY)     (REARRANGED)
!           IGPDT    - BGPDT DATA (INTEGER ARRAY)  (REARRANGED)
!           GPTH     - GRID POINT THICKNESS DATA   (REARRANGED)
!           GPTEMP   - GRID POINT TEMPERATURE DATA (REARRANGED)
!           BGPDM    - BGPDT DATA SAVED IN ORIGINAL FORMAT
!           EGPDT    - BGPDT DATA IN ELEMENT COORD. SYSTEM
!           DGPTH    - GRID POINT THICKNESS DATA
!           GPNORM   - GRID POINT NORMALS
!           EPNORM   - GRID POINT NORMALS IN ELEMENT COORD. SYSTEM
!           NNODE    - THE NO. OF NODES PRESENT IN THE ELEMENT
!           MMN      - ARRAY OF MISSING MIDSIDE NODES
!           NSIL     - INTERNALLY ORDERED SIL ARRAY
!           IORDER   - ARRAY OF ORDER INDICATORS FOR REARRANGED DATA
!           IORDRN   - ARRAY OF ORDER INDICATORS FOR TRIA
!           TEB      - TRANSFORMATION FROM ELEMENT TO BASIC COORD.SYSTEM
!           TUB      - TRANSFORMATION FROM USER TO BASIC COORD. SYSTEM
!           CENTE    - LOCATION OF THE CENTER OF THE ELEMENT
!           AVGTHK   - AVERAGE THICKNESS OF THE ELEMENT
 
 
 INTEGER, INTENT(IN)                      :: mm
 INTEGER, INTENT(IN OUT)                  :: sil(8)
 REAL, INTENT(IN OUT)                     :: bgpdt(4,8)
 INTEGER, INTENT(IN OUT)                  :: igpdt(4,8)
 REAL, INTENT(IN OUT)                     :: gpth(8)
 REAL, INTENT(IN)                         :: elth
 REAL, INTENT(IN OUT)                     :: gptemp(8)
 REAL, INTENT(OUT)                        :: bgpdm(3,8)
 DOUBLE PRECISION, INTENT(OUT)            :: egpdt(4,8)
 DOUBLE PRECISION, INTENT(OUT)            :: dgpth(8)
 DOUBLE PRECISION, INTENT(OUT)            :: gpnorm(4,8)
 DOUBLE PRECISION, INTENT(OUT)            :: epnorm(4,8)
 INTEGER, INTENT(OUT)                     :: nnode
 INTEGER, INTENT(OUT)                     :: mmn(8)
 INTEGER, INTENT(OUT)                     :: nsil(8)
 INTEGER, INTENT(OUT)                     :: iorder(8)
 INTEGER, INTENT(OUT)                     :: iordrn(8)
 DOUBLE PRECISION, INTENT(IN)             :: teb(9)
 DOUBLE PRECISION, INTENT(IN)             :: tub(9)
 DOUBLE PRECISION, INTENT(OUT)            :: cente(3)
 DOUBLE PRECISION, INTENT(OUT)            :: avgthk
 DOUBLE PRECISION, INTENT(IN OUT)         :: tce(63)
 INTEGER, INTENT(IN OUT)                  :: elid
 LOGICAL :: quad
 INTEGER :: ksil(8),kcid(8)
 REAL :: temtem(8), tgrid(4,8),  tmpthk(8)
 DOUBLE PRECISION :: cent(3), ggu(9),ggn(9),  &
     teu(9),smax,smin,sl(3),gge(9), cc,  x31,y31,x42,y42,exi,exj,  &
     aa,bb,ugpdm(3,8)
 
 
 IF (mm /= 3 .AND. mm /= 4 .AND. mm /= 6 .AND. mm /= 8) GO TO 700
!           TRIA3         QUAD4         TRIA6         QUAD8
 
 quad = mm == 8 .OR. mm == 4
 mmx  = 3
 IF (quad) mmx = 4
 nnode = mm
 DO  i = 1,mm
   mmn(i) = sil(i)
   ksil(i)= sil(i)
   IF (sil(i) > 0) CYCLE
   nnode = nnode - 1
 END DO
 
!     FILL IN ARRAY GGU WITH THE COORDINATES OF GRID POINTS 1,2 AND 4
!     (3 FOR TRIA). THIS ARRAY WILL BE USED LATER TO DEFINE THE USER
!     COORDINATE SYSTEM WHILE CALCULATING TRANSFORMATIONS INVOLVING
!     THIS COORDINATE SYSTEM.
 
 DO  i = 1,3
   ii = (i-1)*3
   ij = i
   IF (quad .AND. ij == 3) ij = 4
   DO  j = 1,3
     jj = j + 1
     ggu(ii+j) = DBLE(bgpdt(jj,ij))
   END DO
 END DO
 CALL betrnd (tub,ggu,0,elid)
 
!     STORE INCOMING BGPDT FOR LUMPED MASS AND ELEMENT COORD. SYSTEM
 
 DO  i = 1,3
   i1 = i + 1
   DO  j = 1,mm
     bgpdm(i,j) = bgpdt(i1,j)
   END DO
 END DO
 
!     TRANSFORM BGPDM FROM BASIC TO USER COORD. SYSTEM
 
 DO  i = 1,3
   ip = (i-1)*3
   DO  j = 1,mm
     ugpdm(i,j) = 0.0D0
     DO  k = 1,3
       kk = ip + k
       ugpdm(i,j) = ugpdm(i,j) + tub(kk)*(DBLE(bgpdm(k,j))-ggu(k))
     END DO
   END DO
 END DO
 
 IF (quad) GO TO 200
 
!     FOR TRIA
!     CALCULATE THE CENTER COORDINATES
 
 cente(1) = (ggu(1)+ggu(4)+ggu(7))/3.0D0
 cente(2) = (ggu(2)+ggu(5)+ggu(8))/3.0D0
 cente(3) = (ggu(3)+ggu(6)+ggu(9))/3.0D0
 
!     ESTABLISH THE INTERNAL COORDINATES:
!     X-AXIS IS ALONG THE MIDDLE-SIZED SIDE AND THE XY-PLANE IS
!     DETERMINED BY IT TOGETHER WITH THE SHORTEST SIDE
 
 cc = (ggu(7)-ggu(4))*(ggu(7)-ggu(4)) + (ggu(8)-ggu(5))*(ggu(8)-ggu(5))  &
     + (ggu(9)-ggu(6))*(ggu(9)-ggu(6))
 IF (cc <= 0.0D0) GO TO 700
 sl(1) = DSQRT(cc)
 cc = (ggu(7)-ggu(1))*(ggu(7)-ggu(1)) + (ggu(8)-ggu(2))*(ggu(8)-ggu(2))  &
     + (ggu(9)-ggu(3))*(ggu(9)-ggu(3))
 IF (cc <= 0.0D0) GO TO 700
 sl(2) = DSQRT(cc)
 cc = (ggu(4)-ggu(1))*(ggu(4)-ggu(1)) + (ggu(5)-ggu(2))*(ggu(5)-ggu(2))  &
     + (ggu(6)-ggu(3))*(ggu(6)-ggu(3))
 IF (cc <= 0.0D0) GO TO 700
 sl(3) = DSQRT(cc)
 smax  = sl(1)
 ismax = 1
 DO  i = 2,3
   IF (sl(i) <= smax) CYCLE
   smax  = sl(i)
   ismax = i
 END DO
 smin  = sl(1)
 ismin = 1
 DO  i = 2,3
   IF (sl(i) >= smin) CYCLE
   smin  = sl(i)
   ismin = i
 END DO
 IF (ismax == ismin) ismin = 3
 middl = IABS(ismax-ismin)
 IF (ismax+ismin == 3) middl = 3
 
!     DETECT THE POSSIBLE REVERSAL OF THE INTERNAL Z-AXIS WITH RESPECT
!     TO THE USER Z-AXIS. IF THAT IS THE CASE, SWITCH ISMAX AND ISMIN
!     TO AVOID THE PROBLEM. THE SIDE WITH MEDIUM LENGTH WILL STILL BE
!     THE X-AXIS.
 
 IF (ismax /= MOD(ismin,3)+1) GO TO 120
 iii    = ismin
 ismin  = ismax
 ismax  = iii
 
 120 is3    = 3*(ismax-1)
 ggn(1) = ggu(is3+1)
 ggn(2) = ggu(is3+2)
 ggn(3) = ggu(is3+3)
 
 is3    = 3*(ismin-1)
 ggn(4) = ggu(is3+1)
 ggn(5) = ggu(is3+2)
 ggn(6) = ggu(is3+3)
 
 is3    = 3*(middl-1)
 ggn(7) = ggu(is3+1)
 ggn(8) = ggu(is3+2)
 ggn(9) = ggu(is3+3)
 
 CALL betrnd (teb,ggn,0,elid)
 GO TO 300
 
!     FOR QUAD
!     THE ORIGIN OF THE ELEMENT COORD.SYSTEM IS IN THE MIDDLE OF THE
!     ELEMENT
 
 200 DO  j = 1,3
   cent(j) = 0.0D0
   DO  i = 1,mm
     cent(j) = cent(j) + ugpdm(j,i)/nnode
   END DO
 END DO
 
!     STORE THE CORNER NODE DIFF. IN THE USER COORD. SYSTEM
 
 x31 = ugpdm(1,3) - ugpdm(1,1)
 y31 = ugpdm(2,3) - ugpdm(2,1)
 x42 = ugpdm(1,4) - ugpdm(1,2)
 y42 = ugpdm(2,4) - ugpdm(2,2)
 aa  = x31*x31 + y31*y31
 IF (aa <= 0.0D0) GO TO 700
 aa  = DSQRT(aa)
 bb  = x42*x42 + y42*y42
 IF (bb <= 0.0D0) GO TO 700
 bb  = DSQRT(bb)
 
!     NORMALIZE XIJ'S
 
 x31 = x31/aa
 y31 = y31/aa
 x42 = x42/bb
 y42 = y42/bb
 exi = x31 - x42
 exj = y31 - y42
 
!     STORE GGE ARRAY, THE OFFSET BETWEEN ELEMENT COORD. SYSTEM AND USER
!     COORD. SYSTEM
 
 gge(1) = cent(1)
 gge(2) = cent(2)
 gge(3) = cent(3)
 
 gge(4) = gge(1) + exi
 gge(5) = gge(2) + exj
 gge(6) = gge(3)
 
 gge(7) = gge(1) - exj
 gge(8) = gge(2) + exi
 gge(9) = gge(3)
 
 CALL betrnd (teu,gge,0,elid)
 CALL gmmatd (teu,3,3,0, tub,3,3,0,  teb)
 CALL gmmatd (tub,3,3,1, cent,3,1,0, cente)
 
!     THE ARRAY IORDER STORES THE ELEMENT NODE ID IN INCREASING SIL
!     ORDER.
 
!     IORDER(1)  = NODE WITH LOWEST  SIL NUMBER
!     IORDER(MM) = NODE WITH HIGHEST SIL NUMBER
 
!     ELEMENT NODE NUMBER IS THE INTEGER FROM THE NODE LIST
!     G1,G2,G3,G4,G5,G6,G7,G8 .  THAT IS, THE "I" PART OF THE "GI" AS
!     THEY ARE LISTED ON THE CONNECTIVITY BULK DATA CARD DESCRIPTION.
 
 300 ksild = 99999995
 DO  i = 1,mm
   iorder(i) = 0
   iordrn(i) = 0
   ksil(i) = sil(i)
   IF (sil(i) /= 0) CYCLE
   ksil(i) = ksild
   ksild   = ksild + 1
 END DO
 DO  i = 1,mm
   itemp = 1
   isil  = ksil(1)
   DO  j = 2,mm
     IF (isil <= ksil(j)) CYCLE
     itemp = j
     isil  = ksil(j)
   END DO
   iorder(i) = itemp
   iordrn(i) = itemp
   ksil(itemp) = 99999999
 END DO
 
!     ADJUST EST DATA
 
!     USE THE POINTERS IN IORDER TO COMPLETELY REORDER THE GEOMETRY DATA
!     INTO INCREASING SIL ORDER.
!     DON'T WORRY!! IORDER ALSO KEEPS TRACK OF WHICH SHAPE FUNCTIONS GO
!     WITH WHICH GEOMETRIC PARAMETERS!
 
 DO  i = 1,mm
   ksil(i)  = sil(i)
   tmpthk(i)= gpth(i)
   IF (mm /= 4) temtem(i) = gptemp(i)
   kcid(i) = igpdt(1,i)
   DO  j = 2,4
     tgrid(j,i) = bgpdt(j,i)
   END DO
 END DO
 DO  i = 1,mm
   ipoint  = iorder(i)
   sil(i)  = ksil(ipoint)
   nsil(i) = ksil(ipoint)
   gpth(i) = tmpthk(ipoint)
   IF (mm /= 4) gptemp(i) = temtem(ipoint)
   igpdt(1,i) = kcid(ipoint)
   DO  j = 2,4
     bgpdt(j,i) = tgrid(j,ipoint)
   END DO
 END DO
 
 IF (quad) GO TO 500
 
!     FOR TRIA
!     CREATE THE INTERNAL ORDER OF THE NODES OF ELEMENT IN CONNECTION
!     WITH THE INTERNAL COORDINATE SYSTEM THEN CALCULATE NORMALS
 
 DO  i = 1,mm
   IF (iorder(i) == ismax) iordrn(i) = 1
   IF (iorder(i) == ismin) iordrn(i) = 2
   IF (iorder(i) == middl) iordrn(i) = 3
   IF (iorder(i) == 4) ind4=i
   IF (iorder(i) == 5) ind5=i
   IF (iorder(i) == 6) ind6=i
 END DO
 IF (mm /= 6) GO TO 410
 IF (ismax+ismin == 3) iordrn(ind4) = 4
 IF (ismax+ismin == 4) iordrn(ind6) = 4
 IF (ismax+ismin == 5) iordrn(ind5) = 4
 IF (ismin+middl == 3) iordrn(ind4) = 5
 IF (ismin+middl == 4) iordrn(ind6) = 5
 IF (ismin+middl == 5) iordrn(ind5) = 5
 IF (middl+ismax == 3) iordrn(ind4) = 6
 IF (middl+ismax == 4) iordrn(ind6) = 6
 IF (middl+ismax == 5) iordrn(ind5) = 6
 
 410 DO  i = 1,3
   ii = i + 1
   ip = (i-1)*3
   DO  j = 1,nnode
     egpdt(ii,j) = 0.0D0
     DO  k = 1,3
       kk = ip + k
       egpdt(ii,j) = egpdt(ii,j) + teb(kk)*(DBLE(bgpdt(k+1,j))-ggn(k))
     END DO
   END DO
 END DO
 
!     USE THE POINTERS IN IORDER AND IORDRN TO REORDER MMN
 
 DO  i = 1,mm
   ipoint = iordrn(i)
   jpoint = iorder(i)
   mmn(ipoint) = ksil(jpoint)
 END DO
 
 IF (mm /= 3) GO TO 520
 DO  ii=1,3
   epnorm(1,ii) = 0.0D0
   epnorm(2,ii) = 0.0D0
   epnorm(3,ii) = 0.0D0
   epnorm(4,ii) = 1.0D0
   gpnorm(1,ii) = 0.0D0
   gpnorm(2,ii) = teb(7)
   gpnorm(3,ii) = teb(8)
   gpnorm(4,ii) = teb(9)
 END DO
 GO TO 520
 
!     FOR QUAD - COMPUTE NODAL NORMALS
!     THE COORDINATES OF THE ELEMENT GRID POINTS HAVE TO BE TRANSFORMED
!     FROM THE BASIC COORD. SYSTEM TO THE ELEMENT COORD. SYSTEM
 
 500 iflag = 0
 IF (mm == 4) CALL q4nrmd (bgpdt,gpnorm,iorder,iflag)
 IF (iflag /= 0) GO TO 700
 
 DO  i = 1,3
   ii = i + 1
   ip = (i-1)*3
   DO  j = 1,nnode
     epnorm(ii,j) = 0.0D0
     egpdt (ii,j) = 0.0D0
     DO  k = 1,3
       kk = ip + k
       k1 = k  + 1
       cc = DBLE(bgpdt(k1,j)) - ggu(k) - cente(k)
       epnorm(ii,j) = epnorm(ii,j) + teb(kk)*gpnorm(k1,j)
       egpdt (ii,j) = egpdt (ii,j) + teb(kk)*cc
     END DO
   END DO
 END DO
 
!     SET AVGTHK TO ZERO
 
 520 avgthk = 0.0D0
 DO  i = 1,nnode
   io = iorder(i)
   IF (io > mmx) CYCLE
   
   IF (gpth(i) < 0.0) THEN
     GO TO   700
   ELSE IF (gpth(i) == 0.0) THEN
     GO TO   530
   ELSE
     GO TO   540
   END IF
   530 IF (elth <= 0.0) GO TO 700
   gpth(i)  = elth
   540 dgpth(i) = DBLE(gpth(i))
   avgthk = avgthk + dgpth(i)/nnode
 END DO
 
 DO  i = 1,nnode
   io = iorder(i)
   IF (io      <= mmx) CYCLE
   IF (gpth(i) > 0.0) GO TO 610
   io1 = io  - mmx
   io2 = io1 + 1
   IF (io2 == mmx+1) io2 = 1
   DO  j = 1,mm
     jo = iorder(j)
     IF (jo == io1) ic1 = j
     IF (jo == io2) ic2 = j
   END DO
   gpth (i) = (gpth(ic1)+gpth(ic2))/2.0
   610 dgpth(i) = DBLE(gpth(i))
   avgthk = avgthk + dgpth(i)/nnode
 END DO
 RETURN
 
 700 RETURN 1
END SUBROUTINE shsetd
