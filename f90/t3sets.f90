SUBROUTINE t3sets (ierr,sil,jgpdt,elth,gpth,dgpth,egpdt,gpnorm,  &
        epnorm,iorder,teb,tub,cent,avgthk,lx,ly,edglen,elid)
     
!     SINGLE PRECISION ROUTINE TO DO THE SET-UP FOR TRIA3 ELEMENTS
 
 
!     INPUT :
!           SIL    - ARRAY OF SIL NUMBERS
!           JGPDT  - BGPDT DATA (INTEGER ARRAY)
!           ELTH   - ELEMENT THICKNESS FROM EPT
!           GPTH   - GRID POINT THICKNESS DATA
!           ELID   - ELEMENT ID
!     OUTPUT:
!           IERR   - ERROR FLAG
!           SIL    - ARRAY OF SIL NUMBERS       (REARRANGED)
!           JGPDT  - BGPDT DATA (INTEGER ARRAY) (REARRANGED)
!           GPTH   - GRID POINT THICKNESS DATA  (REARRANGED)
!           DGPTH  - GRID POINT THICKNESS DATA  (HIGH PREC)
!           EGPDT  - BGPDT DATA IN ELEMENT COORD. SYSTEM
!           GPNORM - GRID POINT NORMALS
!           EPNORM - GRID POINT NORMALS IN ELEMENT COORD .SYSTEM
!           IORDER - ARRAY OF ORDER INDICATORS FOR REARRANGED DATA
!           TEB    - TRANSFORMATION FROM ELEMENT TO BASIC COORD. SYSTEM
!           TUB    - TRANSFORMATION FROM USER TO BASIC COORD. SYSTEM
!           CENT   - LOCATION OF THE CENTER OF THE ELEMENT
!           AVGTHK - AVERAGE THICKNESS OF THE ELEMENT
!           LX     - DIMENSION OF ELEMENT ALONG X-AXIS
!           LY     - DIMENSION OF ELEMENT ALONG Y-AXIS
!           EDGLEN - EDGE LENGTHS
 
 
 
 INTEGER, INTENT(OUT)                     :: ierr
 INTEGER, INTENT(IN OUT)                  :: sil(3)
 INTEGER, INTENT(IN OUT)                  :: jgpdt(4,3)
 REAL, INTENT(IN)                         :: elth
 REAL, INTENT(IN OUT)                     :: gpth(3)
 REAL, INTENT(OUT)                        :: dgpth(3)
 REAL, INTENT(OUT)                        :: egpdt(4,3)
 REAL, INTENT(OUT)                        :: gpnorm(4,3)
 REAL, INTENT(OUT)                        :: epnorm(4,3)
 INTEGER, INTENT(OUT)                     :: iorder(3)
 REAL, INTENT(OUT)                        :: teb(9)
 REAL, INTENT(IN OUT)                     :: tub(9)
 REAL, INTENT(OUT)                        :: cent(3)
 REAL, INTENT(OUT)                        :: avgthk
 REAL, INTENT(OUT)                        :: lx
 REAL, INTENT(OUT)                        :: ly
 REAL, INTENT(OUT)                        :: edglen(3)
 INTEGER, INTENT(IN OUT)                  :: elid
 INTEGER :: igpdt(4,3), igrid(4,3),  ksil(3)
 REAL :: bgpdt(4,3), tmpthk(3)
 REAL :: ggu(9), cc,   &
     area2,length,small,x(3),y(3),z(3),edg12(3), edg23(3),edg13(3), axis(3,3)
 EQUIVALENCE      (igpdt(1,1),bgpdt(1,1))
 
 
!     INITIALIZE
 
 ierr  = 0
 nnode = 3
 
 DO  i = 1,nnode
   DO  j = 1,4
     igpdt(j,i) = jgpdt(j,i)
   END DO
 END DO
 
!     SET UP THE USER COORDINATE SYSTEM
 
 DO  i = 1,3
   ii = (i-1)*3
   DO  j = 1,3
     ggu(ii+j) = bgpdt(j+1,i)
   END DO
 END DO
 CALL betrns (tub,ggu,0,elid)
 
!     SET UP THE ELEMENT COORDINATE SYSTEM
 
!     1. SET UP THE EDGE VECTORS AND THEIR LENGTHS
 
 DO  i = 1,nnode
   x(i) = bgpdt(2,i)
   y(i) = bgpdt(3,i)
   z(i) = bgpdt(4,i)
 END DO
 
 cent(1) = (x(1)+x(2)+x(3))/3.0
 cent(2) = (y(1)+y(2)+y(3))/3.0
 cent(3) = (z(1)+z(2)+z(3))/3.0
 
 edg12(1) = x(2) - x(1)
 edg12(2) = y(2) - y(1)
 edg12(3) = z(2) - z(1)
 edglen(1)= edg12(1)**2 + edg12(2)**2 + edg12(3)**2
 IF (edglen(1) == 0.0) GO TO 380
 edglen(1) = SQRT(edglen(1))
 
 edg23(1) = x(3) - x(2)
 edg23(2) = y(3) - y(2)
 edg23(3) = z(3) - z(2)
 edglen(2)= edg23(1)**2 + edg23(2)**2 + edg23(3)**2
 IF (edglen(2) == 0.0) GO TO 380
 edglen(2) = SQRT(edglen(2))
 
 edg13(1) = x(3) - x(1)
 edg13(2) = y(3) - y(1)
 edg13(3) = z(3) - z(1)
 edglen(3)= edg13(1)**2 + edg13(2)**2 + edg13(3)**2
 IF (edglen(3) == 0.0) GO TO 380
 edglen(3) = SQRT(edglen(3))
 
!     2. FIND THE SMALLEST EDGE LENGTH
 
 small = edglen(1)
 nodei = 3
 nodej = 1
 nodek = 2
 
 IF (edglen(2) >= small) GO TO 160
 small = edglen(2)
 nodei = 1
 nodej = 2
 nodek = 3
 160 IF (edglen(3) >= small) GO TO 180
 small = edglen(3)
 nodei = 2
 nodej = 1
 nodek = 3
 
!     3. ESTABLISH AXIS 3 AND NORMALIZE IT
 
 180 CALL saxb (edg12,edg13,axis(1,3))
 
 length = SQRT(axis(1,3)**2 + axis(2,3)**2 + axis(3,3)**2)
 axis(1,3) = axis(1,3)/length
 axis(2,3) = axis(2,3)/length
 axis(3,3) = axis(3,3)/length
 area2     = length
 
!     4. ESTABLISH AXES 1 AND 2 AND NORMALIZE THEM
 
 axis(1,1) = (x(nodej)+x(nodek))/2.0 - x(nodei)
 axis(2,1) = (y(nodej)+y(nodek))/2.0 - y(nodei)
 axis(3,1) = (z(nodej)+z(nodek))/2.0 - z(nodei)
 
 length = SQRT(axis(1,1)**2 + axis(2,1)**2 + axis(3,1)**2)
 axis(1,1) = axis(1,1)/length
 axis(2,1) = axis(2,1)/length
 axis(3,1) = axis(3,1)/length
 
 CALL saxb (axis(1,3),axis(1,1),axis(1,2))
 
 DO  i = 1,3
   teb(i  ) = axis(i,1)
   teb(i+3) = axis(i,2)
   teb(i+6) = axis(i,3)
 END DO
 
 lx = length
 ly = area2/lx
 
 
!     THE ELEMENT COORDINATE SYSTEM IS NOW READY
 
!     THE ARRAY IORDER STORES THE ELEMENT NODE ID IN INCREASING SIL
!     ORDER.
 
!     IORDER(1) = NODE WITH LOWEST  SIL NUMBER
!     IORDER(3) = NODE WITH HIGHEST SIL NUMBER
 
!     ELEMENT NODE NUMBER IS THE INTEGER FROM THE NODE LIST  G1,G2,....
!     THAT IS, THE 'I' PART OF THE 'GI' AS THEY ARE LISTED ON THE
!     CONNECTION BULK DATA CARD DESCRIPTION.
 
 DO  i = 1,nnode
   ksil(i) = sil(i)
 END DO
 
 DO  i = 1,nnode
   itemp = 1
   isil = ksil(1)
   DO  j = 2,nnode
     IF (isil <= ksil(j)) CYCLE
     itemp = j
     isil = ksil(j)
   END DO
   iorder(i) = itemp
   ksil(itemp) = 99999999
 END DO
 
!     USE THE POINTERS IN IORDER TO COMPLETELY REORDER THE GEOMETRY DATA
!     INTO INCREASING SIL ORDER.
 
 DO  i = 1,nnode
   ksil(i) = sil(i)
   tmpthk(i) = gpth(i)
   DO  j = 1,4
     igrid(j,i) = igpdt(j,i)
   END DO
 END DO
 DO  i = 1,nnode
   ipoint = iorder(i)
   sil(i) = ksil(ipoint)
   gpth(i)= tmpthk(ipoint)
   DO  j = 1,4
     igpdt(j,i) = igrid(j,ipoint)
     jgpdt(j,i) = igpdt(j,i)
   END DO
 END DO
 
!     THE COORDINATES OF THE ELEMENT GRID POINTS MUST BE TRANSFORMED
!     FROM THE BASIC COORD. SYSTEM TO THE ELEMENT COORD. SYSTEM
 
 DO  i = 1,3
   ip = (i-1)*3
   DO  j = 1,nnode
     egpdt(i+1,j) = 0.0
     DO  k = 1,3
       cc = bgpdt((k+1),j) - cent(k)
       egpdt(i+1,j) = egpdt(i+1,j) + teb(ip+k)*cc
     END DO
   END DO
 END DO
 
!     SET NODAL NORMALS
 
 DO  i = 1,nnode
   epnorm(1,i) = 0.0
   epnorm(2,i) = 0.0
   epnorm(3,i) = 0.0
   epnorm(4,i) = 1.0
   gpnorm(1,i) = 0.0
   gpnorm(2,i) = teb(7)
   gpnorm(3,i) = teb(8)
   gpnorm(4,i) = teb(9)
 END DO
 
!     SET NODAL THICKNESSES
 
 avgthk = 0.0
 DO  i = 1,nnode
   IF (gpth(i) < 0.0) THEN
     GO TO   380
   ELSE IF (gpth(i) == 0.0) THEN
     GO TO   350
   ELSE
     GO TO   360
   END IF
   350 IF (elth <= 0.0) GO TO 380
   gpth(i) = elth
   360 dgpth(i) = gpth(i)
   avgthk = avgthk + dgpth(i)/nnode
 END DO
 GO TO 400
 
 380 ierr = 1
 400 RETURN
END SUBROUTINE t3sets
