SUBROUTINE flface (TYPE,ect,elt,grid)
     
!     LOCATES THE FLUID GRID POINTS DEFINING THE FACE OF A FLUID
!     ELEMENT.  THE FACE MAY BE SPECIFIED IN TWO MANNERS.
 
!     1) FACE NUMBER - ELT(2) LESS THEN ZERO AND FACE = ELT(3)
!     2) STRUCTURAL ELEMENT WHICH COINCIDES WITH FACE -
!                         ELT(2) = ELEMENT ID AND ELT(3)-ELT(6) = GRIDS
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 INTEGER, INTENT(IN)                      :: ect(10)
 INTEGER, INTENT(IN OUT)                  :: elt(7)
 INTEGER, INTENT(OUT)                     :: grid(4)
 LOGICAL :: error
 INTEGER :: gf(10)   , gs1      ,gs2      ,gs3      ,gs4      ,gf1      ,  &
     gf2      ,gf3      ,gf4      ,gridf(4) ,nface(4) ,  &
     faceid   ,hex1(4,6),hex2(4,6),tetra(4,4)         ,  &
     wedge(4,5)         ,face(4,6,4)
 REAL :: mag      ,r1(3)    ,r2(3)    ,r3(3)    ,ks(3)    ,  &
     kf(3)    ,cs(3)    ,angle(6) ,z        ,heigth(6)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /flbptr/ error    ,icore    ,lcore    ,ibgpdt   ,nbgpdt
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf   ,nout
 EQUIVALENCE     (hex1(1,1),face(1,1,1)),  (hex2(1,1),face(1,1,2)),  &
     (tetra(1,1),face(1,1,3)), (wedge(1,1),face(1,1,4))
 
!     DATA DEFINING FACES OF THE FLUID ELEMENTS
 
!     NUMBER OF GRID POINTS PER ELEMENT
 
!                   FHEX1     FHEX2     FTETRA    FWEDGE
 
 DATA    gridf / 8        ,8        ,4        ,6       /
 
!     NUMBER OF FACES ON EACH ELEMENT
 
 DATA    nface / 6        ,6        ,4        ,5       /
 
!     GRID POINTS WHICH DEFINE FACES FOR FHEX1 ELEMENTS
 
 DATA    hex1  / 1        ,4        ,3        ,2        ,  &
     1        ,2        ,6        ,5        ,  &
     2        ,3        ,7        ,6        ,  &
     3        ,4        ,8        ,7        ,  &
     4        ,1        ,5        ,8        ,  &
     5        ,6        ,7        ,8        /
 
!     GRID POINTS WHICH DEFINE FACES FOR FHEX2 ELEMENTS
 
 DATA    hex2  / 1        ,4        ,3        ,2        ,  &
     1        ,2        ,6        ,5        ,  &
     2        ,3        ,7        ,6        ,  &
     3        ,4        ,8        ,7        ,  &
     4        ,1        ,5        ,8        ,  &
     5        ,6        ,7        ,8        /
 
!     GRID POINTS WHICH DEFINE FACES FOR FTETRA ELEMENTS
 
 DATA    tetra / 1        ,3        ,2        ,-1       ,  &
     1        ,2        ,4        ,-1       ,  &
     2        ,3        ,4        ,-1       ,  &
     3        ,1        ,4        ,-1       /
 
!     GRID POINTS WHICH DEFINE FACES FOR FWEDGE ELEMENTS
 
 DATA    wedge / 1        ,3        ,2        ,-1       ,  &
     1        ,2        ,5        ,4        ,  &
     2        ,3        ,6        ,5        ,  &
     3        ,1        ,4        ,6        ,  &
     4        ,5        ,6        ,-1       /
 
 
!     DETERMINE HOW THE FACE IS SPECIFIED
 
!     SUBTRACT IFP CARD NUMBER OF ELEMENT JUST BEFORE CFHEX1 FROM TYPE
 intype = TYPE - 332
 
 nf = nface(intype)
 IF (elt(2) < 0) GO TO 200
 
!     THE FACE IS DEFINED BY STRUCTURAL GRIDS
 
!     INITIALIZE POINTERS TO GRID POINT DATA
 
 ngrids = 4
 IF (elt(6) < 0) ngrids = 3
 gs1 = ibgpdt + (elt(3)-1)*4
 gs2 = ibgpdt + (elt(4)-1)*4
 gs3 = ibgpdt + (elt(5)-1)*4
 gs4 = -1
 IF (ngrids == 4) gs4 = ibgpdt + (elt(6)-1)*4
 
 ngridf = gridf(intype)
 DO  i = 1,ngridf
   gf(i) = ibgpdt + (ect(i+2)-1)*4
 END DO
 
!     FIND NORMAL VECTOR TO STRUCTURAL ELEMENT FACE
 
 DO  i = 1,3
   r1(i) = z(gs2+i) - z(gs1+i)
   r2(i) = z(gs3+i) - z(gs1+i)
 END DO
 
 ks(1) = r1(2)*r2(3) - r1(3)*r2(2)
 ks(2) = r1(3)*r2(1) - r1(1)*r2(3)
 ks(3) = r1(1)*r2(2) - r1(2)*r2(1)
 
 mag = SQRT(ks(1)**2 + ks(2)**2 + ks(3)**2)
 IF (mag < 1.0E-7) GO TO 8005
 DO  i = 1,3
   ks(i) = ks(i)/mag
 END DO
 
!     FIND AREA OF STRUCTURE FACE AND TOLERANCE USED IN CHECKING
!     SEPERATIOON
 
 area = mag
 IF (gs4 < 0) area = mag/2.0
 tol = .2*SQRT(area)
 
!     FIND CENTROID OF STRUCTURAL FACE
 
 DO  i = 1,3
   cs(i) = z(gs1+i) + z(gs2+i) + z(gs3+i)
   IF (ngrids == 4) cs(i) = cs(i) + z(gs4+i)
   cs(i) = cs(i)/FLOAT(ngrids)
 END DO
 
!     PROCESS EACH FACE OF THE FLUID ELEMENT - FIRST GET GRID POINTERS
!     POINTERS
 
 DO  IF = 1,nf
   i   = face(1,IF,intype)
   gf1 = gf(i)
   i   = face(2,IF,intype)
   gf2 = gf(i)
   i   = face(3,IF,intype)
   gf3 = gf(i)
   i   = face(4,IF,intype)
   gf4 = -1
   IF (i > 0) gf4 = gf(i)
   
!     FIND NORMAL TO FLUID FACE
   
   DO  i = 1,3
     r2(i) = z(gf2+i) - z(gf1+i)
     r3(i) = z(gf3+i) - z(gf1+i)
   END DO
   
   kf(1) = r2(2)*r3(3) - r2(3)*r3(2)
   kf(2) = r2(3)*r3(1) - r2(1)*r3(3)
   kf(3) = r2(1)*r3(2) - r2(2)*r3(1)
   
   mag = SQRT(kf(1)**2 + kf(2)**2 + kf(3)**2)
   IF (mag < 1.0E-7) GO TO 8006
   DO  i = 1,3
     kf(i) = kf(i)/mag
   END DO
   
!     DETERMINE ANGLE BETWEEN FACES
   
   angle(IF) = kf(1)*ks(1) + kf(2)*ks(2) + kf(3)*ks(3)
   IF (ABS(angle(IF)) <= .866) CYCLE
   
!     FIND DISTANCE FROM THE CENTROID OF THE STRUCTURE TO THE FLUID
!     FACE.  THE DISTANCE IS MEASURED ALONG THE NORMAL TO THE
!     FLUID FACE
   
   DO  i = 1,3
     r2(i) = cs(i) - z(gf1+i)
   END DO
   
   heigth(IF) = ABS(kf(1)*r2(1) + kf(2)*r2(2) + kf(3)*r2(3))
   
 END DO
 
!     CHOSE THE FACE OF THE FLUID WITH THE SMALLEST DISTANCE TO THE
!     STRUCTURAL ELEMENT AND WITH THE ANGLE BETWEEN THE TWO FACES LESS
!     THEN 30 DEGREES
 
 dist   = 1.0E+10
 faceid = 0
 DO  IF = 1,nf
   IF (ABS(angle(IF)) <= .866) CYCLE
   IF (heigth(IF) >= dist) CYCLE
   dist   = heigth(IF)
   faceid = IF
 END DO
 IF (faceid == 0) GO TO 8007
 
!     VERIFY THAT THE FACE IS WITHIN PROPER TOLERENCE
 
 IF (dist > tol) GO TO 8008
 
!     IF ANGLE WAS COMPUTED NEGATIVE - SWICTH STRUCTURAL GRIDS AROUND
!     IN ELEMENT TABLE RECORD FOR LATER USE
 
 IF (angle(faceid) >= 0.0) GO TO 300
 IF (ngrids == 3) GO TO 120
 i      = elt(3)
 elt(3) = elt(6)
 elt(6) = i
 i      = elt(4)
 elt(4) = elt(5)
 elt(5) = i
 GO TO 300
 
 120 i      = elt(3)
 elt(3) = elt(5)
 elt(5) = i
 GO TO 300
 
!     THE FACE IS DEFINED BY A FACE ID
 
 200 faceid = elt(3)
 IF (faceid < 1 .OR. faceid > nf) GO TO 8009
 
!     USING THE FACE SPECIFIES OR FOUND - RETURN THE PROPER
!     FLUID GRID POINTS
 
 300 DO  i = 1,4
   j = face(i,faceid,intype)
   IF (j > 0) GO TO 305
   grid(i) = -1
   CYCLE
   305 grid(i) = ect(j+2)
 END DO
 
 RETURN
 
!     ERROR CONDITIONS
 
!     BAD GEOMETRY FOR STRUCTURAL ELEMENT
 
 8005 WRITE (nout,9005) ufm,elt(2)
 GO TO 9000
 
!     BAD GEOMETRY FOR FLUID ELEMENT
 
 8006 WRITE (nout,9006) ufm,IF,ect(1)
 GO TO 9000
 
!     NO FACE WITHIN 30 DEGREES FO STRUCTURAL ELEMENT FACE
 
 8007 WRITE (nout,9007) ufm,ect(1),elt(2)
 GO TO 9000
 
!     FLUID ELEMENT IS NOT WITHIN TOLERENCE RANGE OF STRUCTURAL ELEMENT
 
 8008 WRITE (nout,9008) ufm,ect(1),elt(2)
 GO TO 9000
 
!     ILLEGAL FACE NUMBER
 
 8009 WRITE (nout,9009) ufm,faceid,ect(1)
 
 9000 error = .true.
 RETURN
 
!     ERROR FORMATS
 
 9005 FORMAT (a23,' 8005. BAD GEOMETRY DEFINED FOR STRUCTURAL ELEMENT', i9)
 9006 FORMAT (a23,' 8006. BAD GEOMETRY DEFINED FOR FACE',i9,  &
     ' OF FLUID ELEMENT',i9)
 9007 FORMAT (a23,' 8007. NO FACE ON FLUID ELEMENT',i9,  &
     ' IS WITHIN 30 DEGREES OF STRUCTURAL ELEMENT',i9)
 9008 FORMAT (a23,' 8008. THE DISTANCE BETWEEN FLUID ELEMENT',i9,  &
     ' AND STRUCTURAL ELEMENT',i9, /30X, 'IS GREATER THAN THE ALLOWED TOLERANCE.')
 9009 FORMAT (a23,' 8009. FACE',i9,' SPECIFIED FOR FLUID ELEMENT',i9,  &
     ' IS AN ILLEGAL VALUE.')
END SUBROUTINE flface
