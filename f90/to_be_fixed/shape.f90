SUBROUTINE shape (*,gplst,x,u,pen,deform,iopt,iptl,liplt,opcor)
     
!     IOPT CONTROLS THIS ROUTINE
!     IOPT .LT. 0
!       THE LINEL ARRAY WAS NOT CREATED.  UNIQUE LINES ARE NOT DRAWN.
!     IOPT .GE. 0
!       THE LIPLT ARRAY HAS CONNECTION DATA TO MAKE UNIQUE LINES. SUPLT
!       WILL CREATE THE LINES.  IPTR IS ONE OF THE PARAMETERS.
 
!     REVISED 10/1990 BY G.CHAN/UNISYS, TO INCLUDE BAR, TRIA3 AND QUAE4
!     OFFSET (PEDGE = 3)
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: pen
 INTEGER, INTENT(IN)                      :: deform
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(IN OUT)                  :: iptl
 INTEGER, INTENT(OUT)                     :: liplt(1)
 INTEGER, INTENT(IN OUT)                  :: opcor
 INTEGER :: etyp,ect,NAME(2),gp, elid, te,br,t3,q4,offset,pedge
 REAL :: off(6)
 COMMON /BLANK / ngp,skp1(9),skp2(2),ect,skp21(7),merr
 COMMON /drwdat/ skp15(15),pedge
 DATA    te,br , t3,q4  / 2HTE, 2HBR, 2HT3, 2HQ4 /
 DATA    NAME  / 4HSHAP,1HE/
 
 CALL line (0,0,0,0,0,-1)
 IF (iopt >= 0) GO TO 120
 10 CALL READ (*140,*130,ect,etyp,1,0,i)
 offset = 0
 IF (etyp == br) offset = 6
 IF (etyp == t3 .OR. etyp == q4) offset = 1
 CALL fread (ect,i,1,0)
 ngpel = IABS(i)
 IF (etyp /= te .AND. ngpel < 5) GO TO 30
 
!     NOT A SIMPLE ELEMENT
 
 20 lgpel = ngpel
 ltyp  = etyp
 CALL linel (liplt,ltyp,opcor,lgpel,x,pen,deform,gplst)
 l = IABS(ltyp)
 IF (l-1 < 0) THEN
   GO TO    10
 ELSE IF (l-1 == 0) THEN
   GO TO   115
 ELSE
   GO TO    70
 END IF
 
 30 l = ngpel + 1
 IF (ngpel <= 2 .OR. i < 0) l = ngpel
 ltyp = 10000
 m = 1
 40 CALL fread (ect,elid,1,0)
 IF (elid <= 0) GO TO 10
 CALL fread (ect,lid,1,0)
 CALL fread (ect,liplt,ngpel,0)
 IF (l /= ngpel) liplt(l) = liplt(1)
 
 IF (offset /= 0) CALL fread (ect,off,offset,0)
 IF (offset == 0 .OR. deform /= 0) GO TO 70
 
!     IF THIS IS A BAR, TRIA3 OR QUAD4 ELEMENTS READ IN THE OFFSET
!     NO SCALE FACTOR APPLIES TO OFFSET HERE
 
 IF (offset /= 6) GO TO 50
 
!     BAR OFFSET
 
 off(1) = 0.707*SQRT(off(1)**2 + off(2)**2 + off(3)**2)
 off(2) = 0.707*SQRT(off(4)**2 + off(5)**2 + off(6)**2)
 off(3) = off(1)
 GO TO 70
 
!     TRIA3 AND QUAD4 OFFSET
 
 50 off(1) = 0.707*off(1)
 DO  i = 2,5
   off(i) = off(1)
 END DO
 
!     WRITE THE LINES.  0 FOR SIL MEANS NO LINES DRAWN
 
 70 j = 0
 DO  i = 1,l
   IF (j == 0) GO TO 80
   x1 = x2
   y1 = y2
   80 gp = liplt(i)
   IF (gp == 0) GO TO 110
   gp = IABS(gplst(gp))
   IF (deform /= 0) GO TO 90
   x2 = x(2,gp)
   y2 = x(3,gp)
   IF (offset == 0) GO TO 100
   
!     IF OFFSET IS PRESENT, ADD ARBITRARY AN OFFSET LENGTH TO X2 AND Y2.
!     SINCE THE OFFSET LENGTH IS SO TINY, ITS TRUE DIRECTION IS NOT OF
!     VITAL CONCERNS. THE IDEA HERE IS THAT BIG OFFSET WILL SHOW IN THE
!     PLOT IF ORIGINAL DATA CONTAINS ERRONEOUS AND BIG OFFSET VALUE(S).
   
!     IF OFFSET IS ADDED IN SAME DIRECTION AS THE PLOTTED LINE, ROTATE
!     THE OFFSET LENGTH BY 90 DEGREE
   
   x2 = x2 + off(i)
   xy = xy + off(i)
   IF (ABS((x2-x1)-(y2-y1)) < 0.01) x2 = x2 - 2.*off(i)
   GO TO 100
   90 x2 = u(1,gp)
   y2 = u(2,gp)
   100 IF (j == 0 .OR. j == gp) GO TO 110
   CALL line (x1,y1,x2,y2,pen,0)
   110 j = gp
 END DO
 
 115 IF (l-ltyp < 0) THEN
   GO TO    40
 ELSE IF (l-ltyp == 0) THEN
   GO TO    10
 ELSE
   GO TO    20
 END IF
 
 
 120 IF (pedge == 3) GO TO 130
 CALL suplt (liplt(1),liplt(iptl+1),x,u,gplst,pen,deform)
 130 CALL line (0,0,0,0,0,1)
 IF (iopt < 0) CALL bckrec (ect)
 GO TO 150
 
!     ILLEGAL EOF
 
 140 CALL mesage (-2,ect,NAME)
 150 IF (pedge == 3) RETURN 1
 RETURN
END SUBROUTINE shape
