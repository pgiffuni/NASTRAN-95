SUBROUTINE ofsplt (*,esym,elid,g,offset,x,deform,gplst)
     
!     CALLED ONLY BY LINEL TO PRCESS ELEMENT OFFSET PLOT
!     THIS ROUTINE DRAW THE CBAR, CTRIA3, AND CQUAD4, WITH OFFSET IN
!     PLACE.
 
!     INPUT:
!         ESYM   = BCD, SHOULD BE 'BR', 'T3', OR 'Q4'                BCD
!         ELID   = ELEMENT ID                                         I
!         G      = SIL LIST                                           I
!         OFFSET = 6 COORDINATES (GLOBAL) FOR CBAR,                   I
!                = 1 OFFSET, NORMAL TO PLATE, FOR CTRIA3 OR CQUAD4
!         X      = GRID POINT COORDINATE, ALREADY CONVERTED TO SCREEN
!                  (X-Y) COORDINATES                                  R
!         DEFORM = 0, FOR UNDEFORM PLOT,  .NE.0 FOR DEFORMED OR BOTH  I
!                  THIS ROUTINE WILL NOT PROCESS DEFORMED-OFFSET PLOT
!         OFFSCL = OFFSET MULTIPLICATION FACTOR                       I
!         PEDGE  = OFFSET PLOT FLAG                                   I
!                = 3, PLOT OFFSET ELEMENTS ONLY, SKIP OTHER ELEMNETS
!            NOT = 3, PLOT OFFSET ELEMENTS, RETURN TO PLOT OTHERS
!         PLABEL = FLAG FOR ELEM ID LABEL                             I
!         PEN    = PEN SELECTION, 1-31.  32-62 FOR COLOR FILL         I
!         OFFLAG = HEADING CONTROL                                    I
!         ELSET  = ECT DATA BLOCK. THIS DATA BLOCK WAS MODIFIED IN    I
!                  COMECT TO INCLUDE OFFSET DATA FOR BAR,TRIA3,QUAD4
!         GPLST  = A SUBSET OF GRID POINTS PERTAININGS TO THOSE GRID  I
!                  POINTS USED ONLY IN THIS PLOT
!     LOCAL:
!         SCALE  = REAL NUMBER OF OFFSCL
!         OFF    = OFFSET VALUES FROM ELEMENT DATA IN ELSET DATA BLOCK
!         PN1    = PEN COLOR FOR OFFSET LEG.
!                  IF PEN.GT.1, PN1 = PEN-1. IF PEN.LE.1, PN1 = PEN+1
!         NL     = NO. OF LINES TO BE DRAWN PER ELEMENT
!         DELX   = SMALL OFFSET FROM MIDDLE OF LINE FOR ELEM ID PRINTING
!         0.707  = AN ABITRARY FACTOR TO PUT OFFSET 45 DEGREE OFF GRID
!                  POINT
 
!     TWO METHODS
!     (1) PEDGE .NE. 3
!         AN OFFSET PLOT WITHOUT CONSIDERING ITS TRUE DIRECTION, OFFSET
!         VALUE(S) MAGNIFIED 20 TIMES
!     (2) PEDGE .EQ. 3
!         PLOT WITH TRUE OFFSET DIRECTIONS, AND PLOT, WITH COLOR OPTION,
!         GRID(A)-OFFSET(A)-OFFSET(B)-GRID(B)
!         OFFSET CAN BE SCALE UP BY USER VIA PLOT OFFSET COMMAND,
!         DEFAULT IS NO SCALE UP. (NEW 93)
 
!     A SYMBOL * IS ADDED AT THE TIP OF EACH OFFSET
!     CURRENTLY THE SYMBOLS KBAR,KT3 AND KQ4 ARE NOT USED
 
!     CURRENTLY ONLY CBAR (OFFSET=6), CTRIA3 AND CQUAD4 (OFFSET=1 BOTH)
!     HAVE OFFSET CAPABILITY
 
!     WRITTEN BY G.CHAN/UNISYS   10/1990
 
!     COMMENTS FORM G.C.  3/93
!     THE LOGIC IN COMPUTING THE TRUE OFFSET INVOLVING COORDINATE
!     TRANSFORMATION AT EACH POINT POINT SEEMS SHAKY. MAKE SURE THAT
!     AXIS AND SIGN DATA (FROM PROCES) ARE TRUELY AVAILBLE. ARE THE
!     GIRD POINT XYZ COORDINATES AT HAND IN GLOBAL ALREADY?
!     THE OFFSET PLOT IS QUESTIONABLE.
 
 
 INTEGER, INTENT(IN OUT)                  :: esym
 INTEGER, INTENT(IN OUT)                  :: elid
 INTEGER, INTENT(IN)                      :: g(3)
 INTEGER, INTENT(IN)                      :: offset
 REAL, INTENT(IN)                         :: x(3,1)
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER, INTENT(IN)                      :: gplst(1)
 INTEGER :: offhdg(5),sym(2),offscl,pen
 REAL :: off(3,2),v(3),cstm,SIGN,x1,x2,x3,y1,y2,y3,  &
     xmax,ymax,ymax1,cntx,cnty,cnty4,scale,delx, ofv(3,2)
 COMMON /BLANK / skp1(12),elset
 COMMON /system/ skp2,nout
 COMMON /rstxxx/ cstm(3,3),skp3(12),axis(3),SIGN(3)
 COMMON /xxparm/ skp4(235),offscl
 COMMON /drwdat/ skp5,plabel,skp6,pen,skp7(11),pedge,offlag
 COMMON /pltdat/ skp8(6),xmax,ymax,skp9(15),cntx,cnty
 DATA    kbar  , kt3,kq4 / 2HBR,2HT3,2HQ4 /, sym / 2,0 /
 
 DATA    offhdg/ 4H off,4HSET ,4HSCAL,4HE = ,4H   x    /
 
 CALL fread (elset,off,offset,0)
 
 IF (deform /= 0 .OR. offscl < 0) GO TO 200
 IF (pedge /= 3 .OR. offlag == 1) GO TO 20
 offlag= 1
 cnty4 = 4.*cnty
 ymax1 = ymax - cnty
 scale = 1.0
 IF (pedge /= 3) scale = 20.0
 IF (pedge == 3) scale = FLOAT(offscl)
 mpen  = MOD(pen,31)
 IF (mpen > 1) pn1  = mpen - 1
 IF (mpen <= 1) pn1  = mpen + 1
 
!     ADD OFFSET HEADER LINE
 
 CALL PRINT  (30.*cntx,ymax,1,offhdg,5,0)
 x1 = 48.
 IF (offscl >= 100) x1 = 47.
 CALL typint (x1*cntx,ymax,1,offscl,1,0)
 
 20 x1 = 0.0
 DO  i = 1,offset
   x1 = x1 + ABS(off(i,1))
   ofv(i,1) = off(i,1)
 END DO
 IF (ABS(x1) < 1.0E-7) GO TO 200
 
 nl = 1
 IF (esym == kt3) nl = 3
 IF (esym == kq4) nl = 4
 IF (pedge  /= 3) GO TO 150
 
 j    = ALOG10(FLOAT(elid)) + 1.0
 delx = (j+.03)*cntx
 
!     COMPUTE THE TRUE OFFSET DIRECTION IF PEDGE = 3,
!     OTHERWISE, JUST PLOT OFFSET AT 45 DEGREE
 
 IF (offset == 1) GO TO 90
 
!     CBAR, OFFSET = 6
!     CONVERT OFFSET FROM GLOBAL TO PLOT COORDINATES
 
!     AXIS AND SIGN DATA FROM SUBROUTINE PROCES
 
 DO  k = 1,2
   DO  i = 1,3
     j    = axis(i)
     v(j) = SIGN(i)*ofv(j,k)
   END DO
   DO  j = 1,3
     l = axis(j)
     x1 = 0.0
     DO  i = 1,3
       x1 = x1 + cstm(l,i)*v(i)
     END DO
     off(j,k) = x1*scale
   END DO
 END DO
 GO TO 110
 
!     CTRIA3 AND CQUAD4, OFFSET = 1
!     COMPUTE UNIT NORMAL TO THE PLATE BY CROSS PRODUCT, THEN
!     THE MAGNITUDE OF OFFSET
 
 90 i = g(1)
 j = g(2)
 k = g(3)
 i = gplst(i)
 j = gplst(j)
 k = gplst(k)
 v(1) = (x(2,j)-x(2,i))*(x(3,k)-x(3,i)) - (x(3,j)-x(3,i))*(x(2,k)-x(2,i))
 v(2) = (x(3,j)-x(3,i))*(x(1,k)-x(1,i)) - (x(1,j)-x(1,i))*(x(3,k)-x(3,i))
 v(3) = (x(1,j)-x(1,i))*(x(2,k)-x(2,i)) - (x(2,j)-x(2,i))*(x(1,k)-x(1,i))
 x1   = 0.5*SQRT(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
 v(2) = v(2)/x1
 v(3) = v(3)/x1
 off(2,1) = ofv(1,1)*v(2)*scale
 off(3,1) = ofv(1,1)*v(3)*scale
 off(2,2) = off(2,1)
 off(3,2) = off(3,1)
 
!     DRAW THE ELEMENT LINES AND ELEMENT ID
!     IF COLOR FILL IS REQUESTED, SET PEN TO ZERO ON THE LAST CLOSING-IN
!     EDGE (2- OR 3-DIMESIONAL ELEMENTS ONLY)
 
 110 DO  l = 1,nl
   i  = g(l  )
   j  = g(l+1)
   i  = gplst(i)
   j  = gplst(j)
   x1 = x(2,i)
   y1 = x(3,i)
   x2 = x(2,i) + off(2,1)
   y2 = x(3,i) + off(3,1)
   IF (x2 <   0.1) x2 = 0.1
   IF (x2 >  xmax) x2 = xmax
   IF (y2 < cnty4) y2 = cnty4
   IF (y2 > ymax1) y2 = ymax1
   CALL line (x1,y1,x2,y2,pn1,0)
   CALL symbol (x2,y2,sym,0)
   x3 = x(2,j) + off(2,2)
   y3 = x(3,j) + off(3,2)
   IF (x3 <   0.1) x3 = 0.1
   IF (x3 >  xmax) x3 = xmax
   IF (y3 < cnty4) y3 = cnty4
   IF (y3 > ymax1) y3 = ymax1
   ipen = pen
   IF (pen > 31 .AND. nl >= 3 .AND. l == nl) ipen = 0
   CALL line (x2,y2,x3,y3,ipen,0)
   
   IF (l > 1) CYCLE
   IF (plabel /= 3 .AND. plabel /= 6) GO TO 120
   IF (x2 >= x1) delx = -delx
   x1 = 0.5*(x3 + x2) + delx
   y1 = 0.5*(y3 + y2)
   CALL typint (x1,y1,1,elid,1,0)
   120 IF (nl > 1) CYCLE
   CALL symbol (x3,y3,sym,0)
   x2 = x(2,j)
   y2 = x(3,j)
   CALL line (x3,y3,x2,y2,pen,0)
 END DO
 GO TO 210
 
!     PLOT OFFSET WITHOUT CONSIDERING ITS TRUE OFFSET DIRECTION IN
!     GENERAL PLOT. (SEE 130 LOOP FOR ELEMENTS WITH COLOR FILL)
 
 150 IF (offset == 1) GO TO 160
 v(1) = off(1,1)*off(1,1) + off(2,1)*off(2,1) + off(3,1)*off(3,1)
 v(2) = off(1,2)*off(1,2) + off(2,2)*off(2,2) + off(3,2)*off(3,2)
 v(1) = 0.707*SQRT(v(1))
 v(2) = 0.707*SQRT(v(2))
 GO TO 170
 
 160 v(1) = 0.707*off(1,1)
 v(2) = v(1)
 
 170       v(1) = v(1)*scale
 v(2) = v(2)*scale
 DO  l = 1,nl
   i  = g(l  )
   j  = g(l+1)
   i  = gplst(i)
   j  = gplst(j)
   x1 = x(2,i) + v(1)
   y1 = x(3,i) + v(1)
   x2 = x(2,j) + v(2)
   y2 = x(3,j) + v(2)
   ipen = pen
   IF (pen > 31 .AND. nl >= 3 .AND. l == nl) ipen = 0
   CALL line (x1,y1,x2,y2,ipen,0)
   CALL symbol (x1,y1,sym,0)
   IF (nl == 1) CALL symbol (x2,y2,sym,0)
 END DO
 GO TO 210
 
 200 IF (pedge /= 3 .OR. offscl < 0) RETURN
 210 RETURN 1
END SUBROUTINE ofsplt
