SUBROUTINE draw (gplst,x,u,s,disp,stereo,opcor,buf1)
     
 
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 REAL, INTENT(OUT)                        :: x(3,1)
 REAL, INTENT(OUT)                        :: u(3,1)
 REAL, INTENT(OUT)                        :: s(2,1)
 LOGICAL, INTENT(IN OUT)                  :: disp
 INTEGER, INTENT(IN OUT)                  :: stereo
 INTEGER, INTENT(IN OUT)                  :: opcor
 INTEGER, INTENT(IN OUT)                  :: buf1
 EXTERNAL        andf
 
 INTEGER :: andf,axis(3), daxis,elset,gp,  &
     pcon,pedge,pen,plabel,pltnum,porig,ppen,prject,  &
     pset,pshape,psymbl,psymm,pvectr, sym(2), sucor,v,vec(3),color,defm
 REAL :: a(3),maxdef, SIGN(3), MIN,MAX
 DOUBLE PRECISION :: dr,sum
 COMMON /BLANK / ngp,skp11(3),pltnum,ngpset,skp12(4),skp21(2),elset
 COMMON /xxparm/ pbufsz,camera(5),nopens,papsiz(2),penpap(27),  &
     scale,objmod,skpscl,maxdef,defmax,axes(3),  &
     daxis(3),view(9),vantx1,r0,vantx2(3),d0,vantx3(2),  &
     prject,vantx4,origx1(14),edge(11,4),xy(11,3),  &
     ncntr(51),icntvl,skp24(24),color
 COMMON /pltdat/ model,ploter,reg(4),axymax(14)
 COMMON /rstxxx/ cstm(3,3),MIN(3),MAX(3)
 COMMON /drwdat/ pset,plabel,porig,ppen,pshape,psymbl(2),psymm(6),  &
     pvectr,pcon,pedge
 
!     /DRWDAT/ CONTROLS THIS ROUTINE
!     PLABEL - LABELING GRIDS, ELEMENTS...
!            -N = NONE
!             0 = GID             3 = EID          6 = EID + GID
!             1 = GID + SPC       4 = EID + PID
!             2 = UNDEFINED.      5 = UNDEFINED
!     PSHAPE - WHICH SHAPE OR OUTLINE OPTION TO DRAW...
!             1 = UNDEFORMED      2 = DEFORMED     3 = BOTH
!     PSYMBL(2) - DRAW SYMBOLS IF PSYMBL(1).NE.0
!     PSYMM (6) - SYMMETRY FLAGS...
!           (1) = X AXIS SIGN CHANGE   (4) = X DEFORMATION SIGN CHANGE
!           (2) = Y                    (5) = Y
!           (3) = Z                    (6) = Z
!     PVECTR - DEFORMATION VECTORS DRAWN (AS INTERPRETED BY INTVEC)...
!             0 = NONE
!             1 = X     4 = Z     7 = XYZ        10 = RY    13 = RXZ
!             2 = Y     5 = ZX    8 = UNDEFINED  11 = RXY   14 = RYZ
!             3 = XY    6 = ZY    9 = RX         12 = RZ    15 = RXYZ
!             THE NEGATIVE OF ABOVE, DO NOT DRAW SHAPE.
!    PCON   - NONZERO MEANS CONTOUR PLOT...
!    PEDGE  - 0 = SHAPE DRAWN,
!             1 = OUTLINE (BORDER) DRAWN ACCORDING TO PSHAPE-S
!             2 = HIDDEN LINE PLOT
!             3 = OFFSET PLOT
!             4 THRU N = SHRINK PLOT, ELEMENT SHRUNK BY THIS PERCENT
!             200 +  N = HIDDEN LINE AND SHRINK PLOT, N.GT.2
!             100 = FILL ?
 
!     OPCOR = NO. OF OPEN CORE WORDS AVAILABLE IN S
!             IT IS NOT A POINTER TO S, NOR A OPEN CORE ARRAY IN S
!     BUF1  = BUFFER AVAILABLE AT END OF CORE W.R.T. GPLST = BUFSIZ+1
 
!     OPEN CORE /ZZPLOT/
!     SETID NSETS NDOF      NGP 3*NGPSET 3*NGPSET  OPCOR   N
!     -----+-----+----+----+---+--------+--------+-------+--+--+-+-+-+-+
!          !          N1   N2  I1 (X)   I2 (U)   I3  (S)   DEFBUF..BUF..
!          !(DEFLST)       !
!                          !(GPLST)                      N=2*NGPSET
 
!     NGP    = TOTAL NO. OF GRID POINTS IN THE STRUCTURE
!     NGPSET = NO. OF GRID POINTS USED IN CURRENT SET OF PLOTTING
!     GPLST  = TABLE OF NGP IN LENGTH,
!              GPLST(I) = 0 IF THIS I-TH GRID POINT IS NOT USED FOR THE
!              CURRENT PLOT. OTHERWISE GPLST(I) IS NON-ZERO.
!     X      = X,Y,Z COORDINATES OF THE GRID POINTS CORRESPONDING TO THE
!              NON-ZERO GRID POINTS IN THE GPLST TABLE
!              TOTALLY, THERE ARE NGPSET GRID POINTS IN X
!     U      = X,Y,Z DISPLACEMENTS, ARRANGED SIMILARLY TO X
!     S      = SCRATCH AREA
 
 scalex = 1.0
 IF (prject == 3) scalex = objmod
 
!     SETUP THE PLOTTER REGION.
 
 IF (psymm(1) < 0 .OR. psymm(2) < 0 .OR. psymm(3) < 0) GO TO 10
 reg(1) = edge(porig,1)*axymax(1)
 reg(2) = edge(porig,2)*axymax(2)
 reg(3) = edge(porig,3)*axymax(1)
 reg(4) = edge(porig,4)*axymax(2)
 GO TO 20
 10 reg(1) = 0.0
 reg(2) = 0.0
 reg(3) = axymax(1)
 reg(4) = axymax(2)
 
!     REDUCE THE GRID POINT CO-ORDINATES TO PLOT SIZE + TRANSLATE TO
!     THE SELECTED ORIGIN.
 
 20 DO  i = 1,3
   MIN(i) = +1.e+20
   MAX(i) = -1.e+20
   IF (psymm(i) >= 0) CYCLE
   DO  gp = 1,ngpset
     x(i,gp) = -x(i,gp)
   END DO
 END DO
 CALL proces (x)
 CALL perpec (x,stereo)
 xorig = xy(porig,1)
 IF (stereo /= 0) xorig = xy(porig,2)
 DO  gp = 1,ngpset
   x(2,gp) = scale*x(2,gp) - xorig
   x(3,gp) = scale*x(3,gp) - xy(porig,3)
 END DO
 
 IF (.NOT.disp .OR. maxdef == 0 .OR. defmax == 0) GO TO 120
 
!     PROCESS THE DEFORMATIONS.
!     EXCHANGE AXES, REDUCE THE MAXIMUM DEFORMATION TO -MAXDEF-.
 
 DO  i = 1,3
   axis(i) = IABS(daxis(i))
   SIGN(i) = 1.
   IF (daxis(i) < 0) SIGN(i) = -1.
 END DO
 i = axis(1)
 j = axis(2)
 k = axis(3)
 d = maxdef/defmax
 DO  gp = 1,ngpset
   IF (psymm(4) < 0) u(1,gp) = -u(1,gp)
   IF (psymm(5) < 0) u(2,gp) = -u(2,gp)
   IF (psymm(6) < 0) u(3,gp) = -u(3,gp)
   a(1) = u(i,gp)
   a(2) = u(j,gp)
   a(3) = u(k,gp)
   u(1,gp) = a(1)*SIGN(1)*d
   u(2,gp) = a(2)*SIGN(2)*d
   u(3,gp) = a(3)*SIGN(3)*d
 END DO
 CALL intvec (pvectr)
 
!     IF PVECTR .LT. 0 NO SHAPE WILL BE DRAWN
!     ATTEMPT TO REMOVE DUPLICATE LINES
 
 120 iopt = -1
 sucor = 2*ngpset + 1
 IF (.NOT.disp) sucor = 1
 
!     FIRST DETERMINE OPTIONS - UNIQUE LINES FOR PSHAPE=3 MAY ONLY BE
!     FOR THE UNDERLAY.  ISHAPE = 0 MEANS DRAW THE SHAPE..
 
 ishape = -1
 later  = 0
 IF (pvectr < 0 .OR. (pedge /= 0 .AND. pedge /= 3)) GO TO 130
 ishape = 0
 IF (opcor < ngpset+ngp+1) GO TO 130
 iopt = 0
 defm = 0
 IF (pshape >= 2) defm = 1
 CALL linel (s(sucor,1),iptl,opcor,iopt,x,ppen,defm,gplst)
 IF (pedge == 3) GO TO 500
 IF (iptl  <= 0) iopt = -1
 CALL bckrec (elset)
 130 IF (pshape == 2 .AND. disp) GO TO 260
 
!     DRAW UNDEFORMED SHAPE (USE PEN1 + SYMBOL 2 IF THE DEFORMED SHAPE
!     OR DEFORMATION VECTORS ARE ALSO TO BE DRAWN).
 
 pen = ppen
 IF (disp .AND. pshape > 2) pen = 1
 IF (ishape == 0)  &
     CALL shape (*500,gplst,x,0,pen,0,iopt,iptl,s(sucor,1),opcor)
 IF (pedge < 2) GO TO 140
 CALL hdsurf (gplst,x,0,pen,0,nmax,maxsf,s(sucor,1),buf1,pedge, opcor)
 IF (pedge /= 2 .AND. pedge < 200) GO TO 140
 CALL hdplot (gplst,nmax,maxsf,opcor,buf1)
 GO TO 220
 140 IF (pcon == 0) GO TO 200
 IF (.NOT.disp .OR. pshape < 3) GO TO 210
 later = pcon
 pcon  = 0
 200 IF (pedge == 0 .OR. pedge >= 2) GO TO 220
 210 iopt = -1
 CALL contor (gplst,x,0,u,s(sucor,1),s(sucor,1),pen,0,buf1,opcor)
 IF (pedge == 1) CALL border (gplst,x,0,s(sucor,1),0,buf1,opcor)
 IF (pedge == 1 .OR. color >= 0) GO TO 220
 CALL gopen (elset,gplst(buf1),0)
 CALL shape (*500,gplst,x,0,1,0,iopt,iptl,s(sucor,1),opcor)
 220 pcon = MAX0(pcon,later)
 IF (ppen > 31) CALL shape (*500,gplst,x,0,1,0,iopt,iptl,s(sucor,1),opcor)
 IF (pshape == 1) pcon = 0
 IF (psymbl(1) == 0) GO TO 250
 IF (disp) GO TO 230
 sym(1) = psymbl(1)
 sym(2) = psymbl(2)
 GO TO 240
 230 sym(1) = 2
 sym(2) = 0
 240 CALL gptsym (gplst,x,0,sym,0)
 250 IF (plabel < 0) GO TO 260
 i = plabel/3
 IF (i /= 1) CALL gptlbl (gplst,x,0,0,buf1)
 IF (i < 1) GO TO 260
 CALL elelbl (gplst,x,0,0,buf1)
 CALL bckrec (elset)
 260 IF (.NOT.disp .OR. maxdef == 0.0 .OR. defmax == 0.0) GO TO 500
 IF (pedge == 3) GO TO 500
 IF (pshape < 2 .AND. later == 0) GO TO 350
 
!     ROTATE THE DEFORMATIONS
 
 DO  gp = 1,ngpset
   DO  j  = 1,3
     sum = cstm(j,1)*u(1,gp) + cstm(j,2)*u(2,gp) + cstm(j,3)*u(3,gp)
     IF (j /= 1) GO TO 270
     IF (prject /= 1) dr = d0/(r0-scalex*(x(1,gp)+sum))
     CYCLE
     270 IF (prject /= 1) sum = scalex*dr*sum
     s(j-1,gp) = x(j,gp) + scale*sum
   END DO
 END DO
 
!     DRAW THE DEFORMED SHAPE
 
 IF (pvectr < 0) GO TO 300
 pen = ppen
 IF (pshape == 2 .AND. pvectr /= 0) pen = 1
 IF (pedge == 0) CALL shape (*500,gplst,x,s,pen,1,iopt,iptl,s(sucor,1),opcor)
 IF (pedge < 2) GO TO 300
 CALL hdsurf (gplst,x,s,pen,1,nmax,maxsf,s(sucor,1),buf1,pedge, opcor)
 IF (pedge == 2 .OR. pedge > 200) CALL hdplot (gplst,nmax,maxsf,opcor,buf1)
 300 IF (pcon == 0 .OR. pedge == 2 .OR. pedge > 200) GO TO 310
 IF (icntvl <= 9 .AND. pshape == 1) GO TO 310
 IF (icntvl > 13 .AND. pshape == 1) GO TO 310
 CALL contor (gplst,x,s,u,s(sucor,1),s(sucor,1),pen,0,buf1,opcor)
 IF (pedge == 1 .OR. color >= 0) GO TO 310
 CALL gopen (elset,gplst(buf1),0)
 CALL shape (*500,gplst,x,0,1,0,iopt,iptl,s(sucor,1),opcor)
 310 IF (pedge ==  1) CALL border (gplst,x,s,s(sucor,1),1,buf1,opcor)
 IF (ppen > 31) CALL shape  (*500,gpl,x,0,1,0,iopt,iptl,s(sucor,1),opcor)
 IF (psymbl(1) == 0) GO TO 340
 IF (pshape == 2 .AND. pvectr /= 0) GO TO 320
 sym(1) = psymbl(1)
 sym(2) = psymbl(2)
 GO TO 330
 320 sym(1) = 2
 sym(2) = 0
 330 CALL gptsym (gplst,x,s,sym,1)
 340 IF (plabel < 0 .OR. pshape /= 2) GO TO 350
 i = plabel/3
 IF (i /= 1) CALL gptlbl (gplst,x,s,1,buf1)
 IF (i < 1) GO TO 350
 CALL elelbl (gplst,x,s,1,buf1)
 350 IF (pvectr == 0) GO TO 500
 pvectr = IABS(pvectr)
 
!     PROCESS THE DEFORMATION VECTORS
 
 IF (pvectr <= 7) GO TO 410
 nv = 1
 vec(1) = 0
 vec(2) = 0
 vec(3) = 0
 DO  v = 1,3
   IF (andf(pvectr,2**(v-1)) == 0) CYCLE
   IF (axis(1) == v) vec(1) = 1
   IF (axis(2) == v) vec(2) = 1
   IF (axis(3) == v) vec(3) = 1
 END DO
 GO TO 420
 410 nv = 3
 420 DO  v = 1,nv
   IF (pvectr > 7) GO TO 440
   IF (andf(pvectr,2**(v-1)) == 0) CYCLE
   DO  i = 1,3
     vec(i) = 0
     IF (axis(i) == v) vec(i) = 1
   END DO
   
!     ROTATE THE DEFORMATIONS (VEC = VECTOR DIRECTION TO BE DRAWN)
   
   440 DO  gp = 1,ngpset
     DO  j  = 1,3
       sum = 0.d0
       DO  i = 1,3
         IF (vec(i) /= 0) sum = sum + cstm(j,i)*u(i,gp)
       END DO
       IF (j /= 1) GO TO 460
       IF (prject /= 1) dr = d0/(r0-scalex*(x(1,gp)+sum))
       CYCLE
       460 IF (prject /= 1) sum = scalex*dr*sum
       s(j-1,gp) = x(j,gp) + scale*sum
     END DO
   END DO
   
!     DRAW THE DEFORMATION VECTOR
   
   CALL dvectr (gplst,x,s,ppen)
   IF (psymbl(1) == 0 .OR. pshape == 3) CYCLE
   j = 0
   IF (pshape == 1) j = 1
   CALL gptsym (gplst,x,s,psymbl,j)
 END DO
 
!     END OF PLOT
 
!     IF NOT CONTOUR PLOT, CALL PCOORD TO DRAW A SMALL X-Y-Z COORDINATE
!     TRIAD AT THE LOWER RIGHT CORNER OF PLOT
 
 500 IF (pedge /= 1) CALL pcoord (pen)
 RETURN
END SUBROUTINE draw
