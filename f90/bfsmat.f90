SUBROUTINE bfsmat (nd,NE,nb,np,ntp,length,ntotal,scr1,jf,jl,nas,  &
        fmach,yb,zb,ys,zs,x,delx,ee,xic,sg,cg,ar,ria,  &
        nbea1,nbea2,nasb,nsaray,ncaray,bfs,avr,cbar,  &
        a0,xis1,xis2,kr,nsbea,nt0)
     
!     NOTE:
!     A JUMP (VIA AN ASSIGN STATEMENT) TO 200 AND A JUMP TO 1100 (ALSO
!     VIA AN ASSIGN STATEMENT), INTO THE MIDDLES OF SOME DO LOOPS, ARE
!     ACCEPTABLE ANSI 77 FORTRAN. HOWEVER, IBM COMPILER MAY COMPLAIN.
!     THIS PROBLEM IS NOW ELIMINATED (BY G.C. 9/89)
 
 
!        ND         SYMMETRY FLAG
!        NE         GROUND EFFECTS FLAG
!        NB         NUMBER OF BODIES
!        NP         NUMBER OF PANELS
!        NTP        NUMBER OF LIFTING SURFACE BOXES
!        NTOTAL     NTP + TOTAL NO. OF Y AND Z ORIENTED BODY ELEMENTS
!        LENGTH     NTOTAL + THE TOTAL NUMBER OF Z- AND Y-ORIENTED
!                   SLENDER BODY ELEMENTS
!        SCR1       FILE FOR OUTPUT
!        JF         ROW FOR FIRST ZY BODY
!        JL         ROW FOR LAST  ZY BODY
!        NAS        ARRAY CONTAINING THE NUMBER OF ASSOCIATED BODIES
!                   FOR EACH PANEL
!        FMACH      MACH NUMBER
!        YB         ARRAY OF -Y- COORDINATES OF THE BODIES
!        ZB         ARRAY OF -Z- COORDINATES OF THE BODIES
!        YS         ARRAY OF -Y- COORDINATES OF STRIPS AND BODIES
!        ZS         ARRAY OF -Z- COORDINATES OF STRIPS AND BODIES
!        X          ARRAY OF 3/4 CHORD LOCATIONS OF BOXES AND
!                            1/2 CHORD FOR BODY ELEMENTS
!        DELX       ARRAY OF LENGTHS OF BOXES AND BODY ELEMENTS
!        EE         ARRAY OF THE SEMI-WITH OF STRIPS
!        XIC        ARRAY OF 1/4 CHORD COORDINATES OF BOXES
!        SG         ARRAY OF SINE   OF STRIP DIHEDRAL ANGLE
!        CG         ARRAY OF COSINE OF STRIP DIHEDRAL ANGLE
!        AR         ARRAY OF RATIO OF MAJOR AXES OF BODIES
!        RIA        ARRAY OF RADII OF BODY ELEMENTS
!        NBEA1      ARRAY OF NUMBER OF BODY ELEMENTS PER BODY
!        NBEA2      ARRAY OF THE BODY ORIENTATION FLAGS PER BODY
!        NASB       ARRAY OF THE BODIES ASSOCIATED WITH PANELS
!        NSARAY     ARRAY OF THE NUMBER OF STRIPS PER PANEL
!        NCARAY     ARRAY OF THE NUMBER OF CHORDWISE DIV. PER PANEL
!        BFS        WORK ARRAY FOR TEMPORARY STORAGE OF THE BFS COLS.
!        AVR        ARRAY OF RADII OF BODIES
!        CBAR       REFERENCE CHORD
!        A0         ARRAY OF SLENDER BODY ELEMENT RADII
!        XIS1       ARRAY OF SLENDER BODY ELEMENT LEADING  EDGE COORD.S
!        XIS2       ARRAY OF SLENDER BODY ELEMENT TRAILING EDGE COORD.S
!        KR         REDUCED FREQUENCY
!        NSBEA      ARRAY OF THE NUMBER OF ELEMENTS PER SLENDER BODY
 
 
 INTEGER, INTENT(IN OUT)                  :: nd
 INTEGER, INTENT(IN OUT)                  :: NE
 INTEGER, INTENT(IN)                      :: nb
 INTEGER, INTENT(IN)                      :: np
 INTEGER, INTENT(IN)                      :: ntp
 INTEGER, INTENT(IN)                      :: length
 INTEGER, INTENT(IN OUT)                  :: ntotal
 INTEGER, INTENT(IN OUT)                  :: scr1
 INTEGER, INTENT(OUT)                     :: jf
 INTEGER, INTENT(OUT)                     :: jl
 INTEGER, INTENT(IN)                      :: nas(1)
 REAL, INTENT(IN)                         :: fmach
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: zb(1)
 REAL, INTENT(IN)                         :: ys(1)
 REAL, INTENT(IN)                         :: zs(1)
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(IN)                         :: delx(1)
 REAL, INTENT(IN)                         :: ee(1)
 REAL, INTENT(IN)                         :: xic(1)
 REAL, INTENT(IN)                         :: sg(1)
 REAL, INTENT(IN)                         :: cg(1)
 REAL, INTENT(IN OUT)                     :: ar(1)
 REAL, INTENT(IN OUT)                     :: ria(1)
 INTEGER, INTENT(IN)                      :: nbea1(1)
 INTEGER, INTENT(IN)                      :: nbea2(1)
 INTEGER, INTENT(IN OUT)                  :: nasb(1)
 INTEGER, INTENT(IN)                      :: nsaray(1)
 INTEGER, INTENT(IN)                      :: ncaray(1)
 COMPLEX, INTENT(OUT)                     :: bfs(length,2)
 REAL, INTENT(IN OUT)                     :: avr(1)
 REAL, INTENT(IN)                         :: cbar
 REAL, INTENT(IN)                         :: a0(1)
 REAL, INTENT(IN)                         :: xis1(1)
 REAL, INTENT(IN)                         :: xis2(1)
 REAL, INTENT(IN)                         :: kr
 INTEGER, INTENT(IN)                      :: nsbea(1)
 INTEGER, INTENT(IN OUT)                  :: nt0
 LOGICAL :: last
 
 
 COMPLEX :: fwz, fwy, eikj1 , eikj2
 
 
 
 beta2 = 1.0 - fmach**2
 icol  = 0
 ksp   = 1
 GO TO 3000
 
 50 CONTINUE
 
!     -Y- ORIENTED BODIES AS SENDING ELEMENTS
 
 sgs   =-1.0
 cgs   = 0.0
 nasd  = 0
 izyflg= 3
 ASSIGN 4010 TO ibody
 GO TO 1000
 100 CONTINUE
 
!     - LIFTING SURF. BOXES AS SENDING ELEMENTS
 
 IF (ntp <= 0) GO TO 800
 j     = 1
 jp1   = j
 ibox  = 0
 istrip= 0
 isn   = 0
 ksp   = 1
 
!     LOOP FOR -PANEL-
 
 DO  isp = 1,np
   ns  = nsaray(isp)
   nc  = ncaray(isp)
   ns  = (ns-isn) / nc
   isn = nsaray(isp)
   nasd= nas(isp)
   
!     LOOP FOR -STRIP-
   
   DO  is = 1,ns
     istrip= istrip + 1
     dys   = ys(istrip)
     dzs   = zs(istrip)
     sgs   = sg(istrip)
     cgs   = cg(istrip)
     width = 2.0 * ee(istrip)
     
!     LOOP FOR -BOX-
     
     DO  ib = 1,nc
       ibox = ibox + 1
       dxs  = xic(ibox)
       
       icol = icol + 1
       
       CALL fwmw (nd,NE,sgs,cgs,irb,dria,ar,dxle,dxte,yb,zb,dxs,dys,  &
           dzs,nasd,nasb(ksp),kr,beta2,cbar,avr,fwz,fwy)
       bfs(icol,1) = fwz * (dxte - dxle)
       bfs(icol,2) = fwy * (dxte - dxle)
       bfs(icol,1) = bfs(icol,1) * scale
       bfs(icol,2) = bfs(icol,2) * scale
       
       area = width * delx(ibox)
       bfs(icol,1) = bfs(icol,1) * area
       bfs(icol,2) = bfs(icol,2) * area
     END DO
   END DO
   ksp = ksp + nasd
 END DO
 
!     -Z-  ORIENTED BODIES AS SENDING ELEMENTS
 
 800 CONTINUE
 sgs   = 0.0
 cgs   = 1.0
 nasd  = 0
 izyflg= 1
 ASSIGN   50 TO ibody
 GO TO  1000
 
 
!     *** LOOP FOR EACH INTERFERENCE BODY SENDING ELEMENT
 
 1000 CONTINUE
 INDEX = ntp
 
!     --ISB-- IS THE SENDING BODY
 
 DO  isb  = 1,nb
   IF (nbea2(isb) ==     2 ) GO TO 1070
   IF (nbea2(isb) /= izyflg) GO TO 1850
   1070 dys   = yb(isb)
   nsbe  = nbea1(isb)
   jp1   = 1
   last  = .false.
   dzs   = zb(isb)
   earg2 = 1.0
   
!     --ISBE-- IS THE ELEMENT OF THE SEND BODY
   
   DO  isbe = 1,nsbe
     earg1 = earg2
     INDEX = INDEX + 1
     dxs   = x (INDEX) - delx(INDEX) /4.0
     earg2 = kr * delx(INDEX) / cbar
     
!     CALCULATE THIS COLUMN
     
     icol  = icol + 1
     eikj1 = CMPLX(COS(earg1),-SIN(earg1))
     eikj2 = CMPLX(COS(earg2), SIN(earg2))
     
     CALL fwmw (nd,NE,sgs,cgs,irb,dria,ar,dxle,dxte,yb,zb,dxs,dys,  &
         dzs,nasd,nasb(ksp),kr,beta2,cbar,avr,fwz,fwy)
     bfs(icol,1) = fwz * (dxte - dxle)
     bfs(icol,2) = fwy * (dxte - dxle)
     bfs(icol,1) = bfs(icol,1) * scale
     bfs(icol,2) = bfs(icol,2) * scale
     
     
!     IS THIS THE FIRST COLUMN, YES  BRANCH
     
     IF (isbe == 1) CYCLE
     bfs(icol-1,1) = bfs(icol-1,1)*eikj1 - bfs(icol,1)*eikj2
     bfs(icol-1,2) = bfs(icol-1,2)*eikj1 - bfs(icol,2)*eikj2
   END DO
   CYCLE
   1850 INDEX = INDEX + nbea1(isb)
 END DO
 
!     RETURN TO CALLING POINT - EITHER Y OR Z SENDING BODY ELEM
 
!     *** GO  EITHER TO THE  Y-ORIENTED INTERFERENCE BODY ELEMENT LOOP
!         OR  TO THE LOOP FOR SLENDER BODY SENDING ELEMENTS
 
 GO TO ibody, (50,4010)
 
 
!     CALCULATE EACH ROW OF THE SENDING COLUMN
 
 3000 CONTINUE
 iy   = 0
 jf   = 0
 nw   = length*2
 irow = 0
 
!     --IRB-- IS THE RECEIVING BODY
 
 irb  = 0
 3050 irb  = irb + 1
 IF (irb > nb) GO TO 3900
 nrbe = nsbea(irb)
 itsb = nbea2(irb)
 
 xyb  = yb(irb)
 xzb  = zb(irb)
 scale= 1.0
 IF (nd /= 0 .AND. xyb == 0.0) scale = .5
 IF (NE /= 0 .AND. xzb == 0.0) scale = scale*.5
 
!     --IRBE-- IS THE ELEM. OF THE REC. BODY
 
 irbe = 0
 3060 irbe = irbe + 1
 IF (irbe > nrbe) GO TO 3800
 iy   = iy   + 1
 irow = irow + 1
 dria = a0  (iy)
 dxle = xis1(iy)
 dxte = xis2(iy)
 xx1  = dxle
 xx2  = dxte
 xaa  = dria
 icol = 0
 GO  TO  100
 
 3100 CONTINUE
 SELECT CASE ( itsb )
   CASE (    1)
     GO TO 3110
   CASE (    2)
     GO TO 3120
   CASE (    3)
     GO TO 3130
 END SELECT
 3110 CALL WRITE (scr1,bfs(1,1),nw,0)
 GO TO 3140
 3120 CALL WRITE (scr1,bfs(1,2),nw,0)
 CALL WRITE (scr1,bfs(1,1),nw,0)
 IF (jf == 0) jf = irow
 irow = irow + 1
 GO TO 3140
 3130 CALL WRITE (scr1,bfs(1,2),nw,0)
 irow = irow - 1
 3140 CONTINUE
 GO TO 3060
 3800 CONTINUE
 GO TO 3050
 3900 CONTINUE
 jl = irow
 RETURN
 
 
 4010 CONTINUE
 
 
!     *** LOOP FOR EACH SLENDER BODY SENDING ELEMENT
 
 izyflg= 1
 sgs   = 0.0
 cgs   = 1.0
 4050 CONTINUE
 lsbe  =  0
 DO  lsb = 1,nb
   
!     --LSB-- IS THE INDEX OF THE SLENDER SENDING BODY
   
   IF (nsbea(lsb) ==      0) CYCLE
   IF (nbea2(lsb) ==      2) GO TO 4070
   IF (nbea2(lsb) /= izyflg) GO TO 4097
   4070 CONTINUE
   xeta  =  yb(lsb)
   xzeta =  zb(lsb)
   scale2=  scale
   msbe  =  nsbea(lsb)
   DO  lsbs = 1,msbe
     lsbe  =  lsbe + 1
     icol  =  icol + 1
     xxij  = .50 * xis1(lsbe) + .50 * xis2(lsbe)
     CALL fwmw (nd,NE,sgs,cgs,irb,dria,ar,dxle,dxte,yb,zb,xxij,xeta,  &
         xzeta,nasd,nasb,kr,beta2,cbar,avr,fwz,fwy)
     bfs(icol,1) =  fwz * (dxte - dxle)
     bfs(icol,2) =  fwy * (dxte - dxle)
     bfs(icol,1) =  bfs(icol,1) * scale2
     bfs(icol,2) =  bfs(icol,2) * scale2
   END DO
   CYCLE
   4097 lsbe = lsbe + nsbea(lsb)
 END DO
 IF (izyflg == 3) GO TO  5010
 izyflg = 3
 sgs =-1.0
 cgs = 0.0
 GO TO 4050
 5010 CONTINUE
 
 GO TO 3100
 
END SUBROUTINE bfsmat
