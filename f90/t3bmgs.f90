SUBROUTINE t3bmgs (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,  &
        detjac,shp,bterms,bmatrx)
     
!     [B] MATRIX GENERATOR ROUTINE FOR TRIA3 ELEMENTS
 
!     SINGLE PRECISION ROUTINE TO GENERATE A 9XNDOF [B] MATRIX AT A
!     GIVEN INTEGRATION POINT, USING THE DERIVATIVES OF THE 2-D SHAPE
!     FUNCTIONS.
!     OPTIONALLY, AN 8XNDOF [B] MATRIX IS CONSTRUCTED AND/OR SHEAR TERMS
!     MAY BE DROPPED ALTOGETHER, YIELDING 6XNDOF MATRIX.
!     FOR STRESS RECOVERY, THE EVALUATION POINTS ARE AT THE ELEMENT
!     INTERIOR POINTS RATHER THAN ON THE EDGES.
!     THE CONTENTS OF /TERMS/ ARE USED TO CONSTRUCT THE [B] MATRIX
!     ACCORDING TO THE BEHAVIORAL REQUIREMENTS OF THE ELEMENT.
 
 
!     INPUT :
!           IPT    - POINTER TO THE CURVILNEAR COORDINATES
!           SHEART - LOGICAL INDICATING THE REQUIREMENT FOR OUT-OF-PLANE
!                    SHEAR TERMS
!           IORDER - ARRAY OF INTERNAL SEQUENCE OF NODES
!           EGPDT  - GRID POINT DATA IN THE ELEMENT COORD. SYSTEM
!           DGPTH  - NODAL THICKNESSES
!           AIC    - TRANSFORMATION TO RELIEVE GEOMETRY BIAS
!     OUTPUT:
!           IERR   - ERROR FLAG
!           TH     - THICKNESS AT THE INTEG. PT.
!           DETJAC - DETERMINANT OF JACOBIAN AT THE INTEG. PT.
!           SHP    - ARRAY OF REORDERED SHAPE FUNCTIONS
!           BTERMS - DERIVATIVES WRT THE PHYSICAL COORDINATES
!           BMATRX - STRAIN-DISPLACEMENT RELATIONSHIP
 
 
 
 INTEGER, INTENT(OUT)                     :: ierr
 LOGICAL, INTENT(IN OUT)                  :: sheart
 INTEGER, INTENT(IN)                      :: ipt
 INTEGER, INTENT(IN)                      :: iorder(3)
 REAL, INTENT(IN)                         :: egpdt(4,1)
 REAL, INTENT(IN)                         :: dgpth(1)
 REAL, INTENT(IN)                         :: aic(1)
 REAL, INTENT(OUT)                        :: th
 REAL, INTENT(OUT)                        :: detjac
 REAL, INTENT(OUT)                        :: shp(3)
 REAL, INTENT(OUT)                        :: bterms(1)
 REAL, INTENT(OUT)                        :: bmatrx(1)
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth
 
 REAL :: dshpx(3),dshpe(3),tshp(3),  &
     tdshpx(3),tdshpe(3),vi(2),vj(2),jacob(4),eps,  &
     ptint(2,7),trc(2,3),xsi,xsii,eta,etai,psi,psii, dnx,dny,shpf
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 DATA    eps   /  1.0E-13 /
 DATA    ptint /  0.5,    0.0,    0.5, 0.5,    0.0,    0.5,  &
     0.333333333333333D0, 0.333333333333333D0,  &
     0.166666666666667D0, 0.166666666666667D0,  &
     0.166666666666667D0, 0.666666666666667D0,  &
     0.666666666666667D0, 0.166666666666667D0/
 DATA    trc   /  0.0,    0.0,    1.0, 0.0,    0.0,    1.0/
 
!     INITIALIZE
 
 ierr = 0
 nnode= 3
 nd1  = nnode*6
 nd2  = nd1*2
 nd3  = nd1*3
 nd4  = nd1*4
 nd5  = nd1*5
 nd6  = nd1*6
 nd7  = nd1*7
 nd8  = nd1*8
 nd9  = nd1*9
 
 DO  i = 1,6
   bterms(i) = 0.0
 END DO
 DO  i = 1,nd9
   bmatrx(i) = 0.0
 END DO
 
!     CALCULATE THE SHAPE FUNCTIONS AND THEIR DERIVATIVES, THEN SORT
!     THEM.
 
 xsi = ptint(1,ipt)
 eta = ptint(2,ipt)
 psi = 1.0 - xsi - eta
 
 DO  i = 1,3
   xsii = trc(1,i)
   etai = trc(2,i)
   psii = 1.0 - xsii - etai
   
   shp(i)   = xsi*xsii + eta*etai + psi*psii
   dshpx(i) = xsii - psii
   dshpe(i) = etai - psii
 END DO
 
 DO  i  = 1,nnode
   tshp(i)  = shp(i)
   tdshpx(i)= dshpx(i)
   tdshpe(i)= dshpe(i)
 END DO
 
 DO  i  = 1,nnode
   kk       = iorder(i)
   shp(i)   = tshp(kk)
   dshpx(i) = tdshpx(kk)
   dshpe(i) = tdshpe(kk)
 END DO
 
!     COMPUTE THE ELEMENT THICKNESS
 
 th = 0.0
 DO  ish = 1,nnode
   th = th + shp(ish)*dgpth(ish)
 END DO
 
!     SET UP THE JACOBIAN
 
 DO  i = 1,2
   vi(i) = 0.0
   vj(i) = 0.0
   ii = i + 1
   DO  j = 1,nnode
     vi(i) = vi(i) + egpdt(ii,j)*dshpx(j)
     vj(i) = vj(i) + egpdt(ii,j)*dshpe(j)
   END DO
 END DO
 
!     INVERT THE JACOBIAN
 
 detjac = vi(1)*vj(2) - vi(2)*vj(1)
 IF (detjac >= eps) GO TO 100
 ierr = 1
 RETURN
 
 100 jacob(1) =  vj(2)/detjac
 jacob(2) = -vi(2)/detjac
 jacob(3) = -vj(1)/detjac
 jacob(4) =  vi(1)/detjac
 
 DO  i = 1,4
   IF (ABS(jacob(i)) < eps) jacob(i) = 0.0
 END DO
 
 ipt1 = ipt*2 - 1
 i71  = ipt1
 i72  = ipt1 + 1
 i81  = ipt1 + 6
 i82  = ipt1 + 7
 i91  = ipt1 + 12
 i92  = ipt1 + 13
 
!     LOOP OVER NODES AND BUILD PARTITIONS OF [B]
 
 ip = 0
 DO  i = 1,nnode
   
!     CALCULATE DERIVATIVES WRT THE PHYSICAL COORDINATES.
   
   dnx  = jacob(1)*dshpx(i) + jacob(2)*dshpe(i)
   dny  = jacob(3)*dshpx(i) + jacob(4)*dshpe(i)
   shpf = shp(i)
   
   bterms(i      ) = dnx
   bterms(i+nnode) = dny
   
   IF (.NOT.membrn) GO TO 120
   
!     ROW 1
   
   bmatrx(ip+1) = dnx
   
!     ROW 2
   
   bmatrx(ip+2+nd1) = dny
   
!     ROW 3
   
   bmatrx(ip+1+nd2) = dny
   bmatrx(ip+2+nd2) = dnx
   
   120 IF (.NOT.bendng) GO TO 150
   
!     ROW 4
   
   bmatrx(ip+5+nd3) = -dnx
   
!     ROW 5
   
   bmatrx(ip+4+nd4) =  dny
   
!     ROW 6
   
   bmatrx(ip+5+nd5) = -dny
   bmatrx(ip+4+nd5) =  dnx
   
   IF (.NOT.sheart) GO TO 150
   IF (ipt <  4) GO TO 130
   
!     8-ROW MATRIX
   
!     ROW 7
   
   bmatrx(ip+3+nd6) =  dny
   bmatrx(ip+4+nd6) = -shpf
   
!     ROW 8
   
   bmatrx(ip+3+nd7) =  dnx
   bmatrx(ip+5+nd7) =  shpf
   GO TO 150
   
!     9-ROW MATRIX
   
!     ROW 7
   
   130 bmatrx(ip+3+nd6) =  aic(i71)*dny + aic(i72)*dnx
   bmatrx(ip+4+nd6) = -aic(i71)*shpf
   bmatrx(ip+5+nd6) =  aic(i72)*shpf
   
!     ROW 8
   
   bmatrx(ip+3+nd7) =  aic(i81)*dny + aic(i82)*dnx
   bmatrx(ip+4+nd7) = -aic(i81)*shpf
   bmatrx(ip+5+nd7) =  aic(i82)*shpf
   
!     ROW 9
   
   bmatrx(ip+3+nd8) =  aic(i91)*dny + aic(i92)*dnx
   bmatrx(ip+4+nd8) = -aic(i91)*shpf
   bmatrx(ip+5+nd8) =  aic(i92)*shpf
   
   150 ip = ip + 6
 END DO
 
 RETURN
END SUBROUTINE t3bmgs
