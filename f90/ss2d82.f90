SUBROUTINE ss2d82 (ieqex,neqex,tgrid)
     
!     PHASE 2 OF STRESS DATA RECOVERY FOR 2-D, 8 GRID POINT
!     ISOPARAMETRIC STRUCTURAL ELEMENT
 
!     PH1OUT CONTAINS THE FOLLOWING
!     ELEMENT ID
!     8 SILS
!     TREF
!     ST ARRAY
!     TRANSFORMATION MATRIX FROM GLOBAL TO ELEMENT COORDINATES
!     COORD SYSTEM ID FOR STRESS OUTPUT
!     G MATRIX
!     DNX,DNY AT EACH GRID POINT -EVALUATED 8 TIMES
 
 
 
 INTEGER, INTENT(IN)                      :: ieqex
 INTEGER, INTENT(IN)                      :: neqex
 REAL, INTENT(IN)                         :: tgrid(8)
 DIMENSION  st(3),ta(48),g(9),b(9),db(72),disp(24),  &
     sig(3),bb(72),dnx(8),dny(8),tb(6),temp(9),  &
     istres(3),nsil(1),nph1(1),dn(8),xi(8),eta(8),  &
     pt(3),ex2d82(32),ex2d83(72),sigs(27),sigt(24), iz(1),stress(43)
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform
 COMMON /sdr2x7/ ph1out(100),str(250),forvec(250)
 COMMON /sdr2x8/ disp,dnx,dny,dnz,b,tb,bb,db,sig,ibase,nstrt,npt, is,idtemp
 EQUIVALENCE     (ph1out(1),nph1(1)),(ph1out(1),id),  &
     (nsil(1),ph1out(2)),(tref,ph1out(10)),  &
     (st(1),ph1out(11)),(ta(1),ph1out(14)), (g(1),ph1out(63)),(ph1out(62),id1),  &
     (istres(1),stress(1)),(ldtemp,eltemp), (z(1),iz(1))
 DATA    ex2d82/  &
     1.86603,-.50000,-.50000, .13397,-.50000, .13397,1.86603,-.50000,  &
     .13397,-.50000,-.50000,1.86603,-.50000,1.86603, .13397,-.50000,  &
     .68301,-.18301, .68301,-.18301,-.18301,-.18301, .68301, .68301,  &
     -.18301, .68301,-.18301, .68301, .68301, .68301,-.18301,-.18301/
 DATA    ex2d83/  &
     2.18694,-.98589, .27778,-.98589, .44444,-.12522, .27778,-.12522,  &
     .03528, .27778,-.12522, .03528,-.98589, .44444,-.12522,2.18694,  &
     -.98589, .27778, .03528,-.12522, .27778,-.12522, .44444,-.98589,  &
     .27778,-.98589,2.18694, .27778,-.98589,2.18694,-.12522, .44444,  &
     -.98589, .03528,-.12522, .27778,-.00000,0.00000,-.00000,1.47883,  &
     -.66667, .18784, .00000,-.00000,-.00000,-.00000, .18784, .00000,  &
     -.00000,-.66667,-.00000, .00000,1.47883, .00000,-.00000,0.00000,  &
     .00000, .18784,-.66667,1.47883, .00000,-.00000, .00000,-.00000,  &
     1.47883, .00000,-.00000,-.66667,-.00000,-.00000, .18784,0.00000/
 DATA    xi    / -1., 1., 1.,-1., 0., 1., 0.,-1./
 DATA    eta   / -1.,-1., 1., 1.,-1., 0., 1., 0./
 
!     SET UP DISPLACEMENTS FOR THIS ELEMENT
 
 is  = 0
 DO  i = 1,8
   nstrt = ivec + nsil(i) - 1
   DO  j = 1,3
     is  = is + 1
     npt = nstrt + j - 1
     disp(is) = z(npt)
   END DO
 END DO
 
!     INITIALIZE SOME MATRICES
 
 DO  i = 1,72
   bb(i) = 0.
 END DO
 
!     SET UP INDICATOR FOR GRID POINT TEMPERATURES
 
 idtemp = 0
 DO  i = 1,8
   IF (tgrid(i) /= 0.) GO TO 50
 END DO
 GO TO 60
 50 idtemp = 1
 
!     START LOOPING FOR STRESSES
 
 60 idn = 4
 IF (id1 == 3) idn = 9
 iii = 0
 pt(1) = -0.57735027
 pt(2) = -pt(1)
 IF (id1 == 2) GO TO 133
 pt(1) = -0.77459667
 pt(2) = 0.
 pt(3) = -pt(1)
 133 DO  jii = 1,id1
   DO  jjj = 1,id1
     iii = iii + 1
     
!     COMPUTE BASE POINTER FOR PICKING UP DERIVATIVES
     
     ibase = 71 + 16*(iii-1)
     
     DO  n = 1,8
       nx = n + ibase
       ny = n + ibase + 8
       dnx(n) = ph1out(nx)
       dny(n) = ph1out(ny)
     END DO
     
     DO  n = 1,8
       
!     SET UP THE B MATRIX
       
       DO  i = 1,9
         temp(i) = 0.
         b(i) = 0.
       END DO
       b(1) = dnx(n)
       b(4) = dny(n)
       b(5) = dny(n)
       b(6) = dnx(n)
       
!     TRANSFORM TO ELEMENT COORDINATES
       
       kk = 6*n - 6
       DO  i = 1,6
         k  = kk + i
         tb(i) = ta(k)
       END DO
       CALL gmmats (b,3,2,0,tb,2,3,0,temp(1))
       n3 = 3*n
       bb(n3- 2) = temp(1)
       bb(n3- 1) = temp(2)
       bb(n3   ) = temp(3)
       bb(n3+22) = temp(4)
       bb(n3+23) = temp(5)
       bb(n3+24) = temp(6)
       bb(n3+46) = temp(7)
       bb(n3+47) = temp(8)
       bb(n3+48) = temp(9)
     END DO
     
!     BRING IN G MATRIX
     
     CALL gmmats (g,3,3,0,bb,3,24,0,db)
     
!     COMPUTE STRESSES
     
     CALL gmmats (db,3,24,0,disp,24,1,0,sig)
     
!     STORE GAUSS POINT STRESSES INTO SIGT
     
     i3 = 3*(iii-1)
     DO  i = 1,3
       isub = i3 + i
       sigs(isub) = sig(i)
     END DO
     
!     COMPUTE GAUSS POINT  TEMPERATURES
     
     IF (ldtemp == -1) CYCLE
     IF (idtemp ==  1) GO TO 229
     rgtemp = eltemp - tref
     GO TO 250
     
!     ALL TEMPERATURES ARE DEFAULT VALUE
     
     229 DO  n = 1,4
       dn(n) = .25*(1.+pt(jii)*xi(n))*(1.+pt(jjj)*eta(n))  &
           *(pt(jii)*xi(n)+pt(jjj)*eta(n)-1.)
     END DO
     DO  n = 5,7,2
       dn(n) = .5*(1.-pt(jii)*pt(jii))*(1.+pt(jjj)*eta(n))
     END DO
     DO  n = 6,8,2
       dn(n) = .5*(1.+pt(jii)*xi(n))*(1.-pt(jjj)*pt(jjj))
     END DO
     gstemp = 0.
     DO  n = 1,8
       gstemp = gstemp + dn(n)*tgrid(n)
     END DO
     rgtemp = gstemp - tref
     250 CONTINUE
     DO  i = 1,3
       isub = i3 + i
       sigs(isub) = sigs(isub) - st(i)*rgtemp
     END DO
     
   END DO
 END DO
 
!     MULTIPLY BY TRANSFORMATION FROM GAUSS POINTS TO GRID POINTS
 
 IF (id1 == 2) CALL gmmats (ex2d82,8,4,0,sigs,4,3,0,sigt)
 IF (id1 == 3) CALL gmmats (ex2d83,8,9,0,sigs,9,3,0,sigt)
 
!     FINISH UP
 
 DO  iii = 1,8
   
!     MOVE A ROW OF SIGT INTO SIG
   
   i3 = 3*(iii-1)
   DO  i = 1,3
     isub = i3 + i
     sig(i) = sigt(isub)
   END DO
   
!     STORE STRESSES
   
   jsub  = 5*(iii-1) + 4
   isub1 = ieqex + 1
   isub2 = ieqex + neqex - 1
   DO  jjj = isub1,isub2,2
     ns = iz(jjj)/10
     IF (ns /= nsil(iii)) CYCLE
     istres(jsub) = iz(jjj-1)
     GO TO 162
   END DO
   CALL mesage (-30,164,iz(jjj))
   162 CONTINUE
   istres(jsub+1) = 0
   DO  i = 1,3
     jjsub = jsub + 1 + i
     stress(jjsub) = sig(i)
   END DO
   
!     LOOP FOR OTHER GRID POINTS
   
 END DO
 
!     FINISH UP
 
!     ELEMENT ID
 
 istres(1) = id
 
!     NUMBER OF GRID POINTS PER ELEMENT
 
 istres(2) = 8
 
!     NUMBER OF STRESSES OUTPUT PER ELEMENT
 
 istres(3) = 3
 
 DO  i = 1,43
   str(i) = stress(i)
 END DO
 
 RETURN
END SUBROUTINE ss2d82
