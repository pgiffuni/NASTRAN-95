SUBROUTINE ssold1(itype)
!*****
 
!  E C P T     TETRA          WEDGE          HEXA
!  -----------------------------------------------
!  ECPT( 1) =  EL ID          EL ID          EL ID
!  ECPT( 2) =  MAT-ID         MAT-ID         MAT-ID
!  ECPT( 3) =  GRID-1         GRID-1         GRID-1
!  ECPT( 4) =  GRID-2         GRID-2         GRID-2
!  ECPT( 5) =  GRID-3         GRID-3         GRID-3
!  ECPT( 6) =  GRID-4         GRID-4         GRID-4
!  ECPT( 7) =  CSID-1         GRID-5         GRID-5
!  ECPT( 8) =  X1             GRID-6         GRID-6
!  ECPT( 9) =  Y1             CSID-1         GRID-7
!  ECPT(10) =  Z1             X1             GRID-8
!  ECPT(11) =  CSID-2         Y1             CSID-1
!  ECPT(12) =  X2             Z1             X1
!  ECPT(13) =  Y2             CSID-2         Y1
!  ECPT(14) =  Z2             X2             Z1
!  ECPT(15) =  CSID-3         Y2             CSID-2
!  ECPT(16) =  X3             Z2             X2
!  ECPT(17) =  Y3             CSID-3         Y2
!  ECPT(18) =  Z3             X3             Z2
!  ECPT(19) =  CSID-4         Y3             CSID-3
!  ECPT(20) =  X4             Z3             X3
!  ECPT(21) =  Y4             CSID-4         Y3
!  ECPT(22) =  Z4             X4             Z3
!  ECPT(23) =  EL-TEM         Y4             CSID-4
!  ECPT(24)                   Z4             X4
!  ECPT(25)                   CSID-5         Y4
!  ECPT(26)                   X5             Z4
!  ECPT(27)                   Y5             CSID-5
!  ECPT(28)                   Z5             X5
!  ECPT(29)                   CSID-6         Y5
!  ECPT(30)                   X6             Z5
!  ECPT(31)                   Y6             CSID-6
!  ECPT(32)                   Z6             X6
!  ECPT(33)                   ELTEMP         Y6
!  ECPT(34)                                  Z6
!  ECPT(35)                                  CSID-7
!  ECPT(36)                                  X7
!  ECPT(37)                                  Y7
!  ECPT(38)
!  ECPT(39)                                  CSID-8
!  ECPT(40)                                  X8
!  ECPT(41)                                  Y8
!  ECPT(42)                                  Z8
!  ECPT(43)                                  EL-TEMP
!*****
 
 INTEGER, INTENT(IN OUT)                  :: itype
 REAL :: nu
 INTEGER :: necpt(100),nphi(170), m(14,4)
 
 COMMON / sdr2x5/ ecpt(100),phiout(170)
 COMMON / sdr2x6 /       cmat(18,8)         ,beta(8)  &
     ,temp(18)           ,elvol ,vol                ,fact  &
     ,npts               ,nel ,mfirst             ,nrow  &
     ,itest              ,floc ,j1                 ,jloc  &
     ,kpt                ,ntemp ,GE(36)             ,h(4,4)  &
     ,r(3,3)             ,ti(9)
 COMMON / matin / nmat,matflg,eltemp
 COMMON / matout/ e,g,nu,rho,alfa,tempo
 
 EQUIVALENCE (nphi(1),phiout(1))
 EQUIVALENCE (necpt(1),ecpt(1))
 
 DATA  m( 1,1),m( 1,2),m( 1,3),m( 1,4)/ 1   ,2   ,3   ,4 /
 
 DATA  m( 2,1),m( 2,2),m( 2,3),m( 2,4)/ 1   ,2   ,3   ,6 /
 DATA  m( 3,1),m( 3,2),m( 3,3),m( 3,4)/ 1   ,2   ,6   ,5 /
 DATA  m( 4,1),m( 4,2),m( 4,3),m( 4,4)/ 1   ,4   ,5   ,6 /
 
 DATA  m( 5,1),m( 5,2),m( 5,3),m( 5,4)/ 1   ,2   ,3   ,6 /
 DATA  m( 6,1),m( 6,2),m( 6,3),m( 6,4)/ 1   ,3   ,4   ,8 /
 DATA  m( 7,1),m( 7,2),m( 7,3),m( 7,4)/ 1   ,3   ,8   ,6 /
 DATA  m( 8,1),m( 8,2),m( 8,3),m( 8,4)/ 1   ,5   ,6   ,8 /
 DATA  m( 9,1),m( 9,2),m( 9,3),m( 9,4)/ 3   ,6   ,7   ,8 /
 DATA  m(10,1),m(10,2),m(10,3),m(10,4)/ 2   ,3   ,4   ,7 /
 DATA  m(11,1),m(11,2),m(11,3),m(11,4)/ 1   ,2   ,4   ,5 /
 DATA  m(12,1),m(12,2),m(12,3),m(12,4)/ 2   ,4   ,5   ,7 /
 DATA  m(13,1),m(13,2),m(13,3),m(13,4)/ 2   ,5   ,6   ,7 /
 DATA  m(14,1),m(14,2),m(14,3),m(14,4)/ 4   ,5   ,7   ,8 /
 
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 100
   CASE (    2)
     GO TO 110
   CASE (    3)
     GO TO 120
   CASE (    4)
     GO TO 130
 END SELECT
!*****
!     THE TYPE OF THE ELEMENT DETERMINES THE FOLLOWING PARAMETERS
!*****
 100 npts  =4
 nel   =1
 mfirst=1
 GO TO 140
 110 npts  =6
 nel   =3
 mfirst=2
 GO TO 140
 120 npts  =8
 nel   =5
 mfirst=5
 GO TO 140
 130 npts  =8
 nel   =10
 mfirst=5
 140 CONTINUE
!*****
!     ZERO OUT ARRAYS
!*****
 elvol =0.0
 DO   j =1,npts
   beta(j) =0.0
   DO   i =1,18
     cmat(i,j) = 0.0
   END DO
 END DO
!*****
!     LOOP ON SUBELEMENTS
!*****
 DO  me =1,nel
   nrow = mfirst +me -1
!*****
!     J  CORRESPONDS TO THE X,Y,AND Z LOCATIONS OF EACH CONNECTED POINT
!*****
   DO  j =1,3
     j1 = m(nrow,1)*4 +npts +j -1
!*****
!     I  CORRESPONDS TO POINTS 2,3,AND 4
!*****
     DO  i= 1,3
       jloc = m(nrow,i+1)*4 + npts + j - 1
!*****
!     ECPT(JLOC) IS THE JTH COMPONENT OF POINT I+1
!*****
       r(i,j) = ecpt(jloc) -ecpt(j1)
     END DO
   END DO
!*****
!     INVERT THE GEOMETRY MATRIX EXPLICITLY USING VECTOR OPERATORS
!*****
   CALL saxb(   r(1,3),r(1,2),temp)
   h(2,1) = temp(1) + temp(2) + temp(3)
   h(2,2) = r(2,2)*r(3,3) -r(3,2)*r(2,3)
   h(2,3) = r(3,2)*r(1,3) -r(1,2)*r(3,3)
   h(2,4) = r(1,2)*r(2,3) -r(2,2)*r(1,3)
   CALL saxb(   r(1,1),r(1,3),temp)
   h(3,1) = temp(1)+ temp(2) +temp(3)
   h(3,2) = r(2,3)*r(3,1) -r(3,3)*r(2,1)
   h(3,3) = r(3,3)*r(1,1) -r(1,3)*r(3,1)
   h(3,4) = r(1,3)*r(2,1) -r(2,3)*r(1,1)
   CALL saxb(   r(1,1),r(1,2),temp)
   h(4,1) =-temp(1) -temp(2)-temp(3)
   h(4,2) = r(2,1)*r(3,2) -r(3,1)*r(2,2)
   h(4,3) = r(3,1)*r(1,2) -r(1,1)*r(3,2)
   h(4,4) = r(1,1)*r(2,2) -r(2,1)*r(1,2)
   vol =   (r(1,3)*temp(1) + r(2,3)*temp(2) +r(3,3)*temp(3))/6.0
   elvol = elvol +vol
   DO  i = 1,4
     kpt = m(nrow,i)
     beta(kpt) = beta(kpt) + vol
     cmat(1, kpt) = h(2,i)/6.0  + cmat(1 ,kpt)
     cmat(5, kpt) = h(3,i)/6.0  + cmat(5 ,kpt)
     cmat(9, kpt) = h(4,i)/6.0  + cmat(9 ,kpt)
     cmat(11,kpt) = h(4,i)/6.0  + cmat(11,kpt)
     cmat(12,kpt) = h(3,i)/6.0  + cmat(12,kpt)
     cmat(13,kpt) = h(4,i)/6.0  + cmat(13,kpt)
     cmat(15,kpt) = h(2,i)/6.0  + cmat(15,kpt)
     cmat(16,kpt) = h(3,i)/6.0  + cmat(16,kpt)
     cmat(17,kpt) = h(2,i)/6.0  + cmat(17,kpt)
   END DO
 END DO
!*****
!     END OF ELEMENT LOOP
!*****
!*****
!     CMAT CONTAINS THE SUM OF THE STRAIN -DISPLACEMENT MATRICES
!                       TIMES THE VOLUME OF THE CONNECTED TETRAHEDRON
 
!     CALL THE MATERIAL  ROUTINE TO OBTAIN PARAMETERS
!*****
 nmat = necpt(2)
 matflg =1
 eltemp = ecpt(5*npts+3)
 CALL mat (ecpt(1))
 fact = e /((1.0+nu)*(1.0-2.0*nu) )
 DO  i = 1,36
   GE(i) = 0.0
 END DO
 GE(1) = fact*(1.0-nu)
 GE(2) = fact* nu
 GE(3) = GE(2)
 GE(7) = GE(2)
 GE(8) = GE(1)
 GE(9) = GE(2)
 GE(13)= GE(2)
 GE(14)= GE(2)
 GE(15)= GE(1)
 GE(22)= g
 GE(29)= g
 GE(36)= g
!*****
!     EACH CMAT MATRIX IS PREEMULTIPLIED BY THE STRESS-STRAIN GE MATRIX
!         AND DIVIDED BY THE SUM OF THE VOLUMES.
!     IF NECESSARY THE MATRIX IS POST-MULTIPLIED BY A GLOBAL TRANSFORM T
 
!     LOOP ON GRID POINTS
!*****
 DO  i =1,npts
   nphi(i+1) = necpt(i+2)
   k= npts+i +8
   phiout(k) =  beta(i)/(4.0*elvol)
   icord = npts +4*i -1
   DO  j =1,18
     cmat(j,i) = cmat(j,i)/elvol
   END DO
   k= npts*2 +18*i -9
   IF(necpt(icord) /= 0) GO TO 1200
   CALL gmmats( GE,6,6,0,cmat(1,i),6,3,0,phiout(k) )
   CYCLE
   1200 CALL transs(necpt(icord), ti )
   CALL gmmats( cmat(1,i),6,3,0,ti,3,3,0,temp)
   CALL gmmats( GE,6,6,0,temp,6,3,0, phiout(k) )
 END DO
 
 nphi(1) = necpt(1)
 phiout(npts+2) = tempo
 temp(1)= alfa
 temp(2)= alfa
 temp(3)= alfa
 temp(4)= 0.0
 temp(5)= 0.0
 temp(6)= 0.0
!*****
!     THE THERMAL EXPANSION VECTOR IS MULTIPLIED BY THE STRESS-STRAIN
!       MATRIX,GE
!*****
 CALL gmmats( GE,6,6,0,temp(1),6,1,0,phiout(npts+3) )
!*****
!     THE OUTPUT ARRAY IS NOW COMPLETE
!*****
 RETURN
!*****
!     PHIOUT CONTAINS THE FOLLOWING WHERE N IS THE NUMBER OF CORNERS
 
!              ELEMENT ID
!              N SILS
!              T SUB 0
!              6 THERMAL STRESS COEFFICIENTS
!              N VOLUME RATIO COEFFICIENTS
!              N 6 BY 3 MATRICES RELATING STRESS TO DISPLACEMENTS
 
!*****
END SUBROUTINE ssold1
