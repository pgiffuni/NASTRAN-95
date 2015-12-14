SUBROUTINE stpax2 (sorc,ti)
     
!     THIS ROUTINE IS PHASE II OF STRESS RECOVERY FOR THE TRAPEZOIDAL
!     CROSS SECTION RING
 
!     OUTPUTS FROM PHASE I ARE THE FOLLOWING..
!     IDEL, IGP(4), TZ, SEL(360), TS(06), AK(144), PHI(14)
!     AKUPH(48), AKPH2(16), SELP1(120), SELP2(180), SELP3(60)
 
!     ANY GROUP OF STATEMENTS PREFACED BY AN IF STATEMENT CONTAINING
!     ...KSYS78 OR LSYS78 ...  INDICATES CODING NECESSARY FOR THIS
!     ELEMENT*S PIEZOELECTRIC CAPABILITY
 
!     KSYS78 = 0   ELASTIC, NON-PIEZOELECTRIC MATERIAL
!     KSYS78 = 1   ELECTRICAL-ELASTIC COUPLED, PIEZOELETRIC MATERIAL
!     KSYS78 = 2   ELASTIC ONLY, PIEZOELECTRIC MATERIAL
!     LSYS78 = .TRUE. IF KSYS78 = 0, OR 2
 
 
 INTEGER, INTENT(IN)                      :: sorc
 REAL, INTENT(IN)                         :: ti(4)
 LOGICAL :: zero,zeron,lsys78
 INTEGER :: iblock(62,14),istres(100),iforce(25),elemid, iclock(62,14)
 REAL :: nphi
 DIMENSION  dum3(225),stres(100),force(25),akuph(48),  &
     akph2(16),selp1(120),selp2(180),selp3(60),d4(4),  &
     d15(15),d30(30),dispp(4),echrg(4),eflux(15)
 
!     SDR2 VARIABLE CORE
 
 COMMON /zzzzzz/ zz(1)
 
!     SDR2 BLOCK FOR POINTERS AND LOADING  TEMPERATURES
 
 COMMON /sdr2x4/ dum1(33),icstm,ncstm,ivec,ivecn,templd,eldefm, dum4(12),ktype
 
!     SCRATCH BLOCK
 
 COMMON /sdr2x8/ disp(12),eforc(12),estres(30),harm,n,sinphi,  &
     conphi,nphi,nangle,elemid,unu(93),nelhar,kangle, klemid
 
!     SDR2 INPUT AND OUTPUT BLOCK
 
 COMMON /sdr2x7/ idel,igp(4),tz,sel(360),ts(6),ak(144),phi(14),  &
     dum2(424),BLOCK(62,14),clock(62,14)
 
 COMMON /system/ ksystm(77),ksys78
 COMMON /sdr2de/ dum5(33), ipart
 COMMON /condas/ consts(5)
 EQUIVALENCE     (iblock(1,1),BLOCK(1,1)),(iclock(1,1),clock(1,1)),  &
     (dum3(1),idel),(dum3(101),stres(1),istres(1)),  &
     (dum3(201),force(1),iforce(1)),(consts(4),degrad),  &
     (ldtemp,templd),(dum2(1),akuph(1)),  &
     (dum2(49),akph2(1)),(dum2(65),selp1(1)),  &
     (dum2(185),selp2(1)),(dum2(365),selp3(1)),  &
     (unu(1),d4(1)),(unu(5),d15(1)),(unu(20),d30(1))
 DATA    zeron / .false. /
 DATA    iosorc/ 0       /
 
 elemid = idel / 1000
 nelhar = idel - elemid*1000
 klemid = elemid
 lsys78 =.false.
 IF (ksys78 == 0 .OR. ksys78 == 2) lsys78 = .true.
 
!     SET BLOCK = 0 IF HARMONIC = 0
 
 n = nelhar - 1
 IF (n /= 0) GO TO 21
 IF (n == 0 .AND. zeron .AND. iosorc /= sorc) GO TO 14
 zeron  = .true.
 iosorc = sorc
 DO  i = 2,62
   DO  j = 1,14
     IF (ktype /= 2 .OR. ipart /= 2) BLOCK(i,j) = 0.0
     clock(i,j) = 0.0
   END DO
 END DO
 
!     SET ANGLES CONTROL FOR SUMMATION
 
 zero = .false.
 j = 0
 DO  i = 1,14
   IF (phi(i) == 0.0) THEN
     GO TO    18
   ELSE
     GO TO    17
   END IF
   18 IF (zero) CYCLE
   zero = .true.
   17 j = j + 1
   BLOCK(1,j) = phi(i)
   clock(1,j) = phi(i)
 END DO
 j = j + 1
 IF (j > 14) GO TO 21
 iblock(1,j) = 1
 iclock(1,j) = 1
 GO TO 21
 14 zeron = .false.
 21 harm  = n
 
!     INITIALIZE LOCAL VARIABLES
 
 ndof  = 3
 numpt = 4
 n     = ndof*numpt
 nsp   = 5
 ncomp = 6
 ns    = nsp*ncomp
 
!     FIND GRID POINTS DISPLACEMENTS
 
 k = 0
 DO  i = 1,numpt
   iloc = ivec + igp(i) - 2
   
   IF (lsys78) GO TO 90
   ilocp = iloc + 4
   dispp(i) = zz(ilocp)
   90 CONTINUE
   
   DO  j = 1,ndof
     iloc = iloc + 1
     k = k + 1
     disp(k) = zz(iloc)
   END DO
 END DO
 
!     COMPUTE THE GRID POINT FORCES
 
 CALL gmmats (ak(1),n,n,0, disp(1),n,1,0, eforc(1))
 
 DO  i = 1,4
   echrg(i) = 0.0
 END DO
 
 IF (lsys78) GO TO 125
 CALL gmmats (akuph(1),n,numpt,0, dispp(1),numpt,1,0, d15(1))
 DO  i = 1,12
   eforc(i) = eforc(i) + d15(i)
 END DO
 
 CALL gmmats (akuph(1),n,numpt,1, disp(1),n,1,0, d4(1))
 CALL gmmats (akph2(1),numpt,numpt,0, dispp(1),numpt,1,0, echrg(1))
 DO  i = 1,4
   echrg(i) = echrg(i) + d4(i)
 END DO
 125 CONTINUE
 
!     COMPUTE THE STRESSES
 
 CALL gmmats (sel(1),ns,n,0, disp(1),n,1,0, estres(1))
 
 DO  i = 1,15
   eflux(i) = 0.0
 END DO
 
 IF (lsys78) GO TO 145
 CALL gmmats (selp1(1),ns,numpt,0, dispp(1),numpt,1,0, d30(1))
 DO  i = 1,30
   estres(i) = estres(i) + d30(i)
 END DO
 
 CALL gmmats (selp2(1),15,n,0, disp(1),n,1,0, eflux(1))
 CALL gmmats (selp3(1),15,numpt,0, dispp(1),numpt,1,0, d15(1))
 DO  i = 1,15
   eflux(i) = eflux(i) + d15(i)
 END DO
 145 CONTINUE
 
!     COMPUTE THERMAL STRESS IF IT IS EXISTS
 
 IF (ldtemp == -1) GO TO 300
 k = 0
 t = tz
 IF (harm > 0.0) t = 0.0
 DO  i = 1,nsp
   dt = ti(i) - t
   IF (i == 5) dt = (ti(1)+ti(2)+ti(3)+ti(4))/4.0 - t
   DO  j = 1,ncomp
     k = k + 1
     estres(k) = estres(k) - dt*ts(j)
   END DO
 END DO
 300 CONTINUE
 
!     BRANCH TO INSERT HARMONIC STRESSES AND FORCES INTO BLOCK OR CLOCK
 
!     KTYPE = 1 - REAL OUTPUT, STORED IN BLOCK, NOTHING IN CLOCK
!     KTYPE = 2 - COMPLEX OUTPUT
!     IPART = 1 - IMAGINARY PART OF COMPLEX OUTPUT, STORED IN BLOCK
!     IPART = 2 - REAL PART OF COMPLEX OUTPUT, STORED IN CLOCK
 
 IF (ktype == 2 .AND. ipart == 2) GO TO 550
 
!     INSERT HARMONIC STRESSES AND FORCES INTO BLOCK
 
 DO  i = 1,14
   IF (iblock(1,i) == 1) GO TO 380
   IF (harm == 0.0) GO TO 350
   nphi   = harm*BLOCK(1,i)*degrad
   sinphi = SIN(nphi)
   conphi = COS(nphi)
   SELECT CASE ( sorc )
     CASE (    1)
       GO TO 330
     CASE (    2)
       GO TO 310
   END SELECT
   
   310 CONTINUE
   DO  ie = 1,5
     ke   = 9*(ie-1)
     kepz = 6*(ie-1)
     BLOCK(2+ke,i) = BLOCK(2+ke,i) + conphi*estres(1+kepz)
     BLOCK(3+ke,i) = BLOCK(3+ke,i) + conphi*estres(2+kepz)
     BLOCK(4+ke,i) = BLOCK(4+ke,i) + conphi*estres(3+kepz)
     BLOCK(5+ke,i) = BLOCK(5+ke,i) + conphi*estres(4+kepz)
     BLOCK(6+ke,i) = BLOCK(6+ke,i) + sinphi*estres(5+kepz)
     BLOCK(7+ke,i) = BLOCK(7+ke,i) + sinphi*estres(6+kepz)
     
     IF (lsys78) CYCLE
     kepz2 = kepz/2
     BLOCK( 8+ke,i) = BLOCK( 8+ke,i) + conphi*eflux (1+kepz2)
     BLOCK( 9+ke,i) = BLOCK( 9+ke,i) + conphi*eflux (2+kepz2)
     BLOCK(10+ke,i) = BLOCK(10+ke,i) + sinphi*eflux (3+kepz2)
   END DO
   
   DO  ir = 1,4
     kr   = 4*(ir-1)
     krpz = 3*(ir-1)
     BLOCK(47+kr,i) = BLOCK(47+kr,i) + conphi*eforc(1+krpz)
     BLOCK(48+kr,i) = BLOCK(48+kr,i) + sinphi*eforc(2+krpz)
     BLOCK(49+kr,i) = BLOCK(49+kr,i) + conphi*eforc(3+krpz)
     kr3 = 1 + krpz/3
     IF(.NOT.lsys78) BLOCK(50+kr,i) = BLOCK(50+kr,i) +conphi*echrg(kr3)
   END DO
   CYCLE
   
   330 CONTINUE
   DO  ie = 1,5
     ke   = 9*(ie-1)
     kepz = 6*(ie-1)
     BLOCK(2+ke,i) = BLOCK(2+ke,i) + sinphi*estres(1+kepz)
     BLOCK(3+ke,i) = BLOCK(3+ke,i) + sinphi*estres(2+kepz)
     BLOCK(4+ke,i) = BLOCK(4+ke,i) + sinphi*estres(3+kepz)
     BLOCK(5+ke,i) = BLOCK(5+ke,i) + sinphi*estres(4+kepz)
     BLOCK(6+ke,i) = BLOCK(6+ke,i) - conphi*estres(5+kepz)
     BLOCK(7+ke,i) = BLOCK(7+ke,i) - conphi*estres(6+kepz)
     
     IF (lsys78) CYCLE
     kepz2 = kepz/2
     BLOCK( 8+ke,i) = BLOCK( 8+ke,i) + sinphi*eflux(1+kepz2)
     BLOCK( 9+ke,i) = BLOCK( 9+ke,i) + sinphi*eflux(2+kepz2)
     BLOCK(10+ke,i) = BLOCK(10+ke,i) - conphi*eflux(3+kepz2)
   END DO
   
   DO  ir = 1,4
     kr   = 4*(ir-1)
     krpz = 3*(ir-1)
     BLOCK(47+kr,i) = BLOCK(47+kr,i) + sinphi*eforc(1+krpz)
     BLOCK(48+kr,i) = BLOCK(48+kr,i) - conphi*eforc(2+krpz)
     BLOCK(49+kr,i) = BLOCK(49+kr,i) + sinphi*eforc(3+krpz)
     kr3 = 1 + krpz/3
     IF(.NOT.lsys78) BLOCK(50+kr,i) = BLOCK(50+kr,i) +sinphi*echrg(kr3)
   END DO
   CYCLE
   
   350 DO  ie = 1,5
     ke   = 9*(ie-1)
     kepz = 6*(ie-1)
     BLOCK(2+ke,i) = estres(1+kepz)
     BLOCK(3+ke,i) = estres(2+kepz)
     BLOCK(4+ke,i) = estres(3+kepz)
     BLOCK(5+ke,i) = estres(4+kepz)
     BLOCK(6+ke,i) = estres(5+kepz)
     BLOCK(7+ke,i) = estres(6+kepz)
     
     IF (lsys78) CYCLE
     kepz2 = kepz/2
     BLOCK( 8+ke,i) = eflux(1+kepz2)
     BLOCK( 9+ke,i) = eflux(2+kepz2)
     BLOCK(10+ke,i) = eflux(3+kepz2)
   END DO
   
   DO  ir = 1,4
     kr   = 4*(ir-1)
     krpz = 3*(ir-1)
     BLOCK(47+kr,i) = eforc(1+krpz)
     BLOCK(48+kr,i) = eforc(2+krpz)
     BLOCK(49+kr,i) = eforc(3+krpz)
     kr3 = 1 + krpz/3
     IF(.NOT.lsys78) BLOCK(50+kr,i) = echrg(kr3)
   END DO
   
 END DO
 
!     COPY STRESSES AND FORCES INTO OUTPUT BLOCKS
 
 380 CONTINUE
 j = 2
 k = 1
 l = 0
 istres (1) = elemid
 istres (2) = nelhar
 DO  i = 1,ns
   j = j + 1
   stres(j) = estres(i)
   
   IF (i/6 /= k) CYCLE
   k = k + 1
   DO  ii = 1,3
     j = j + 1
     l = l + 1
     stres(j) = eflux(l)
   END DO
   
 END DO
 k = 0
 j = 2
 l = 1
 iforce(1) = elemid
 iforce(2) = nelhar
 DO  i  = 1,numpt
   DO  kk = 1,ndof
     j = j + 1
     k = k + 1
     force(j) = eforc(k)
     
     IF (k/3 /= l) CYCLE
     j = j + 1
     force(j) = echrg(l)
     l = l + 1
     
   END DO
 END DO
 
 IF (ktype == 1 .OR. (ktype == 2 .AND. ipart == 1)) GO TO 1001
 550 CONTINUE
 
!     INSERT HARMONIC STRESSES AND FORCES INTO CLOCK
 
 DO  i = 1,14
   IF (iclock(1,i) == 1) GO TO 700
   IF (harm == 0.0) GO TO 660
   nphi   = harm*clock(1,i)*degrad
   sinphi = SIN(nphi)
   conphi = COS(nphi)
   SELECT CASE ( sorc )
     CASE (    1)
       GO TO 630
     CASE (    2)
       GO TO 600
   END SELECT
   600 CONTINUE
   
   DO  ie = 1,5
     ke   = 9*(ie-1)
     kepz = 6*(ie-1)
     clock(2+ke,i) = clock(2+ke,i) + conphi*estres(1+kepz)
     clock(3+ke,i) = clock(3+ke,i) + conphi*estres(2+kepz)
     clock(4+ke,i) = clock(4+ke,i) + conphi*estres(3+kepz)
     clock(5+ke,i) = clock(5+ke,i) + conphi*estres(4+kepz)
     clock(6+ke,i) = clock(6+ke,i) + sinphi*estres(5+kepz)
     clock(7+ke,i) = clock(7+ke,i) + sinphi*estres(6+kepz)
     
     IF (lsys78) CYCLE
     kepz2 = kepz/2
     clock( 8+ke,i) = clock( 8+ke,i) + conphi*eflux (1+kepz2)
     clock( 9+ke,i) = clock( 9+ke,i) + conphi*eflux (2+kepz2)
     clock(10+ke,i) = clock(10+ke,i) + sinphi*eflux (3+kepz2)
   END DO
   
   DO  ir = 1,4
     kr   = 4*(ir-1)
     krpz = 3*(ir-1)
     clock(47+kr,i) = clock(47+kr,i) + conphi*eforc(1+krpz)
     clock(48+kr,i) = clock(48+kr,i) + sinphi*eforc(2+krpz)
     clock(49+kr,i) = clock(49+kr,i) + conphi*eforc(3+krpz)
     kr3 = 1 + krpz/3
     IF(.NOT.lsys78) clock(50+kr,i) = clock(50+kr,i) +conphi*echrg(kr3)
   END DO
   CYCLE
   
   630 CONTINUE
   DO  ie = 1,5
     ke   = 9*(ie-1)
     kepz = 6*(ie-1)
     clock(2+ke,i) = clock(2+ke,i) + sinphi*estres(1+kepz)
     clock(3+ke,i) = clock(3+ke,i) + sinphi*estres(2+kepz)
     clock(4+ke,i) = clock(4+ke,i) + sinphi*estres(3+kepz)
     clock(5+ke,i) = clock(5+ke,i) + sinphi*estres(4+kepz)
     clock(6+ke,i) = clock(6+ke,i) - conphi*estres(5+kepz)
     clock(7+ke,i) = clock(7+ke,i) - conphi*estres(6+kepz)
     
     IF (lsys78) CYCLE
     kepz2 = kepz/2
     clock( 8+ke,i) = clock( 8+ke,i) + sinphi*eflux(1+kepz2)
     clock( 9+ke,i) = clock( 9+ke,i) + sinphi*eflux(2+kepz2)
     clock(10+ke,i) = clock(10+ke,i) - conphi*eflux(3+kepz2)
   END DO
   
   DO  ir = 1,4
     kr   = 4*(ir-1)
     krpz = 3*(ir-1)
     clock(47+kr,i) = clock(47+kr,i) + sinphi*eforc(1+krpz)
     clock(48+kr,i) = clock(48+kr,i) - conphi*eforc(2+krpz)
     clock(49+kr,i) = clock(49+kr,i) + sinphi*eforc(3+krpz)
     kr3 = 1 + krpz/3
     IF(.NOT.lsys78) clock(50+kr,i) = clock(50+kr,i) +sinphi*echrg(kr3)
   END DO
   CYCLE
   
   660 DO  ie = 1,5
     ke   = 9*(ie-1)
     kepz = 6*(ie-1)
     clock(2+ke,i) = estres(1+kepz)
     clock(3+ke,i) = estres(2+kepz)
     clock(4+ke,i) = estres(3+kepz)
     clock(5+ke,i) = estres(4+kepz)
     clock(6+ke,i) = estres(5+kepz)
     clock(7+ke,i) = estres(6+kepz)
     
     IF (lsys78) CYCLE
     kepz2 = kepz/2
     clock( 8+ke,i) = eflux(1+kepz2)
     clock( 9+ke,i) = eflux(2+kepz2)
     clock(10+ke,i) = eflux(3+kepz2)
   END DO
   
   DO  ir = 1,4
     kr   = 4*(ir-1)
     krpz = 3*(ir-1)
     clock(47+kr,i) = eforc(1+krpz)
     clock(48+kr,i) = eforc(2+krpz)
     clock(49+kr,i) = eforc(3+krpz)
     kr3 = 1 + krpz/3
     IF(.NOT.lsys78) clock(50+kr,i) = echrg(kr3)
   END DO
   
 END DO
 
!     COPY STRESSES AND FORCES INTO OUTPUT BLOCKS
 
 700 CONTINUE
 j = 2
 k = 1
 l = 0
 istres (1) = elemid
 istres (2) = nelhar
 DO  i = 1,ns
   j = j + 1
   stres(j) = estres(i)
   
   IF (i/6 /= k) CYCLE
   k = k + 1
   DO  ii = 1,3
     j = j + 1
     l = l + 1
     stres(j) = eflux(l)
   END DO
 END DO
 
 k = 0
 j = 2
 l = 1
 iforce(1) = elemid
 iforce(2) = nelhar
 DO  i  = 1,numpt
   DO  kk = 1,ndof
     j = j + 1
     k = k + 1
     force(j) = eforc(k)
     
     IF (k/3 /= l) CYCLE
     j = j + 1
     force(j) = echrg(l)
     l = l + 1
   END DO
 END DO
 
 1001 CONTINUE
 
 RETURN
END SUBROUTINE stpax2
