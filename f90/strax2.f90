SUBROUTINE strax2 (sorc,ti)
     
!     THIS ROUTINE IS PHASE II OF STRESS DATA FOR THE TRIANGULAR
!     CROSS SECTION RING
 
!     OUTPUTS FROM PHASE I ARE THE FOLLOWING -
!     IDEL IGP(3) TZ SEL(54) TS(4) AK(81) PHI(14)
!     AKUPH(27) AKPH2(9) SELP1(18) SELP2(27) SELP3(9)
 
!     ANY GROUP OF STATEMENTS PREFACED BY AN IF STATEMENT CONTAINING
!     ...KSYS78 OR LSYS78 ...  INDICATES CODING NECESSARY FOR THIS
!     ELEMENT*S PIEZOELECTRIC CAPABILITY
 
!     KSYS78 = 0   ELASTIC, NON-PIEZOELECTRIC MATERIAL
!     KSYS78 = 1   ELECTRICAL-ELASTIC COUPLED, PIEZOELETRIC MATERIAL
!     KSYS78 = 2   ELASTIC ONLY, PIEZOELECTRIC MATERIAL
!     LSYS78 = .TRUE. IF KSYS78 = 0, OR 2
 
 
 INTEGER, INTENT(IN)                      :: sorc
 REAL, INTENT(IN)                         :: ti(3)
 LOGICAL :: zero,zeron,lsys78
 INTEGER :: iblock(22,14),istres(100),iforce(25),elemid, iclock(22,14)
 REAL :: nphi
 DIMENSION  dum3(225),stres(100),force(25),akuph(27),  &
     akph2(9),selp1(18),selp2(27),selp3(9),d3(3),d6(6),  &
     d9(9),dispp(3),echrg(3),eflux(3)
 
!     SDR2 VARIABLE CORE
 
 COMMON /zzzzzz/ zz(1)
 
!     SDR2 BLOCK FOR POINTERS AND LOADING  TEMPERATURES
 
 COMMON /sdr2x4/ dum1(33),icstm,ncstm,ivec,ivecn,templd,eldefm, dum4(12),ktype
 
!     SCRATCH BLOCK
 
 COMMON /sdr2x8/ disp(9),eforc(9),estres(9),harm,n,sinphi,conphi,  &
     nphi,nangle,elemid,unu(123),nelhar
 
!     SDR2 INPUT AND OUTPUT BLOCK
 
 COMMON /sdr2x7/ idel,igp(3),tz,sel(54),ts(6),ak(81),phi(14),  &
     dum2(90),BLOCK(22,14),clock(22,14)
 
 COMMON /sdr2de/ dum5(33), ipart
 COMMON /condas/ consts(5)
 COMMON /system/ ksystm(77),ksys78
 EQUIVALENCE     (iblock(1,1),BLOCK(1,1)),(iclock(1,1),clock(1,1)),  &
     (dum3(1),idel),(ldtemp,templd), (dum3(109),stres(9),istres(9),eflux(1)),  &
     (dum3(201),force(1),iforce(1)), (dum2(1),selp1(1)),(dum2(19),akph2(1)),  &
     (dum2(28),akuph(1)),(dum2(55),selp2(1)),  &
     (dum2(82),selp3(1)),(consts(4),degrad),  &
     (unu(1),d3(1)),(unu(4),d6(1)),(unu(10),d9(1))
 DATA    zeron / .false. /
 DATA    iosorc/ 0       /
 
 lsys78 = .false.
 IF (ksys78 == 0 .OR. ksys78 == 2) lsys78 = .true.
 
 elemid = idel/1000
 nelhar = idel - elemid*1000
 
!     SET BLOCK = 0 IF HARMONIC = 0
 
 n = nelhar - 1
 IF (n /= 0) GO TO 21
 IF (n == 0 .AND. zeron .AND. iosorc /= sorc) GO TO 14
 zeron  = .true.
 iosorc = sorc
 DO  i = 2,22
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
 
!     INITIALIZE LOCAT VARIABLES
 
 ndof  = 3
 numpt = 3
 n     = ndof*numpt
 nsp   = 1
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
     k    = k + 1
     disp(k) = zz(iloc)
   END DO
 END DO
 
!     COMPUTE THE GRID POINT FORCES
 
 CALL gmmats (ak(1),n,n,0, disp(1),n,1,0, eforc(1))
 
 DO  i = 1,3
   echrg(i) = 0.0
 END DO
 
 IF (lsys78) GO TO 125
 CALL gmmats (akuph(1),n,numpt,0, dispp(1),numpt,1,0, d9(1))
 DO  i = 1,9
   eforc(i) = eforc(i) + d9(i)
 END DO
 
 CALL gmmats (akuph(1),n,numpt,1, disp(1),n,1,0, d3(1))
 CALL gmmats (akph2(1),numpt,numpt,0, dispp(1),numpt,1,0, echrg(1))
 DO  i = 1,3
   echrg(i) = echrg(i) + d3(i)
 END DO
 
!     COMPUTE THE STRESSES
 
 125 CALL gmmats (sel(1),ns,n,0, disp(1),n,1,0, estres(1))
 
 DO  i = 1,3
   eflux(i) = 0.0
 END DO
 
 IF (lsys78) GO TO 145
 CALL gmmats (selp1(1),ns,numpt,0, dispp(1),numpt,1,0, d6(1))
 DO  i = 1,6
   estres(i) = estres(i) + d6(i)
 END DO
 
 CALL gmmats (selp2(1),numpt,n,0, disp(1),n,1,0, eflux(1))
 CALL gmmats (selp3(1),numpt,numpt,0, dispp(1),numpt,1,0, d3(1))
 
 DO  i = 1,3
   eflux(i) = eflux(i) + d3(i)
 END DO
 
!     COMPUTE THERMAL STRESS IF IT IS EXISTS
 
 145 IF (ldtemp == -1) GO TO 300
 dt = tz
 IF (harm > 0.0) dt = 0.0
 dt = (ti(1)+ti(2)+ti(3))/3.0 - dt
 DO  i = 1,ns
   estres(i) = estres(i) - dt*ts(i)
 END DO
 
!     BRANCH TO INSERT HARMONIC STRESSES AND FORCES INTO BLOCK OR CLOCK
 
!     KTYPE = 1 - REAL OUTPUT, STORED IN BLOCK, NOTHING IN CLOCK
!     KTYPE = 2 - COMPLEX OUTPUT
!     IPART = 1 - IMAGINARY PART OF COMPLEX OUTPUT, STORED IN BLOCK
!     IPART = 2 - REAL PART OF COMPLEX OUTPUT, STORED IN CLOCK
 
 300 IF (ktype == 2 .AND. ipart == 2) GO TO 505
 
!     INSERT HARMONIC STRESSES AND FORCES INTO BLOCK
 
 DO  i = 1,14
   IF (iblock(1,i) == 1) GO TO 390
   IF (harm /= 0.0) GO TO 330
   DO  iwa = 1,6
     BLOCK(iwa+1,i) = estres(iwa)
     BLOCK(iwa+7,i) = eforc (iwa)
   END DO
   BLOCK(14,i) = eforc(7)
   BLOCK(15,i) = eforc(8)
   BLOCK(16,i) = eforc(9)
   
   IF (lsys78) GO TO 320
   BLOCK(17,i) = eflux(1)
   BLOCK(18,i) = eflux(2)
   BLOCK(19,i) = eflux(3)
   BLOCK(20,i) = echrg(1)
   BLOCK(21,i) = echrg(2)
   BLOCK(22,i) = echrg(3)
   320 CONTINUE
   CYCLE
   330 CONTINUE
   nphi   = harm*BLOCK(1,i)*degrad
   sinphi = SIN(nphi)
   conphi = COS(nphi)
   
   SELECT CASE ( sorc )
     CASE (    1)
       GO TO 360
     CASE (    2)
       GO TO 340
   END SELECT
   
   340 BLOCK( 2,i) = BLOCK( 2,i) + conphi*estres(1)
   BLOCK( 3,i) = BLOCK( 3,i) + conphi*estres(2)
   BLOCK( 4,i) = BLOCK( 4,i) + conphi*estres(3)
   BLOCK( 5,i) = BLOCK( 5,i) + conphi*estres(4)
   BLOCK( 6,i) = BLOCK( 6,i) + sinphi*estres(5)
   BLOCK( 7,i) = BLOCK( 7,i) + sinphi*estres(6)
   BLOCK( 8,i) = BLOCK( 8,i) + conphi*eforc(1)
   BLOCK( 9,i) = BLOCK( 9,i) + sinphi*eforc(2)
   BLOCK(10,i) = BLOCK(10,i) + conphi*eforc(3)
   BLOCK(11,i) = BLOCK(11,i) + conphi*eforc(4)
   BLOCK(12,i) = BLOCK(12,i) + sinphi*eforc(5)
   BLOCK(13,i) = BLOCK(13,i) + conphi*eforc(6)
   BLOCK(14,i) = BLOCK(14,i) + conphi*eforc(7)
   BLOCK(15,i) = BLOCK(15,i) + sinphi*eforc(8)
   BLOCK(16,i) = BLOCK(16,i) + conphi*eforc(9)
   IF (lsys78) GO TO 350
   BLOCK(17,i) = BLOCK(17,i) + conphi*eflux(1)
   BLOCK(18,i) = BLOCK(18,i) + conphi*eflux(2)
   BLOCK(19,i) = BLOCK(19,i) + sinphi*eflux(3)
   BLOCK(20,i) = BLOCK(20,i) + conphi*echrg(1)
   BLOCK(21,i) = BLOCK(21,i) + conphi*echrg(2)
   BLOCK(22,i) = BLOCK(22,i) + conphi*echrg(3)
   350 CONTINUE
   CYCLE
   360 BLOCK( 2,i) = BLOCK( 2,i) + sinphi*estres(1)
   BLOCK( 3,i) = BLOCK( 3,i) + sinphi*estres(2)
   BLOCK( 4,i) = BLOCK( 4,i) + sinphi*estres(3)
   BLOCK( 5,i) = BLOCK( 5,i) + sinphi*estres(4)
   BLOCK( 6,i) = BLOCK( 6,i) - conphi*estres(5)
   BLOCK( 7,i) = BLOCK( 7,i) - conphi*estres(6)
   BLOCK( 8,i) = BLOCK( 8,i) + sinphi*eforc(1)
   BLOCK( 9,i) = BLOCK( 9,i) - conphi*eforc(2)
   BLOCK(10,i) = BLOCK(10,i) + sinphi*eforc(3)
   BLOCK(11,i) = BLOCK(11,i) + sinphi*eforc(4)
   BLOCK(12,i) = BLOCK(12,i) - conphi*eforc(5)
   BLOCK(13,i) = BLOCK(13,i) + sinphi*eforc(6)
   BLOCK(14,i) = BLOCK(14,i) + sinphi*eforc(7)
   BLOCK(15,i) = BLOCK(15,i) - conphi*eforc(8)
   BLOCK(16,i) = BLOCK(16,i) - sinphi*eforc(9)
   IF (lsys78) GO TO 370
   BLOCK(17,i) = BLOCK(17,i) + sinphi*eflux(1)
   BLOCK(18,i) = BLOCK(18,i) + sinphi*eflux(2)
   BLOCK(19,i) = BLOCK(19,i) - conphi*eflux(3)
   BLOCK(20,i) = BLOCK(20,i) + sinphi*echrg(1)
   BLOCK(21,i) = BLOCK(21,i) + sinphi*echrg(2)
   BLOCK(22,i) = BLOCK(22,i) + sinphi*echrg(3)
   370 CONTINUE
 END DO
 
!     COPY STRESSES AND FORCES INTO OUTPUT BLOCKS
!     FLUXES ARE EQUIVALENCED INTO STRES(J)
!     CHARGES ARE WRITTEN INTO FORCE(J)
 
 390 j = 2
 istres (1) = elemid
 istres (2) = nelhar
 DO  i = 1,ncomp
   j = j + 1
   stres(j) = estres(i)
 END DO
 k = 0
 j = 2
 iforce(1) = elemid
 iforce(2) = nelhar
 DO  i  = 1,numpt
   DO  kk = 1,ndof
     j = j + 1
     k = k + 1
     force(j) = eforc(k)
     
     IF (k /= 3 .AND. k /= 6 .AND. k /= 9) CYCLE
     j  = j + 1
     k3 = k/3
     force(j) = echrg(k3)
   END DO
 END DO
 
 IF (ktype == 1 .OR. (ktype == 2 .AND. ipart == 1)) GO TO 1000
 
!     INSERT HARMONIC STRESSES AND FORCES INTO CLOCK
 
 505 DO  i = 1,14
   IF (iclock(1,i) == 1) GO TO 600
   IF (harm /= 0.0) GO TO 530
   DO  iwa = 1,6
     clock(iwa+1,i) = estres(iwa)
     clock(iwa+7,i) = eforc (iwa)
   END DO
   clock(14,i) = eforc(7)
   clock(15,i) = eforc(8)
   clock(16,i) = eforc(9)
   
   IF (lsys78) GO TO 520
   clock(17,i) = eflux(1)
   clock(18,i) = eflux(2)
   clock(19,i) = eflux(3)
   clock(20,i) = echrg(1)
   clock(21,i) = echrg(2)
   clock(22,i) = echrg(3)
   520 CONTINUE
   CYCLE
   530 CONTINUE
   nphi   = harm*clock(1,i)*degrad
   sinphi = SIN(nphi)
   conphi = COS(nphi)
   
   SELECT CASE ( sorc )
     CASE (    1)
       GO TO 560
     CASE (    2)
       GO TO 540
   END SELECT
   
   540 clock( 2,i) = clock( 2,i) + conphi*estres(1)
   clock( 3,i) = clock( 3,i) + conphi*estres(2)
   clock( 4,i) = clock( 4,i) + conphi*estres(3)
   clock( 5,i) = clock( 5,i) + conphi*estres(4)
   clock( 6,i) = clock( 6,i) + sinphi*estres(5)
   clock( 7,i) = clock( 7,i) + sinphi*estres(6)
   clock( 8,i) = clock( 8,i) + conphi*eforc(1)
   clock( 9,i) = clock( 9,i) + sinphi*eforc(2)
   clock(10,i) = clock(10,i) + conphi*eforc(3)
   clock(11,i) = clock(11,i) + conphi*eforc(4)
   clock(12,i) = clock(12,i) + sinphi*eforc(5)
   clock(13,i) = clock(13,i) + conphi*eforc(6)
   clock(14,i) = clock(14,i) + conphi*eforc(7)
   clock(15,i) = clock(15,i) + sinphi*eforc(8)
   clock(16,i) = clock(16,i) + conphi*eforc(9)
   
   IF (lsys78) GO TO 550
   clock(17,i) = clock(17,i) + conphi*eflux(1)
   clock(18,i) = clock(18,i) + conphi*eflux(2)
   clock(19,i) = clock(19,i) + sinphi*eflux(3)
   clock(20,i) = clock(20,i) + conphi*echrg(1)
   clock(21,i) = clock(21,i) + conphi*echrg(2)
   clock(22,i) = clock(22,i) + conphi*echrg(3)
   550 CONTINUE
   CYCLE
   
   560 clock( 2,i) = clock( 2,i) + sinphi*estres(1)
   clock( 3,i) = clock( 3,i) + sinphi*estres(2)
   clock( 4,i) = clock( 4,i) + sinphi*estres(3)
   clock( 5,i) = clock( 5,i) + sinphi*estres(4)
   clock( 6,i) = clock( 6,i) - conphi*estres(5)
   clock( 7,i) = clock( 7,i) - conphi*estres(6)
   clock( 8,i) = clock( 8,i) + sinphi*eforc(1)
   clock( 9,i) = clock( 9,i) - conphi*eforc(2)
   clock(10,i) = clock(10,i) + sinphi*eforc(3)
   clock(11,i) = clock(11,i) + sinphi*eforc(4)
   clock(12,i) = clock(12,i) - conphi*eforc(5)
   clock(13,i) = clock(13,i) + sinphi*eforc(6)
   clock(14,i) = clock(14,i) + sinphi*eforc(7)
   clock(15,i) = clock(15,i) - conphi*eforc(8)
   clock(16,i) = clock(16,i) + sinphi*eforc(9)
   IF (lsys78) GO TO 570
   clock(17,i) = clock(17,i) + sinphi*eflux(1)
   clock(18,i) = clock(18,i) + sinphi*eflux(2)
   clock(19,i) = clock(19,i) - conphi*eflux(3)
   clock(20,i) = clock(20,i) + sinphi*echrg(1)
   clock(21,i) = clock(21,i) + sinphi*echrg(2)
   clock(22,i) = clock(22,i) + sinphi*echrg(3)
   570 CONTINUE
 END DO
 
!     COPY STRESSES AND FORCES INTO OUTPUT BLOCKS
!     FLUXES ARE EQUIVALENCED INTO STRES(J)
!     CHARGES ARE WRITTEN INTO FORCE(J)
 
 600 j = 2
 istres (1) = elemid
 istres (2) = nelhar
 DO  i = 1,ncomp
   j = j + 1
   stres(j) = estres(i)
 END DO
 k = 0
 j = 2
 iforce(1) = elemid
 iforce(2) = nelhar
 DO  i  = 1,numpt
   DO  kk = 1,ndof
     j = j + 1
     k = k + 1
     force(j) = eforc(k)
     
     IF (k /= 3 .AND. k /= 6 .AND. k /= 9) CYCLE
     j  = j + 1
     k3 = k/3
     force(j) = echrg(k3)
   END DO
 END DO
 
 1000 RETURN
END SUBROUTINE strax2
