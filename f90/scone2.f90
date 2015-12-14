SUBROUTINE scone2 (sorc)
     
!     PHASE II OF STRESS DATA RECOVERY
 
!     OUTPUTS FROM PHASE I ARE THE FOLLOWING (TOTAL OF 118 WORDS) -
!     1) ELEMENT ID
!     2 AND 3) SILS A AND B
!     4) S SUB T
!     5) N
!     6) I
!     7) Z1
!     8) Z2
!     9 THRU 22) PHI-S
!     23 THRU 118) TWO 8X6 S MATRICES
 
 
 INTEGER, INTENT(IN OUT)                  :: sorc
 LOGICAL :: zero
 INTEGER :: sil(2),iforce(8),istres(100),elemid, iblock(9,14)
 REAL :: nphi,phi(14),force(7),s(96),stress(18),zoff(2),iii
 COMMON /condas/ consts(5)
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,dum11(3)
 COMMON /sdr2x7/ commun(225)
 COMMON /sdr2x8/ vec(8),sum(8),sig(3),sig1,sig2,sig12,temp,  &
     delta,theta,npoint,zoveri,ipt,BLOCK(9,14),  &
     nelhar,elemid,harm,n,sinphi,conphi,nphi,nangle
 EQUIVALENCE     (consts(4),degra),(nelem,commun(1)),  &
     (sil(1),commun(2)),(iii,commun(6)),  &
     (zoff(1),commun(7)),(phi(1),commun(9)),  &
     (s(1),commun(23)),(iblock(1,1),BLOCK(1,1)),  &
     (stress(1),commun(101),istres(1)), (force (1),commun(201),iforce(1))
 DATA    nelold/ -1 /
 
 DO  i = 1,8
   sum(i) = 0.0
 END DO
 
 elemid = nelem/1000
 nelhar = nelem - elemid*1000
 
!     ZERO OUT BLOCK IF THIS IS FIRST CALL WITH HARMONIC = 0 FOR THIS
!     ELEMENT
 
 n = nelhar - 1
 IF (n /= 0) GO TO 21
 IF (elemid == nelold) GO TO 21
 nelold = elemid
 DO  i = 2,9
   DO  j = 1,14
     BLOCK(i,j) = 0.0
   END DO
 END DO
 
!     INSERT ANGLES FOR OUTPUT INTO FIRST ROW OF BLOCK
 
 zero = .false.
 j = 0
 DO  i = 1,14
   IF (phi(i) == 0.0) THEN
     GO TO    15
   ELSE
     GO TO    17
   END IF
   15 IF (zero) CYCLE
   zero = .true.
   17 j = j + 1
   BLOCK(1,j) = phi(i)
 END DO
 j = j + 1
 IF (j <= 14) iblock(1,j) = 1
 21 harm = n
 
 DO  i = 1,2
   
!     DISPLACEMENT VECTOR POINTER
   
   npoint = ivec + sil(i) - 1
   
   CALL gmmats (s(48*i-47),8,6,0, z(npoint),6,1,0, vec(1))
   
   DO  j = 1,8
     sum(j) = sum(j) + vec(j)
   END DO
 END DO
 
!     INSERT HARMONIC STRESSES AND FORCES INTO BLOCK FOR THIS HARMONIC
 
 DO  i = 1,14
   IF (iblock(1,i) == 1) GO TO 50
   nphi   = harm*BLOCK(1,i)*degra
   sinphi = SIN(nphi)
   conphi = COS(nphi)
   SELECT CASE ( sorc )
     CASE (    1)
       GO TO 35
     CASE (    2)
       GO TO 36
   END SELECT
   35 BLOCK(2,i) = BLOCK(2,i) + sinphi*sum(1)
   BLOCK(3,i) = BLOCK(3,i) + sinphi*sum(2)
   BLOCK(4,i) = BLOCK(4,i) - conphi*sum(3)
   BLOCK(5,i) = BLOCK(5,i) + sinphi*sum(4)
   BLOCK(6,i) = BLOCK(6,i) + sinphi*sum(5)
   BLOCK(7,i) = BLOCK(7,i) - conphi*sum(6)
   BLOCK(8,i) = BLOCK(8,i) + sinphi*sum(7)
   BLOCK(9,i) = BLOCK(9,i) - conphi*sum(8)
   CYCLE
   36 BLOCK(2,i) = BLOCK(2,i) + conphi*sum(1)
   BLOCK(3,i) = BLOCK(3,i) + conphi*sum(2)
   BLOCK(4,i) = BLOCK(4,i) + sinphi*sum(3)
   BLOCK(5,i) = BLOCK(5,i) + conphi*sum(4)
   BLOCK(6,i) = BLOCK(6,i) + conphi*sum(5)
   BLOCK(7,i) = BLOCK(7,i) + sinphi*sum(6)
   BLOCK(8,i) = BLOCK(8,i) + conphi*sum(7)
   BLOCK(9,i) = BLOCK(9,i) + sinphi*sum(8)
 END DO
 
!     COPY FORCES INTO FORCE OUTPUT BLOCK
 
 50 iforce(1) = elemid
 iforce(2) = nelhar
 force (3) = sum(4)
 force (4) = sum(5)
 force (5) = sum(6)
 force (6) = sum(7)
 force (7) = sum(8)
 
!     COMPUTE STRESSES AT Z1 AND Z2
 
 istres(1) = elemid
 istres(2) = nelhar
 
 DO  i = 1,2
   zoveri = 0.0
   IF (iii /= 0.0) zoveri = zoff(i)/iii
   
   DO  j = 1,3
     sig(j) = sum(j) + sum(j+3)*zoveri
   END DO
   
   ipt = 8*i - 6
   stress(ipt+1) = zoff(i)
   stress(ipt+2) = sig(1)
   stress(ipt+3) = sig(2)
   stress(ipt+4) = sig(3)
   istres(ipt+5) = 1
   istres(ipt+6) = 1
   istres(ipt+7) = 1
   istres(ipt+8) = 1
 END DO
 
 RETURN
END SUBROUTINE scone2
