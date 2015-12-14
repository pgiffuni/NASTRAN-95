SUBROUTINE pidck (pfile,geom2,nopid,z)
     
!     THIS ROUTINE CHECKS THE UNIQUNESS OF PROPERTY IDS FOR ALL ELEMENTS
!     THAT HAVE PID FIELDS
 
!     IT SHOULD BE CALLED ONLY ONCE BY IFP
!     IT DOES NOT OPEN NOR CLOSE ANY GINO FILE.
 
!     DESIGN REQUIREMENT -
 
!     IF PID IS REFERENCED BY AN ELEMENT, THE PID MUST RESIDE ON THE
!     THIRD FIELD OF THE ELEMENT INPUT CARD.
!     INPUT FILES - GEOM2 AND PROPERTY FILE (EPT).
 
!     THIS VERSION INCLUDES SPECIAL HANDLING OF THE CQUAD4 AND CTRIA3
!     ELEMENTS WHICH USE AND SHARE MORE THAN ONE STANDARD PROPERTY CARD.
!     THE PROPERTY TYPE IDS OF THE PSHELL, PCOMP, PCOMP1 AND PCOMP2
!     MUST NOT BE INTERRUPTED BY ANOTHER PROPERTY TYPE. (I.E. NO OTHER
!     PROPERTY TYPE SHOULD HAVE AN ID PLACED IN BETWEEN 5502 THRU 5802).
!     NOTICE THAT THE PSHELL CARD HAS FIXED LENGTH WHILE THE 3 PCOMPI
!     CARDS HAVE VARIABLE LENGTH.
 
!     WRITTEN BY G.CHAN/UNISYS, SEPT. 1983
 
 
 INTEGER, INTENT(IN OUT)                  :: pfile
 INTEGER, INTENT(IN OUT)                  :: geom2
 INTEGER, INTENT(OUT)                     :: nopid
 INTEGER, INTENT(OUT)                     :: z(1)
 LOGICAL :: abort
 INTEGER :: NAME(2), flag,     x(3),     quad4,    pshell,  &
     pcomp(3)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,      uwm,      uim
 COMMON /system/ ibuf,     nout,     abort,    skip(42), kdum(9)
 COMMON /gpta1 / nelem,    last,     incr,     NE(1)
 DATA    quad4 , pshell,   pcomp                      /  &
     5408 , 5802,     5502,     5602,     5702   /
 DATA    NAME  / 4HPIDC,   4HK       /
 
!     UPDATE /GPTA1/ IF DUMMY ELEMENTS ARE PRESENT
 
 DO  i = 1,9
   IF (kdum(i) == 0) CYCLE
   k  = kdum(i)
   ng = k/10000000
   nc = (k-ng*10000000)/10000
   np = (k-ng*10000000 - nc*10000)/10
   k  = (51+i)*incr
   NE(k+ 6) = 2 + ng + nc
   NE(k+ 9) = 2 + np
   NE(k+10) = ng
 END DO
 
!     CREATE A PROPERTY ID TABLE IN Z FROM /GPTA1/ DATA BLOCK FOR THOSE
!     ELEMENTS THAT HAVE PROPERTY CARDS
!     4 WORDS PER ENTRY
!       WORD 1, PROPERTY TYPE CODE  (EPT-ID)
!       WORD 2, LENGTH OF PROPERTY CARD  (EPTWDS)
!       WORD 3, ELEMENT TYPE CODE   (ECT-ID)
!       WORD 4, LENGTH OF ELEMENT CARD (ECTWDS), PLUS POINTER TO GPTA1
 
 ii = 0
 DO  i = 1,last,incr
   IF (NE(i+6) == 0) CYCLE
   z(ii+1) = NE(i+6)
   z(ii+2) =-NE(i+8)
   z(ii+3) = NE(i+3)
   z(ii+4) = NE(i+5)*10000 + i
   ii = ii + 4
 END DO
 
!     ADD 3 MORE PROPERTY CARDS (PCOMP, PCOMP1, PCOMP2) FOR CQUAD4 (64)
!     AND CTRIA3
!     NOTE - THESE THREE ARE OPEN-ENDED, AND WE SET WORD 2 TO -8888
!          - WE GIVE THEM LOCALLY NEW QUAD4 IDS IN THE 3RD WORD, SO THAT
!            ELEMENT CQUAD4 AND ELEMENT CTRIA3 WILL PICK THEM UP VIA
!            THE PSEHLL DATA LATER.
 
 i = (64-1)*incr + 1
 IF (NE(i+3) /= quad4) CALL mesage (-37,0,NAME)
 DO  j = 1,3
   z(ii+1) = pcomp(j)
   z(ii+2) = -8888
   z(ii+3) = quad4 - j
   z(ii+4) = NE(i+5)*10000 + i
   ii = ii + 4
 END DO
 
!     SORT THIS 4-ENTRY Z-TABLE BY THE FIRST WORD.
!     SET WORD 2 TO -9999 IF ELEMENT USES THE SAME PROPERTY CARD AS THE
!     PREVIOUS ELEMENT.
 
 i4 = ii/4
 CALL sort (0,0,4,1,z,ii)
 DO  i = 5,ii,4
   IF (z(i) == z(i-4)) z(i+1) = -9999
 END DO
 
!     READ FROM PFILE ALL PID INTO REMAINING CORE. REPLACE WORD 2 IN THE
!     Z-TABLE BY PID BEGIN-ENDING POINTERS
 
 jj = ii + 1
 IF (nopid == 1) GO TO 210
 CALL REWIND (pfile)
 120  CALL fwdrec (*360,pfile)
 130  CALL READ (*190,*190,pfile,x,3,0,flag)
!     2147483647  = 2**31-1
 IF (x(1) == 2147483647) GO TO 190
 CALL bisloc (*120,x(1),z,4,i4,k)
 140  kp1 = k + 1
 IF (z(kp1) /= -9999) GO TO 150
 k = k - 4
 GO TO 140
 150  nwds = -z(kp1)
 IF (nwds <= 0) GO TO 120
 komp = 0
 IF (nwds /= 8888) GO TO 155
 komp = 1
 nwds = 8
 155  z(kp1) = (jj*10000) + (jj-1)
 jb = jj
 160  CALL READ (*360,*130,pfile,z(jj),nwds,0,flag)
 IF (komp == 0) GO TO 167
 165  CALL READ (*360,*130,pfile,j,1,0,flag)
 IF (j /= -1) GO TO 165
 167  je = MOD(z(kp1),10000)
 IF (je < jb) GO TO 180
 DO  j = jb,je
   IF (z(jj) == z(j)) GO TO 160
 END DO
 180  z(kp1) = z(kp1) + 1
 jj = jj + 1
 GO TO 160
 190  CALL REWIND (pfile)
 jj = jj - 1
 IF (jj <= ii) nopid = -1
 
!     RESET THE PSHELL POINTERS TO INCLUDE THE PCOMP GROUP IDS.
!     MAKE SURE THIS GROUP ARE ALL TOGETHER, NOT SEPERATED BY OTHER
!     PROPERTY CARD
!     THERE ARE 2 PSHELL CARDS, ONE FROM CQUAD4 AND ONE FROM CTRIA3,
!     MAKE SURE THE FIRST PSHELL POINTER IS USED
 
 CALL bisloc (*210,pshell,z,4,i4,kp1)
 IF (z(kp1+1) == -9999) kp1 = kp1 - 4
 IF (z(kp1- 4) /= pcomp(3) .OR. z(kp1-8) /= pcomp(2) .OR.  &
     z(kp1-12) /= pcomp(1)) GO TO 380
 j = z(kp1+1)
 IF (j <= 0) j = 0
 jb = j/10000
 je = MOD(j,10000)
 IF (jb == 0) jb = 9999999
 DO  i = 1,3
   CALL bisloc (*370,pcomp(i),z,4,i4,k)
   IF (z(k+1) <= 0) CYCLE
   j = z(k+1)/10000
   k = MOD(z(k+1),10000)
   IF (j < jb) jb = j
   IF (k > je) je = k
 END DO
 IF (jb /= 9999999) z(kp1+1) = (jb*10000) + je
 
!     RESET POINTERS FOR THOSE PROPERTY ID COMMON TO MORE THAN ONE TYPE
!     OF ELEMENTS, AND
!     MOVE THE THIRD ENTRY IN THE Z-TABLE TO FIRST, FOR ELEMENT SORT
 
 210  DO  i = 1,ii,4
   z(i) = z(i+2)
   j = i + 1
   IF (z(j) > 0) CYCLE
   IF (z(j) == -9999) z(j) = z(j-4)
 END DO
 CALL sort (0,0,4,1,z,ii)
 
!     READ IN CONNECTING ELEMENTS, ONE BY ONE, FROM GEOM2 FILE, AND
!     CHECK THE EXISTENCE OF THE PROPERTY ID IF IT IS SPECIFIED.
 
 kk = jj + 1
 CALL REWIND (geom2)
 230  CALL fwdrec (*360,geom2)
 240  CALL READ (*300,*300,geom2,x,3,0,flag)
 CALL bisloc (*230,x(1),z,4,i4,k)
 nwds = z(k+3)/10000
 IF (nwds <= 0) GO TO 230
 j = z(k+1)
 IF (j <= 0) GO TO 270
 jb = j/10000
 je = MOD(j,10000)
 250  CALL READ (*360,*240,geom2,z(kk),nwds,0,flag)
 jj1 = z(kk+1)
 DO  j = jb,je
   iz = IABS(z(j))
   IF (jj1 /= iz) CYCLE
   z(j) = -iz
   GO TO 250
 END DO
 CALL mesage (30,10,z(kk))
 abort = .true.
 GO TO 250
 270  j = MOD(z(k+3),10000)
 CALL mesage (30,11,NE(j))
 abort = .true.
 GO TO 230
 300  CALL REWIND (geom2)
 IF (abort .OR. nopid /= 0) GO TO 350
 
!     PREPARE AN ACTIVE PROPERTY ID LIST FOR SUBROUTINE MATCK
 
 j  = ii + 1
 ii = 1
 DO  i = j,jj
   IF (z(i) >= 0) GO TO 310
   ii = ii + 1
   z(ii) = -z(i)
   CYCLE
   310  z(kk) = z(i)
   kk = kk + 1
 END DO
 z(1) = ii
 
!     Z(2,...II) CONTAINS A LIST OF ACTIVE PROPERTY IDS, UN-SORTED,
!     REFERENCED BY ELEMENTS IN GEOM2 FILE.  Z(1) = LENGTH OF THIS LIST
 
 jj1 = jj + 1
 kk  = kk - 1
 IF (kk < jj1) RETURN
 WRITE  (nout,330) uim
 330  FORMAT (a29,', THE FOLLOWING PROPERTY IDS ARE PRESENT BUT NOT ',  &
     'USED -')
 WRITE  (nout,340) (z(j),j=jj1,kk)
 340  FORMAT (/5X,12I9)
 RETURN
 
!     SET Z(1) TO ZERO IF NO ACTIVE PROPERTY LIST EXISTS.
 
 350  z(1) = 0
 RETURN
 
 360  j = -2
 GO TO 400
 370  WRITE  (nout,375)
 375  FORMAT ('0*** CAN NOT LOCATE PSHELL OR PCOMP DATA IN /GPTA1/')
 GO TO 390
 380  WRITE  (nout,385) z(kp1),pshell,z(kp1-4),pcomp(3),  &
     z(kp1-8),pcomp(2),z(kp1-12),pcomp(1)
 385  FORMAT ('0*** ERROR IN /GPTA1/ PCOMP ARRANGEMENT',(/3X,2I7))
 390  j = -37
 400  CALL mesage (j,0,NAME)
 RETURN
END SUBROUTINE pidck
