SUBROUTINE comect (ele,MAX)
     
!     REVISED  10/1990 BY G.CHAN/UNISYS
!              TO INCLUDE OFFSET DATA FOR CBAR, CTRIA3 AND CQUAD4 IN
!              THE ECT2 DATA BLOCK
!              (6 COORDINATE VALUES FOR THE BAR, AND 1 OFFSET VALUE
!              FOR EACH OF THE TWO PLATES, ARE ADDED AFTER THE GRID
!              DATA)
 
 
 INTEGER, INTENT(OUT)                     :: ele(1)
 INTEGER, INTENT(IN)                      :: MAX
 INTEGER :: idrec(3), elid(2),TYPE,ihx2(20),ihx3(32),  &
     ect1,ect2,bufsiz,b1,b2,gp(32),outrew,rew,m1(18),  &
     NAME(2),ERR(5),ept,pid,ix(1),pcomp(12)
 REAL :: offset(1)
 COMMON /BLANK / skp1(12),ect1,skp2(7),merr,skp3(10),ect2
 COMMON /system/ bufsiz
 COMMON /zzzzzz/ x(1)
 COMMON /gpta1 / nel,last,incr,NE(1)
 EQUIVALENCE     (offset(1),gp(1)), (ix(1),x(1))
 DATA    NAME  / 4H com,4HECT /,  outrew,rew,inrew / 1, 1, 0 /
 DATA    pcomp / 5502,25,2, 5602,14,2, 5702,13,2, 5802,17,17 /
!                     PCOMP      PCOMP1     PCOMP2     PSHELL
 DATA    nm1   / 18    /,  &
     m1    / 4H(33X, 4H,2A4, 4H,18H, 4HIGNO, 4HRING, 4H ele,  &
     4HMENT, 4H (2A, 4H4,32, 4HH) w, 4HITH , 4HMORE,  &
     4H tha, 4HN 32, 4H con, 4HNECT, 4HIONS, 4H.)  /, ilxx  / 2HXX /
 DATA    ihx2  / 1,1,3,3,5,5,7,7,1,3,5,7,13,13,15,15,17,17,19,19/
 DATA    ihx3  / 1,1,4,4,4,7,7,7,10,10,10,1,1,4,7,10,21,24,27,30,  &
     21,21,24,24,24,27,27,27,30,30,30,21            /
 
 b1 = korsz(x) - (3*bufsiz+2)
 b2 = b1 + bufsiz + 3
 ERR(1) = 4
 ERR(2) = NAME(1)
 ERR(3) = NAME(2)
 
!     IF EPT FILE IS PRESENT, AND ANY OF THE PSHELL, PCOMP, PCOMP1 AND
!     PCOMP2 CARDS IS ALSO PRESENT, CREATE A TABLE OF PROPERTY ID AND
!     OFFSET DATA, TO BE USE LATER BY CTRIA3 OR CQURD4 ELEMENTS
 
 jcomp = b1
 ept = 104
 CALL OPEN (*40,ept,x(b1),inrew)
 CALL READ (*30,*30,ept,ix,2,1,m)
 CALL CLOSE (ept,rew)
 CALL preloc (*40,x(b1),ept)
 n = 1
 DO  i = 1,12,3
   idrec(1) = pcomp(i)
   idrec(2) = idrec(1)/100
   CALL locate (*20,x(b1),idrec,j)
   k = pcomp(i+1)
   j = pcomp(i+2)
   10 CALL READ (*20,*20,ept,x,k,0,m)
   IF (x(j) == 0.0) GO TO 10
   jcomp = jcomp - 2
   ix(jcomp ) = ix(1)
   x(jcomp+1) =  x(j)
   GO TO 10
   n = n + 1
 END DO
 30 CALL CLOSE (ept,rew)
 kcomp = b1 - 1
 
!     CONSTRUCT A LIST OF INDICES IN THE ECT FOR USE WITH GPECT IN THE
!     PLOT MODULE BY CONTOUR PLOTTING
 
 40 CALL gopen (ect1,x(b1),inrew)
 DO  j = 1,MAX
   ele(j) = 0
 END DO
 i = 1
 60 CALL READ (*130,*80,ect1,idrec,3,0,m)
 DO  j = 1,nel
   idx = (j-1)*incr
   IF (NE(idx+4) == idrec(1)) GO TO 100
 END DO
 CALL skprec (ect1,1)
 GO TO 60
 80 CALL mesage (-3,ect1,NAME)
 90 CALL mesage (-2,ect1,NAME)
 100 lect = NE(idx+6) - 1
 110 CALL READ (*90,*60,ect1,ele(i),1,0,m)
 CALL fread (ect1,0,-lect,0)
 i = i + 1
 IF (i > MAX) CALL mesage (-8,0,NAME)
 GO TO 110
 
 120 CALL mesage (-1,ect1,NAME)
 
 130 lele = i - 1
 CALL CLOSE (ect1,rew)
 
 CALL preloc (*120,x(b1),ect1)
 CALL gopen (ect2,x(b2),outrew)
 DO  n = 1,nel
   idx = (n-1)*incr
   
!     IF SCALAR CONNECTION POSSIBLE FOR ELEMENT THEN SKIP IT
   
   IF (NE(idx+11) /= 0) CYCLE
   
!     SKIP DUMMY ELEMENTS AND POINT ELEMENTS
   
   IF (NE(idx+10)-1  <= 0) CYCLE
   IF (NE(idx+16) == ilxx) CYCLE
   CALL locate (*290,x(b1),NE(idx+4),i)
   ngpel = NE(idx+10)
   IF (ngpel > 32) GO TO 270
   
   CALL WRITE (ect2,n,1,0)
   CALL WRITE (ect2,ngpel,1,0)
   140 CALL READ (*280,*280,ect1,elid,1,0,i)
   
!     FIND THIS ELEMENTS POINTER IN THE ECT
   
   DO  i = 1,lele
     IF (ele(i) == elid(1)) GO TO 160
   END DO
   CALL mesage (-37,0,NAME)
   160 elid(2) = i
   
!     DETERMINE NUMBER ENTRIES FOR SKIPPING TO GRID ENTRIES
   
   i = NE(idx+13) - 2
   IF (n == 52) GO TO 190
!              CHBDY
   
   IF (n == 64 .OR. n == 83) GO TO 240
!           CQUAD4        CTRIA3
   IF (i < 0) THEN
     GO TO   120
   ELSE IF (i == 0) THEN
     GO TO   180
   END IF
   
   170 CALL fread (ect1,0,-i,0)
   180 CALL fread (ect1,gp,ngpel,0)
   IF (n == 34) GO TO 230
!               CBAR
   
   CALL fread (ect1,0,-(NE(idx+6)-ngpel-i-1),0)
   GO TO 200
   
!     SPECIAL HANDLING FOR CHBDY
!     IF TYPE IS NEGATIVE, SAVE TYPE FLAG AFTER GRIDS.
   
   190 CALL fread (ect1,0,-1,0)
   CALL fread (ect1,TYPE,1,0)
   CALL fread (ect1,gp,8,0)
   CALL fread (ect1,0,-(NE(idx+6)-ngpel-i-1),0)
   IF (TYPE < 0) GO TO 140
   IF (TYPE == 6) TYPE = 3
   gp(9) = TYPE
   GO TO 220
   
!     SPCIAL HANDLING OF IHEX2 AND IHEX3 WITH ZERO GRIDS
   
   200 IF (n /= 66 .AND. n /= 67) GO TO 220
   DO  j = 1,ngpel
     IF (gp(j) /= 0) CYCLE
     k = ihx3(j)
     IF (n == 66) k = ihx2(j)
     gp(j) = gp(k)
   END DO
   
   220 CALL WRITE (ect2,elid,2,0)
   CALL WRITE (ect2,gp,ngpel,0)
   GO TO 140
   
!     SPECIAL HANDLING OF THOSE ELEMENTS HAVING GRID OFFSET.
!     ADD THESE OFFSET DATA AFTER THE GRID POINTS
   
!     (1) CBAR ELEMENT, 2 OFFSET VECTORS (6 VALUES)
   
   230 CALL WRITE (ect2,elid,2,0)
   CALL WRITE (ect2,gp,ngpel,0)
   CALL fread (ect1,0,-6,0)
   CALL fread (ect1,offset,6,0)
   CALL WRITE (ect2,offset,6,0)
   GO TO 140
   
!     (2) CTRIA3 AND CQUAD4 ELEMENTS, ONE OFFSET DATA NORMAL TO PLATE.
!         OFFSET DATA COULD BE ON ELEMENT CARD OR ON PSHELL OR PCOMPI
!         CARDS
   
   240 CALL fread (ect1,pid,1,0)
   CALL fread (ect1,gp,ngpel,0)
   CALL WRITE (ect2,elid,2,0)
   CALL WRITE (ect2,gp,ngpel,0)
   j = 5
   IF (n == 64) j = 6
   CALL fread (ect1,0,-j,0)
   CALL fread (ect1,offset,1,0)
   IF (offset(1) /= 0.0) GO TO 260
   IF (jcomp     ==  b1) GO TO 260
   DO  i = jcomp,kcomp,2
     IF (ix(i) /= pid) CYCLE
     offset(1) = x(i+1)
     EXIT
   END DO
   260 CALL WRITE (ect2,offset,1,0)
   GO TO 140
   
!     ELEMENT TYPE WITH MORE THAN 32 CONNECTIONS
   
   270 ERR(4) = NE(idx+1)
   ERR(5) = NE(idx+2)
   CALL wrtprt (merr,ERR,m1,nm1)
   CALL skprec (ect1,1)
   CYCLE
   
   280 CALL WRITE (ect2,0,0,1)
 END DO
 
 CALL clstab (ect2,rew)
 CALL CLOSE  (ect1,rew)
 RETURN
END SUBROUTINE comect
