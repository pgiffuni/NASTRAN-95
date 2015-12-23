SUBROUTINE linkup (*,NAME)
 
 IMPLICIT INTEGER (a-z)
 INTEGER, INTENT(IN)                      :: NAME(2)
 EXTERNAL        lshift,rshift,andf,orf,complf
 
 COMMON /machin/ machx
 COMMON /lnklst/ itop,ibot,isn,kind,itype,mask1,mask2,mask3
 COMMON /zzzzzz/ z(1)
 
!     HASH INTO TABLE
 
 SELECT CASE ( machx )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 10
   CASE (    3)
     GO TO 10
   CASE (    4)
     GO TO 20
   CASE (    5)
     GO TO 30
   CASE (    6)
     GO TO 30
   CASE (    7)
     GO TO 10
   CASE (    8)
     GO TO 30
   CASE (    9)
     GO TO 30
   CASE (   10)
     GO TO 30
   CASE (   11)
     GO TO    30
   CASE (   12)
     GO TO 40
   CASE (   13)
     GO TO 10
   CASE (   14)
     GO TO 40
   CASE (   15)
     GO TO 40
   CASE (   16)
     GO TO 30
   CASE (   17)
     GO TO 30
   CASE (   18)
     GO TO 30
   CASE (   19)
     GO TO 30
   CASE (   20)
     GO TO 30
   CASE (   21)
     GO TO    30
   CASE (   22)
     GO TO 30
 END SELECT
 
!     IBM AND UNIVAC
 
 10 itotal = NAME(1) + NAME(2)
 GO TO 50
 
!     60-BIT MACHINE
 
 20 itotal = rshift(NAME(1),18) + rshift(NAME(2),18)
 GO TO 50
 
!     32-BIT MACHINES
 
 30 itotal = rshift(NAME(1), 1) + rshift(NAME(2), 1)
 GO TO 50
 
!     64-BIT MACHINES
 
 40 itotal = rshift(NAME(1),32) + rshift(NAME(2),32)
 
 50 ihash = 4*IABS(MOD(itotal,250)) + 4
 k = andf(z(ihash),mask1)
 IF (k /= 0) GO TO 60
 
!     NO HASH CHAIN FOUND - CREATE CHAIN
 
 z(ihash) = z(ihash) + itop
 GO TO 90
 
!     HASH CHAIN FOUND - CHECK PRESENCE OF NAME
 
 60 IF (z(k) /= NAME(1) .OR. z(k+1) /= NAME(2)) GO TO 70
 ikind = rshift(andf(z(k+3),mask3),28)
 IF ((ikind+1)/2 == (kind+1)/2) GO TO 100
 70 l = andf(z(k+3),mask2)
 IF (l == 0) GO TO 80
 k = rshift(l,14)
 GO TO 60
 80 z(k+3) = z(k+3) + lshift(itop,14)
 
!     NO ENTRY FOUND - CREATE ENTRY
 
 90 z(itop  ) = NAME(1)
 z(itop+1) = NAME(2)
 z(itop+2) = lshift(itype,28)
 z(itop+3) = z(itop+3) + lshift(IABS(kind),28)
 itop = itop + 4
 IF (itop >= ibot) RETURN 1
 IF (kind <    0) RETURN
 k = itop - 4
 
!     ADD STATEMENT NUMBER TO LIST
 
 100 l = andf(z(k+2),mask1)
 IF (l /= 0) GO TO 110
 
!     LIST IS EMPTY - START LIST
 
 z(k+2) = z(k+2) + ibot
 GO TO 120
 
!     CHAIN ENTRY ON LIST
 
 110 l    = rshift(andf(z(k+2),mask2),14)
 z(l) = andf(z(l),complf(mask2))
 z(l) = orf(z(l),lshift(ibot,14))
 
!     ADD ENTRY TO LIST
 
 120 z(ibot)= orf(lshift(kind,28),isn)
 z(k+2) = andf(z(k+2),complf(mask2))
 z(k+2) = z(k+2) + lshift(ibot,14)
 ibot   = ibot - 1
 IF (itop >= ibot) RETURN 1
 RETURN
END SUBROUTINE linkup
