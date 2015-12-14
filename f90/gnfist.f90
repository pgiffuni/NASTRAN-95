SUBROUTINE gnfist (filenm,fistnm,modno)
     
 
 INTEGER, INTENT(IN)                      :: filenm(2)
 INTEGER, INTENT(IN)                      :: fistnm
 INTEGER, INTENT(OUT)                     :: modno
 EXTERNAL        andf
 INTEGER :: andf,fiat, fist, fistx,oscar
 COMMON /xfist / fist(2)
 COMMON /xfiat / fiat(3)
 COMMON /xdpl  / idpl(3)
 COMMON /oscent/ oscar(7)
 COMMON /ipurge/ ipval(5)
 COMMON /isosgn/ isval(34)
 COMMON /ixsfa / ixval(5)
 COMMON /system/ skip(23),icfiat
 DATA    mask1 / 65535 /,  mask  / 32767 /
!             MASK1 = O177777   MASK  = O77777
 
 DO  k = 1,5
   ipval(k) = 0
   ixval(k) = 0
 END DO
 DO  k = 5,34
   isval(k) = 0
 END DO
 
 isval(1) = 3
 isval(2) = 3
 isval(3) = 1
 isval(4) = 2
 
 ixval(3) = 10
 
 IF (filenm(1) == 0 .AND. filenm(2) == 0) RETURN
 
!     SEARCH FIAT FOR MATCHING FILE
 
 lfiat = fiat(3)
 k = 5
 DO  j = 1,lfiat
   IF (filenm(1) == fiat(k) .AND. filenm(2) == fiat(k+1)) GO TO 30
   k = k + icfiat
 END DO
 
!     FILE NOT IN FIAT - IF INPUT FILE ASSUME PURGED
 
 IF (fistnm > 100 .AND. fistnm < 200) GO TO 40
 
!     MUST CALL IN FILE ALLOCATOR
 
 20 CALL xsfa (modno)
 modno = -modno
 RETURN
 
!     IF FILE POINTER = 77777 NO ENTRY IS MADE IN FIST
 
 30 IF (andf(fiat(k-1),mask) == mask) RETURN
 IF (fistnm <= 100 .OR. fistnm >= 300) GO TO 170
 IF (fistnm >= 200) GO TO 120
 
 
!     INPUT FILE
!     ==========
 
!     SEE IF IT EXISTS
 
 IF (fiat(k+2) /= 0 .OR. fiat(k+3) /= 0 .OR. fiat(k+4) /= 0) GO TO 170
 IF (icfiat == 11 .AND. (fiat(k+7) /= 0 .OR. fiat(k+8) /= 0 .OR.  &
     fiat(k+9) /= 0)) GO TO 170
 
!     INPUT FILE NOT GENERATED ACCORDING TO FIAT - CHECK DPL
 
 40 i1 = oscar(7)*3 + 5
 j1 = idpl(3) *3 + 1
 l  = fiat(3) *icfiat - 2
 DO  j = 4,j1,3
   IF (idpl(j) == filenm(1) .AND. idpl(j+1) == filenm(2)) GO TO 60
 END DO
 RETURN
 
!     FILE IN DPL - ZERO FIAT ENTRY SO FILE ALLOCATOR WILL UNPOOL IT.
!     DO THIS FOR OTHER LIKE I/P FILES IN OSCAR ENTRY.
 
 60 DO  i = 8,i1,3
   IF (oscar(i) == 0) CYCLE
   
!     SEARCH FIAT
   
   DO  k = 4,l,icfiat
     IF (oscar(i) == fiat(k+1) .AND. oscar(i+1) == fiat(k+2)) GO TO 80
   END DO
   
!     FILE NOT IN FIAT - CHECK NEXT INPUT FILE
   
   CYCLE
   
!     FILE IN FIAT - CHECK DPL IF FIAT TRAILER IS ZERO
   
   80 IF (fiat(k+3) /= 0 .OR. fiat(k+4) /= 0 .OR. fiat(k+5) /= 0 .OR.  &
       andf(mask,fiat(k)) == mask) CYCLE
   IF (icfiat == 11 .AND. (fiat(k+8) /= 0 .OR. fiat(k+9) /= 0 .OR.  &
       fiat(k+10) /= 0)) CYCLE
   DO  j = 4,j1,3
     IF (idpl(j) == fiat(k+1) .AND. idpl(j+1) == fiat(k+2)) GO TO 100
   END DO
   CYCLE
   
!     FILE IS IN DPL - ZERO OUT FIAT ENTRY
   
   100 fiat(k) = andf(mask1,fiat(k))
   IF (andf(mask,fiat(k)) == mask) fiat(k) = 0
   fiat(k+1) = 0
   fiat(k+2) = 0
 END DO
 
!     CALL FILE ALLOCATOR AND UNPOOL FILES
 
 GO TO 20
 
 
!     OUTPUT FILE
!     ===========
 
!     SEARCH DPL FOR FILE NAME
 
 120 j1 = idpl(3)*3 + 1
 DO  m = 4,j1,3
   IF (idpl(m) == filenm(1) .AND. idpl(m+1) == filenm(2)) GO TO 140
 END DO
 GO TO 170
 
!     FILE NAME IS IN DPL - PURGE IT AND ALL EQUIV FILE FROM DPL
 
 140 idpl(m  ) = 0
 idpl(m+1) = 0
 l = idpl(m+2)
 DO  j = 4,j1,3
   IF (j == m .OR. l /= idpl(j+2)) CYCLE
   idpl(j  ) = 0
   idpl(j+1) = 0
   idpl(j+2) = 0
 END DO
 
!     IF THIS IS LAST FILE ON POOL TAPE, DECREASE FILE COUNT IN DPL
 
 IF (andf(l,mask) /= idpl(1)-1) GO TO 160
 idpl(  1) = idpl(1) - 1
 idpl(m+2) = 0
 
!     IF DELETED FILES ARE AT END OF DPL, DECREMENT ENTRY COUNT
 
 160 IF (idpl(j1) /= 0 .OR. idpl(j1+1) /= 0 .OR. idpl(j1+2) /= 0) GO TO 170
 idpl(3) = idpl(3) - 1
 j1 = idpl(3)*3 + 1
 GO TO 160
 
!     CHECK FOR FIST TABLE OVERFLOW
 
 170 IF (fist(1) <= fist(2)) CALL mesage (-20,IABS(modno),filenm)
 fist(2) = fist(2)   + 1
 fistx   = fist(2)*2 + 1
 fist(fistx  ) = fistnm
 fist(fistx+1) = k - 2
 IF (fistnm < 300) RETURN
 
!     ZERO TRAILER FOR SCRATCH FILE
 
 fiat(k+2) = 0
 fiat(k+3) = 0
 fiat(k+4) = 0
 IF (icfiat == 8) GO TO 180
 fiat(k+7) = 0
 fiat(k+8) = 0
 fiat(k+9) = 0
 180 RETURN
END SUBROUTINE gnfist
