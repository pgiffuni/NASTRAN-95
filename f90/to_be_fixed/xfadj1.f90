SUBROUTINE xfadj1 (bf,shift,sd)
     
!     XFADJ1 ADJUSTS 4 CHARACTER FIELDS, LEFT OR RIGHT, 2 OR 4 FIELDS
!     AT A TIME
 
!     BF    = ADDR OF LEFT MOST FIELD
!     SHIFT = LSHIFT OR RSHIFT
!     SD   = 0 SINGLE (2 FIELDS), 1 DOUBLE (4 FIELDS)
!     RIGHT SHIFTING CAUSES INSERTION OF LEADING ZEROS
 
 
 INTEGER, INTENT(IN OUT)                  :: bf(1)
 INTEGER, INTENT(IN)                      :: shift
 INTEGER, INTENT(IN OUT)                  :: sd
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf,shift
 LOGICAL :: dec
 DIMENSION  bk(6),mk(6),sft(3)
 COMMON /machin/ mach
 COMMON /xsrtcm/ bimsk1(6),bimsk2(5),bimsk3(4),bimsk4(4),bimsk5(2),  &
     bimsk6,bkmsk1(8),bkmsk2,shifts(4),icon1,icon2,  &
     star,plus,dollar,starl,slash,sftm,mask,BLANK,mka, is,mbit4
 EQUIVALENCE     (bk(1) ,bkmsk1(2)),(mk(1),bimsk1(1)),  &
     (sft(1),shifts(2)),(blks, bkmsk1(8)), (bkx   ,bkmsk1(1))
 
!     DATA     BK   / 4H0000,4H0000,4H0000,4H000 ,4H00  ,4H0   /
!     DATA     (MK(I),I=1,6) /O777777000000,O777700000000,O770000000000,
!    1                        O000000770000,O000077770000,O007777770000/
!     DATA     (SFT(I),I=1,3)/6,12,18/
!     DATA     BLKS / 4H    /,    BKX/4H0000/
 
 
!     INITIALIZE ROUTINES
 
 dec = mach == 5 .OR. mach == 6 .OR. mach == 21
 IF (shift(mk(3),sft(1)) /= 0) GO TO 10
 
!     LEFT SHIFT REQUESTED
 
 blk= blks
 i1 = 1
 i2 = 2
 i3 = 3
 i4 = 4
 j  = 3
 GO TO 30
 
!     RIGHT SHIFT REQUESTED
 
 10 blk= bkx
 j  = 4
 IF (sd == 0) GO TO 20
 
!     DOUBLE FIELD
 
 i1 = 4
 i2 = 3
 i3 = 2
 i4 = 1
 GO TO 30
 
!     SINGLE FIELD
 
 20 i1 = 2
 i2 = 1
 
!     TOTAL FIELD SHIFTS
 
 30 n = 0
 40 IF (j == 4 .AND. bf(i1) /= blks) GO TO 60
 IF (j == 3 .AND. bf(i1) /= blks .AND. bf(i1) /= bkx) GO TO 60
 bf(i1) = bf(i2)
 bf(i2) = blk
 IF (sd == 0) GO TO 50
 n = n + 1
 bf(i2) = bf(i3)
 bf(i3) = bf(i4)
 bf(i4) = blk
 IF (n /= 3) GO TO 40
 50 IF (bf(i1) == blks) RETURN
 
!     CHARACTER SHIFTS BETWEEN FIELDS
 
 60 n = 0
 IF (j == 3) GO TO 150
 
!     RIGHT
 
 ii = i1
 70 IF (bf(ii) /= blks) GO TO 80
 bf(ii) = bkx
 GO TO 110
 80 IF (bf(ii) == bkx) GO TO 110
 90 IF (.NOT.dec) ihld = rshift(andf(mk(3),bf(ii)),1)
 IF (     dec) ihld = khrfn4(rshift(khrfn4(khrfn1(bkmsk2,1, bf(ii),1)),1))
 IF (ihld /= icon1) GO TO 100
 n = n + 1
 IF (.NOT.dec) bf(ii) = lshift(bf(ii),sft(1))
 IF (     dec) bf(ii) = khrfn3( bkmsk2,bf(ii),1,1)
 IF (n < 3) GO TO 90
 GO TO 120
 100 IF (n /= 0) GO TO 120
 110 ii = ii - 1
 IF (ii == 0) GO TO 130
 GO TO 70
 120 n2 = 4 - n
 IF (.NOT.dec) bf(ii) = orf(rshift(bf(ii),sft(n)),bk(n2))
 IF (     dec) bf(ii) = khrfn3( bk(n2),bf(ii),n,0)
 n = 0
 GO TO 110
 130 n = 0
 
!     RIGHT
 
 140 IF (dec) GO TO 141
 IF (andf(mk(4),bf(i1)) /= andf(mk(4),blks)) GO TO 170
 GO TO 160
 141 IF (khrfn1(mk(4),4,bf(i1),4) /= khrfn1(mk(4),4,blks,4)) GO TO 170
 GO TO 160
 
!     LEFT
 
 150 IF (.NOT.dec) ihld = rshift(andf(mk(3),bf(i1)),1)
 IF (     dec) ihld = khrfn4(rshift(khrfn4(khrfn1(bkmsk2,1,bf(i1), 1)),1))
 IF (ihld /= icon1 .AND. ihld /= icon2) GO TO 170
 160 n = n + 1
 IF (.NOT.dec) bf(i1) = shift(bf(i1),sft(1))
 IF (     dec) bf(i1) = khrfn3( bkmsk2,bf(i1),1,4-j)
 IF (n >= 3) GO TO 180
 IF (j == 3) GO TO 150
 GO TO 140
 170 IF (n == 0) RETURN
 180 IF (j == 4) GO TO 190
 
!     LEFT SHIFTS
 
 n1 = n
 n2 = n + 3
 GO TO 200
 
!     RIGHT SHIFTS
 
 190 n1 = 7 - n
 n2 = 4 - n
 200 n3 = 4 - n
 IF (.NOT.dec) bf(i1) = orf(andf(mk(n1),bf(i1)),andf(mk(n2),  &
     isft(bf(i2),sft(n3),j)))
 IF (     dec) bf(i1) = khrfn3(bf(i1),bf(i2),n3,j-3)
 bf (i1) = orf(bf(i1),bkmsk2)
 IF (.NOT.dec) bf(i2) = orf(andf(mk(n1),shift(bf(i2),sft(n))), bk(n2))
 IF (     dec) bf(i2) = khrfn3( bk(n2),bf(i2),n,4-j)
 IF (sd == 0) RETURN
 
 IF (.NOT.dec) bf(i2) = orf(andf(mk(n1),bf(i2)),andf(mk(n2),  &
     isft(bf(i3),sft(n3),j)))
 IF (     dec) bf(i2) = khrfn3( bf(i2),bf(i3),n3,j-3 )
 bf(i2) = orf(bf(i2),bkmsk2)
 IF (bf(i2) == blk) RETURN
 
 IF (.NOT.dec) bf(i3) = orf(andf(mk(n1),shift(bf(i3),sft(n))), bk(n2))
 IF (     dec) bf(i3) = khrfn3(bk(n2),bf(i3),n,4-j)
 IF (.NOT.dec) bf(i3) = orf(andf(mk(n1),bf(i3)),andf(mk(n2),  &
     isft(bf(i4),sft(n3),j)))
 IF (     dec) bf(i3) = khrfn3(bf(i3),bf(i4),n3,j-3)
 bf(i3) = orf(bf(i3),bkmsk2)
 IF (bf(i3) == blk) RETURN
 
 IF (.NOT.dec) bf(i4) = orf(andf(mk(n1),shift(bf(i4),sft(n))), bk(n2))
 IF (     dec) bf(i4) = khrfn3(bk(n2),bf(i4),n,4-j)
 RETURN
END SUBROUTINE xfadj1
