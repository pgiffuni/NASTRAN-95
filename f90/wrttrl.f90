SUBROUTINE wrttrl (filblk)
     
 
 INTEGER, INTENT(IN OUT)                  :: filblk(7)
 EXTERNAL        lshift,rshift,andf,orf
 INTEGER :: fiat,fist,NAME(2),orf,rshift,andf, filbk(7),lb(2)
 REAL :: words(4)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /machin/ mach
 COMMON /xfiat / fiat(3)
 COMMON /xfist / fist(2)
 COMMON /xsortx/ isav(6)
 COMMON /l15_l8/ l15,l8
 COMMON /system/ system(175)
 COMMON /logout/ lout
 EQUIVALENCE     (system(2),iout), (system(24),icfiat), (system(40),nbpw)
 DATA    mbit  / 0 /   ,words / 1.0, 2.0, 2.0, 4.0 /
 DATA    NAME  / 4HWRTT,4HRL  /
 DATA    mask  / 65535 /
 
 
!     IF ICFIAT= 8, WRTTRL WILL PACK SIX SIXTEEN BIT POSITIVE INTEGERS
!     INTO THREE THIRTY-TWO BIT WORDS AND STORE THEM IN THE FIAT
!     NO SUCH PACKING IF ICFIAT=11
 
 
!     SEARCH FIST FOR THE FILE
 
!     WRTTRL WILL NOT CHANGE TRAILER FOR 100 SERIES FILES
 
 IF (filblk(1) > 100 .AND. filblk(1) < 199) CALL mesage (-40, filblk(1),NAME)
 
!     ONLY MATGEN, OPTION 10, SENDS FILE 199 OVER HERE
 
 IF (filblk(1) == 199) filblk(1) = 101
 
!     THIS 'NEXT TO SIGN' MBIT IS SET BY SDCOMP AND SDCMPS
 
 mbit = lshift(1,nbpw-2 - (nbpw-32))
 nout = iout
 nout = lout
 
!     VERIFY SQUARE AND SYMM. MATRICES
 
 IF (l8 == 0 .OR. l15 == 0) GO TO 20
 IF (filblk(7)  <  mbit) GO TO 20
 IF (filblk(4) /= 1 .AND. filblk(4) /= 6) GO TO 20
 IF (filblk(2) == filblk(3)) GO TO 20
 CALL fname (filblk(1),lb(1))
 WRITE  (iout,10) swm,lb(1),lb(2),filblk(2),filblk(3),filblk(4)
 10 FORMAT (a27,', DATA BLOCK ',2A4,1H,,i9,3H by,i8,', IS MIS-LABLED',  &
     ' SQUARE OR SYMM.  (FORM=',i3,1H), /5X,  &
     'TURN DIAGS 1, 8 AND 15 ON FOR ERROR TRACEBACK')
 CALL sswtch (1,n)
!WKBD IF (N .NE. 0) CALL ERRTRC ('WRTTRL  ',10)
 
 20 CONTINUE
 n = fist(2)*2 + 1
 DO  i = 3,n,2
   IF (fist(i) /= filblk(1)) CYCLE
   INDEX = fist(i+1) + 1
   GO TO 40
 END DO
 CALL mesage (-11,filblk(1),NAME)
 
!     IF (1) BIT 'NEXT TO SIGN BIT' IS ON IN FILBLK(7), (2) FILBLK(2)
!     AND FILBLK(3), WHICH ARE COLUMN AND ROW, ARE NON ZEROS, AND
!     FILBLK(5), WHICH IS TYPE, IS 1,2,3 OR 4, THE INCOMING TRAILER IS
!     A MATRIX TRAILER. IN THIS CASE FILBLK(7) IS CONVERTED TO A DENSITY
!     PERCENTAGE BEFORE STORING IN THE FIAT.
 
 40 IF (filblk(7) < mbit) GO TO 50
 count = filblk(7) - mbit
 i = filblk(5)
 IF (filblk(2) == 0 .OR. filblk(3) == 0 .OR. i < 1 .OR. i > 4) GO TO 50
 fn = filblk(2)
 fm = filblk(3)
 filblk(7) = (count/(fn*fm*words(i)))*1.e4 + 1.0E-3
 IF (filblk(7) == 0 .AND. filblk(6) /= 0) filblk(7) = 1
 50 CONTINUE
 
 IF (l8 == 0) GO TO 100
 WRITE  (nout,60,ERR=70) fiat(INDEX+1),fiat(INDEX+2), (filblk(i),i=2,7)
 60 FORMAT (' *** DIAG 8, MESSAGE -- TRAILER FOR DATA BLOCK ',2A4, 2H =,6I10)
 GO TO 100
 70 CALL sswtch (1,n)
 IF (n == 0) GO TO 100
 WRITE  (nout,80,ERR=90) (filblk(i),i=2,7)
!IBMR 6/93 80 FORMAT (3H  (,6O20,1H))
 80 FORMAT (3H  (,6I8,1H))
!WKBR   90 CALL ERRTRC ('WRTTRL  ',70)
 90 CONTINUE
 
!     IF ICFIAT IS 8, PACK THE TRAILER INFORMATION IN THE FIAT.
!     BEFORE PACKING MAKE SURE NUMBERS ARE POSITIVE AND .LE. 16 BITS.
 
!     IF ICFIAT IS 11, 6 TRAILER WORDS ARE STORED DIRECTLY INTO 4TH,
!     5TH, 6TH, 9TH, 10TH AND 11TH WORD OF A FIAT ENTRY
 
 100 IF (icfiat == 11) GO TO 120
 DO  i = 2,7
   filblk(i) = andf(mask,IABS(filblk(i)))
 END DO
 fiat(INDEX+ 3) = orf(filblk(3),lshift(filblk(2),16))
 fiat(INDEX+ 4) = orf(filblk(5),lshift(filblk(4),16))
 fiat(INDEX+ 5) = orf(filblk(7),lshift(filblk(6),16))
 GO TO 130
 120 fiat(INDEX+ 3) = filblk(2)
 fiat(INDEX+ 4) = filblk(3)
 fiat(INDEX+ 5) = filblk(4)
 fiat(INDEX+ 8) = filblk(5)
 fiat(INDEX+ 9) = filblk(6)
 fiat(INDEX+10) = filblk(7)
 130 IF (fiat(INDEX) >= 0) GO TO 150
 
!     FIND EQUIVALENCED FILES IN FIAT AND WRITE TRAILER ON THEM
 
 iucb  = andf(fiat(INDEX),mask)
 iendf = fiat(3)*icfiat - 2
 DO  i = 4,iendf,icfiat
   IF (fiat(i) >= 0) CYCLE
   
!     PICK UP UNIT CONTROL BLOCK
   
   itucb = andf(fiat(i),mask)
   IF (itucb /= iucb) CYCLE
   
!     FOUND FILE
   
   fiat(i+ 3) = fiat(INDEX+ 3)
   fiat(i+ 4) = fiat(INDEX+ 4)
   fiat(i+ 5) = fiat(INDEX+ 5)
   IF (icfiat == 8) CYCLE
   fiat(i+ 8) = fiat(INDEX+ 8)
   fiat(i+ 9) = fiat(INDEX+ 9)
   fiat(i+10) = fiat(INDEX+10)
 END DO
 
!     SAVE THE TRAILER IN ISAV IF FILE IS SCRATCH 1
!     (SAVED FOR GINOFILE MODULE, SUBROUTINE GINOFL)
 
 150 IF (filblk(1) /= 301) RETURN
 isav(1) = fiat(INDEX+ 3)
 isav(2) = fiat(INDEX+ 4)
 isav(3) = fiat(INDEX+ 5)
 IF (icfiat == 8) GO TO 160
 isav(4) = fiat(INDEX+ 8)
 isav(5) = fiat(INDEX+ 9)
 isav(6) = fiat(INDEX+10)
 160 RETURN
 
 
 ENTRY rdtrl (filbk)
!     ===================
 
!     RDTRL WILL UNPACK THE THREE WORDS STORED IN THE FIAT AND RETURN
!     THE SIX WORDS OF TRAILER INFORMATION
 
 
!     SEARCH THE FIST FOR THE FILE
 
 n = fist(2)*2 + 1
 DO  i = 3,n,2
   IF (fist(i) /= filbk(1)) CYCLE
   INDEX = fist(i+1) + 1
   GO TO 210
 END DO
 
!     FILE WAS NOT FOUND, SET THE FILE NAME NEGATIVE
 
 filbk(1) = -IABS(filbk(1))
 RETURN
 
!     CHECK FIAT ENTRY 8 OR 11 WORDS PER ENTRY
 
 210 IF (icfiat == 11) GO TO 220
 
!     8 WORD ENTRY, UNPACK THE TRAILER INFORMATION
 
 filbk(2) = rshift(fiat(INDEX+3),16)
 filbk(3) = andf(fiat(INDEX+3),mask)
 filbk(4) = rshift(fiat(INDEX+4),16)
 filbk(5) = andf(fiat(INDEX+4),mask)
 filbk(6) = rshift(fiat(INDEX+5),16)
 filbk(7) = andf(fiat(INDEX+5),mask)
 GO TO 230
 
!     11 WORD ENTRY, TRAILER NOT PACKED
 
 220 filbk(2) = fiat(INDEX+ 3)
 filbk(3) = fiat(INDEX+ 4)
 filbk(4) = fiat(INDEX+ 5)
 filbk(5) = fiat(INDEX+ 8)
 filbk(6) = fiat(INDEX+ 9)
 filbk(7) = fiat(INDEX+10)
 
 230 RETURN
END SUBROUTINE wrttrl
