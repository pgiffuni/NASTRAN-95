SUBROUTINE push (in,bcd,icol,nchar,flag)
     
!     THIS ROUTINE IS USED TO PLACE BCD CHARACTERS OR INTEGERS FROM II
!     ARRAY INTO THE BCD STRING . IF FLAG = 1 AN INTEGER IS INPUT
 
 
 INTEGER, INTENT(IN)                      :: in(1)
 INTEGER, INTENT(OUT)                     :: bcd(1)
 INTEGER, INTENT(IN)                      :: icol
 INTEGER, INTENT(IN)                      :: nchar
 INTEGER, INTENT(IN OUT)                  :: flag
 EXTERNAL        orf
 LOGICAL :: first
 INTEGER :: orf, cperwd, ii(18),BLANK, digit,numbs(10)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ isys,iout,nogo,idum(35),nbpc,nbpw,ncpw
 DATA    numbs / 1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9 /
 DATA    cperwd/ 4 /, first / .true. /, BLANK /1H  /, minus / 4H-   /
 
 IF (.NOT.first) GO TO 15
 first = .false.
 
!     REMOVE BLANKS FROM NUMBERS, AND ZERO FILL
 
 ish   = ncpw - 1
 DO  i = 1,10
   isave = krshft(numbs(i),ish)
   numbs(i) = klshft(isave,ish)
 END DO
 isave = krshft(minus,ish)
 minus = klshft(isave,ish)
 nx    = ncpw - cperwd
 ixtra = nx*nbpc
 ibl   = 0
 IF (nx == 0) GO TO 15
 ib1   = krshft(BLANK,ish)
 DO  i = 1,nx
   ibl   = orf(ibl,klshft(ib1,i-1))
 END DO
 
 15 IF (nchar <= 0) RETURN
 IF (nchar+icol > 128) GO TO 70
 nin   = (nchar-1)/cperwd + 1
 DO  i = 1,nin
   ii(i) = in(i)
 END DO
 IF (flag /= 1) GO TO 50
 
!     INTEGER HAS BEEN INPUT - 1 WORD ONLY
 
!     FIND POWER OF 10 = NUMBER OF CHARACTERS
 
 ix    = IABS(in(1))
 DO  i = 1,8
   ix    = ix/10
   IF (ix == 0) GO TO 30
 END DO
 GO TO 80
 30 ic    = i
 IF (in(1) <  0) ic = ic + 1
 IF (ic > nchar) GO TO 80
 ii(2) = BLANK
 ix    = IABS(in(1))
 IF (ic <= cperwd) GO TO 40
 
!     NUMBER TAKES TWO WORDS
 
 m     = ic - cperwd
 ii(2) = krshft(BLANK,m)
 DO  i = 1,m
   ij    = ix/10
   digit = IABS(ix-10*ij) + 1
   ix    = ij
   iadd  = numbs(digit)
   ii(2) = orf(ii(2),krshft(iadd,m-i))
 END DO
 
 ic    = cperwd
 
!     FIRST WORD SET HERE FOR BOTH CASES
 
 40 ii(1) = krshft(BLANK,ic)
 DO  i = 1,ic
   IF (i == ic .AND. in(1) < 0) CYCLE
   ij    = ix/10
   digit = IABS(ix-10*ij) + 1
   ix    = ij
   iadd  = numbs(digit)
   ii(1) = orf(ii(1),krshft(iadd,ic-i))
 END DO
 IF (in(1) < 0) ii(1) = orf(ii(1),minus)
 
 50 iwrd  = (icol-1)/cperwd + 1
 icl   = icol - (iwrd-1)*cperwd
 lwrd  = (icol+nchar-2)/cperwd + 1
 lcol  = icol + nchar - (lwrd-1)*cperwd - 1
 m1    = (icl-1)*nbpc
 m2    = cperwd*nbpc - m1
 m3    = m2 + (ncpw-cperwd)*nbpc
 
!     M1 IS THE NUMBER OF BITS FOR THE  FIRST SET OF CHARACTERS
!     M2 IS THE NUMBER OF BITS FOR THE SECOND SET OF CHARACTERS
!     M3 IS THE NUMBER OF BITS FOR THE RIGHT HALF OF THE WORD
 
!     IADD IS THE CURRENT WORKING WORD, IADD1 IS THE SPILL
 
 isave = krshft(bcd(iwrd),m3/nbpc)
 iadd1 = klshft(isave,m3/nbpc)
 k = 0
 DO  i = iwrd,lwrd
   k = k + 1
   
!     SPLIT INPUT WORD INTO TWO SETS
   
!     MOVE LEFT HALF TO RIGHT SIDE OF IADD AND ADD IADD1
   
   isave = krshft(ii(k),(m1+ixtra)/nbpc)
   iadd  = orf(klshft(isave,ixtra/nbpc),iadd1)
   
!     IF THIS ISNT THE LAST WORD MOVE THE RIGHT HALF TO IADD1 AND INSERT
   
   IF (i >= lwrd) CYCLE
   isave = krshft(ii(k),ixtra/nbpc)
   iadd1 = klshft(isave,m3/nbpc)
   
   bcd(i) = orf(iadd,ibl)
   
 END DO
 
!     LAST WORD PROCESSED HERE, REMOVE EXTRA CHARACTERS
 
 ish   = ncpw - lcol
 isave = krshft(iadd ,ish)
 iadd  = klshft(isave,ish)
 isave = klshft(bcd(lwrd),lcol)
 bcd(lwrd) = orf(iadd,krshft(isave,lcol))
 RETURN
 
 70 WRITE  (iout,75) ufm,nchar,in
 75 FORMAT (a23,' 6015. TOO MANY CHARACTERS TO BE INSERTED IN A DMAP',  &
     ' LINE', /6X,4H n = , i8 ,6X, 6HWORD =,a4)
 GO TO 90
 80 WRITE  (iout,85) ufm,in
 85 FORMAT (a23,' 6016. TOO MANY DIGITS TO BE INSERTED IN DMAP.',  &
     2X,'VALUE =',i12)
 
 90 nogo = 1
 RETURN
END SUBROUTINE push
