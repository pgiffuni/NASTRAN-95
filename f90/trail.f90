SUBROUTINE trail
     
!     MODULE TO INTERROGATE OR ALTER ANY VALUE OF A 6 WORD MATRIX
!     OR TABLE TRAILER
 
!     DMAP CALL
 
!        TRAILER  DB / /*OPT*/WORD/S,N,VALUE $
 
!     INPUT DATA BLOCKS
 
!        DB - DATA BLOCK FOR WHICH TRAILER IS TO BE ALTERED OR READ
 
!     PARAMETERS
 
!        OPT   - BCD,INPUT.
!                RETURN - VALUE OF SPECIFIED TRAILER WORD IS TO
!                         BE RETURNED
!                STORE  - VALUE OF SPECIFIED TRAILER WORD IS TO
!                         CHANGED
!        WORD  - INTEGER,INPUT. DESIRED WORD OF TRAILER
!        VALUE - INTEGER,INPUT OR OUTPUT. LOCATION WHERE VALUED WILL
!                RETURNED OR FROM WHICH REPLACEMENT VALUE WILL BE
!                TAKEN.  RETURNED NEGATIVE IF DB IS PURGED.
 
!     FOR MATRIX DATA BLOCKS, THE TRAILER POSITIONS ARE AS FOLLOWS
 
!        WORD 1 - NUMBER OF COLUMNS
!        WORD 2 - MUNBER OF ROWS
!        WORD 3 - MATRIX FORM
!        WORD 4 - TYPE OF ELEMENTS
!        WORD 5 - MAXIMUM NUMBER OF NON-ZERO WORDS IN ANY ONE COLUMN
!        WORD 6 - MATRIX DENSITY * 100
 
 EXTERNAL        lshift  ,andf    ,orf
 INTEGER :: db      ,opt     ,word     ,value    ,store(2) ,  &
     mcb(7)  ,fiat    ,fist     ,RETURN(2),modnam(2), orf     ,andf
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /xfiat / fiat(2)
 COMMON /xfist / fist(2)
 COMMON /BLANK / opt(2)  ,word    ,value
 COMMON /system/ ibuf    ,nout    ,dum21(21),icfiat
 DATA    store / 4HSTOR  ,4HE    /
 DATA    RETURN/ 4HRETU  ,4HRN   /
 DATA    modnam/ 4HTRAI  ,4HLER  /
 
!     GET TRAILER
 
 db = 101
 mcb(1) = db
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) GO TO 70
 
!     TEST ILLEGAL PARAMETER VALUES AND BRANCH ON OPT
 
 IF (word < 1 .OR. word > 6) GO TO 100
 IF (opt(1) == RETURN(1) .AND. opt(2) == RETURN(2)) GO TO 10
 IF (opt(1) == store(1)  .AND. opt(2) == store(2) ) GO TO 20
 GO TO 300
 
!     RETURN OPTION
 
 10 value = mcb(word+1)
 RETURN
 
!     STORE OPTION
 
!     SEARCH FIST FOR THE FILE
 
 20 n = fist(2)*2 + 1
 DO  i = 3,n,2
   IF (fist(i) /= db) CYCLE
   INDEX = fist(i+1) + 1
   GO TO 40
 END DO
 GO TO 70
 
!     PACK THE TRAILER INFORMATION INTO THE REQUESTED WORD.
!     MAKE SURE THE NUMBER IS POSITIVE AND .LE. 16 BITS IF ICFIAT=8
 
 40 IF (value  <  0) GO TO 200
 IF (icfiat == 11) GO TO 60
 IF (value > 65535) GO TO 200
 iw = (word+1)/2 + 2
 IF (word == (word/2*2)) GO TO 50
 
!     WORD IS ODD
 
 mask = 65535
 fiat(INDEX+iw) = orf(andf(fiat(INDEX+iw),mask),lshift(value,16))
 RETURN
 
!     WORD IS EVEN
 
 50 mask = lshift(65535,16)
 fiat(INDEX+iw) = orf(andf(fiat(INDEX+iw),mask),value)
 RETURN
 
!     ICFIAT = 11, TRAILER WORDS ARE NOT PACKED
 
 60 iw = 2
 IF (word >= 4) iw = 4
 fiat(INDEX+iw+word) = value
 RETURN
 
!     PURGED DATA BLOCK
 
 70 value = -1
 RETURN
 
!     ERROR CONDITIONS
 
 100 WRITE  (nout,110) ufm,word
 110 FORMAT (a23,' 2202.  PARAMETER, WORD, HAS ILLEGAL VALUE OF',i9)
 GO TO 500
 
 200 WRITE  (6,210) ufm,value
 210 FORMAT (a23,' 2202.  PARAMETER, VALUE, HAS ILLEGAL VALUE OF',i9)
 GO TO 500
 
 300 WRITE  (nout,310) ufm,opt
 310 FORMAT (a23,' 2202.  PARAMETER, OPT, HAS ILLEGAL VALUE OF ',2A4)
 
 500 CALL mesage (-37,0,modnam)
 RETURN
END SUBROUTINE trail
