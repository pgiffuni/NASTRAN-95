SUBROUTINE xrgdtb (lu)
     
!     XRGDTB PROCESSES THE CARD AND FILE NAME RESTART TABLES
!     THIS SUBROUTINE IS CALLED ONLY BY XRGDFM
 
!     WRITTEN BY  RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!       LU            FORTRAN UNIT NUMBER FOR THE RIGID FORMAT FILE
!     /SYSTEM/
!       OPTAPE        OUTPUT UNIT NUMBER FOR THE PRINT FILE
!     /XRGDXX/
!       ICHAR         ARRAY IN 80A1 FORMAT CONTAINING CARD IMAGE
!       ISCR          FILE NUMBER ON WHICH TABLES ARE WRITTEN
!       ITYPE         TYPE OF TABLE BEING PROCESSED-('CARD'OR'FILE')
!       LIMIT         LOWER/UPPER LIMITS FOR VALUES IN THE TABLE
!       RECORD        CARD IMAGE IN 20A4 FORMAT
 
!     OUTPUT
!     /XRGDXX/
!       ICOL          COLUMN WITHIN CARD BEING PROCESSED
!       ICOUNT        NUMBER OF ALPHA CHARACTERS WITHIN A NAME
!       IERROR        ERROR FLAG - NON-ZERO IF ERROR OCCURRED
!       NAME          NAME OF THE SUBROUTINE
!       NUMBER        VALUE OF NUMBER RETURNED BY XRGNUM
 
!     LOCAL VARIABLES
!       ASTRSK          CONTAINS THE VALUE 1H*
!       BLANK           CONTAINS THE VALUE 1H
!       COMENT          CONTAINS THE VALUE OF 4H$$$$
!       DOLLAR          CONTAINS THE VALUE OF 1H$
!       ICOLUM          COLUMN NUMBER OF THE NEXT CHARACTER WITHIN
!                       A NAME
 
!     FUNCTIONS
!       1. CALLS READ AND XDCODE FOR EACH CARD WITHIN THE TABLE.
!       2. CALLS XRGNUM TO PROCESS ALL NUMBERS
!       3. CALLS XECODE TO PROCESS ALL NAMES
!       4. ALL ENTRIES READ ARE EXPECTED TO BE IN THE FOLLOWING
!          FORMAT:    NNNN    NAME  NAME  NAME  NAME  NAME ...
!          WHERE NNNN IS ANY NUMBER.
 
!     SUBROUTINES CALLED - XRGNUM,XECODE,READ,WRITE
 
!     ERRORS  MESSAGES 8028,8034,8029,8036 MAY BE ISSUED
 
 
 INTEGER, INTENT(IN OUT)                  :: lu
 INTEGER :: record, BLANK, dollar, astrsk, optape, coment
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ksystm(100)
 COMMON /xrgdxx/ irestr, nsubst, iphase, icol  , NUMBER, itype ,  &
     istate, ierror, num(2), ind   , nument        ,  &
     record(20)    , ICHAR(80)     , limit(2)      ,  &
     icount, idmap , iscr  , NAME(2), member(2)    , ignore
 EQUIVALENCE     (ksystm( 2),optape), (ksystm(39), nbpc),  &
     (ksystm(40), nbpw ), (ksystm(41), ncpw)
 DATA    BLANK / 1H /,  dollar / 1H$ /, astrsk / 1H* /
 DATA    coment/ 4H$$$$ /
 
 100  NUMBER  = 0
 NAME(1) = 0
 READ (lu,150,ERR=710,END=710) record
 150  FORMAT (20A4)
 CALL xdcode
 IF (record(1) == coment) GO TO 100
 IF (ICHAR(1) == dollar .AND. ICHAR(2) == astrsk) GO TO 800
 icol = 1
 200  IF (ICHAR(icol) == BLANK .OR. icol > 80) GO TO 500
 IF (NUMBER /= 0) GO TO 300
 CALL xrgnum
 IF (NUMBER == 0) GO TO 720
 IF (NUMBER >= limit(1) .AND. NUMBER <= limit(2)) GO TO 200
 GO TO 730
 300  icount = 1
 350  icolum = icol + icount
 IF (ICHAR(icolum) == BLANK .OR. icolum > 80) GO TO 400
 icount = icount + 1
 IF (icount <= 8) GO TO 350
 GO TO 740
 400  IF (icount == 0) GO TO 350
 CALL xecode
 CALL WRITE (iscr,NAME,2,0)
 CALL WRITE (iscr,NUMBER,1,0)
 icol = icol + icount
 GO TO 200
 500  IF (icol >= 80) GO TO 600
 icol = icol + 1
 GO TO 200
 600  IF (NUMBER == 0 .OR. NAME(1) == 0) GO TO 750
 GO TO 100
 
!     ERRORS
 
 710  WRITE  (optape,715) ufm,member
 715  FORMAT (a23,' 8027, UNEXPECTED EOF ENCOUNTERED ON FILE ',2A4,  &
     ' IN SUBROUTINE XRGDTB.')
 GO TO 770
 720  WRITE  (optape,725) ufm,record
 725  FORMAT (a23,' 8028, EXPECTED TO FIND AN INTEGER IN THE FIRST ',  &
     'FIELD OF THE FOLLOWING CARD', //20X,20A4)
 GO TO 760
 730  WRITE  (optape,735) ufm,NUMBER,record,limit,itype
 735  FORMAT (a23,' 8029, THE VALUE',i4,' GIVEN IN THE FIRST FIELD OF',  &
     ' THE FOLLOWING CARD', //20X,20A4, //5X,'IS OUTSIDE THE ',  &
     'RANGE OF',i5,1H-,i4,6H for ',A4,8H' cards.)
 GO TO 760
 740  WRITE  (optape,745) ufm,record
 745  FORMAT (a23,' 8029, THE FOLLOWING CARD CONTAINS NAMES THAT ARE' ,  &
     'COMPRISED OF MORE THAN 8 CHARACTERS', //20X,20A4)
 GO TO 760
 750  WRITE  (optape,755) ufm,record
 755  FORMAT (a23,' 8036, MISSING FIELDS ON THE FOLLOWING CARD', /20X, 20A4)
 760  ierror = 1
 GO TO 100
 770  ierror = 1
 800  CALL WRITE (iscr,0,0,1)
 RETURN
END SUBROUTINE xrgdtb
