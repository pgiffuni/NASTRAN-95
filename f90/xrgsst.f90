SUBROUTINE xrgsst (newsol)
     
!     PURPOSE - XRGSST PROCESSES SUBSTRUCTURE CONTROLS CARDS IN
!               A RIGID FORMAT (I.E., THE ****PHS- CARDS)
 
!     AUTHOR  - RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!       ARGUMENTS
!           NEWSOL     SOLUTION NUMBER
!       OTHER
!         /SYSTEM/
!           OPTAPE     UNIT NUMBER CONTAINING THE PRINT FILE
!         /XRGDXX/
!           ICHAR      CONTAINS THE CARD IMAGE IN 80A1 FORMAT
!           IDMAP      CURRENT DMAP SEQUENCE NUMBER
!           IPHASE     PHASE NUMBER
!           RECORD     CARD IMAGE IN 20A4 FORMAT
 
!     OUTPUT
!         /XRGDXX/
!           ICOL       COLUMN NUMBER LAST PROCESSED
!           IERROR     ERROR FLAG - NON-ZERO IF AN ERROR OCCURRED
 
!     LOCAL VARIABLES
!         BEGIN        CONTAINS THE VALUE 1HB
!         BLANK        CONTAINS THE VALUE 1H
!         DELETE       CONTAINS THE VALUE 1HD
!         EIGHT        CONTAINS THE VALUE 1H8
!         END          CONTAINS THE VALUE 1HE
!         FIVE         CONTAINS THE VALUE 1H5
!         ICFLAG       FLAG TO DISTINGUISH WHICH COMMON BLOCK IS
!                      BEING PROCESSED
!                      =1, /PHAS11/ ; =2, /PHAS25/ ; =3, /PHAS28/
!                      =4, /PHAS31/ ; =5, /PHAS37/
!         IFLAG        FLAG FOR THE KIND OF COMMAND BEING PROCESSED
!                      =1, FOR INSERT; =2, FOR DELETE;
!                      =3, FOR DELETE BEGIN; =4 FOR DELETE END
!         IMAP         2 WORD ARRAY FOR DMAP NUMBERS
!         IND11        INDEX FOR COMMON /PHAS11/
!         IND25        INDEX FOR COMMON /PHAS25/
!         IND28        INDEX FOR COMMON /PHAS28/
!         IND31        INDEX FOR COMMON /PHAS31/
!         IND37        INDEX FOR COMMON /PHAS37/
!         LFLAG        ARRAY USED FOR THE LAST FLAG (I.E., IFLAG)
!                      THAT WAS APPLIED TO A GIVEN COMMON - THIS
!                      IS USED TO CHECK FOR MATCHING 'DB' AND 'DE'
!                      SUBCOMMANDS
!         NMAP         NUMBER OF DMAP NUMBERS IN THE ARRAY IMAP
!         ONE          CONTAINS THE VALUE 1H1
!         SEVEN        CONTAINS THE VALUE 1H7
 
!     FUNCTIONS
!       XRGSST PROCESSES SUBSTRUCTURE CONTROL COMMANDS WITHIN THE
!       RIGID FORMAT.  THE COMMANDS ARE OF THE FOLLOWING FORMAT;
!       ****PHS- I=   (OR INSTEAD OF I=; D=, DB= OR DE= ) WHERE
!       THE '-' OF PHS IS THE PHASE NUMBER AND = REFERS TO THE
!       APPROPIATE ASCM== SUBROUTINE.  FOR THE I= SUBCOMMAND,
!       TWO NUMBERS ( N AND 0 ) ARE ADDED TO THE APPROPIATE
!       COMMON.  FOR THE D= SUBCOMMAND, TWO NUMBERS ( N1 AND N1 )
!       ARE ADDED TO THE APPROPIATE COMMON.  FOR THE DB=
!       SUBCOMMAND, ONE NUMBER IS ADDED TO THE COMMON AND
!       FOR THE DE= SUBCOMMAND, ONE NUMBER IS ADDED TO THE COMMON.
!       THE NUMBER THAT IS ADDED TO THE COMMONS
!       IS THE CURRENT DMAP SEQUENCE NUMBER AS FOUND IN THE
!       VARIABLE IDMAP.
!       THE I= COMMAND CORRESPONDS TO A DMAP ALTER INSERT
!       OF THE FORM ALTER N,0.  THE D= SUBCOMMAND CORRESPONDS
!       TO THE DMAP DELETE COMMAND  ALTER N1,N1.  THE DB=
!       SUBCOMMAND GIVES THE FIRST OF A RANGE OF DMAP NUMBERS
!       STATEMENTS TO BE DELETED AND THE DE= GIVES THE LAST
!       VALUE OF THE RANGE OF DMAP STATEMENTS TO BE DELETED.
!       THE COMMONS ARE NAMED PHAS== WHERE THE FIRST = REFERS
!       TO THE PHASE NUMBER AND THE SECOND = REFERS TO THE
!       APPROPIATE ASCM== SUBROUTINE.
 
!     SUBROUTINES CALLED - XDCODE
 
!     CALLING SUBROUTINES - XRGDFM
 
!     ERRORS
!       ERROR MESSAGES 8031,8032,8033,8035 ARE ISSUED
 
 
 INTEGER, INTENT(IN OUT)                  :: newsol
 INTEGER :: record,DELETE,BLANK,begin,END,icflag,optape,  &
     lflag(5),imap(2),one,five,seven,eight
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /xrgdxx/ irestr,nsubst,iphase,icol,NUMBER,itype,istate,  &
     ierror,num(2),ind,nument,record(20),ICHAR(80),  &
     limit(2),icount,idmap,iscr,NAME(2),member(2), ignore
 COMMON /phas11/ ipas11( 8)
 COMMON /phas25/ ipas25(14)
 COMMON /phas28/ ipas28(14)
 COMMON /phas31/ ipas31( 2)
 COMMON /phas37/ ipas37( 6)
 COMMON /system/ isysbf,optape,dum(98)
 DATA    blank / 1H  /, DELETE/ 1HD /, begin / 1HB /
 DATA    end   / 1HE /, lflag / 5*0 /
 DATA    one   / 1H1 /, five  / 1H5 /, seven / 1H7 /, eight / 1H8 /
 DATA    ind11 / 0   /, ind25 / 0   /, ind28 / 0   /, ind31 / 0   /
 DATA    ind37 / 0   /, insert/ 1HI /

CALL xdcode
icol  = 9
10   IF (ICHAR(icol) == BLANK ) GO TO 50
IF (ICHAR(icol) /= insert) GO TO 20
iflag = 1
nmap  = 2
imap(1) = idmap
imap(2) = 0
GO TO 100
20   IF (ICHAR(icol) /= DELETE) GO TO 710
icol  = icol + 1
IF (ICHAR(icol) == begin) GO TO 30
IF (ICHAR(icol) == END  ) GO TO 40
iflag = 2
nmap  = 2
imap(1) = idmap
imap(2) = idmap
GO TO 110
30   iflag = 3
nmap  = 1
imap(1) = idmap
GO TO 100
40   iflag = 4
nmap  = 1
imap(1) = idmap
GO TO 100
50   IF (icol >= 80) GO TO 800
icol  = icol + 1
GO TO 10
100  icol  = icol + 1
110  IF (iphase /= 1) GO TO 120
IF (ICHAR(icol) /= one) GO TO 710
icflag = 1
GO TO 200
120  IF (iphase /= 2) GO TO 140
IF (ICHAR(icol) /= five) GO TO 130
icflag = 2
GO TO 200
130  IF (ICHAR(icol) /= eight) GO TO 710
icflag = 3
GO TO 200
140  IF (ICHAR(icol) /= one) GO TO 150
icflag = 4
GO TO 200
150  IF (ICHAR(icol) /= seven) GO TO 710
icflag = 5
200  IF (iflag == 4 .AND. lflag(icflag) /= 3) GO TO 730
IF (iflag == 3 .AND. lflag(icflag) == 3) GO TO 740
IF (iflag <= 2 .AND. lflag(icflag) == 3) GO TO 740
lflag(icflag) = iflag
icol = icol + 1
SELECT CASE ( icflag )
  CASE (    1)
    GO TO 210
  CASE (    2)
    GO TO 220
  CASE (    3)
    GO TO 230
  CASE (    4)
    GO TO 240
  CASE (    5)
    GO TO 250
END SELECT
210  IF (ind11+nmap > 8) GO TO 720
DO  k = 1,nmap
  ind11 = ind11 + 1
  ipas11(ind11) = imap(k)
END DO
GO TO 800
220  IF (ind25+nmap > 14) GO TO 720
DO  k = 1,nmap
  ind25 = ind25 + 1
  ipas25(ind25) = imap(k)
END DO
GO TO 800
230  IF (ind28+nmap > 14) GO TO 720
DO  k = 1,nmap
  ind28 = ind28 + 1
  ipas28(ind28) = imap(k)
END DO
GO TO 800
240  IF (ind31+nmap > 2) GO TO 720
DO  k = 1,nmap
  ind31 = ind31 + 1
  ipas31(ind31) = imap(k)
END DO
GO TO 800
250  IF (ind37+nmap > 6) GO TO 720
DO  k = 1,nmap
  ind37 = ind37 + 1
  ipas37(ind37) = imap(k)
END DO
GO TO 800

!     ERRORS

710  j = 0
k = 1
WRITE  (optape,715) ufm,icol,record,j,(i,i=1,8),k,(j,i=1,8)
715  FORMAT (a23,' 8031, INVALID PARAMETER NEAR COLUMN ',i3,  &
    ' IN THE FOLLOWING CARD', //20X,20A4, /,(20X,i1,i9,7I10))
ierror = 1
GO TO 770
720  WRITE  (optape,725) ufm,iphase,record
725  FORMAT (a23,' 8032, ',19H' TOO MANY '****phs,i1, 9H' ENTRIES,  &
    ' error occurred on card', //20X,20A4)
GO TO 770
730  WRITE  (optape,735) ufm,record
735  FORMAT (a23,' 8033, ',34H a 'DE' ENTRY has no matching 'DB',  &
    ' ENTRY - ERROR ON CARD', //20X,20A4)
GO TO 770
740  WRITE  (optape,745) ufm,record
745  FORMAT (a23,' 8035, ', 41H attemp TO nest 'DB's OR no matching 'DE',  &
    ' - ERROR OCCURRED ON THE FOLLOWING CARD', /20X,20A4)
770  ierror = 1
800  RETURN
END SUBROUTINE xrgsst
