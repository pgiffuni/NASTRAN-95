SUBROUTINE xrgsub (irestb,subset)
!****
!     PURPOSE - XRGSUB PROCESSES THE ****SBST CONTROL CARD IN
!               RIGID FORMAT DATA BASE
 
!     AUTHOR  - RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!       ARGUMENTS
!         SUBSET        SUBSET NUMBERS GIVEN BY THE USER
!       OTHER
!         /XRGDXX/
!           NSUBST      NUMBER OF SUBSET NUMBERS GIVEN BY USER
!           NUM         2 WORD ARRAY CONTAINING A RANGE OF NUMBERS
!                       FROM THE LIST OF NUMBERS ON THE ****SBST
!                       CONTROL CARD
!         NUMENT        NUMBER OF WORDS PER ENTRY IN THE MODULE
!                       EXECUTION DECISION TABLE
 
!     OUTPUT
!       ARGUMENTS
!         IRESTB        MODULE EXECUTION DECISION TABLE ENTRY FOR
!                       CURRENT DMAP STATEMENT
!       OTHER
!         /XRGDXX/
!           ICOL      COLUMN NUMBER BEING PROCESSED ON THE CARD
!           IERROR    ERROR FLAG - NON-ZERO IF AN ERROR OCCURRED
!           IGNORE    IGNORE FLAG SET TO NON-ZERO IF THE DMAP
!                     STATEMENT IS TO BE DELETED BY THE SUBSET
!           LIMIT     LOWER/UPPER LIMITS OF VALUES WITHIN AN
!                     ENTRY ON THE CARD
 
!     LOCAL VARIABLES
!       IEND          VALUE OF NUM(2)
!       ISTR          VALUE OF NUM(1)
 
!     FUNCTIONS
!        XRGSUB CALLS XRGDEV TO EXTRAPOLATE THE NUMBER FROM THE
!        THE CARD AND THEN IT COMPARES THE NUMBER(S) WITH THOSE
!        SUPPLIED BY THE USER AS SUBSETS.  IF A MATCH IS FOUND,
!        IGNORE IS SET AND THE MODULE EXECUTION DECISION TABLE
!        ENTRY IS SET TO ZERO.  CHECKS CONTINUE UNTIL ALL VALUES
!        GIVEN ON ****SBST CARD HAD BEEN CHECK OR UNTIL A MATCH
!        IS FOUND
 
!     SUBROUTINES CALLED - XDCODE,XRGDEV
 
!     CALLING SUBROUTINES - XRGDFM
 
!     ERRORS - NONE
 
!****
 
 INTEGER, INTENT(OUT)                     :: irestb(7)
 INTEGER, INTENT(IN)                      :: subset(12)
 INTEGER :: record
 COMMON /xrgdxx/ irestr,nsubst,iphase,icol,NUMBER,itype,  &
     istate,ierror,num(2),ind,nument, record(20),ICHAR(80),limit(2),  &
     icount,idmap,iscr,NAME(2),member(2),ignore
 
 icol   = 9
 ierror = 0
 CALL xdcode
 limit(1) = 1
 limit(2) = 12
 200   CALL xrgdev
 IF (ierror /= 0 .OR. icol > 80) GO TO 700
 istr = num(1)
 iend = num(2)
 DO  k  = istr,iend
   DO  kk = 1,nsubst
     IF (k == subset(kk)) GO TO 500
   END DO
 END DO
 icol = icol + 1
 GO TO 200
 500   DO  k = 1,nument
   irestb(k) = 0
 END DO
 ignore = 1
 700   RETURN
END SUBROUTINE xrgsub
