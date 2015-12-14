SUBROUTINE xrgdcf (irestb)
     
!     PURPOSE - XRGDCF PROCESSES THE '****CARD', '****FILE' AND
!               '****RFMT' CONTROL CARDS WITHIN THE RIGID DMAP
!               DATA BASE.
 
!     AUTHOR  - RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!       /SYSTEM/
!         NOUT         UNIT NUMBER FOR THE OUTPUT PRINT FILE
!       /XRGDXX/
!         NUM          VALUE OF THE NUMBER OR RANGE OF NUMBERS
!                      IN THE CURRENT FIELD BEING PROCESSED
 
!     OUTPUT
!       ARGUMENTS
!         IRESTB       THE MODULE EXECUTION DECISION TABLE
!       OTHER
!         /XRGDXX/
!           ICOL       CURRENT COLUMN NUMBER BEING PROCESSED IN
!                      THE CARD
!           IERROR     ERROR FLAG - NON-ZERO IF AN ERROR OCCURRS
 
!     LOCAL VARIABLES
!       IBIT           BIT NUMBER FOR FLAG IN THE MODULE EXEC. DEC.
!                      TABLE
!       IEND           LAST NUMBER OF RANGE OF NUMBERS READ FROM
!                      THE CURRENT FIELD
!       ISTR           SAME AS IEND EXCEPT FIRST NUMBER
!       IWORD          SAME AS IBIT BUT REFERS TO THE WORD NUMBER
 
!     FUNCTIONS
!       XRGDCF PROCESSES THE ABOVE TYPES OF CARDS WHICH ALL HAVE
!       FORMATS AS FOLLOWS:  '****XXXX   M1,M2,..'
!       WHERE M- IS IN ANY OF THE FOLLOWING FORMS ( NNN  OR NNN-NNN).
!       NNN IS AN INTEGER NUMBER AND THE '-' REFERS TO A RANGE
!       WHERE THE RANGE MUST BE IN ASCENDING ORDER.
!       XRGDCF CALLS XDCODE TO CONVERT THE CARD IMAGE TO 80A1 AND
!       CALLS XRGDEV TO VALIDATE THE SYNTAX AND TO GET A M-
!       ENTRY FROM THE CARD.  BASED ON THE VALUE(S) RETURNED IN
!       NUM, THE CORRESPONDING BITS ARE TURNED ON IN THE MODULE
!       EXECUTION DECISION TABLE.  PROCESSING CONTINUES UNTIL ALL
!       FIELDS OF THE CARD HAVE BEEN PROCESSED.
 
 
!     SUBROUTINES CALLED - XDCODE, XRGDEV
 
!     CALLING SUBROUTINES - XRGRFM
 
!     ERRORS - NONE
 
 
 INTEGER, INTENT(OUT)                     :: irestb(7)
 EXTERNAL        orf
 INTEGER :: record, orf
 COMMON /system/ isysbf, nout  , dum(98)
 COMMON /xrgdxx/ irestr, nsubst, iphase, icol   , NUMBER, itype ,  &
     istate, ierror, num(2), ind    , nument        ,  &
     record(20)    , ICHAR(80)      , limit(2)      ,  &
     icount, idmap , iscr  , NAME(2), member(2)     , ignore
 
 ierror = 0
 icol   = 9
 CALL xdcode
 10   CALL xrgdev
 IF (ierror /= 0 .OR. icol > 80) GO TO 30
 istr   = num(1)
 iend   = num(2)
 DO  k = istr,iend
   iword = (k-1)/31
   ibit  = 2**(31*iword + 31 - k)
   irestb(iword+1) = orf(irestb(iword+1),ibit)
 END DO
 icol = icol + 1
 GO TO 10
 30   CONTINUE
 RETURN
END SUBROUTINE xrgdcf
