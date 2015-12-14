SUBROUTINE xrgnum
     
!     XRGNUM PROCESSES THE NUMBER ON A CARD OR FILE NAME TABLE ENTRY
!     THIS ROUTINE IS CALLED ONLY BY XRGDTB
 
!     WRITTEN BY  RPK CORPORATION; DECEMBER, 1983
 
!     INPUT
!       /SYSTEM/
!         OPTAPE       UNIT NUMBER FOR THE OUTPUT PRINT FILE
!       /XRGDXX/
!         ICHAR        CONTAINS THE CARD IMAGE IN 80A1 FORMAT
!         ICOL         CURRENT COLUMN BEING PROCESSED
!         RECORD       CONTAINS THE CARD IMAGE IN 20A4 FORMAT
 
!     OUTPUT
!       /XRGDXX/
!         ICOL         CURRENT COLUMN BEING PROCESSED
!         NUMBER       VALUE OF THE NUMBER IN INTEGER FORMAT
 
!     LOCAL VARIABLES
!         BLANK          CONTAINS THE VALUE 1H
!         IFRCOL         FIRST COLUMN TO BE EXAMINED BY XRGNUM
!         NEWNUM         INTEGER VALUE OF THE CHARACTER IN THE CURRENT
!                        COLUMN
!         NUMS           CONTAINS THE ALPHA VALUES 1,2,...0
 
!     THE CARD IS SCANED TO FIND THE VALUE OF THE NUMBER IN THE FIRST
!     FIELD OF THE CARD
 
!     MESSAGE 8030 MAY BE ISSUED
 
 INTEGER :: record, optape, BLANK , nums(10)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /xrgdxx/ irestr, nsubst, iphase, icol  , NUMBER, itype ,  &
     istate, ierror, num(2), ind   , nument,  &
     record(20)    , ICHAR(80)     , limit(2),  &
     icount, idmap , iscr  ,NAME(2), member(2), ignore
 COMMON /system/ isysbf, optape, dum(98)
 DATA    nums  / 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9, 1H0/
 DATA    BLANK / 1H  /
 
 ifrcol = icol
 NUMBER = 0
 50   IF (icol >= 80) GO TO 350
 IF (ICHAR(icol) == BLANK) GO TO 200
 DO  k = 1,10
   IF (ICHAR(icol) /= nums(k)) CYCLE
   newnum = MOD(k,10)
   NUMBER = NUMBER*10 + newnum
   GO TO 150
 END DO
 GO TO 250
 150  icol = icol + 1
 GO TO 50
 200  icol = icol + 1
 IF (NUMBER == 0) GO TO 50
 GO TO 350
 250  NUMBER = 0
 j = 0
 k = 1
 WRITE  (optape,300) ufm,ifrcol,record,j,(i,i=1,8),k,(j,i=1,8)
 300  FORMAT (a23,' 8030, EXPECTED AN INTEGER NEAR COLUMN',i3,  &
     ' IN THE FOLLOWING CARD', //20X,20A4, /,(20X,i1,i9,7I10))
 350  RETURN
END SUBROUTINE xrgnum
