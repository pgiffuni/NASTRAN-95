SUBROUTINE dbmalb ( lenreq, INDEX )
!********************************************************************
!     DBMALB - ALLOCATES A MEMORY BLOCK OF LENGTH "LENREQ"
!              FROM THE FREE CHAIN AND RETURNS THE
!              POINTER IN MEMORY FOR THE BLOCK IN "INDEX" (RELATIVE TO.
!              /DBM/.
 
!     EACH FREE BLOCK IN MEMORY HAS THE FOLLOWING FORMAT:
!         WORD 1  POINTER TO PREVIOUS FREE BLOCK (=0, IF FIRST)
!         WORD 2  POINTER TO NEXT FREE BLOCK (=0, IF END OF CHAIN)
!         WORD 3  NUMBER OF WORDS AVAILABLE IN THIS FREE BLOCK
 
!     NOTE:  IDBFRE POINTS TO THE FIRST FREE BLOCK OF THE CHAIN
!********************************************************************
 INCLUDE 'DSIOF.COM'
 INCLUDE 'ZZZZZZ.COM'
 next  = idbfre
 IF ( idbfre == 0 ) GO TO 701
!  OBTAIN THE LENGTH OF THE FREE BLOCK
 10    lenavl = mem( next + 2 )
 IF ( lenavl >= lenreq ) GO TO 100
!  MEMORY NOT AVAILABLE IN THIS BLOCK, CHECK FOR OTHER BLOCKS
 next  = mem( next + 1 )
!  IF NO MORE FREE BLOCKS, RETURN WITH INDEX SET TO -1
 IF ( next == 0 ) GO TO 701
 GO TO 10
!  RETURN POINTER FOR THIS BLOCK
 100   INDEX = next
 IF ( lenreq /= lenavl ) GO TO 200
!  COME HERE WHEN REQUESTED BLOCK SAME SIZE AS FREE BLOCK
 next  = mem( next+1 )
 iprev = mem(INDEX   )
 IF ( iprev == 0 ) GO TO 110
 IF ( next  == 0 ) GO TO 120
!  CONNECT THE PREVIOUS FREE BLOCK WITH THE NEXT FREE BLOCK
 mem( iprev+1 ) = next
 mem( next    ) = iprev
 GO TO 700
 110   IF ( next == 0 ) GO TO 130
!  NO PREVIOUS BLOCK, SET IDBFRE TO POINT TO NEW FIRST FREE BLOCK
 idbfre = next
 mem( next ) = 0
 GO TO 700
!  PREVIOUS BLOCK EXITS BUT BLOCK ALLOCATED WAS LAST IN CHAIN
 120   mem( iprev+1 ) = 0
 GO TO 700
!  NO MORE FREE BLOCKS EXIST, SET IDBFRE TO ZERO
 130   idbfre = 0
 GO TO 700
!  COME HERE WHEN FREE BLOCK HAS MORE SPACE THEN REQUESTED
 200   newind = INDEX + lenreq + 4
 iprev  = mem( INDEX  )
 next   = mem( INDEX+1)
!  CHECK TO DETERMINE IF ANY SPACE REMAINS
 IF ( ( lenavl-lenreq-4 ) <= 0 ) GO TO 240
!  RECOMPUTE FREE SPACE AND SET UP CHAIN WORDS
 mem( newind+2 ) = lenavl - lenreq - 4
 IF ( iprev == 0 ) GO TO 210
 IF ( next  == 0 ) GO TO 220
!  CONNECT TO PREVIOUS AND NEXT FREE BLOCK
 mem( newind  ) = iprev
 mem( newind+1) = next
 mem( iprev+1 ) = newind
 mem( next    ) = newind
 GO TO 700
 210   IF ( next == 0 ) GO TO 230
!  NO PREVIOUS BLOCK, NEWLY CREATED BLOCK BECOMES THE FIRST FREE BLOCK
 idbfre         = newind
 mem( newind  ) = 0
 mem( newind+1) = next
 mem( next    ) = newind
 GO TO 700
!  PREVIOUS BLOCK EXISTS BUT THE NEWLY CREATED BLOCK IS LAST
 220   mem( iprev+1 ) = newind
 mem( newind  ) = iprev
 mem( newind+1) = 0
 GO TO 700
!  NEW BLOCK IS THE ONLY FREE BLOCK
 230   idbfre  = newind
 mem( newind  ) = 0
 mem( newind+1) = 0
 GO TO 700
!  FREE CHAIN IS EXHAUSTED
 240   idbfre = 0
 701   INDEX = -1
 RETURN
 700   CONTINUE
 RETURN
END SUBROUTINE dbmalb
