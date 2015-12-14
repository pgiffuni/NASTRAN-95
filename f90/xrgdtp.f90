SUBROUTINE xrgdtp
!****
!    PURPOSE - XRGDTP DETERMINES A TYPE CODE FOR A CHARACTER
 
!    AUTHOR  - RPK CORPORATION; DECEMBER, 1983
 
!    INPUT
!      /XRGDXX/
!        ICHAR       AN ARRAY IN 80A1 FORMAT
!        ICOL        CURRENT ELEMENT IN THE ARRAY ICHAR
 
!    OUTPUT
!      /XRGDXX/
!        ITYPE       TYPE CODE ASSOCIATED WITH THE CHARACTER
!                    =1, IF CHARACTER IS A NUMBER
!                    =2, IF CHARACTER IS A ','
!                    =3, IF CHARACTER IS A '-'
!                    =4, IF CHARACTER IS A BLANK
!                    =5, OTHERWISE
!        NUMBER      INTEGER VALUE FOR CHARACTER OF ITYPE=1
 
!    LOCAL VARIABLES
!      DELIM         3 WORD ARRAY CONTAINING A COMMA, DASH, AND BLANK
!      NUMS          10 WORD ARRAY OF ALPHA NUMBERS 1,2..0
!      K             K DO LOOP INDEX TO SEARCH DELIM ARRAY
 
!   SUBROUTINES CALLED - NONE
 
!   CALLING SUBROUTINES - XRGDEV
 
!    FUNCTIONS - XRGDTP EXAMINES THE CHARACTER IN ICHAR(ICOL)
!                TO DETERMINE ITS TYPE CODE.
 
!    ERRORS - NONE
 
!****
 INTEGER :: record
 INTEGER :: nums( 10 ), delim( 3 )
 COMMON / xrgdxx / irestr, nsubst, iphase, icol  , NUMBER, itype  &
     ,                 istate, ierror, num(2), ind   , nument  &
     ,                 record(20)    , ICHAR(80)     , limit(2)  &
     ,                 icount, idmap , iscr  , NAME(2), member(2)  &
     ,                 ignore
 DATA nums / 1H1, 1H2, 1H3, 1H4, 1H5, 1H6, 1H7, 1H8, 1H9, 1H0 /
 DATA delim/ 1H,, 1H-, 1H                                     /
 
 DO  k = 1,3
   IF ( ICHAR( icol ) /= delim( k ) ) CYCLE
   itype = k + 1
   GO TO 30
 END DO
 DO  k = 1, 10
   IF ( ICHAR( icol ) /= nums( k ) ) CYCLE
   itype = 1
   NUMBER = MOD( k,10 )
   GO TO 30
 END DO
 itype = 5
 30   RETURN
END SUBROUTINE xrgdtp
