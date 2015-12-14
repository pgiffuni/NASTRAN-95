SUBROUTINE smcrtr ( zr, zd )
     
!  THIS SUBROUTINE MOVES DATA FROM STRINGS TO OPEN CORE AND PERFORMS
!  ANY TYPE CONVERSIONS REQUIRED.  KTYPE IS THE TYPE THAT THE
!  DECOMPOSITION IS TO BE DONE.  MTYPE IS THE TYPE OF INPUT DATA ON
!  THE MATRIX TO BE DECOMPOSED.  ISKIP IS THE NUMBER OF TERMS AT THE
!  BEGINNING OF THE STRING TO SKIP OVER.
 
 
 REAL, INTENT(OUT)                        :: zr(10)
 DOUBLE PRECISION, INTENT(OUT)            :: zd(10)
 
 DOUBLE PRECISION :: xnd(10)
 INCLUDE 'SMCOMX.COM'
 COMMON  /zzzzzz/  xns(10)
 EQUIVALENCE       ( xns, xnd )
 EQUIVALENCE       ( mblk(6), mterms ), (mblk(5), mstr )
 EQUIVALENCE       ( mblk(4), mrow   ), (mblk(2), mtype)
 SELECT CASE ( ktype )
   CASE (    1)
     GO TO  1000
   CASE (    2)
     GO TO  2000
   CASE (    3)
     GO TO  3000
   CASE (    4)
     GO TO  4000
 END SELECT
 1000  SELECT CASE ( mtype )
   CASE (    1)
     GO TO  1010
   CASE (    2)
     GO TO  1020
   CASE (    3)
     GO TO  1030
   CASE (    4)
     GO TO  1040
 END SELECT
 
! INPUT IS RS AND DECOMPOSITION TO BE DONE IN RS
 
 1010  CONTINUE
 
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zr( indexv + k - 1 ) = xns( istr + k - 1 )
 END DO
 indexv = indexv + num
 GO TO 7000
 
! INPUT IS RD AND DECOMPOSITION TO BE DONE IN RS
 
 1020  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zr( indexv + k - 1 ) = xnd( istr + k - 1 )
 END DO
 indexv  = indexv  + num
 GO TO 7000
 
! INPUT IS CS AND DECOMPOSITION TO BE DONE IN RS
 
 1030  CONTINUE
 istr   = mstr   + iskip*2
 num    = mterms - iskip
 DO  k = 1, num
   zr( indexv + k - 1 ) = xns( istr + (k-1)*2 )
 END DO
 indexv = indexv + num
 GO TO 7000
 
! INPUT IS CD AND DECOMPOSITION TO BE DONE IN RS
 
 1040  CONTINUE
 istr   = mstr   + iskip*2
 num    = mterms - iskip
 DO  k = 1, num
   zr( indexv + k - 1 ) = xnd( istr + (k-1)*2 )
 END DO
 indexv  = indexv  + num
 GO TO 7000
 2000  SELECT CASE ( mtype )
   CASE (    1)
     GO TO  2010
   CASE (    2)
     GO TO  2020
   CASE (    3)
     GO TO  2030
   CASE (    4)
     GO TO  2040
 END SELECT
 
! INPUT IS RS AND DECOMPOSITION TO BE DONE IN RD
 
 2010  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zd( indexvd + k - 1 ) = xns( istr + k - 1 )
 END DO
 indexv  = indexv  + 2*num
 indexvd = indexvd + num
 GO TO 7000
 
! INPUT IS RD AND DECOMPOSITION TO BE DONE IN RD
 
 2020  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zd( indexvd + k - 1 ) = xnd( istr + k - 1 )
 END DO
 indexvd = indexvd + num
 indexv  = indexv  + 2*num
 GO TO 7000
 
! INPUT IS CS AND DECOMPOSITION TO BE DONE IN RD
 
 2030  CONTINUE
 istr   = mstr + iskip*2
 num    = mterms - iskip
 DO  k = 1, num
   zd( indexvd + k - 1 ) = xns( istr + (k-1)*2 )
 END DO
 indexv  = indexv  + 2*num
 indexvd = indexvd + num
 GO TO 7000
 
! INPUT IS CD AND DECOMPOSITION TO BE DONE IN RD
 
 2040  CONTINUE
 istr   = mstr + iskip*2
 num    = mterms - iskip
 DO  k = 1, num
   zd( indexvd + k - 1 ) = xnd( istr + (k-1)*2 )
 END DO
 indexvd = indexvd + num
 indexv  = indexv  + 2*num
 GO TO 7000
 3000  SELECT CASE ( mtype )
   CASE (    1)
     GO TO  3010
   CASE (    2)
     GO TO  3020
   CASE (    3)
     GO TO  3030
   CASE (    4)
     GO TO  3040
 END SELECT
 
! INPUT IS RS AND DECOMPOSITION TO BE DONE IN CS
 
 3010  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zr( indexv + (k-1)*2 )   = xns( istr + k - 1 )
   zr( indexv + (k-1)*2+1 ) = 0.0
 END DO
 indexv = indexv + 2*num
 GO TO 7000
 
! INPUT IS RD AND DECOMPOSITION TO BE DONE IN CS
 
 3020  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zr( indexv + (k-1)*2 )     = xnd( istr + k - 1 )
   zr( indexv + (k-1)*2 + 1 ) = 0.0D0
 END DO
 indexv  = indexv  + 2*num
 GO TO 7000
 
! INPUT IS CS AND DECOMPOSITION TO BE DONE IN CS
 
 3030  CONTINUE
 istr   = mstr + iskip*2
 num    = ( mterms - iskip ) * 2
 DO  k = 1, num
   zr( indexv + k - 1 ) = xns( istr + k - 1 )
 END DO
 indexv = indexv + num
 GO TO 7000
 
! INPUT IS CD AND DECOMPOSITION TO BE DONE IN CS
 
 3040  CONTINUE
 istr   = mstr + iskip*2
 num    = ( mterms - iskip ) * 2
 DO  k = 1, num
   zr( indexv + k - 1 ) = xnd( istr + k - 1 )
 END DO
 indexv  = indexv  + num
 GO TO 7000
 4000  SELECT CASE ( mtype )
   CASE (    1)
     GO TO  4010
   CASE (    2)
     GO TO  4020
   CASE (    3)
     GO TO  4030
   CASE (    4)
     GO TO  4040
 END SELECT
 
! INPUT IS RS AND DECOMPOSITION TO BE DONE IN CD
 
 4010  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zd( indexvd + (k-1)*2 )   = xns( istr + k - 1 )
   zd( indexvd + (k-1)*2+1 ) = 0.0
 END DO
 indexv = indexv   + 4*num
 indexvd = indexvd + 2*num
 GO TO 7000
 
! INPUT IS RD AND DECOMPOSITION TO BE DONE IN CD
 
 4020  CONTINUE
 istr   = mstr   + iskip
 num    = mterms - iskip
 DO  k = 1, num
   zd( indexvd + (k-1)*2 )     = xnd( istr + k - 1 )
   zd( indexvd + (k-1)*2 + 1 ) = 0.0D0
 END DO
 indexv  = indexv  + 4*num
 indexvd = indexvd + 2*num
 GO TO 7000
 
! INPUT IS CS AND DECOMPOSITION TO BE DONE IN CD
 
 4030  CONTINUE
 istr   = mstr + iskip*2
 num    = ( mterms - iskip ) * 2
 DO  k = 1, num
   zd( indexvd + k - 1 ) = xns( istr + k - 1 )
 END DO
 indexv  = indexv + 2*num
 indexvd = indexvd + num
 GO TO 7000
 
! INPUT IS CD AND DECOMPOSITION TO BE DONE IN CD
 
 4040  CONTINUE
 istr   = mstr + iskip*2
 num    = ( mterms - iskip ) * 2
 DO  k = 1, num
   zd( indexvd + k - 1 ) = xnd( istr + k - 1 )
 END DO
 indexv  = indexv  + 2*num
 indexvd = indexvd + num
 GO TO 7000
 7000  RETURN
END SUBROUTINE smcrtr
