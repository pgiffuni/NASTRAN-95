SUBROUTINE smcdmp (  zi, zr, zd )
     
! SMCDMP DUMPS THE CONTENTS OF THE COLUMN DATA AS STORED IN MEMORY
! AND POINTED TO BY THE DIRECTORY CREATED BY SUBROUTINE SMCPH1
 
 
 INTEGER, INTENT(IN)                      :: zi(10)
 REAL, INTENT(IN OUT)                     :: zr(10)
 DOUBLE PRECISION, INTENT(IN OUT)         :: zd(10)
 
 
 
 INCLUDE   'SMCOMX.COM'
 
 DO  icol = 1, ncol
   ind = ( icol-1 ) * 4 + 1
   indexr = zi( ind )
   IF ( indexr == 0 ) GO TO 700
   icolum  = zi( indexr   )
   lenrows = zi( indexr+1 )
   nwords  = zi( indexr+2 )
   lenvals = zi( indexr+3 )
   WRITE( nout,901) icolum
   901   FORMAT(/,' -----------------------------COLUMN NUMBER   =',i10)
   WRITE ( nout, 900 )  zi(ind), zi(ind+1), zi(ind+2), zi(ind+3)  &
       ,                    lenrows, lenvals, nwords
   900   FORMAT(20X,' COLUMN DIRECTORY: INDEX                    =',i10  &
       ,/    ,20X,' FIRST COLUMN DATA NEEDED FOR THIS PIVOT    =',i10  &
       ,/    ,20X,' LAST PIVOT COLUMN TO USE THIS COLUMN       =',i10  &
       ,/    ,20X,' SAVPOS                                     =',2X,z8  &
       ,/    ,20X,' NUMBER OF WORDS DEFINING ROW NUMBERS       =',i10  &
       ,/    ,20X,' NUMBER OF NON-ZERO TERMS IN THIS COLUMN    =',i10  &
       ,/    ,20X,' TOTAL NUMBER OF WORDS ALLOCATED FOR COLUMN =',i10 )
   indexr = indexr + 4
   indexv = indexr+lenrows
   indexvd= indexv / 2 + 1
   iend   = indexv
   50    irow   = zi( indexr )
   nterms = zi( indexr+1 )
   WRITE( nout,905) irow, nterms, ktype
   905   FORMAT(' ROW, TERMS, TYPE=',3I10)
   SELECT CASE ( ktype )
     CASE (    1)
       GO TO  100
     CASE (    2)
       GO TO  200
     CASE (    3)
       GO TO  300
     CASE (    4)
       GO TO  400
   END SELECT
   100   WRITE( nout,906) (zr(indexv+k-1),k=1,nterms)
   906   FORMAT( 5E16.8 )
   indexv = indexv + nterms
   indexr = indexr + 2
   GO TO 500
   200   WRITE( nout,907) (zd(indexvd+k-1),k=1,nterms)
   907   FORMAT( 5D16.8 )
   indexv = indexv  + 2*nterms
   indexvd= indexvd + nterms
   indexr = indexr + 2
   GO TO 500
   300   WRITE( nout,906) (zr(indexv+k-1),k=1,2*nterms)
   indexv = indexv + 2*nterms
   indexr = indexr + 2
   GO TO 500
   400   WRITE( nout,907) (zd(indexvd+k-1),k=1,2*nterms)
   indexv = indexv  + 4*nterms
   indexvd= indexvd + 2*nterms
   indexr = indexr + 2
   GO TO 500
   500   IF ( indexr >= iend ) CYCLE
   GO TO 50
   700   CONTINUE
   WRITE( nout,901) icol
   WRITE ( nout, 902 )  zi(ind), zi(ind+1), zi(ind+2), zi(ind+3)
   902   FORMAT(20X,' COLUMN DIRECTORY: INDEX                    =',i10  &
       ,/    ,20X,' FIRST COLUMN DATA NEEDED FOR THIS PIVOT    =',i10  &
       ,/    ,20X,' LAST PIVOT COLUMN TO USE THIS COLUMN       =',i10  &
       ,/    ,20X,' SAVPOS                                     =',2X,z8  &
       ,/    ,20X,' NUMBER OF WORDS DEFINING ROW NUMBERS       = N/A'  &
       ,/    ,20X,' NUMBER OF NON-ZERO TERMS IN THIS COLUMN    = N/A'  &
       ,/    ,20X,' TOTAL NUMBER OF WORDS ALLOCATED FOR COLUMN = N/A' )
   WRITE( nout,908)
   908   FORMAT(20X,' ------COLUMN HAS BEEN PUT TO SPILL FILE-----')
 END DO
 7777  RETURN
END SUBROUTINE smcdmp
