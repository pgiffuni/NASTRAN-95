SUBROUTINE ferlts (ifile,dz,dy,zm)
     
!   FEER MATRIX TRANSPOSE MULTIPLY  (SINGLE PRECISION)
!   SEE SUBROUTINE FERRDM FOR CONTENTS OF SMAPOS AND HOW THE MATRIX
!   DATA IS STORED IN MEMORY.
 
 
 INTEGER, INTENT(IN)                      :: ifile(7)
 REAL, INTENT(IN)                         :: dz(1)
 REAL, INTENT(OUT)                        :: dy(1)
 REAL, INTENT(IN)                         :: zm(1)
 REAL :: dsum
 REAL :: dcore(1)
 INTEGER :: smapos
 COMMON  /feerim/ nidsma    ,nidlt     ,nidorv    ,nltli  &
     ,                nsmali    ,ibfsma    ,ibflt  &
     ,                ibforv    ,smapos(7) ,ltpos(7)
 COMMON  /unpakx/ iprc      ,ip        ,np        ,incr
 COMMON  /zzzzzz/ icore(1)
 EQUIVALENCE      ( dcore(1),icore(1) )
 
 n     = ifile(2)
 iccol = 1
 IF ( nidsma == 0 ) GO TO 1005
 mem   = nidsma
 ilcol = smapos( 1 )
 DO  i = 1,n
   iccol = i
! CHECK TO SEE IF REMAINING DATA IS ON THE FILE AND NOT IN MEMORY
   IF ( iccol > ilcol ) GO TO 1000
   dy(i) = 0.
   dsum  = 0.
   5 icol  = icore(mem)
   IF( icol /= i ) CYCLE
   ntms  = icore(mem+1)
   ip    = icore(mem+2+ntms)
   np    = ip+ntms-1
   indx  = mem+1
   ii    =  0
   DO  j = ip,np
     ii = ii +1
     dsum  = dsum + dcore(indx+ii) * dz(j)
   END DO
   dy(i) = dsum
   mem   = mem+4+ntms
   GO TO 5
 END DO
 GO TO 7000
 1000 CONTINUE
 CALL dsspos ( ifile, smapos(2), smapos(3), smapos(4) )
 GO TO 1008
 1005 CALL REWIND ( ifile )
 CALL skprec ( ifile, 1 )
 1008 CONTINUE
 incr  = 1
 iprc  = ifile(5)
 DO  i = iccol, n
   dy(i) = 0.
   ip    = 0
   CALL unpack(*1020,ifile,zm(1))
   ii    = 0
   dsum  = 0.0
   DO  j = ip,np
     ii    = ii +1
     dsum  = dsum + zm(ii) * dz(j)
   END DO
   dy(i) = dsum
   1020 CONTINUE
 END DO
 7000 CONTINUE
 RETURN
END SUBROUTINE ferlts
