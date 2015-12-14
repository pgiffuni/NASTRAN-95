SUBROUTINE sma3b(iflag,izk)
     
!     THIS ROUTINE PROCESSES A GENERAL ELEMENT FROM GEI
 
!     IT  PRODUCES A ZE MATRIX OR A ZINVS MATRIX AND A SE MATRIX
 
!     ASSUMES GEI SITS AT BEGINNING OF UI SET AND IS OPEN TO READ
 
 
 INTEGER, INTENT(OUT)                     :: iflag
 INTEGER, INTENT(IN OUT)                  :: izk
 DOUBLE PRECISION :: d11
 INTEGER :: sysbuf,ze,se,gei,se1,ze1,zinvs
 DIMENSION ze(7),se(7)
 
 COMMON  /zzzzzz/ z(1)
 COMMON /system/ sysbuf
 COMMON /zblpkx/d11(2),idx
 COMMON/genely/gei,dum1(2),ze1,se1,id1,zinvs,dum2(22),id(7), dum4(35),m,n
 
!     COMPUTE LENGTH OF VARIABLE CORE
 
 nz  = korsz(z)-sysbuf
 iflag=-1
 nz =nz -sysbuf
 
! SKIP M+N WORDS ON GEI FILE
 
 CALL fread(gei,z,m+n,0)
 
! READ FLAG VARIABLE FOR Z OR K MATRIX
 
 CALL fread(gei,izk,1,0)
 
! IF Z MATRIX INPUT,WRITE ON ZE1 FILE
! IF K MATRIX INPUT,WRITE ON ZINVS FILE
 
 CALL makmcb(ze,ze1,m,6,2)
 IF (izk == 2) ze(1)=zinvs
 CALL makmcb(se,se1,n,2,2)
 
! READY FOR PACKING MATRICES
 
 
! OPEN ZE MATRIX
 
 CALL gopen(ze,z(nz+1),1)
 
!     LOOP ON M COLUMNS OF ZE
 
 DO  i=1,m
   CALL bldpk(2,2,ze(1),0,0)
   DO  j=1,m
     CALL fread(gei,z,1,0)
     d11(1) = z(1)
     idx = j
     CALL zblpki
   END DO
   CALL bldpkn(ze(1),0,ze)
 END DO
 CALL CLOSE( ze(1),1)
 CALL wrttrl( ze )
 IF(n == 0) GO TO 50
 iflag =1
 
!     NOW BUILD SE TRANSPOSE
 
 
!     OPEN AND WRITE HEADER
 
 CALL gopen(se,z(nz+1),1)
 
!     LOOP ON N COLUMNS OF SE
!     LOOP ON M COLUMNS OF SE  TRANSPOSE
 
 DO  i=1,m
   CALL bldpk(2,2,se(1),0,0)
   DO  j=1,n
     CALL fread(gei,z,1,0)
     d11(1) = -z(1)
     idx = j
     CALL zblpki
   END DO
   CALL bldpkn(se(1),0,se)
 END DO
 CALL CLOSE(se(1),1)
 CALL wrttrl(se)
 
!     BACKSPACE GEI SO UD AND UI AVAILABLE LATER
 
 50 CALL bckrec(gei)
 CALL CLOSE(gei,2)
 id(1) = id1
 id(2)=m
 id(3)=m
 id(4)=8
 id(5)=2
 id(6)=1
 id(7)=0
 RETURN
END SUBROUTINE sma3b
