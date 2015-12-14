SUBROUTINE rombdk (b,precis,itdone,fintg,k,x)
     
!     INTEGRATE F(X) FROM 0.0 TO X = B
 
!         B      = UPPER LIMIT
!         NOSIG  = NUMBER OF CORRECT SIGNIFICANT DIGITS DESIRED
!                  (NOT MORE THAN 7) = 5
!         PRECIS = 0.0  UPON RETURN, PRECIS = ACTUAL NUMBER
!                  OF SIGNIFICANT DIGITS ATTAINED
!         NUM    = MAXIMUM NUMBER OF HALVINGS OF B-A TO BE MADE
!                  (NOT MORE THAN 99) = 15
!     UPON RETURN FROM ROMBER, THE VALUE OF THE INTEGRAL WILL BE FOUND
!     IN FINTG
 
!     IT IS CUSTOMARY TO MEASURE THE PRECISION OF LARGE NUMBERS IN
!     TERMS OF NUMBER OF SIGNIFICANT DIGITS AND THE ACCURACY OF SMALL
!     NUMBERS IN TERMS OF NUMBER OF SIGNIFICANT DECIMALS.  TO CONFORM
!     TO THIS PRACTICE, THE SUBROUTINE TERMINATES WHEN EITHER OF THESE
!     CONDITIONS IS MET.
 
 
 DOUBLE PRECISION, INTENT(IN)             :: b
 DOUBLE PRECISION, INTENT(OUT)            :: precis
 INTEGER, INTENT(OUT)                     :: itdone
 DOUBLE PRECISION, INTENT(OUT)            :: fintg
 INTEGER, INTENT(IN OUT)                  :: k
 DOUBLE PRECISION, INTENT(OUT)            :: x(6)
 DOUBLE PRECISION :: den,const, faaaa,faaab,  &
     faaac,faaad,f,faaae,faaaf,faaag,faaah,diff
 DIMENSION  faaaa(20),faaab(20)
 
 faaac = 0.00001D0
 iaaaa = 1
 faaad = b
 x(1)  = 0.0D0
 ASSIGN 100 TO iret
 SELECT CASE ( k )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
 END SELECT
 100 CONTINUE
 faaae = f
 x(1)  = b
 ASSIGN 200 TO iret
 SELECT CASE ( k )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
 END SELECT
 200 CONTINUE
 faaae = faaae+f
 faaaa(1) = 0.5*faaad*faaae
 9988 faaad = 0.5*faaad
 iaaac = 2**(iaaaa-1)
 faaae = 0.0
 
 iaaad = 0
 9986 iaaad = iaaad+1
 faaaf = iaaad
 x(1) = (2.0*faaaf-1.0)*faaad
 ASSIGN 300 TO iret
 SELECT CASE ( k )
   CASE (    1)
     GO TO 1000
   CASE (    2)
     GO TO 2000
   CASE (    3)
     GO TO 3000
 END SELECT
 300 CONTINUE
 faaae = faaae + f
 IF (iaaad < iaaac) GO TO 9986
 
 faaab(1) = 0.5*faaaa(1) + faaad*faaae
 iaaaa = iaaaa+1
 DO  iaaad = 2,iaaaa
   faaag = 4.0**(iaaad-1)
   faaah = faaag-1.0
   iaaaf = iaaad-1
   faaab(iaaad) = (faaag*faaab(iaaaf) - faaaa(iaaaf))/faaah
 END DO
 iaaac = 2*iaaac+1
 diff  = faaab(iaaaa) - faaaa(iaaaa-1)
 IF (DABS(diff)-DABS(faaac*faaab(iaaaa)) < 0.0) THEN
   GO TO  9979
 END IF
 9981 DO  iaaad = 1,iaaaa
   faaaa(iaaad)  = faaab(iaaad)
 END DO
 IF (iaaaa < 15) GO TO 9988
 9979 precis = diff
 itdone = iaaaa - 1
 fintg  = faaab(iaaaa)
 RETURN
 
!     THIS CODE REPLACES D4K
 
 1000 CONTINUE
 IF (x(1) == 0.0D0) GO TO 1010
 den = x(3) - x(2)*x(5)*(1.d0-DCOS(x(1))) + x(2)*x(4)*DSIN(x(1))
 f   = x(1)**(x(6)-1.d0)*DSIN(x(1))**2/den
 GO TO 1020
 1010 f = 0.0D0
 1020 GO TO iret, (100,200,300)
 
!     THIS CODE REPLACES D5K
 
 2000 CONTINUE
 IF (x(1) == 0.0D0) GO TO 2010
 den = x(3) - x(2)*x(5)*(1.d0-DCOS(x(1))) + x(2)*x(4)*DSIN(x(1))
 f   = x(1)**(x(6)-1.d0)*2.*DSIN(x(1))*DCOS(x(1))/den
 GO TO 2020
 2010 f=0.0D0
 2020 GO TO iret, (100,200,300)
 
!     THIS CODE REPLACES D6K
 
 3000 CONTINUE
 den = x(3) - x(2)*x(5)*(1.d0-DCOS(x(1))) + x(2)*x(4)*DSIN(x(1))
 IF (x(6) == 1.d0) const = 1.d0
 IF (x(6) /= 1.d0) const = x(1)**(x(6)-1.d0)
 IF (den == 0.0D0) GO TO 3010
 f = const*DCOS(x(1))**2/den
 GO TO 3020
 3010 f = 0.0D0
 3020 GO TO iret, (100,200,300)
END SUBROUTINE rombdk
