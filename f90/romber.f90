SUBROUTINE romber (b,precis,itdone,fintg,k,x)
     
!     TO INTEGRATE F(X) FROM X = A TO X = B
 
!         A = LOWER LIMIT = 0.0
!         B = UPPER LIMIT
!         NOSIG = NUMBER OF CORRECT SIGNIFICANT DIGITS DESIRED
!                 (NOT MORE THAN 7) = 5
!         PRECIS = 0.0  UPON RETURN, PRECIS = ACTUAL NUMBER
!                       OF SIGNIFICANT DIGITS ATTAINED
!         NUM = MAXIMUM NUMBER OF HALVINGS OF B-A TO BE MADE
!               (NOT MORE THAN 99) = 15
!         PROGRAMMER MUST ALSO PROGRAM SUBROUTINE FUNCT(F,X)
!         GIVING F AS THE DESIRED FUNCTION OF X, UPON
!         DEMAND OF ROMBER
!         UPON RETURN FROM ROMBER, THE VALUE OF THE INTEGRAL
!         WILL BE FOUND IN FINTG
!         IT IS CUSTOMARY TO MEASURE THE PRECISION OF LARGE NUMBERS IN
!     TERMS OF NUMBER OF SIGNIFICANT DIGITS AND THE ACCURACY OF SMALL
!     NUMBERS IN TERMS OF NUMBER OF SIGNIFICANT DECIMALS.  TO CONFORM
!     TO THIS PRACTICE, THE SUBROUTINE TERMINATES WHEN EITHER OF THESE
!     CONDITIONS IS MET.
 
 
 REAL, INTENT(IN)                         :: b
 REAL, INTENT(OUT)                        :: precis
 INTEGER, INTENT(OUT)                     :: itdone
 REAL, INTENT(OUT)                        :: fintg
 INTEGER, INTENT(IN OUT)                  :: k
 REAL, INTENT(OUT)                        :: x(6)
 
 DIMENSION   faaaa(20),faaab(20)
 
 faaac = 0.00001
 iaaaa = 1
 faaad = b
 x(1)  = 0.0
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
 9986 iaaad = iaaad + 1
 faaaf = iaaad
 x(1)  = (2.0*faaaf-1.0)*faaad
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
 faaae = faaae+f
 IF (iaaad < iaaac) GO TO 9986
 faaab(1) = 0.5*faaaa(1)+faaad*faaae
 iaaaa = iaaaa+1
 DO  iaaad = 2,iaaaa
   faaag = 4.0**(iaaad-1)
   faaah = faaag-1.0
   iaaaf = iaaad-1
   faaab(iaaad) = (faaag*faaab(iaaaf)-faaaa(iaaaf))/faaah
 END DO
 iaaac = 2*iaaac+1
 diff  = faaab(iaaaa)- faaaa(iaaaa-1)
 IF (ABS(diff) - ABS(faaac*faaab(iaaaa)) < 0.0) THEN
   GO TO  9979
 END IF
 9981 DO  iaaad = 1,iaaaa
   faaaa(iaaad) = faaab(iaaad)
 END DO
 IF (iaaaa < 15) GO TO 9988
 9979 precis = diff
 itdone = iaaaa - 1
 fintg  = faaab(iaaaa)
 RETURN
 
!     THIS CODE REPLACES F4
 
 1000 CONTINUE
 IF (x(1) == 0.0) GO TO 1020
 den = x(3) - x(2)*x(5) + x(2)*x(5)*COS(x(1))+x(2)*x(4)*SIN(x(1))
 f   = x(1)**(x(6)-1.0)* SIN(x(1))**2 / den
 1010 GO TO iret, (100,200,300)
 1020 f = 0.0
 GO TO 1010
 
!     THIS CODE REPLACES F5
 
 2000 CONTINUE
 IF (x(1) == 0.0) GO TO 2020
 den = x(3) - x(2)*x(5) + x(2)*x(5)*COS(x(1))+x(2)*x(4)*SIN(x(1))
 f   = x(1)**(x(6)-1.0)*2. *SIN(x(1))*COS(x(1))/den
 2010 GO TO iret, (100,200,300)
 2020 f = 0.0
 GO TO 2010
 
!     THIS CODE REPLACES F6
 
 3000 CONTINUE
 den = x(3) - x(2)*x(5) + x(2)*x(5)*COS(x(1))+x(2)*x(4)*SIN(x(1))
 IF (x(6) == 1.0) const = 1.0
 IF (x(6) /= 1.0) const = x(1)**(x(6) - 1.0)
 IF (den  == 0.0) GO TO 3020
 f = const*COS(x(1))**2 / den
 3010 GO TO iret, (100,200,300)
 3020 f = 0.0
 GO TO 3010
END SUBROUTINE romber
