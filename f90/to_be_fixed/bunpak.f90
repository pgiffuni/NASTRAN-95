SUBROUTINE bunpak (ig,i,nj,jg)
     
!     THIS ROUTINE WORKS SIMILARLY AS BUNPK EXCEPT IT UNPACKS A WHOLE
!     (I-TH) ROW OF GRID NUMBERS (1 THRU NJ) FROM IG TABLE, AND SAVES
!     THE UNPACKED DATA IN JG ARRAY.
!     (BUNPK UNPACKS ONLY AN ELEMENT OF GRID NUMBER IN IG TABLE)
 
!     THIS ROUTINE GREATLY INCREASES BANDIT INTERNAL EFFICIENCY
!     WRITTEN BY G.CHAN/UNISYS,    MAY 1988
 
 
 INTEGER*2, INTENT(IN)                    :: ig(1)
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(IN)                      :: nj
 INTEGER, INTENT(OUT)                     :: jg(1)
 IMPLICIT INTEGER (a-z)
 
!DC   NEXT 2 LINES FOR CDC AND UNIVAC ONLY
!     EXTERNAL         ANDF,     RSHIFT
!     INTEGER          ANDF,     RSHIFT  ,IG(1)
 
!     NEXT LINE FOR IBM, VAX, AND MACHINES THAT HAVE INTEGER*2
 
 
 INTEGER :: nam(2)
 COMMON /system/  ibuf,     nout
 COMMON /bandb /  nbit,     dum3b(3), ipass,    nw
 COMMON /bands /  dum4s(4), ii1,      maxdeg,   dum4(4),  mask
 DATA    nam   /  4HUNPA  , 4HK       /
 
 IF (nj <= maxdeg) GO TO 20
 WRITE  (nout,10) nj,maxdeg
 10   FORMAT ('0 *** BUNPAK .GT. MAXDEG',2I7)
 CALL errtrc (nam)
 
 20   ipass = ipass+nj
 n1 = i
 
!     ********************************************
!     UNIVAC AND CDC MACHINES
!     ********************************************
 
!     DO 40 N=1,NJ,NW
!     N2 = IG(N1)
!     N3 = N+NW-1
!     DO 30 M=1,NW
!     JG(N3) = ANDF(N2,MASK)
!     IF (M .EQ. NW) GO TO 40
!     N2 = RSHIFT(N2,NBIT)
!  30 N3 = N3-1
!  40 N1 = N1+II1
!     RETURN
 
!     ********************************************
!     IBM AND VAX MACHINES
!     ********************************************
 
 DO  n=1,nj
   jg(n) = ig(n1)
   n1 = n1+ii1
 END DO
 RETURN
END SUBROUTINE bunpak
