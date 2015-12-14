SUBROUTINE ifp1g (itype,case,isub1)
     
!     MAKE SURE THIS VERSION ALSO WORKS IN UNIVAC, IBM, CDC AND 64-BIT
!     MACHINES
!     ================================================================
!     IZZZBB = 0 (ALL BITS ZERO)
!     IBEN   = FIRST BYTE BLANK, REST IS ZERO FILL
!     EQUAL  = FIRST BYTE EQUAL, REST IS ZERO FILL
 
 
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER, INTENT(OUT)                     :: case(200,2)
 INTEGER, INTENT(IN OUT)                  :: isub1
 INTEGER :: CHAR,core(1),corey(401),equal,title
 COMMON /output/ title(32)
 COMMON /zzzzzz/ corex(1)
 COMMON /ifp1a / skip1(3),nwpc,ncpw4,skip2(4),izzzbb,istr,skip3(2), iben,equal
 EQUIVALENCE     (corex(1),corey(1)), (core(1),corey(401))
 
!     FIND EQUAL SIGN AND COPY REMAINING DATA ON CARD
 
!     OR FIND THE FIRST BLANK CHARACTER AFTER THE FIRST NON-BLANK WORD
!     (USED ONLY FOR ITYPE = 8, PTITLE, AXIS TITLE ETC. WHERE EQUAL SIGN
!     IS OPTIONAL AND NOT MANDATORY)
 
 k  = -1
 i2 = nwpc - 2
 DO  i = 1,i2
   DO  j = 1,ncpw4
     CHAR = khrfn1(izzzbb,1,core(i),j)
     IF (CHAR == equal) GO TO 170
     IF (CHAR /= iben .AND. k == -1) k = 0
     IF (CHAR == iben .AND. k == 0) k = i*100 + j
   END DO
 END DO
 IF (itype /= 8) GO TO 170
 i  = k/100
 j  = MOD(k,100)
 170 k  = (itype-1)*32
 k1 = k + 38
 IF (itype == 8) k1 = 0
 IF (j  /= ncpw4) GO TO 180
 i = i + 1
 j = 0
 180 j = j + 1
 ipos  = 1
 its   = k + 1
 isave = izzzbb
 DO  ii = i,i2
   190 isave = khrfn1(isave,ipos,core(ii),j)
   ipos  = ipos + 1
   IF (ipos > 4) GO TO 210
   200 j = j + 1
   IF (j <= ncpw4) GO TO 190
   j = 1
   CYCLE
   210 ipos = 1
   IF (itype == 7) GO TO 220
   IF (istr-1 == 0) THEN
     GO TO   230
   END IF
   220 title(its) = isave
   GO TO 240
   230 case(k1+1,isub1) = isave
   k1 = k1 + 1
   240 isave = izzzbb
   its = its + 1
   GO TO 200
 END DO
 DO  i = ipos,4
   isave = khrfn1(isave,i,iben,1)
 END DO
 IF (itype == 7) GO TO 270
 IF (istr-1 == 0) THEN
   GO TO   280
 END IF
 270 title(its) = isave
 GO TO 290
 280 case(k1+1,isub1) = isave
 290 RETURN
END SUBROUTINE ifp1g
