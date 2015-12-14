SUBROUTINE optp1a (elt,elop,ele,dtyp)
     
 
 INTEGER, INTENT(OUT)                     :: elt(1)
 INTEGER, INTENT(OUT)                     :: elop(2,1)
 REAL, INTENT(OUT)                        :: ele(1)
 INTEGER, INTENT(IN)                      :: dtyp(1)
 INTEGER :: iwd(28),count,est,sysbuf,outtap,  &
     ycor,ecor,b1p1,ie(1),ipt(21),imat(1),NAME(2),
 
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / skp1(2),count,skp2(2),ycor,b1p1,npow,  &
     nelw,nwdse,nprw,nwdsp,skp3, mpt,skp4(3),est,skp5(2),neltyp,itype(21)
 COMMON /optpw1/ ecor,e(1)
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 COMMON /system/ sysbuf,outtap
 COMMON /matin / matid,inflag,temp,pla,sinth,costh
 COMMON /matout/ omat(20)
 COMMON /names / nrd,noeor,nwrt,nweor
 EQUIVALENCE     (e(1),ie(1)), (omat(1),imat(1)), (k1,c1)
 DATA    NAME  / 4H opt,4HP1A  /
 
!     POINTER TO IPT ARRAY - ZERO MEANS ELEMENT NOT USED.
!     UPDATE IPT DIMENSIONS AS NEW ELEMENTS ARE ADDED
 
 DATA    ipt   / 15,17,21,11,23,  25,11,13,11, 1,  5,9,7,9,19,  &
     9, 9, 3,27,27,   0/
 
!     WORD POINTER TO EST AND MATERIAL STRESS LIMITS
!     WORD 1 = 100*WORD TO OPTIMIZE (EST - IF.NE.0) + ALTERNATE
!     WORD 2 = 100*WORD FOR STRESS LIMIT + ALTERNATE
!              WHERE 1 = SHEAR
!                    2 = TENSION/COMPRESSION
!                    3 = ANY/ALL NONZERO
 
 DATA    iwd   / 506,201 , 500,200 , 700,100 , 709,303 , 700,300 ,
!                      11        13        15        17        19  &
 800,300 , 810,303 ,1718,202 , 910,202 ,1011,303 ,
!                      21        23        25       27  &
 1300,303 , 800,303 , 800,303 ,1400,300 /
 
 nelw   = 0
 sinth  = 0.0
 costh  = 1.0
 pla    = 0.0
 inflag = 2
 
!     COPY POINTER ARRAY INTO CORE
 
 DO  i = 1,ntypes
   elt(i) = dtyp(i)
 END DO
 
!     ZERO OUT POINTER ARRAY
 
 i1 = 2*(npow+1)
 DO  i = 2,i1
   elop(i,1) = 0
 END DO
 elop(1,1) = 1
 
!     READ IN ELEMENT TYPE
 
 30 CALL READ (*120,*170,est,ietyp,1,noeor,i)
 IF (ietyp > ntypes) GO TO 110
 intyp = dtyp(ietyp)
 IF (intyp <= 0) GO TO 110
 
!     DECODE LIMITS NEEDED
 
 i  = ipt(intyp)
 j2 = iwd(i)
 j1 = j2/100
 j2 = j2 - j1*100
 i2 = iwd(i+1)
 i1 = i2/100
 i2 = i2 - i1*100
 nest = (ietyp-1)*incr + 12
 nest = NE(nest)
 IF (nest > ecor) GO TO 200
 40 CALL READ (*160,*115,est,e,nest,noeor,k1)
 matid = ie(j1-1)
 IF (matid == 0) matid = ie(j2-1)
 temp = e(nest)
 CALL mat (ie(1))
 
!     TEST IF PERTINENT STRESS LIMITS ARE ZERO
 
 k1 = 0
 k2 = 0
 IF (i1 == 2 .AND. i2 == 2) GO TO 50
 
!     SHEAR
 
 IF (omat(15) /= 0.0 ) GO TO 50
 IF (i1 /= 2) k1 = 1
 IF (i2 == 1  .OR. i2 == 3) k2 = 1
 50 IF (i1 == 1 .AND. i2 <= 1) GO TO 70
 
!     TENSION
 
 IF (omat(13) /= 0.0) GO TO 60
 IF (i1 > 1) k1 = k1 + 1
 IF (i2 > 1) k2 = k2 + 1
 
!     COMPRESSION
 
 60 IF (omat(14) /= 0.0) GO TO 70
 IF (i1 > 1) k1 = k1 + 1
 IF (i2 > 1) k2 = k2 + 1
 
 70 IF (k1 >= i1 .AND. k2 >= i2) GO TO 40
 
!     CHECK IF PROPERTY IS NONZERO AND STORE INFO IN PID POINTER
 
 IF (e(j1) /= 0.0) GO TO 90
 80 IF (e(j2) == 0.0) GO TO 40
 
 IF (k2 >= i2) GO TO 40
 
!     ALTERNATE PROPERTY USED
 
 k1 = j2*100 + i2
 GO TO 100
 
 90 IF (k1 >= i1) GO TO 80
 
!     PRIMARY PROPERTY USED
 
 k1 = j1*100 + i1
 100 IF (nelw+5 > ycor) GO TO 180
 ele(nelw+1) = e(1)
 ele(nelw+2) = omat(13)
 ele(nelw+3) = omat(14)
 ele(nelw+4) = omat(15)
 
!     NOTE, K1 = C1
 
 ele(nelw+5) = c1
 nelw = nelw + nwdse
 GO TO 40
 
!     NEW ELEMENT TYPE
 
 110 CONTINUE
 CALL fread (est,0,0,nweor)
 IF (ietyp > ntypes) GO TO 120
 IF (intyp <= 0) GO TO 30
 115 elop(1,intyp+1) = nelw + 1
 GO TO 30
 
!     EOF
 
 120 i1 = npow + 1
 DO  i = 2,i1
   IF (elop(1,i) > 0) CYCLE
   elop(1,i) = elop(1,i-1)
 END DO
 IF (nelw /= 0) GO TO 150
 CALL page2 (-2)
 WRITE  (outtap,140) ufm
 140 FORMAT (a23,' 2295, NO ELEMENTS EXIST FOR OPTIMIZATION.')
 count = -1
 150 RETURN
 
!     ILLEGAL EOF
 
 160 CALL mesage (-2,est,NAME)
 
!     ILLEGAL EOR
 
 170 CALL mesage (-3,est,NAME)
 
!     INSUFFICIENT CORE
 
 180 CALL page2 (-2)
 WRITE  (outtap,190) ufm,NAME,b1p1,ie(1)
 190 FORMAT (a23,' 2296, INSUFFICIENT CORE ',2A4,1H(,i10,' ), ELEMENT', i9)
 nelw = 0
 GO TO 150
 200 CALL page2 (-2)
 WRITE (outtap,190) NAME,ecor,ietyp
 nelw = 0
 GO TO 150
END SUBROUTINE optp1a
