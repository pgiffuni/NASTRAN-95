SUBROUTINE optp1c (elt,elop,pr)
     
 
 INTEGER, INTENT(IN OUT)                  :: elt(1)
 INTEGER, INTENT(IN)                      :: elop(2,2)
 INTEGER, INTENT(IN OUT)                  :: pr(1)
 INTEGER :: count, ept, sysbuf,outtap,  &
     ycor,prcor,prc,NAME(2),card(2),dtyp(21),b1p1,ENTRY
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / skp1(2),count,skp2(2),ycor,b1p1,npow,  &
     skp3(2),nprw,nwdsp,skp4, skp5,ept,skp6(6),ENTRY(21)
 COMMON /optpw1/ prcor,prc(1)
 COMMON /zzzzzz/ x(1)
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 COMMON /system/ sysbuf,outtap
 COMMON /names / nrd,noeor,nwrt,nweor
 EQUIVALENCE     (m1,rm1)
 DATA    NAME  / 4H opt,4HP1C  /,  rm1 / -1.0 /
 
!      PROPERTY CORRELATOR TO EST DESIGN VARIABLE (100*EST LOCATION).
!      THIS VALUE ADDS/SUBTRACTS FROM EST ENTRY TO GET EPT LOCATION.
!      ENTRY IS MADE BY THE ELT ARRAY (SEQUENTIAL LIST OF NUMBERS WITH
!      ZEROS FOR ELEMENTS NOT USED).
 
 DATA    dtyp  / -14, -6, -10, -5, -5, -5, -5, -5, -5, -2, &
!              BR  EB   IS  QM  M1  M2  QP  Q1  Q2  RD  &
 -4, -4,  -4, -4, -7, -4, -4, -2, -5, -5, &
!              SH  TB   T1  T2  T6  TM  TP  TU  Q4  T3  &
 0/
 
 jetyp = 1
 idps  = elop(2,1)
 idpe  = elop(2,2) - 1
 
 DO  ietyp = 1,ntypes
   IF (elt(ietyp) <= 0) CYCLE
   npr = (idpe+1-idps)/nwdsp
   IF (npr < 0) THEN
     GO TO    30
   ELSE IF (npr == 0) THEN
     GO TO    70
   END IF
   
   10 idx = ENTRY(jetyp)
   idx = incr*(idx-1)
   idp = idx + 7
   card(1) = NE(idp  )
   card(2) = NE(idp+1)
   IF (NE(idp+2) > prcor) GO TO 130
   
   CALL locate (*110,x(b1p1),card,i)
   icpr = pr(idps)
   icpt = idps
   
   20 CALL READ (*150,*160,ept,prc,NE(idp+2),noeor,i)
   
!     SEQUENTIAL PROPERTY SEARCH.  PROPERTIES THAT ARE UNSORTED ON EPT
!     WILL FAIL.  THIS MAY OCCUR FOR 2 PID/CARD (E.G., QDMEM, QUAD2,
!     SHEAR, TRIA2, TRMEM).
   
   IF (prc(1)-icpr < 0.0) THEN
     GO TO    20
   ELSE IF (prc(1)-icpr == 0.0) THEN
     GO TO    50
   END IF
   
!     LOGIC OR UNSORTED FILE ERROR
   
   30 CALL page2 (-2)
   WRITE  (outtap,40) sfm,ietyp,prc(1),NAME
   40 FORMAT (a25,' 2299, INCORRECT LOGIC FOR ELEMENT TYPE',i4,  &
       ', PROPERTY',i9,2H (,2A4,2H).)
   GO TO 100
   
!     PROPERTY IN CORE LOCATED.
   
   50 npr = npr - 1
   pr(icpt+5) = 0
   pr(icpt+4) = m1
   
!     LOCATE VARIABLE AS SET BY OPTP1A
   
   j1 = pr(icpt+1)/100
   j2 = j1+dtyp(jetyp)
   pr(icpt+3) = prc(j2)
   pr(icpt+2) = prc(j2)
   
!     ICPT+0, +1 SET BY OPTP1A
   
   icpt = icpt + nwdsp
   IF (icpt > idpe) GO TO 60
   icpr = pr(icpt)
   GO TO 20
   
!     NEW ELEMENT TYPE COMING
   
   60 IF (npr > 0) GO TO 30
   CALL fread (ept,0,0,nweor)
   70 idps  = idpe  + 1
   jetyp = jetyp + 1
   IF (jetyp > npow) GO TO 90
   idpe = elop(2,jetyp+1) - 1
 END DO
 
 
 90 RETURN
 
!     ERRORS
 
 100 count = -1
 GO TO 90
 
!     UNABLE TO LOCATE SORTED PID
 
 110 WRITE  (outtap,120) sfm,NAME,prc(1)
 120 FORMAT (a25,' 2300, ',2A4,'UNABLE TO LOCATE PROPERTY',i10,  &
     ' ON EPT OR IN CORE.')
 GO TO 100
 
!     INSUFFICIENT CORE /OPTPW1/
 
 130 CALL page2 (-2)
 WRITE  (outtap,140) ufm,NAME,prcor,ietyp
 140 FORMAT (a23,' 2296. INSUFFICIENT CORE ',2A4,1H(,i10,' ), ELEMENT', i9)
 GO TO 100
 
!     ILLEGAL EOF
 
 150 CALL mesage (-2,ept,NAME)
 
!     ILLEGAL EOR
 
 160 CALL mesage (-3,ept,NAME)
 GO TO 100
END SUBROUTINE optp1c
