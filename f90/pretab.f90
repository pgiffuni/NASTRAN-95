SUBROUTINE pretab (ditf,rz,inz,buf,lcrgvn,lcused,tabnol,list)
     
!     SUBROUTINE PRETAB READS TABLES INTO OPEN CORE, SETS UP TABLE
!     DICTIONARIES WHICH ARE LATER USED WHEN THE CALLING ROUTINE
!     REQUESTS A FUNCTIONAL VALUE FROM A TABLE VIA A CALL TO THE ENTRY
!     POINT TAB.
 
!     REVISED  7/92, BY G.CHAN/UNISYS.
!     1. NEW REFERENCE TO THE OPEN CORE ARRAY SUCH THAT THE SOURCE CODE
!        IS UP TO ANSI FORTRAN 77 STANDARD
!     2. LOGARITHMIC SCALE ENHANCEMENT
 
!     ARGUMENT LIST -
 
!     DITF     THE GINO NAME OF THE FILE ON WHICH THE TABLES RESIDE.
!     RZ       THE OPEN CORE ARRAY. RZ IS USED AS REAL BY THIS ROUTINE.
!     INZ      SAME ADDRESS AS RZ.  USED AS INTEGER IN THIS ROUTINE.
!     BUF      A BUFFER TO BE USED BY SUBROUTINE PRELOC.
!     LCRGVN   THE LENGTH OF OPEN CORE GIVEN TO PRETAB.
!     LCUSED   THE AMOUNT OF CORE USED BY PRETAB.
!     TABNOL   LIST OF TABLE NUMBERS THAT THE USER WILL BE REFERENCING.
!              TABNOL(1) = N IS THE NUMBER OF TABLES TO BE REFERENCED.
!              TABNOL(2),...,TABNOL(N+1) CONTAIN THE TABLE NUMBERS. NOTE
!              THAT 0 IS AN ADMISSIBLE TABLE NUMBER IN THE TABLE NO.
!              LIST.  TABLE NO. 0 DEFINES A FUNCTION WHICH IS IDENTICAL-
!              LY = 0 FOR ALL VALUES OF THE INDEPENDENT VARIABLE.
!     LIST     ARRAY OF CONTROL WORDS FOR SUBROUTINE LOCATE AND TABLE
!              TYPES.
!              LIST(1) = M IS THE NO. OF TRIPLES WHICH FOLLOW IN LIST.
!              THE FIRST TWO WORDS OF EACH TRIPLE ARE THE SUBROUTINE
!              LOCATE CONTROL WORDS AND THE THIRD WORD IS THE TABLE TYPE
!              = 1,2,3,4, OR 5.
!     LNTH     = 12 WORDS PER TABLE ENTRY
 
 
 INTEGER, INTENT(IN)                      :: ditf
 REAL, INTENT(IN OUT)                     :: rz(1)
 INTEGER, INTENT(IN OUT)                  :: inz(1)
 REAL, INTENT(IN OUT)                     :: buf(1)
 INTEGER, INTENT(IN)                      :: lcrgvn
 INTEGER, INTENT(OUT)                     :: lcused
 INTEGER, INTENT(IN OUT)                  :: tabnol(1)
 INTEGER, INTENT(IN)                      :: list(1)
 LOGICAL :: part1
 INTEGER :: dit     , iary(8),tabno  ,tabtyp ,tabido ,NAME(2) ,  &
     clsrw  ,tabid  ,offset ,sctyp
 REAL :: y(2)   , z(1)   , px(2,2)
 COMPLEX :: sum    ,a      ,b      ,term
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /condas/ pi     ,twopi  ,radeg  ,degra  ,s4pisq
 COMMON /system/ ibuf   ,nout
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA     clsrw, neor  ,NAME         ,px              , lnth  /  &
     1    , 0     ,4HPRET,4HAB  ,3.,2.,1.339,1.0 , 12    /
 
!     INITIALIZE
 
 offset = locfx(inz(1)) - locfx(iz(1))
 IF (offset < 0) CALL errtrc ('pretab  ',5)
 dit  = ditf
 idic = 0 + offset
 part1= .true.
 lim  = tabnol(1)
 icrq = lnth*lim - lcrgvn
 IF (icrq >= 0) GO TO 1080
 
!     SET UP TABLE NUMBERS IN DICTIONARY
 
!     FOR EACH TABLE THE DICTIONARY ENTRY IS AS FOLLOWS -
 
!       LOC.  1      TABLE NUMBER
!       LOC.  2      TABLE TYPE(1,2,3, 4, OR 5)
!       LOC.  3      POINTER TO 1ST  ENTRY IN TABLE.
!       LOC.  4      POINTER TO LAST ENTRY IN TABLE.
!       LOC.  5      *
!       LOC.  6      *
!       LOC.  7      *
!       LOC.  8      * LOCATIONS 5 THRU 11 CONTAIN TABLE PARAMETERS.
!       LOC.  9      *
!       LOC. 10      *
!       LOC. 11      *
!       LOC. 12      SCALE TYPE - LINEAR-LINER(0), LOG-LOG(1), LINEAR-
!                                 LOG(2), LOG-LINEAR(3)
 
 DO  i = 1,lim
   iz(idic+1) = tabnol(i+1)
   jlow = idic + 2
   jlim = idic + lnth
   DO  j = jlow,jlim
     iz(j) = 0
   END DO
   idic  = idic + lnth
 END DO
 idicl = 1 + offset
 idich = idic
 
!     READ THE CARDS REFERENCED VIA THE TABNOL AND LIST ARRAY.
 
 itable = idic
 CALL preloc (*1010,buf,dit)
 limjj = tabnol(1)
 lim   = list(1)
 jj    = 1
 30 jj3   = 3*jj - 1
 CALL locate (*110,buf,list(jj3),flag)
 
!     READ 8 WORDS INTO THE ARRAY IARY
 
 40 CALL READ (*1020,*110,dit,iary,8,neor,flag)
 tabno = iary(1)
 sctyp = iary(8)
 
!     DETERMINE IF THIS TABLE NUMBER IS IN THE USER SUPPLIED LIST OF
!     TABLE NUMBERS
 
 DO  j = 1,limjj
   IF (tabno == IABS(tabnol(j+1))) IF (tabno-tabnol(j+1)) 60,70,60
 END DO
 
!     THIS TABLE IS NOT CALLED FOR.  READ THE TABLE SERIALLY UNTIL AN
!     END OF TABLE INDICATOR (TWO MINUS ONES FOR TABLE TYPES 1,2,3 AND
!     ONE MINUS ONE FOR TABLE TYPE 4
 
 nwds = 2
 IF (list(3*jj+1) == 4) nwds = 1
 55 CALL READ (*1020,*1040,dit,iary(2),nwds,neor,iflag)
 IF (iary(2) == -1) GO TO 40
 GO TO 55
 
!     THERE ARE TWO DIFFERENT TABLES WITH THE SAME NUMBER -- FATAL ERROR
 
 60 iary(1) = tabno
 iary(2) = list(3*jj-1)
 CALL mesage (-30,88,iary)
 
!     THIS IS A NEW TABLE.  SET TABLE NUMBER NEGATIVE AND DEFINE WORDS
!     2 AND 3 OF THE PROPER DICTIONARY ENTRY.
 
 70 tabtyp = list(3*jj+1)
 tabnol(j+1) = -tabnol(j+1)
 INDEX  = lnth*(j-1) + offset
 iz(INDEX+2) = tabtyp
 iz(INDEX+3) = itable + 1
 
!     READ THE TABLE INTO CORE.
 
 nwdsrd = 2
 IF (tabtyp == 4) nwdsrd = 1
 ii = itable + 1
 80 CALL READ (*1020,*1040,dit,z(ii),nwdsrd,neor,flag)
 IF (iz(ii) == -1) GO TO 90
 ii   = ii + nwdsrd
 icrq = ii - lcrgvn - offset
 IF (icrq >= 0) GO TO 1080
 GO TO 80
 
!     STORE THE LAST LOCATION OF THE TABLE IN IZ(INDEX+4)
 
 90 iz(INDEX+4) = ii - nwdsrd
 
!     STORE THE PARAMETERS ON THE TABLE CARD IN WORDS 5 THRU 11 OF THE
!     PROPER DICTIONARY ENTRY.
 
 lx = INDEX + 4
 DO  k = 2,8
   lx = lx + 1
   iz(lx) = iary(k)
 END DO
 iz(lx+1) = sctyp
 
!     STORE THE CORRECT 0TH ADDRESS OF THE NEXT TABLE IN ITABLE
 
 itable = iz(INDEX+4)
 
!     IF THE TABLE IS A POLYNOMIAL EVALUATE THE END POINTS.
 
 IF (tabtyp /= 4) GO TO 108
 l  = INDEX + 1
 xx = (z(l+6) - z(l+4))/z(l+5)
 ASSIGN 470 TO igoto
 GO TO  440
 102 ASSIGN 480 TO igoto
 xx = (z(l+7) - z(l+4))/z(l+5)
 GO TO 440
 108 itable = itable + 1
 GO TO 40
 
!     TEST TO SEE IF ALL OF THE REQUESTED TABLES HAVE BEEN FOUND. IF
!     ALL TABLES HAVE NOT BEEN FOUND, GO TO NEXT TRIPLE IN LIST ARRAY
 
 110 IF (jj >= lim) GO TO 120
 DO  i = 1,limjj
   IF (tabnol(i+1) > 0) GO TO 117
 END DO
 GO TO 120
 117 jj = jj + 1
 GO TO 30
 
!     SET ALL ENTRIES IN TABNOL BACK TO THEIR ORIGINAL POSITIVE STATUS.
!     IF AN ENTRY IS STILL POSITIVE, THIS IMPLIES THE TABLE WAS NOT
!     FOUND IN THE DIT AND A FATAL ERROR CONDITION EXISTS.
 
 120 iflag = 0
 DO  i = 1,limjj
   IF (tabnol(i+1) <= 0) GO TO 130
   CALL mesage (30,89,tabnol(i+1))
   iflag = 1
   CYCLE
   130 tabnol(i+1) = -tabnol(i+1)
 END DO
 IF (iflag /= 0) CALL mesage (-37,0,NAME)
 
!     WRAP-UP PRETAB
 
 CALL CLOSE (dit,clsrw)
 part1  = .false.
 tabido = -1
 xo     = -10.0E+37
 lcused = itable + 1 - offset
 icheck = 123456789
 RETURN
 
!     ENTRY TAB COMPUTES THE FUNCTIONAL VALUE Y AT THE ABSCISSA X FOR
!     THE FUNCTION DEFINED BY THE TABLE WHOSE NUMBER IS TABID
 
 
 ENTRY tab (tabid,x,y)
!     =====================
 
 IF (icheck /= 123456789) CALL errtrc ('pretab  ',200)
 ASSIGN 251 TO ihop
 
 IF (tabid == tabido .AND. x == xo) GO TO 210
 tabido = tabid
 xo     = x
 GO TO 220
 210 y(1) = yo
 RETURN
 220 IF (tabid /= 0) GO TO 230
 y(1) = 0.0
 yo   = 0.0
 RETURN
 
!     SEARCH THE TABLE DICTIONARY TO FIND THE TABLE NUMBER
 
 230 DO  ii = idicl,idich,lnth
   IF (tabid == iz(ii)) GO TO 250
 END DO
 
!     TABID COULD NOT BE FOUND IN THE DICTIONARY - FATAL ERROR
 
 CALL mesage (-30,90,tabid)
 250 l = ii
 itype = iz(l+ 1)
 sctyp = iz(l+11) + 1
 GO TO ihop, (251,501)
 251 CONTINUE
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 260
   CASE (    2)
     GO TO 270
   CASE (    3)
     GO TO 280
   CASE (    4)
     GO TO 290
   CASE (    5)
     GO TO 295
 END SELECT
 
!     TABLE TYPE = 1
 
!     A  RGUMENT = X
 
 260 xx = x
 GO TO 300
 
!     TABLE TYPE = 2
 
!     ARGUMENT = (X - X1)
 
 270 xx = x - z(l+4)
 GO TO 300
 
!     TABLE TYPE = 3
 
!     ARGUMENT = (X - X1)/X2
 
 280 xx = (x - z(l+4))/z(l+5)
 GO TO 300
 
!     TABLE TYPE = 4
 
!     ARGUMENT = (X - X1)/X2
 
 290 xx = (x - z(l+4))/z(l+5)
 GO TO 400
 
!     TABLE TYPE = 5
 
!     TABRNDG CARD FUNTION ONLY
 
 295 CONTINUE
 
!     PICK UP TYPE
 
 lx = iz(l+4)
 
!     P US ONE OVER TERM IN PX TABLE BASED ONL TYPE
 
 p = 1./px(lx,1)
 
!     CONPUTE K SQUARED FROM PX TABLE
 
 xksq = px(lx,2)*px(lx,2)
 
!     RETRIEVE LU (L/U) FROM TABLE PARAMS
 
 xlu = z(l+5)
 xx  = 2.*z(l+6)**2*xlu
 xlu = xlu*xlu
 wsq = s4pisq*xo*xo
 tr  = xksq*xlu*wsq
 prop= xx*(1.+2.*(p+1.)*tr)/(1.+tr)**(p+1.5)
 GO TO 500
 
!     ROUTINE TO PERFORM LINEAR INTERPOLATION FOR FUNCTION IN A TABLE.
!     L POINTS TO THE ENTRY IN THE TABLE DICTIONARY WHICH DEFINES THE
!     TABLE. THE ARGUMENT IS XX. THE FUNCTIONAL VALUE IS STORED IN PROP.
!     EXTRAPOLATION IS MADE IF XX IS OUTSIDE THE LIMITS OF THE TABLE.
!     HENCE THERE ARE NO ERROR RETURNS.
!     HOWEVER, IF FUNCTION OVERFLOWED ON EXTRAPOLATION OUTSIDE TABLE
!     LIMITS, A FATAL MESSAGE IS ISSUED.
 
 300 itabl = iz(l+2)
 ntabl = iz(l+3)
 up    = 1.0
 IF (z(itabl) > z(itabl+2)) up = -1.0
 kxx1  = itabl
 IF ((xx - z(itabl))*up <= 0.0) GO TO 350
 kxx1  = ntabl - 2
 IF ((xx - z(ntabl))*up >= 0.0) GO TO 350
 klo = 1
 khi = (ntabl - itabl)/2  +  1
 310 kx  = (klo + khi + 1)/2
 kxx = (kx - 1)*2 + itabl
 IF ((xx - z(kxx))*up < 0.0) THEN
   GO TO   320
 ELSE IF ((xx - z(kxx))*up == 0.0) THEN
   GO TO   370
 ELSE
   GO TO   330
 END IF
 320 khi = kx
 GO TO 340
 330 klo = kx
 340 IF (khi-klo /= 1) GO TO 310
 kxx1 = (klo - 1)*2 + itabl
 IF (kxx ==      kxx1) GO TO 350
 IF (xx  == z(kxx1+2)) GO TO 360
 350 SELECT CASE ( sctyp )
   CASE (    1)
     GO TO 355
   CASE (    2)
     GO TO 351
   CASE (    3)
     GO TO 352
   CASE (    4)
     GO TO 353
 END SELECT
 351 CALL loglog (z(kxx1),z(kxx1+1),z(kxx1+2),z(kxx1+3),xx,prop)
 GO TO 500
 352 CALL smilog (z(kxx1),z(kxx1+1),z(kxx1+2),z(kxx1+3),xx,prop)
 GO TO 500
 353 CALL logsmi (z(kxx1),z(kxx1+1),z(kxx1+2),z(kxx1+3),xx,prop)
 GO TO 500
 355 prop = (xx - z(kxx1))*(z(kxx1+3) - z(kxx1+1))/(z(kxx1+2)  &
     - z(kxx1)) + z(kxx1+1)
 IF (ABS(prop) < 1.0E-36) prop = 0.0
 IF (ABS(prop) < 1.0E+36) GO TO 500
 IF (up > 0. .AND. (xx < z(itabl) .OR. xx > z(ntabl)))GO TO 1050
 IF (up < 0. .AND. (xx > z(itabl) .OR. xx < z(ntabl)))GO TO 1050
 GO TO 500
 360 kxx = kxx1 + 2
 370 IF (xx == z(kxx-2)) GO TO 380
 IF (xx == z(kxx+2)) GO TO 390
 prop = z(kxx+1)
 GO TO 500
 380 prop = (z(kxx-1) + z(kxx+1))/2.0
 GO TO 500
 390 prop = (z(kxx+1) + z(kxx+3))/2.0
 GO TO 500
 
!     POLYNOMIAL EVALUATION
 
 400 IF (xx - (z(l+6) - z(l+4))/z(l+5) > 0.0) THEN
   GO TO   420
 END IF
 410 prop = z(l+8)
 GO TO 500
 420 IF (xx - (z(l+7) - z(l+4))/z(l+5) < 0.0) THEN
   GO TO   440
 END IF
 430 prop = z(l+9)
 GO TO 500
 440 nn   = iz(l+3)
 prop = z(nn)
 450 IF (nn <= iz(l+2)) GO TO 460
 prop = prop*xx + z(nn-1)
 nn   = nn - 1
 GO TO 450
 460 IF (part1) GO TO igoto, (470,480)
 GO TO 500
 470 z(l+8) = prop
 GO TO 102
 480 z(l+9) = prop
 GO TO 40
 
!     TAB WRAP-UP
 
 500 y(1) = prop
 yo   = y(1)
 RETURN
 
 
 ENTRY tab1 (tabid,x,y)
!     ======================
 
!     ENRTY FOR TABLE TRANSFORM
 
 ASSIGN 501 TO ihop
 GO TO 220
 501 CONTINUE
 
!     L  POINTS  TO TABLE
!     ITYPE IS THE TABLE TYPE
 
 itabl = iz(l+2)
 ntabl = iz(l+3)
 omega = twopi*x
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 510
   CASE (    2)
     GO TO 520
   CASE (    3)
     GO TO 530
   CASE (    4)
     GO TO 540
 END SELECT
 
!     TABLED1
 
 510 CONTINUE
 x1 = 0.0
 x2 = 1.0
 GO TO 550
 
!     TABLED2
 
 520 CONTINUE
 x1 = z(l+4)
 x2 = 1.0
 GO TO 550
 
!     TABLED3
 
 530 CONTINUE
 x1 = z(l+4)
 x2 = z(l+5)
 GO TO 550
 
!     TABLED4
 
 540 CONTINUE
 
!     EVALUATE SUM
 
 550 CONTINUE
 sum = CMPLX(0.0,0.0)
 k   = itabl
 551 CONTINUE
 yi   = z(k+1)
 xi   = z(k)
 yip1 = z(k+3)
 xip1 = z(k+2)
 omegax = omega*x2*(xip1-xi)
 CALL ifte2 (omegax,rp,cp)
 p    =-omega*(x1 + x2*xip1)
 a    = CMPLX(0.,p)
 b    = CMPLX(rp,cp)
 term = CEXP(a)*b*yip1
 p    =-omega*(x1 + x2*xi)
 a    = CMPLX(0.,p)
 b    = CMPLX(rp,-cp)
 term = term + CEXP(a)*b*yi
 term = term*(xip1- xi)*.5
 sum  = sum  + term
 k    = k + 2
 IF (k < ntabl) GO TO 551
 
!     FINISH FUNCTION
 
 sum  = sum*x2
 y(1) = REAL(sum)
 y(2) = AIMAG(sum)
 RETURN
 
!     FATAL ERROR MESSAGES
 
 1010 mn = -1
 GO TO 1100
 1020 mn = -2
 GO TO 1100
 1040 mn = -3
 GO TO 1100
 1050 WRITE  (nout,1055) ufm,iz(l)
 1055 FORMAT (a23,' 3308, TABLE',i9,' INTERPOLATION ERROR', /5X,  &
     'FUNCTION OVERFLOWED WHEN EXTRAPOLATION WAS MADE OUTSIDE ',  &
     'TABLE GIVEN RANGE.')
 mn = -37
 GO TO 1100
 1080 mn = -8
 dit= icrq
 1100 CALL mesage (mn,dit,NAME)
 RETURN
END SUBROUTINE pretab
