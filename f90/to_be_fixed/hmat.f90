SUBROUTINE hmat (id)
     
!     MAT ROUTINE FOR USE IN -HEAT- FORMULATIONS ONLY.
 
!         CALL PREHMA (Z)  SETUP CALL MADE BY SMA1A, EMGTAB, ETC.
 
!         CALL HMAT (ELID)  ELEMENT ROUTINE CALLS
 
 
!     REVISED BY G.CHAN/UNISYS
!     5/90 - THE THERMAL CONDUCTIVITY OR CONVECTIVE FILM COEFFICIENT K,
!            IS TIME DEPENDENT IF MATT4 REFERS TO TABLEM5. TIME STEP IS
!            DEFINED VIA TSTEP IN /HMATDD/. IF TIME STEP IS NOT USED,
!            TSTEP SHOULD BE -999.
!            (TSTEP IS INITIALIZED TO -999. WHEN PREHMA IS CALLED)
!     7/92 - NEW REFERENCE TO OPEN CORE ARRAY SUCH THAT THE SOURCE CODE
!            IS UP TO ANSI FORTRAN 77 STANDARD.
 
 
 INTEGER, INTENT(IN OUT)                  :: id
 LOGICAL :: any4    ,any5    ,anyt4   ,anyt5   ,linear  , anytab
 INTEGER :: NAME(2) ,sysbuf  ,outpt   ,flag    ,core    ,  &
     dit     ,oldmid  ,oldflg  ,clsrew  ,cls     ,  &
     TYPE    ,mat4(2) ,mat5(2) ,matt4(2),matt5(2), tset    ,offset  ,tablst(16)
 REAL :: card(10),rz(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm     ,uwm     ,uim     ,sfm
 COMMON /zzzzzz/ z(1)
 COMMON /matin / matid   ,inflag  ,eltemp  ,dum(1)  ,s       ,c
 COMMON /hmtout/ buf(7)
 COMMON /names / rd      ,rdrew   ,wrt     ,wrtrew  ,clsrew  ,cls
 COMMON /hmatdd/ ihmatx  ,nhmatx  ,mpt     ,dit     ,linear  , anytab  ,tstep
 COMMON /system/ ksystm(65)
 EQUIVALENCE     (ksystm( 1),sysbuf) ,(ksystm( 2),outpt ) ,  &
     (ksystm(10),tset  ) ,(ksystm(56),itherm) , (f4,n4)    ,(f5,n5)
 DATA    NAME  / 4HHMAT,4H    /, noeor  / 0 /
 DATA    mat4  / 2103  ,21    /
 DATA    mat5  / 2203  ,22    /
 DATA    matt4 / 2303  ,23    /
 DATA    matt5 / 2403  ,24    /
 DATA    tablst/ 5, 105,1,1,  205,2,2, 305,3,3, 405,4,4, 505,5,5 /
 
!     CALL BUG (4HHMTI,0,MATID,6)
 
 IF (icheck /= 123456789) CALL errtrc ('hmat    ',0)
 GO TO 200
 
 
 ENTRY prehma (rz)
!     =================
 
 IF (itherm > 0) THEN
   GO TO    10
 ELSE
   GO TO   500
 END IF
 10 icheck = 123456789
 offset = locfx(rz(1)) - locfx(z(1))
 IF (offset < 0) CALL errtrc ('hmat    ',10)
 tstep  = -999.
 ihmat  = ihmatx + offset
 nhmat  = nhmatx + offset
 lbuf   = nhmat  - sysbuf
 core   = lbuf   - ihmat
 IF (core < 10) CALL mesage (-8,0,NAME)
 CALL preloc (*125,z(lbuf),mpt)
 
!     LOCATE MAT4 CARDS AND BLAST THEM INTO CORE.
 
 anytab =.false.
 any4   =.false.
 imat4  = ihmat + 1
 nmat4  = ihmat
 CALL locate (*40,z(lbuf),mat4,flag)
 CALL READ (*480,*30,mpt,z(imat4),core,noeor,iwords)
 CALL mesage (-8,0,NAME)
 30 nmat4 =  nmat4 + iwords
 40 mat4s = (nmat4 - imat4 + 1)/3
 IF (mat4s > 0) any4 = .true.
 
!     LOCATE MATT4 CARDS AND BLAST THEM INTO CORE IF THERE WERE ANY MAT4
 
 anyt4  =.false.
 imatt4 = nmat4 + 1
 nmatt4 = nmat4
 IF (.NOT.any4 .OR. (tset == 0 .AND. tstep < 0.)) GO TO 60
 CALL locate (*60,z(lbuf),matt4,flag)
 CALL READ (*480,*50,mpt,z(imatt4),core,noeor,iwords)
 CALL mesage (-8,0,NAME)
 50 nmatt4 =  nmatt4 + iwords
 60 matt4s = (nmatt4 - imatt4 + 1)/2
 IF (matt4s > 0) anyt4 = .true.
 
!     LOCATE MAT5 CARDS AND BLAST THEM INTO CORE.
 
 any5  =.false.
 imat5 = nmatt4 + 1
 nmat5 = nmatt4
 CALL locate (*80,z(lbuf),mat5,flag)
 CALL READ (*480,*70,mpt,z(imat5),core,noeor,iwords)
 CALL mesage (-8,0,NAME)
 70 nmat5 =  nmat5 + iwords
 80 mat5s = (nmat5 - imat5 + 1)/8
 IF (mat5s > 0) any5 = .true.
 
!     LOCATE MATT5 CARDS AND BLAST THEM INTO CORE IF THERE WERE ANY MAT5
 
 anyt5  =.false.
 imatt5 = nmat5 + 1
 nmatt5 = nmat5
 IF (.NOT.any5 .OR. (tset == 0 .AND. tstep < 0.)) GO TO 100
 CALL locate (*100,z(lbuf),matt5,flag)
 CALL READ (*480,*90,mpt,z(imatt5),core,noeor,iwords)
 CALL mesage (-8,0,NAME)
 90 nmatt5 =  nmatt5 + iwords
 100 matt5s = (nmatt5 - imatt5 + 1)/7
 IF (matt5s > 0) anyt5 = .true.
 CALL CLOSE (mpt,clsrew)
 
!     IF A TEMPERATURE SET IS SPECIFIED -DIT- IS NOW READ INTO CORE,
!     PROVIDING ANY MATT4 OR MATT5 CARDS WERE PLACED INTO CORE.
 
 IF ((tset == 0 .AND. tstep < 0.) .OR.  &
     (.NOT.anyt4 .AND. .NOT.anyt5)) GO TO 130
 
!     BUILD LIST OF TABLE NUMBERS POSSIBLE FOR REFERENCE
 
 kk = 0
 itabno = nmatt5 + 1
 ntabno = itabno
 
 IF (matt4s <= 0) GO TO 110
 loop108:  DO  i = imatt4,nmatt4,2
   f4 = z(i+1)
   IF (n4 > 0) THEN
     GO TO   102
   ELSE
     GO TO   108
   END IF
   102 IF (kk > 0) THEN
     GO TO   103
   ELSE
     GO TO   107
   END IF
   103 DO  j = itabno,ntabno
     f5 = z(j)
     IF (n4 == n5) CYCLE loop108
   END DO
   
!     ADD NEW TABLE ID TO LIST
   
   107 ntabno = ntabno + 1
   z(ntabno) = z(i+1)
   kk = 1
 END DO loop108
 
 110 IF (matt5s <= 0) GO TO 120
 DO  i = imatt5,nmatt5,7
   j1 = i + 1
   j2 = i + 6
   loop117:  DO  j = j1,j2
     f4 = z(j)
     IF (n4 > 0) THEN
       GO TO   111
     ELSE
       GO TO   117
     END IF
     111 IF (kk > 0) THEN
       GO TO   113
     ELSE
       GO TO   115
     END IF
     113 DO  k = itabno,ntabno
       f5 = z(k)
       IF (n4 == n5) CYCLE loop117
     END DO
     
!     ADD NEW TABLE ID TO LIST
     
     115 ntabno = ntabno + 1
     z(ntabno) = z(j)
     kk = 1
   END DO loop117
 END DO
 
 120 n4 = ntabno - itabno
 z(itabno) = f4
 
!     CALL BUG (4HTABL,120,Z(ITABNO),NTABNO-ITABNO+1)
 
 IF (n4 > 0) THEN
   GO TO   122
 ELSE
   GO TO   130
 END IF
 122 CALL sort (0,0,1,1,z(itabno+1),n4)
 
!     OK READ IN DIRECT-INPUT-TABLE (DIT)
 
 idit  = ntabno + 1
 igbuf = nhmat - sysbuf - 2
 lz    = igbuf - idit  - 1
 IF (lz < 10) CALL mesage (-8,0,NAME)
 CALL pretab (dit,z(idit),z(idit),z(igbuf),lz,lused,z(itabno), tablst)
 ndit  = idit + lused
 nhmat = ndit + 1
 
!     CALL BUG (4HDITS,123,Z(IDIT),NDIT-IDIT+1)
 
 GO TO 140
 
!     WRAP UP THE PRE-HMAT SECTION
 
 125 nhmat = ihmat - 1
 any4  =.false.
 any5  =.false.
 GO TO 140
 130 nhmat  = nmatt5
 140 oldmid = 0
 oldflg = 0
 oldsin = 0.0
 oldcos = 0.0
 oldtem = 0.0
 oldstp = 0.0
 s      = 0.0
 c      = 0.0
 dum(1) = 0.0
 eltemp = 0.0
 nhmatx = nhmat - offset
 
!     CHECK FOR DUPLICATE MATID-S ON BOTH MAT4 AND MAT5 CARDS.
 
 IF (.NOT.any4 .OR. .NOT.any5) GO TO 490
 j4 = imat4
 j5 = imat5
 f4 = z(j4)
 f5 = z(j5)
 150 IF (n4 - n5 < 0) THEN
   GO TO   160
 ELSE IF (n4 - n5 == 0) THEN
   GO TO   180
 ELSE
   GO TO   170
 END IF
 
!     MAT4 ID IS LESS THAN MAT5 ID
 
 160 j4 = j4 + 3
 IF (j4 > nmat4) GO TO 490
 f4 = z(j4)
 GO TO 150
 
!     MAT5 ID IS LESS THAN MAT4 ID.
 
 170 j5 = j5 + 8
 IF (j5 > nmat5) GO TO 490
 f5 = z(j5)
 GO TO 150
 
!     ID OF MAT4 IS SAME AS THAT OF MAT5
 
 180 WRITE  (outpt,190) uwm,n4
 190 FORMAT (a25,' 2155, MAT4 AND MAT5 MATERIAL DATA CARDS HAVE SAME ',  &
     'ID =',i14, /5X,'MAT4 DATA WILL BE SUPPLIED WHEN CALLED ', 'FOR THIS ID.')
 GO TO 170
 
!                 DATA RETURNED IF MAT-ID       DATA RETURNED IF MAT-ID
!     INFLAG      IS ON A MAT4 CARD             IS ON A MAT5 CARD.
!     =================================================================
 
!       1           1- K                               1- KXX
!                   2- CP                              2- CP
 
!       2           1- K                               1- KXXB
!                   2- 0.0                             2- KXYB
!                   3- K                               3- KYYB
!                   4- CP                              4- CP
 
!       3           1- K                               1- KXX
!                   2- 0.0                             2- KXY
!                   3- 0.0                             3- KXZ
!                   4- K                               4- KYY
!                   5- 0.0                             5- KYZ
!                   6- K                               6- KZZ
!                   7- CP                              7- CP
 
!       4           1- CP                              1- CP
 
 
 
 
!     DATA LOOK UP SECTION.  FIND MAT-ID IN CARD IMAGES.
 
 
 200 IF (inflag - oldflg == 0) THEN
   GO TO   210
 ELSE
   GO TO   260
 END IF
 210 IF (matid  - oldmid == 0) THEN
   GO TO   220
 ELSE
   GO TO   260
 END IF
 220 IF (eltemp - oldtem == 0.0) THEN
   GO TO   225
 ELSE
   GO TO   260
 END IF
 225 IF (tstep  - oldstp == 0.0) THEN
   GO TO   230
 ELSE
   GO TO   260
 END IF
 230 IF (TYPE  ==    4) GO TO 250
 IF (s      - oldsin == 0.0) THEN
   GO TO   240
 ELSE
   GO TO   260
 END IF
 240 IF (c      - oldcos == 0.0) THEN
   GO TO   250
 ELSE
   GO TO   260
 END IF
 
!     ALL INPUTS SEEM TO BE SAME THUS RETURN IS MADE.
 
 250 GO TO 490
 
!     FIND POINTER TO SECOND WORD OF CARD IMAGE WITH MAT-ID DESIRED.
!     AMONG EITHER MAT4S OR MAT5S.
 
 260 oldflg = inflag
 oldmid = matid
 oldcos = c
 oldsin = s
 oldtem = eltemp
 oldstp = tstep
 linear = .true.
 IF (.NOT.any4) GO TO 270
 CALL bisloc (*270,matid,z(imat4),3,mat4s,jpoint)
 j = imat4 + jpoint
 TYPE = 4
 GO TO 280
 270 IF (.NOT. any5) GO TO 460
 CALL bisloc (*460,matid,z(imat5),8,mat5s,jpoint)
 j = imat5 + jpoint
 TYPE = 5
 
!     IF A THERMAL SET IS REQUESTED (TSET.NE.0) THEN A FACTOR, WHICH IS
!     A FUNCTION OF THE AVERAGE ELEMENT TEMPERATURE, ELTEMP, (OR TIME
!     STEP, TSTEP) AND THE TABULATED VALUE IN TABLEMI, IS USED AS A
!     MULTIPLIER TO THE K-TERMS IN MAT4 OR MATT5
 
!     IF THE MATERIAL ID IS FOUND ON A -MAT4- AN ATTEMPT IS MADE TO FIND
!     A CORRESPONDING -MATT4- CARD.  LIKEWISE THIS IS DONE IF THE
!     MATERIAL ID IS FOUND ON A -MAT5- CARD WITH RESPECT TO A -MATT5-
!     CARD. IF THE -MAT4- OR -MAT5- HAS A RESPECTIVE -MATT4- OR -MATT5-
!     CARD, THEN THE THERMAL CONDUCTIVITY OR THE CONVECTIVE FILM COEFF.
!     K, IS TEMPERATURE DEPENDENT IF TABLEM1, TABLEM2, TABLEM3 AND
!     TABLEM4 ARE REFERENECED. K IS TIME DEPENDENT IF TABLEM5 IS USED.
!     THE K-TERMS OF THE -MAT4- OR -MAT5- CARDS WILL BE MODIFIED BY
!     USING -ELTEMP- AND THE -DIT- AS SPECIFIED IN THE RESPECTIVE FIELDS
!     OF THE RESPECTIVE -MATT4- OR -MATT5- CARD.  A ZERO T(K) IN A
!     PARTICULAR FIELD OF THE RESPECTIVE -MATT4- OR -MATT5- CARD IMPLIES
!     NO TEMPERATURE DEPENDENCE FOR THAT RESPECTIVE K VALUE.
!     -DIT- TABLES TABLEM1, TABLEM2, TABLEM3, TABLEM4 AND TABLEM5 MAY BE
!     USED.
 
 
!     MOVE MAT CARD INTO SPECIAL BUF WHERE IT CAN BE MODIFIED IF
!     NECESSARY
 
 280 DO  i = 1,10
   card(i) = z(j)
   j = j + 1
 END DO
 
!     CHECK FOR EXISTENCE OF A THERMAL SET REQUEST OR TIME STEP.
 
 IF ((tset == 0 .AND. tstep < 0.) .OR. inflag == 4) GO TO 350
 
!     IF -MAT4- CARD, FIND THE -MATT4- CARD
!     (IF NO MATT4 ASSUME NO TEMPERATURE OR TIME DEPENDENCE)
 
 IF (TYPE == 5) GO TO 300
 iwords = 2
 imat   = imatt4
 mats   = matt4s
 GO TO 310
 300 iwords = 7
 imat   = imatt5
 mats   = matt5s
 310 IF (mats > 0) THEN
   GO TO   315
 ELSE
   GO TO   350
 END IF
 315 CALL bisloc (*350,matid,z(imat),iwords,mats,jpoint)
 itemp =  imat + jpoint
 ntemp = itemp + iwords - 2
 
!     Z(I) FIELDS SPECIFYING A NON-ZERO TABLE IMPLY TEMPERATURE (OR
!     TIME) DEPENDENCE ON CORRESPONDING FIELDS OF THE MAT4 OR MAT5
!     STORED IN THE ARRAY -CARD-.
 
 kk = 0
 DO  i = itemp,ntemp
   kk = kk + 1
   f4 = z(i)
   IF (n4 > 0) THEN
     GO TO   320
   ELSE
     GO TO   340
   END IF
   
!     OK TEMPERATURE (OR TIME) DEPENDENCE.
   
   320 IF (tset  >  0) x = eltemp
   IF (tstep >= 0.) x = tstep
   CALL tab (n4,x,factor)
   card(kk) = card(kk)*factor
   linear = .false.
 END DO
 
!     BRANCH ON INFLAG.
 
 350 IF (inflag < 1 .OR. inflag > 4) GO TO 440
 SELECT CASE ( inflag )
   CASE (    1)
     GO TO 360
   CASE (    2)
     GO TO 380
   CASE (    3)
     GO TO 400
   CASE (    4)
     GO TO 420
 END SELECT
 
!     INFLAG = 1
 
 360 IF (TYPE == 5) GO TO 370
 buf(1) = card(1)
 buf(2) = card(2)
 GO TO 490
 370 buf(1) = card(1)
 buf(2) = card(7)
 GO TO 490
 
!     INFLAG = 2
 
 380 IF (TYPE == 5) GO TO 390
 buf(1) = card(1)
 buf(2) = 0.0
 buf(3) = buf(1)
 buf(4) = card(2)
 GO TO 490
 390 csq = c*c
 ssq = s*s
 cs  = c*s
 cs2kxy = cs *2.0*card(2)
 buf(1) = csq* card(1) - cs2kxy   +  ssq*card(4)
 buf(2) = cs *(card(1) - card(4)) + (csq - ssq)*card(2)
 buf(3) = ssq* card(1) + cs2kxy   +  csq*card(4)
 buf(4) = card(7)
 GO TO 490
 
!     INFLAG = 3
 
 400 IF (TYPE == 5) GO TO 410
 buf(1) = card(1)
 buf(2) = 0.0
 buf(3) = 0.0
 buf(4) = buf(1)
 buf(5) = 0.0
 buf(6) = buf(1)
 buf(7) = card(2)
 GO TO 490
 410 buf(1) = card(1)
 buf(2) = card(2)
 buf(3) = card(3)
 buf(4) = card(4)
 buf(5) = card(5)
 buf(6) = card(6)
 buf(7) = card(7)
 GO TO 490
 
!     INFLAG = 4.  RETURN ONLY CP.
 
 420 IF (TYPE == 5) GO TO 430
 buf(1) = card(2)
 GO TO 490
 430 buf(1) = card(7)
 GO TO 490
 
!     ERROR CONDITIONS
 
 440 WRITE  (outpt,450) sfm,inflag
 450 FORMAT (a25,' 2156, ILLEGAL INFLAG =',i14,' RECEIVED BY HMAT.')
 GO TO 520
 460 WRITE  (outpt,470) ufm,matid
 470 FORMAT (a23,' 2157, MATERIAL ID =',i14,  &
     ' DOES NOT APPEAR ON ANY MAT4 OR MAT5 MATERIAL DATA CARD.')
 GO TO 520
 480 CALL mesage (-2,mpt,NAME)
 GO TO 520
 
!     RETURN LOGIC
 
 490 CONTINUE
 
!     CALL BUG (4HHMAT,490,BUF,7)
 
 RETURN
 
!     ERROR - HMAT CALLED IN NON-THERMAL PROBLEM.
 
 500 WRITE  (outpt,510) sfm
 510 FORMAT (a25,' 3062, HMAT MATERIAL ROUTINE CALLED IN A NON-HEAT-',  &
     'TRANSFER PROBLEM.')
 520 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE hmat
