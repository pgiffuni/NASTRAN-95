SUBROUTINE emgcng
     
!     THIS ROUTINE OF THE -EMG- MODULE READS -CNGRNT- CARD
!     IMAGES, IF ANY, FROM GEOM2 AND BUILDS A PAIRED LIST.
 
!     ON EACH -CNGRNT- DATA CARD THE  FIRST ID (NEED NOT BE THE SMALLEST
!     ID) BECOMES THE PRIMARY ID.  THIS ID WILL BE PAIRED WITH A ZERO
!     NOW AND A NEGATIVE DICTIONARY-TABLE  ADDRESS LATER.  AS SOME OF
!     THE ID-S APPEARING ON THE -CNGRNT- DATA CARD MAY NOT EVEN BE IN
!     THE PROBLEM, THE FIRST ID OF A CONGRUENT GROUP REFERENCED WILL
!     RESULT IN THE ELEMENT COMPUTATIONS AND THE SETTING OF A DICTIONARY
!     FILE TABLE ADDRESS WITH THE PRIMARY ID.
 
 LOGICAL :: anycon, error, heat
 INTEGER :: z, geom2, sysbuf, buf, subr(2), cngrnt(2), est,  &
     cstm, dit, dictn, rd, wrt, wrtrew, rdrew, cls, clsrew, precis, flag, flags
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm, uwm
 COMMON /system/ ksystm(65)
 COMMON /names / rd, rdrew, wrt, wrtrew, clsrew, cls
 COMMON /emgfil/ est, cstm, mpt, dit, geom2, mats(3), dictn(3)
 COMMON /emgprm/ icore, jcore, ncore, icstm, ncstm, imat, nmat,  &
     ihmat, nhmat, idit, ndit, icong, ncong, lcong,  &
     anycon, flags(3), precis, error, heat, icmbar, lcstm, lmat, lhmat
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm(1), sysbuf), (ksystm(2), nout)
 DATA    subr  / 4HEMGC,  4HNG  /, noeor / 0 /, cngrnt / 5008,50 /
 
 buf = ncore - sysbuf - 2
 IF (buf <= jcore) CALL mesage (-8,jcore-buf,subr)
 anycon= .false.
 icong = jcore
 ncong = jcore - 1
 lcong = 0
 
!     LOCATE -CNGRNT- BULK DATA CARDS IF ANY.
 
 CALL preloc (*90,z(buf),geom2)
 CALL locate (*80,z(buf),cngrnt,flag)
 
!     PROCESS ONE DATA CARD
 
 10 IF (ncong+2 >= buf) GO TO 35
 CALL READ (*40,*40,geom2,z(ncong+1),1,noeor,iwords)
 z(ncong+2) = 0
 idprim = z(ncong+1)
 ncong  = ncong + 2
 
!     READ ANY SECONDARY IDS.
 
 20 IF (ncong+2 >= buf) GO TO 35
 CALL READ (*40,*40,geom2,z(ncong+1),1,noeor,iwords)
 
!     CHECK FOR THE FOLLOWING CONDITION
 
!     CONDITION 1
!     ------------
 
!     A SECONDARY ID ON THIS CARD IS THE SAME AS THE PRIMARY ID
!     ON THIS CARD.  THE SECONDARY ID IS IGNORED AND THE CONDITION
!     IS INDICATED BY A USER INFORMATION MESSAGE.
 
 IF (z(ncong+1) /= idprim) GO TO 25
 
!     THE ABOVE CONDITION EXISTS
 
 CALL page2 (3)
 WRITE (nout,2010) uwm,idprim
 GO TO 20
 
 25 IF (z(ncong+1) < 0.0) THEN
   GO TO    10
 ELSE IF (z(ncong+1) == 0.0) THEN
   GO TO    20
 END IF
 30 z(ncong+2) = idprim
 ncong = ncong + 2
 GO TO 20
 
!     INSUFFICIENT CORE TO PROCESS ALL -CNGRNT- CARDS
 
 35 icrq = ncong + 2 - buf
 CALL page2 (2)
 WRITE (nout,2050) uwm,icrq
 
!     NO MORE -CNGRNT- CARDS
 
 40 lcong = ncong - icong + 1
 IF (lcong <= 0) GO TO 80
 CALL sort (0,0,2,1,z(icong),lcong)
 
!     CHECK FOR THE FOLLOWING ADDITIONAL CONDITIONS
 
!     CONDITION 2
!     -----------
 
!     A PRIMARY ID ON A CNGRNT CARD IS ALSO USED AS A SECONDARY
!     ID ON ANOTHER CNGRNT CARD.  THIS RESULTS IN A USER FATAL
!     MESSAGE.
 
!     CONDITION 3
!     -----------
 
!     A SECONDARY ID IS SPECIFIED AS CONGRUENT TO MORE THAN ONE
!     PRIMARY ID.  THIS ALSO RESULTS IN A USER FATAL MESSAGE.
 
!     CONDITION 4
!     -----------
 
!     A SECONDARY ID IS REDUNDANTLY SPECIFIED.  THE REDUNDANCIES ARE
!     IGNORED AND THE CONDITION IS INDICATED BY A USER INFORMATION
!     MESSAGE.
 
 nogo   = 0
 ncong1 = ncong - 2
 DO  i = icong,ncong1,2
   IF (z(i  ) /= z(i+2)) CYCLE
   IF (z(i+1) == z(i+3)) CYCLE
   nogo = 1
   IF (z(i+1) /= 0 .AND. z(i+3) /= 0) GO TO 420
   
!     THIS IS CONDITION 2 DESCRIBED ABOVE
   
   WRITE (nout,2020) ufm,z(i)
   CYCLE
   
!     THIS IS CONDITION 3 DESCRIBED ABOVE
   
   420 WRITE (nout,2030) ufm,z(i)
   
 END DO
 IF (nogo == 1) CALL mesage (-37,0,subr)
 ncong2 = ncong1
 DO  i = icong,ncong1,2
   IF (z(i) <      0) CYCLE
   IF (z(i) /= z(i+2)) CYCLE
   j = i + 2
   450 DO  k = j,ncong2,2
     z(k  ) = z(k+2)
     z(k+1) = z(k+3)
   END DO
   lcong = lcong - 2
   ncong = ncong - 2
   z(ncong2-1) = -1
   ncong2 = ncong2 - 2
   IF (z(j) == z(i)) GO TO 450
   IF (z(i+1) ==  0) CYCLE
   
!     THIS IS CONDITION 4 DESCRIBED ABOVE
   
   CALL page2 (2)
   WRITE (nout,2040) uwm,z(i)
   
 END DO
 
!     REPLACE PRIMARY ID ASSOCIATED WITH EACH SECONDARY ID
!     WITH LOCATION OF PRIMARY ID IN TABLE.
 
 lnum   = lcong / 2
 icongz = icong - 1
 DO  i = icong,ncong,2
   IF (z(i+1) == 0.0) THEN
     GO TO    60
   END IF
   50 kid = z(i+1)
   CALL bisloc (*60,kid,z(icong),2,lnum,j)
   z(i+1) = icongz + j
   60 CONTINUE
 END DO
 
!     TABLE IS COMPLETE
 
 80 CALL CLOSE (geom2,clsrew)
 IF (ncong > icong) anycon = .true.
 jcore = ncong + 1
 90 RETURN
 
 
 2010 FORMAT (a25,' 3169, PRIMARY ID',i9,' ON A CNGRNT CARD ALSO USED ',  &
     'AS A SECONDARY ID ON THE SAME CARD.', /5X, 'SECONDARY ID IGNORED.')
 2020 FORMAT (a23,' 3170, PRIMARY ID',i9,' ON A CNGRNT CARD ALSO USED ',  &
     'AS A SECONDARY ID ON ANOTHER CNGRNT CARD.')
 2030 FORMAT (a23,' 3171, SECONDARY ID',i9,  &
     ' SPECIFIED AS CONGRUENT TO MORE THAN ONE PRIMARY ID.')
 2040 FORMAT (a25,' 3172, SECONDARY ID',i9,' REDUNDANTLY SPECIFIED ON ',  &
     'CNGRNT CARDS.  REDUNDANCY IGNORED.')
 2050 FORMAT (a25,' 3182, INSUFFICIENT CORE TO PROCESS ALL CNGRNT ',  &
     'CARDS.  ADDITIONAL CORE NEEDED =',i8,7H words.)
 
END SUBROUTINE emgcng
