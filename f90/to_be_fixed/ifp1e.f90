SUBROUTINE ifp1e (isubc,symseq,nwdsc,i81,icaste)
     
!     IFP1E WRITES CASECC OUT FROM CASE
 
 
 INTEGER, INTENT(IN OUT)                  :: isubc(5)
 INTEGER, INTENT(IN OUT)                  :: symseq(1)
 INTEGER, INTENT(IN)                      :: nwdsc
 INTEGER, INTENT(OUT)                     :: i81
 INTEGER, INTENT(OUT)                     :: icaste
 LOGICAL :: bit64
 INTEGER :: case(200,2),BLANK,casecc, core(1),corey(401)
 COMMON /zzzzzz/ corex(1)
 COMMON /xifp1 / BLANK,bit64
 COMMON /ifp1a / scr1,casecc,is,nwpc,ncpw4,nmodes,icc,nset,  &
     nsym,zzzzbb,istr,isub,lencc,iben,equal,IEOR
 EQUIVALENCE    (corex(1),corey(1),case(1,1)),(core(1),corey(401))
 DATA    NONE  / 4HNONE/
 
!     INITIALIZE
 
 
!     SKIP FILTER INTO SUBCASES FOR SYM SUBCASES
 
 DO  i = 1,16
   IF (case(i,2) == 0) case(i,2) = case(i,1)
 END DO
 IF (case(38,2) == 0) case(38,2) = case(38,1)
 IF (nsym > 1 .AND. case(16,2) == 0) GO TO 1140
 DO  i = 1,7
   ik = (i-1)*3 + 17
   IF (case(ik,2) /= 0) GO TO 1125
   DO  j = 1,3
     ii = ik + j - 1
     case(ii,2) = case(ii,1)
   END DO
   1125 iword = case(ik,2)
   IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
   IF (iword == NONE) case(ik,2) = 0
 END DO
 loop1170:  1140 DO  j = 1,3
   DO  i = 1,32
     k = 32*j + i + 6
     iword = case(k,2)
     IF (bit64) CALL mvbits (BLANK,0,32,iword,0)
     IF (iword /= BLANK) CYCLE loop1170
   END DO
   DO  i = 1,32
     k = 32*j + i + 6
     case(k,2) = case(k,1)
   END DO
 END DO loop1170
 j = 129
 DO  i = 1,5
   case(j,2) = isubc(i)
   j = j + 1
 END DO
 DO  i = 135,lencc
   IF (case(i,2) == 0) case(i,2) = case(i,1)
 END DO
!     IMOV = CASE(136,2)*100000000  !! VAX/IBM INTGER OVERFLOW FOR ANOMA
 imov = case(136,2)
 IF (imov < 0) imov = 0
 imov = imov*100000000
 case(136,2) = IABS(case(136,2))
 case(2,2) = case(2,2) + imov
 case(3,2) = case(3,2) + imov
 IF (case(7,2) /= 0) case(7,2) = case(7,2) + imov
 IF (case(8,2) /= 0) case(8,2) = case(8,2) + imov
 icaste = case(8,2)
 DO  iloop = 1,nmodes
   IF (case(1,2) > 99999999) CALL ifp1d (-625)
   
!     CHECK FOR METHOD AND LOAD IN SAME SUBCASE
   
   IF (case(5,2) /= 0 .AND. case(4,2)+case(6,2)+case(7,2) /= 0)  &
       CALL ifp1d (-627)
   IF (case(4,2) == case(6,2) .AND. case(4,2) /= 0 .OR.  &
       case(6,2) == case(7,2) .AND. case(6,2) /= 0 .OR.  &
       case(4,2) == case(7,2) .AND. case(4,2) /= 0) CALL ifp1d (-628)
   CALL WRITE (casecc,case(1,2),lencc,0)
   case(1,2) = case(1,2) + 1
   IF (case(16,2) <= 0) GO TO 1200
   ido = case(lencc,2)
   CALL WRITE (casecc,symseq(1),ido,0)
   1200 IF (nset == 0) GO TO 1220
   ip = nwdsc + 1
   DO  i = 1,nset
     nwor = core(ip)
     CALL WRITE (casecc,core(ip-1),   2,0)
     CALL WRITE (casecc,core(ip+2),nwor,0)
     ip = ip + nwor + 3
   END DO
   1220 CALL WRITE (casecc,core(1),0,1)
 END DO
 nmodes = 1
 IF (nset == 0) GO TO 1270
 
!     REMOVE ALL SETS REFERING TO SUBCASE ONLY
 
 iup  = nwdsc
 ip   = nwdsc
 nset1= nset
 imov = 0
 DO  i = 1,nset
   IF (core(ip+2) /= 1) GO TO 1250
   IF (imov       == 0) GO TO 1240
   ido = core(ip+1) + 3
   DO  j = 1,ido
     ii = iup + j - 1
     ik = ip  + j - 1
     core(ii) = core(ik)
   END DO
   1240 iup = iup+core(ip+1) + 3
   ip  = ip +core(ip+1) + 3
   CYCLE
   1250 imov = 1
   nset1= nset1 - 1
   ip   = ip + core(ip+1) + 3
 END DO
 nset = nset1
 i81  = iup
 1270 CONTINUE
 DO  i = 1,lencc
   case(i,2) = 0
   IF (i > 38 .AND. i < 135) case(i,2) = BLANK
 END DO
 CALL ifp1f (*1281,iword,i2)
 DO  i = 1,5
   isubc(i) = core(i2)
   i2 = i2 + 1
 END DO
 1281 RETURN
END SUBROUTINE ifp1e
