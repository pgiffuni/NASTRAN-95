SUBROUTINE ssg1
     
 INTEGER :: pg,slt,bgpdt,cstm,sil,ecpt,mpt,gptt,edt,casecc,  &
     core(166),sysbuf,iword(4),mcb(7),subnam(2)
 DIMENSION       pg(7),ilist(360),ary(1),defml(2),idefml(2),  &
     iary(1),gvect(1080)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /loadx / lc,slt,bgpdt,OLD,cstm,sil,isil,ecpt,mpt,gptt,edt,  &
     n(3),lodc,mass,nobld,idit,dum6(6)
 COMMON /BLANK / nrowsp,loadnn
 COMMON /system/ sysbuf,iotpe,dum53(53),itherm
 COMMON /loads / nload,iptr
 COMMON /zzzzzz/ icore(1)
 EQUIVALENCE     (core(1),icore(1),iary(1),ary(1)), (defml(1),idefml(1))
 DATA    iword / 4,6,7,162/
 DATA    subnam/ 4HSSG1,4H    /
 
!     MODIFY OPEN CORE POINTER IPTR FOR MAGNETICS PROBLEM
 
 iptr  = MAX0(nrowsp,166)
 mcb(1)= 105
 CALL rdtrl (mcb(1))
 IF (mcb(1) > 0) iptr = MAX0(3*nrowsp,3*mcb(2),166)
 
!     INITIALIZE.
 
 lc    = korsz(icore(1))
 nllst = lc - 2*sysbuf
 slt   = 101
 bgpdt = 102
 cstm  = 103
 sil   = 104
 ecpt  = 105
 mpt   = 106
 gptt  = 107
 edt   = 108
 mass  = 109
 casecc= 110
 idit  = 111
 lodc  = 201
!             205 = NEWSLT (THERMAL)
 pg(1) = 301
 icr2  = 302
 icr3  = 303
 DO  i = 2,7
   pg(i) = 0
 END DO
 pg(3) = nrowsp
 pg(4) = 2
 pg(5) = 1
 
!     AVOID CALCULATING UNUSED LOADS
 
!     NEDT  = NUMBER OF ELEMENT  DEFORMATIONS
!     NTEMP = NUMBER OF THERMAL LOADS
!     NCENT = NUMBER OF CENTRIFUGAL LOADS
 
 CALL ssg1a (n1,ilist(1),nedt,ntemp,ncent,casecc,iharm)
 n1a = n1 + 1
 lc  = lc - sysbuf
 CALL OPEN (*310,pg(1),icore(lc+1),1)
 CALL WRITE (pg(1),pg(1),2,1)
 ngrav = 0
 nex   = n1 + ntemp + nedt + ncent
 IF (n1 == 0) GO TO 21
 
!     MODIFY SLT -QVOL-, -QBDY1-, -QBDY2-, AND -QVECT- CARDS.
 
 newslt = icr3
 IF (itherm /= 0) newslt = 205
 islt = slt
 CALL ssgslt (slt,newslt,ecpt)
 slt = newslt
 CALL extern (nex,ngrav,gvect(1),ilist(1),pg(1),n1,iharm)
 
!     RESET -SLT- TO ORIGINAL SLT DATA BLOCK
 
 slt = islt
 n1  = n1 - ngrav
 21 IF (ntemp == 0) THEN
   GO TO    40
 END IF
 30 CALL templ (ntemp,ilist(n1+1),pg(1))
 n1 = n1 + ntemp
 40 IF (nedt == 0) THEN
   GO TO    60
 END IF
 50 CALL edtl (nedt,ilist(n1+1),pg(1))
 n1 = n1 + nedt
 60 CALL CLOSE (pg,1)
 CALL wrttrl (pg(1))
 IF (ngrav == 0) THEN
   GO TO   100
 END IF
 90 CONTINUE
 
!     CHECK TO SEE IF THE MASS MATRIX IS PURGED
 
 mcb(1) = mass
 CALL rdtrl (mcb(1))
 IF (mcb(1) <= 0) CALL mesage (-56,0,iword)
 CALL gravl1 (ngrav,gvect(1),icr2,iharm)
 
!     USE LOAD FILE AS SCRATCH NOTHING ON IT NOW
 
 CALL ssg2b (mass,icr2,0,icr3,0,1,1,lodc)
 CALL gravl2 (ngrav,icr3,pg(1))
 n1 = n1 + ngrav
 100 ipont1 = iptr + 2
 ipont  = iptr + 1
 nload  = 0
 DO  i = 1,nllst
   iary(i) = 0
 END DO
 CALL OPEN (*320,casecc,icore(lc+1),0)
 lc1  = lc - sysbuf
 islt = 0
 CALL OPEN (*130,slt,icore(lc1+1),0)
 islt = 1
 DO  i = 1,n1a
   CALL fwdrec (*270,slt)
 END DO
 130 DO  i = 1,loadnn
   CALL fwdrec (*320,casecc)
 END DO
 ifrst = 0
 150 CALL READ (*250,*250,casecc,core(1),166,1,flag)
 IF (ifrst /= 0) GO TO 151
 ifrst = 1
 ispcn = core(3)
 mpcn  = core(2)
 151 CONTINUE
 
!     TEST FOR SYMMETRY, BUCKLING OR DIFFERENTIAL STIFFNESS.
 
 IF (core(16) /= 0 .OR. core(5) /= 0 .OR. core(138) /= 0) GO TO 150
 IF (core(3) /= ispcn .OR. core(2) /= mpcn) GO TO 250
 inull = 0
 DO  k = 1,4
   i = iword(k)
   IF (itherm /= 0 .AND. i == 7) CYCLE
   IF (core(i) == 0) CYCLE
   DO  j = 1,n1
     IF (core(i) == ilist(j)) GO TO 220
   END DO
   
!     COMBINATION CARD
   
   inull = 1
   170 CALL READ (*270,*330,slt,idefml(1),2,0,iflag)
   IF (core(i) == idefml(1)) GO TO 190
   180 CALL READ (*270,*330,slt,idefml(1),2,0,iflag)
   IF (idefml(2) == -1) GO TO 170
   GO TO 180
   190 a = defml(2)
   200 CALL READ (*270,*330,slt,idefml(1),2,0,iflag)
   IF (idefml(2) ==  -1) GO TO 210
   IF (ipont+1 > nllst) GO TO 340
   iary(ipont  ) = iary(ipont) + 1
   iary(ipont1 ) = idefml(2)
   ary(ipont1+1) = a*defml(1)
   ipont1 = ipont1 + 2
   GO TO 200
   210 CALL bckrec (slt)
   CYCLE
   220 iary(ipont) = iary(ipont) + 1
   IF (ipont+1 > nllst) GO TO 340
   iary(ipont1 ) = core(i)
   ary(ipont1+1) = 1.0
   ipont1 = ipont1 + 2
   inull  = 1
 END DO
 IF (inull == 0) GO TO 260
 240 ipont = ipont + iary(ipont)*2 + 1
 nload = nload + 1
 ipont1= ipont1+ 1
 GO TO 150
 250 CALL CLOSE (casecc,1)
 IF (islt == 1) CALL CLOSE (slt,1)
 CALL combin (pg(1),ilist(1),n1)
 RETURN
 
 260 iary(ipont) = 1
 IF (ipont+1 > nllst) GO TO 340
 iary(ipont1 ) =-1
 ary(ipont1+1) = 1.0
 ipont1 = ipont1 + 2
 GO TO 240
 
 270 ip1 = slt
 280 ip2 =-1
 290 CALL mesage (ip2,ip1,subnam)
 ip1 = casecc
 GO TO 280
 310 ip1 = pg(1)
 GO TO 280
 320 ip1 = casecc
 GO TO 280
 330 ip2 =-2
 ip1 = slt
 GO TO 290
 
 340 i = icore(i)
 nwds = 0
 350 CALL READ (*330,*360,slt,core(1),lc,0,iflag)
 nwds = nwds + lc
 GO TO 350
 360 nwds = nwds + iflag
 WRITE  (iotpe,370) ufm,i,nllst,nwds
 370 FORMAT (a23,' 3176, INSUFFICIENT OPEN CORE AVAILABLE TO PROCESS ',  &
     'ALL LOAD CARD COMBINATIONS IN MODULE SSG1.',  &
     /32X,'CURRENT LOAD ID BEING PROCESSED IS',i9,1H.,  &
     /32X,'OPEN CORE AVAILABLE IS',i9,' WORDS.',  &
     /32X,'ADDITIONAL OPEN CORE REQUIRED IS',i9,' WORDS.')
 ip1 = 0
 ip2 =-61
 GO TO 290
END SUBROUTINE ssg1
