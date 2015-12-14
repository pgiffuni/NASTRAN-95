SUBROUTINE exlvl (nos,md,NAME,z,nwds)
     
!     EXLVL ADDS A SUBSTRUCTURE TO THE RESIDENT SOF FOR THE SOFIN
!     OPERATION.  IT USES THE DIT AND MDI DATA WRITTEN ON THE EXTERNAL
!     FILE BY SOFOUT TO RESTORE THE HL, CS, AND LL POINTERS IN THE MDI.
 
 
 INTEGER, INTENT(IN)                      :: nos
 INTEGER, INTENT(IN)                      :: md(4,1)
 INTEGER, INTENT(IN)                      :: NAME(2)
 INTEGER, INTENT(OUT)                     :: z(2)
 INTEGER, INTENT(IN OUT)                  :: nwds
 EXTERNAL lshift   ,rshift   ,andf     ,orf
 LOGICAL :: mdiup
 INTEGER :: buf      ,andf     ,rshift    ,  &
     ps       ,cs       ,hl       ,tp       , subr(2)  ,orf
 COMMON  /zzzzzz/   buf(1)
 COMMON  /system/   sysbuf   ,nout     ,x1(6)    ,nlpp      , x2(2)    ,line
 COMMON  /sof   /   x3(34)   ,mdiup
 DATA     subr  /   4HEXLV   ,4HL      /
 
 
!     ADD THE NEW SUBSTRUCTURE TO THE RESIDENT DIT.
 
 CALL fdsub (NAME,i)
 IF (i /= -1) GO TO 6104
 CALL crsub (NAME,i)
 IF (nos <= 0) GO TO 200
 z(1) = NAME(1)
 z(2) = NAME(2)
 nss  = 1
 iss  = 1
 
!     DECODE THE OLD MDI ENTRY
 
 5 DO  i = 1,nos
   IF (md(1,i) /= z(2*iss-1) .OR. md(2,i) /= z(2*iss)) CYCLE
   ps = andf(md(3,i),1023)
   tp = andf(rshift(md(3,i),20),1023)
   ll = rshift(md(4,i),20)
   cs = andf(rshift(md(4,i),10),1023)
   hl = andf(md(4,i),1023)
   iold = i
   GO TO 15
 END DO
 
!     SET NEW MDI POINTERS FOR HL, CS, AND LL IF THE SUBSTRUCTURES OF
!     THE ORIGINATING SOF WHICH ARE INDICATED THEREBY EXIST.
 
 
!     HIGHER LEVEL (HL)
 
 15 m = 0
 IF (hl == 0) GO TO 30
 CALL fdsub (md(1,hl),i)
 IF (i > 0) GO TO 20
 CALL crsub (md(1,hl),i)
 nss = nss + 1
 IF (2*nss > nwds) GO TO 9008
 z(2*nss-1) = md(1,hl)
 z(2*nss  ) = md(2,hl)
 20 m  = i
 hl = i
 
!     COMBINED SUBSTRUCTURE (CS)
 
 30 IF (cs == 0) GO TO 60
 CALL fdsub (md(1,cs),j)
 IF (j > 0) GO TO 50
 CALL crsub (md(1,cs),j)
 nss = nss + 1
 IF (2*nss > nwds) GO TO 9008
 z(2*nss-1) = md(1,cs)
 z(2*nss  ) = md(2,cs)
 50 m  = orf(m,lshift(j,10))
 cs = j
 
!     LOWER LEVEL (LL)
 
 60 IF (ll == 0) GO TO 90
 CALL fdsub (md(1,ll),j)
 IF (j > 0) GO TO 80
 CALL crsub (md(1,ll),j)
 nss = nss + 1
 IF (2*nss > nwds) GO TO 9008
 z(2*nss-1) = md(1,ll)
 z(2*nss  ) = md(2,ll)
 80 m  = orf(m,lshift(j,20))
 ll = j
 
!     UPDATE THE MDI
 
 90 CALL fdsub (z(2*iss-1),j)
 CALL fmdi (j,i)
 buf(i+1) = lshift(tp,20)
 buf(i+2) = m
 mdiup    =.true.
 
!     WRITE USER MESSAGES
 
 nl = 2
 IF (ll /= 0) nl = nl + 1
 IF (cs /= 0) nl = nl + 1
 IF (hl /= 0) nl = nl + 1
 IF (ps /= 0) nl = nl + 3
 IF (line+nl > nlpp) CALL page
 line = line + nl
 WRITE (nout,63470) z(2*iss-1),z(2*iss)
 IF (hl == 0) GO TO 100
 CALL fdit (hl,i)
 WRITE (nout,63471) buf(i),buf(i+1)
 100 IF (cs == 0) GO TO 130
 CALL fdit (cs,i)
 WRITE (nout,63472) buf(i),buf(i+1)
 130 IF (ll == 0) GO TO 160
 CALL fdit (ll,i)
 WRITE (nout,63473) buf(i),buf(i+1)
 160 IF (ps == 0) GO TO 170
 WRITE (nout,63590) z(2*iss-1),z(2*iss)
 170 iss = iss + 1
 IF (iss-nss > 0) THEN
   GO TO   210
 ELSE
   GO TO     5
 END IF
 
!     SUBSTRUCTURE ADDED TO SOF SUCCESSFULLY
 
 200 WRITE (nout,63470) NAME
 210 RETURN
 
!     SUBSTRUCTURE NAME WAS DUPLICATED
 
 6104 CALL smsg (4,0,NAME)
 RETURN
 
!     INSUFFICIENT CORE
 
 9008 CALL mesage (-8,0,subr)
 RETURN
 
!     MESSAGE TEXT
 
 63470 FORMAT (49H0*** user information message 6347, substructure ,  &
     2A4,18H added TO the sof.)
 63471 FORMAT (5X, 25HHIGHER level substructure,2X,2A4)
 63472 FORMAT (5X, 25HCOMBINED substructure    ,6(2X,2A4))
 63473 FORMAT (5X, 25HLOWER level substructure ,7(2X,2A4))
 63590 FORMAT (49H0*** user information message 6359, substructure ,  &
     2A4,41H was originally a secondary substructure./36X,  &
     42HON this sof, it is a primary substructure.)
END SUBROUTINE exlvl
