SUBROUTINE apd12
     
 EXTERNAL        orf
 LOGICAL :: ls,lc,dlb
 INTEGER :: auset(6,2),orf,pspa,uk,usa,iz(1),nam(2),  &
     eid,pid,cp,cidbx,acsid,silb,scr1,scr2,scr3,scr4,  &
     scr5,ecta,bgpa,gpla,useta,sila,cstma,acpt,buf10, buf11,buf12,iax(20),  &
     ca2s,ca2e,ca3s,ca3e,ca4s,ca4e, pa2s,pa2e,pa3s,pa3e,pa4s,pa4e
 COMMON /system/ sysbuf,iut
 COMMON /apd1c / eid,pid,cp,nspan,nchord,lspan,lchord,igid,  &
     x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm,  &
     ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd,scr1,  &
     scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,sila,  &
     cstma,acpt,buf10,buf11,buf12,next,left,isiln,  &
     ncam,naef1,naef2,nca1,nca2,ca2s,ca2e,ca3s,ca3e,  &
     ca4s,ca4e,npa1,npa2,pa2s,pa2e,pa3s,pa3e,pa4s,pa4e
 COMMON /apd12c/ key(5),auset,usa,uk,ncam2,nasb,ippc
 COMMON /bitpos/ ibit(64)
 COMMON /two   / itwo(32)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE    (z(1),iz(1)),(eid,iax(1))
 DATA    nam   /4HAPD1,4H2   /
 
 i17 = ibit(17)
 i18 = ibit(18)
 i19 = ibit(19)
 i20 = ibit(20)
 pspa= orf(itwo(i17),itwo(i20))
 usa = orf(pspa,itwo(i18))
 uk  = orf(itwo(i19),itwo(i20))
 DO  j = 1,2
   DO  i = 1,6
     auset(i,j) = usa
   END DO
 END DO
 auset(3,2) = uk
 auset(5,2) = uk
 ncam = ((nca2-nca1)+1)/16
 IF (nca1 == 0) ncam = 0
 ncam2 = ((ca2e-ca2s)+1)/16
 IF (ca2s == 0) ncam2 = 0
 lca = 16
 
!     CREATE IGID SEQUENCE ARRAY
 
 nigid1 = next
 igid2  = next
 nigid  = next
 nx = nca1
 j  = nigid1
 IF (ncam == 0) GO TO 15
 DO  i = 1,ncam
   iz(j) = iz(nx+7)
   j = j + 1
   iz(j) = nx
   nx = nx + lca
   j = j + 1
 END DO
 
!     SORT IGID ARRAY ON IGID
 
 CALL sort (0,0,2,1,iz(nigid1),2*ncam)
 15 IF (ncam2 == 0) GO TO 30
 nx = ca2s
 nigid2 = j
 DO  i = 1,ncam2
   iz(j) = iz(nx+7)
   j = j + 1
   iz(j) = nx
   nx = nx + lca
   j = j + 1
 END DO
 CALL sort (0,0,2,1,iz(nigid2),2*ncam2)
 igid2 = nigid2
 30 nextc = j
 IF (ncam == 0) GO TO 500
 nigid = nigid1
 
!     OUTTER LOOP PROCESSES CAERO1 CARDS
 
 DO  i = 1,ncam
   
!     SET APD1 INPUT COMMON BLOCK
   
   nc = iz(nigid+1) - 1
   
!     MOVE CAERO TO COMMON
   
   DO  j = 1,16
     n1 = j + nc
     iax(j) = iz(n1)
   END DO
   mcstm  = mcstm + 1
   iz(nc+2) = mcstm
   
!     FIND PAERO1 CARD
   
   IF (npa1 == 0) GO TO 890
   DO  j = npa1,npa2,8
     ippc = j
     IF (pid == iz(j)) GO TO 270
   END DO
   GO TO 890
   270 xop  = .25
   x1p  = .75
   alzo = 0.0
   
!     FIND AEFACT ARRAYS IF PRESENT
   
   jspan  = nspan
   jchord = nchord
   IF (lspan == 0) GO TO 280
   CALL apdoe (lspan,iz,naef1,naef2,ispan,jspan)
   IF (ispan == 0) GO TO 850
   ispan = ispan + 1
   jspan = jspan - 1
   280 IF (lchord == 0) GO TO 350
   CALL apdoe (lchord,iz,naef1,naef2,ichord,jchord)
   IF (ichord == 0) GO TO 860
   ichord = ichord + 1
   jchord = jchord - 1
   350 CONTINUE
   
!     CHECK IF FIRST OR LAST ENTRY IN IGID SET
   
   ls = .false.
   IF (i == 1) GO TO 370
   IF (iz(nigid) == iz(nigid-2)) GO TO 380
   370 ls = .true.
   380 lc = .false.
   dlb= .false.
   IF (i == ncam) GO TO 390
   IF (iz(nigid) == iz(nigid+2)) GO TO 400
   390 lc = .true.
   
!     CHECK FOR CAERO2 ELEMENT
   
   IF (ncam2 == 0) GO TO 400
   IF (nigid2 > nextc) GO TO 50
   40 IF (iz(nigid2) > iz(nigid)) GO TO 50
   IF (iz(nigid) == iz(nigid2)) dlb = .true.
   IF (dlb) GO TO 50
   nigid2 = nigid2 + 2
   IF (nigid2 > nextc) GO TO 50
   GO TO 40
   50 CONTINUE
   IF (dlb) lc = .false.
   
!     CALL APD1 TO MANUFACTURE BOXES
   
   400 CALL apd1 (z(ispan),jspan,z(ichord),jchord,ls,lc)
   nchord = jchord
   nspan  = jspan
   iz(nc+4) = nspan
   iz(nc+5) = nchord
   iz(nc+8) = 1
   IF (.NOT.dlb) GO TO 410
   
!     PROCESS CAERO2 WITH CAERO1
   
   CALL apd2 (1,iz(next),iz(igid2 ),nextc,iz(nigid))
   410 nigid = nigid + 2
 END DO
 
!     PROCESS CAERO2 CARDS NOT PROCESSED YET
 
 500 IF (ncam2 == 0) GO TO 1000
 CALL apd2 (0,iz(next),iz(igid2 ),nextc,iz(nigid))
 1000 RETURN
 
!     ERROR MESSAGES
 
 812 CALL mesage (-61,0,nam)
 850 CALL emsg (0,2326,1,2,0)
 WRITE  (iut,851) eid,lspan
 851 FORMAT (10X,19HCAERO1 element no. ,i8,28H references aefact card n  &
     o. ,i8,22H which does NOT EXIST.)
 GO TO 812
 860 CALL emsg (0,2327,1,2,0)
 WRITE (iut,851) eid,lchord
 GO TO 812
 890 CALL emsg (0,2323,1,2,0)
 WRITE  (iut,891) pid,eid
 891 FORMAT (10X,16HPAERO1 card no. , i8,31H referenced by caero1 card  &
     no. ,i8,20H but does NOT EXIST.)
 GO TO 812
END SUBROUTINE apd12
