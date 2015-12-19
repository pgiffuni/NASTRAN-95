SUBROUTINE apd3
     
 EXTERNAL        orf
 LOGICAL :: cntrl1,cntrl2,crank1,crank2
 INTEGER :: nam(2),iz(1),back,pspa,ret,ic(16),eid,pid,cidbx,  &
     silb,scr1,ecta,bgpa,gpla,useta,sila,acpt,buf10,  &
     ca3s,ca3e,pa3s,pa3e,auset(6,2),silc,orf,usa,uk,  &
     eidb,sildx(2),acsix(4),cid(5),necta(2)
 REAL :: vx1(3),vx2(3),acpl(3,3),rb1(3)
 DIMENSION       ihead(10),bnd(24)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nk,nj,luseta
 COMMON /apd1c / eid,pid,cp,nspan,nchord,lspan,lchord,igid,  &
     x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm,  &
     ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd,scr1,  &
     scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,sila,  &
     cstma,acpt,buf10,buf11,buf12,next,left,isiln, ncam,naef1,naef2,  &
     nca1,nca2,ca2s,ca2e,ca3s,ca3e,ca4s,ca4e,  &
     npa1,npa2,pa2s,pa2e,pa3s,pa3e,pa4s,pa4e
 COMMON /apd1d / icpl(14),yp4,s1,c1,xp2,xp3,xp4,ra1(3)
 COMMON /zzzzzz/ z(1)
 COMMON /bitpos/ ibit(64)
 COMMON /two   / itwo(32)
 COMMON /system/ sysbuf,NOT
 EQUIVALENCE     (icpl(3),rb1(1)),(icpl(6),acpl(1,1)),  &
     (necta(1),eidb),(necta(2),cid(1)), (acsix(2),vx2(1)),(sildx(1),icid),  &
     (z(1),iz(1)),(eid,ic(1)), (crank1,ihead(3)),(crank2,ihead(4)),  &
     (cntrl1,ihead(5)),(cntrl2,ihead(6))
 DATA     nam  / 4HAPD3,4H    /
 
 nogo = 0
 lca  = 16
 nc3  = ((ca3e-ca3s)+1)/lca
 ncam = ncam+nc3
 
!     INITIAL SETUP
 
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
 ihead(1) = 3
 silc = silb
 
!     LOOP ON NC3 MOVING CAERO3 CARD TO COMMON
 
 DO  i = 1,nc3
   n = (i-1)*lca - 1
   DO  j = 1,lca
     ic(j) = iz(ca3s+n+j)
   END DO
   mcstm = mcstm + 1
   iz(ca3s+n+2) = mcstm
   iz(ca3s+n+8) = 3
   acsix(1) = mcstm
   
!     GET POINTS IN PROPER COORD SYSTEM
   
   CALL apdcs
   
!     FIND PAERO3 CARD
   
   j = pa3s
   20 IF (j >= pa3e) GO TO 999
   IF (iz(j) == pid) GO TO 30
   j = j + 4 + iz(j+3)
   GO TO 20
   30 ipid = j
   ihead(7) = iz(ipid+1)
   crank1 = .false.
   crank2 = .false.
   IF (z(ipid+5) > 0.0) crank1 = .true.
   IF (z(ipid+7) > 0.0) crank2 = .true.
   cntrl1 = .false.
   cntrl2 = .false.
   IF (iz(ipid+2) > 0 ) cntrl1 = .true.
   IF (iz(ipid+2) == 2 ) cntrl2 = .true.
   
!     GENERATE AERO POINTS FOR CAERO3  PUT POINTS 1-4 IN BGPDT
   
   DO  j = 13,24
     bnd(j) = 0.0
   END DO
   vx1(3) = 0.0
   kk = 1
   ASSIGN 50 TO back
   ibs = ncrd + 1
   vx1(1) = 0.0
   vx1(2) = 0.0
   bnd(1) = 0.0
   bnd(2) = 0.0
   GO TO 300
   50 ASSIGN 60 TO back
   vx1(1) = x12
   vx1(2) = 0.0
   bnd(7) = x12
   bnd(8) = 0.0
   GO TO 300
   60 vx1(1) = xp4
   vx1(2) = yp4
   bnd(5) = xp4
   bnd(6) = yp4
   ASSIGN 70 TO back
   GO TO 300
   70 ASSIGN 80 TO back
   vx1(1) = xp4 + x43
   vx1(2) = yp4
   bnd(11) = vx1(1)
   bnd(12) = vx1(2)
   GO TO 300
   
!     ADD POINTS 5 AND 6 IF THEY EXIST
   
   80 bnd(3) = bnd(5)
   bnd(4) = bnd(6)
   IF (.NOT.crank1) GO TO 90
   vx1(1) = z(ipid+4)
   vx1(2) = z(ipid+5)
   bnd(3) = vx1(1)
   bnd(4) = vx1(2)
   ASSIGN 90 TO back
   GO TO 300
   90 bnd(9) = bnd(11)
   bnd(10) = bnd(12)
   IF (.NOT.crank2) GO TO 100
   vx1(1) = z(ipid+6)
   vx1(2) = z(ipid+7)
   bnd(9) = vx1(1)
   bnd(10)= vx1(2)
   ASSIGN 100 TO back
   GO TO 300
   
!      ADD CONTROLS
   
   100 IF (.NOT.cntrl1) GO TO 120
   ASSIGN 101 TO back
   vx1(1)  = z(ipid+8)
   vx1(2)  = z(ipid+9)
   bnd(15) = vx1(1)
   bnd(16) = vx1(2)
   GO TO 300
   101 ASSIGN 102 TO back
   vx1(1)  = z(ipid+10)
   vx1(2)  = z(ipid+11)
   bnd(13) = vx1(1)
   bnd(14) = vx1(2)
   GO TO 300
   102 ASSIGN 103 TO back
   vx1(1)  = z(ipid+12)
   vx1(2)  = z(ipid+13)
   bnd(17) = vx1(1)
   bnd(18) = vx1(2)
   GO TO 300
   103 ASSIGN 104 TO back
   vx1(1)  = z(ipid+14)
   vx1(2)  = z(ipid+15)
   bnd(21) = vx1(1)
   bnd(22) = vx1(2)
   GO TO 300
   104 IF (.NOT.cntrl2) GO TO 120
   ASSIGN 105 TO back
   vx1(1)  = z(ipid+16)
   vx1(2)  = z(ipid+17)
   bnd(19) = vx1(1)
   bnd(20) = vx1(2)
   GO TO 300
   105 ASSIGN 120 TO back
   vx1(1)  = z(ipid+18)
   vx1(2)  = z(ipid+19)
   bnd(23) = vx1(1)
   bnd(24) = vx1(2)
   GO TO 300
   
!     CONNECT POINT TO BOXES FOR ECTA
   
   120 eidb   = eid
   cid(1) = ibs
   cid(2) = ibs + 1
   cid(5) = ibs
   IF (crank1) GO TO 121
   IF (crank2) GO TO 122
   cid(3) = ibs + 3
   cid(4) = ibs + 2
   GO TO 124
   121 IF (crank2) GO TO 123
   cid(3) = ibs + 3
   cid(4) = ibs + 4
   GO TO 124
   122 cid(3) = ibs + 4
   cid(4) = ibs + 2
   GO TO 124
   123 cid(3) = ibs + 5
   cid(4) = ibs + 4
   124 CONTINUE
   CALL WRITE (ecta,necta,6,0)
   eidb   = eidb + 1
   cid(1) = ibs + 2
   cid(2) = ibs + 3
   cid(5) = ibs + 2
   ibs    = ibs + 4
   IF (.NOT.crank1 .AND. .NOT.crank2) GO TO 130
   IF (crank1 .AND. crank2) GO TO 125
   cid(3) = ibs
   cid(4) = cid(5)
   ibs = ibs + 1
   GO TO 129
   125 cid(3) = ibs + 1
   cid(4) = ibs
   ibs    = ibs + 2
   129 CALL WRITE (ecta,necta,6,0)
   eidb   = eidb + 1
   130 IF (.NOT.cntrl1) GO TO 135
   cid(1) = ibs + 2
   cid(2) = ibs + 3
   cid(3) = ibs + 1
   cid(4) = ibs
   cid(5) = ibs + 2
   CALL WRITE (ecta,necta,6,0)
   eidb   = eidb + 1
   IF (.NOT.cntrl2) GO TO 135
   cid(3) = ibs+5
   cid(4) = ibs+4
   CALL WRITE (ecta,necta,6,0)
   
!     FIND CONTROL POINTS FOR ELEMENT
   
   135 CALL apdoe (nspan,iz,naef1,naef2,ilw,nww)
   IF (ilw == 0) GO TO 998
   IF (nww < 6) GO TO 998
   IF (MOD(nww,2) /= 0) GO TO 998
   ilc1 = 0
   ilc2 = 0
   nwc1 = 0
   nwc2 = 0
   IF (.NOT.cntrl1) GO TO 140
   CALL apdoe (nchord,iz,naef1,naef2,ilc1,nwc1)
   IF (ilc1 == 0) GO TO 997
   IF (nwc1 < 6) GO TO 997
   IF (MOD(nwc1,2) /= 0) GO TO 997
   IF (.NOT.cntrl2) GO TO 140
   CALL apdoe (lspan,iz,naef1,naef2,ilc2,nwc2)
   IF (ilc2 == 0) GO TO 996
   IF (nwc2 < 6) GO TO 996
   IF (MOD(nwc2,2) /= 0) GO TO 996
   140 ihead( 8) = nww/2
   ihead( 9) = nwc1/2
   ihead(10) = nwc2/2
   ihead( 2) = ihead(8)+ihead(9)+ihead(10)
   nk = nk + ihead(2)
   nj = nj + ihead(2)
   iz(ca3s+n+4) = ihead(2)
   iz(ca3s+n+5) = 1
   
!     START THE ACPT AND ADD THE CONTROL POINTS IN A LOOP
   
   CALL WRITE (acpt,ihead,10,0)
   CALL WRITE (acpt,bnd,24,0)
   eidb = eid - 1
   kk  = 2
   nn  = nww
   kkk = ilw - 1
   ASSIGN 150 TO ret
   GO TO 190
   150 IF (ihead(9) == 0) GO TO 180
   ASSIGN 160 TO ret
   nn  = nwc1
   kkk = ilc1 - 1
   GO TO 190
   160 IF (ihead(10) == 0) GO TO 180
   ASSIGN 180 TO ret
   nn  = nwc2
   kkk = ilc2 - 1
   GO TO 190
   180 CALL WRITE (acpt,0,0,1)
   
!     GEOMETRY CHECKS
   
   nm = 0
   IF (bnd(1) >  bnd(3)) nm = 1
   IF (bnd(3) >  bnd(5)) nm = 1
   IF (bnd(15) > bnd(17)) nm = 1
   IF (cntrl2 .AND. bnd(17) > bnd(19)) nm = 1
   IF (bnd(16) < bnd(14)) nm = 1
   IF (bnd(18) < bnd(22)) nm = 1
   IF (bnd(20) < bnd(24)) nm = 1
   IF (nm == 1) nogo = 1
   IF (nm == 1) WRITE (NOT,1851) ufm,eid
   1851 FORMAT (a23,' 2278, PLANFORM GEOMETRY FOR CAERO3 ID',i9,  &
       ' IS IN ERROR', /5X,'CHECK SWEEP  ANGLE FOR LEADING EDGE ',  &
       'OR CONTROL SURFACE HINGE LINE.')
   CYCLE
   
!     PUT CONTROL POINTS IN TABLE
   
   190 j = 2
   195 CONTINUE
   vx1(1) = z(kkk+j  )
   vx1(2) = z(kkk+j+1)
   CALL WRITE (acpt,vx1,2,0)
   ASSIGN 200 TO back
   GO TO 300
   200 CONTINUE
   j = j + 2
   IF (j <= nn) GO TO 195
   GO TO ret, (150,160,180)
   
!     BGPA  GPL  USET
   
   300 CALL gmmats (acpl,3,3,0,vx1,3,1,0,vx2)
   DO  k = 1,3
     vx2(k) = vx2(k) + rb1(k)
   END DO
   CALL WRITE (bgpa,acsix,4,0)
   IF (kk == 2) GO TO 320
   cidbx = cidbx + 1
   icid  = cidbx
   GO TO 330
   320 eidb  = eidb + 1
   icid  = eidb
   330 CALL WRITE (gpla,icid,1,0)
   CALL WRITE (useta,auset(1,kk),6,0)
   
!     SIL AND EQEXIN
   
   ncrd  = ncrd + 1
   silc  = silc + 6
   isiln = isiln + 6
   luseta= silc
   sildx(2) = 10*silc + 1
   CALL WRITE (sila,silc,1,0)
   CALL WRITE (scr2,isiln,1,0)
   CALL WRITE (scr2,silc,1,0)
   CALL WRITE (scr1,icid,2,0)
   GO TO back, (50,60,70,80,90,100,101,102,103,104,105,120,200)
 END DO
 silb = silc
 IF (nogo == 1) GO TO 1001
 1000 RETURN
 
!     ERROR MESSAGES
 
 996 i = lspan
 GO TO 9941
 997 i = nchord
 GO TO 9941
 998 i = nspan
 9941 WRITE  (NOT,9942) ufm,i,eid
 9942 FORMAT (a23,' 2429, WRONG NUMBER OF WORDS OR CARD NOT FOUND FOR ',  &
     'CARD ID',i9, /29X,'ASSOCIATED WITH CAERO3 ID',i9)
 GO TO 1001
 999 CALL emsg (0,2323,1,2,0)
 WRITE  (NOT,891) pid,eid
 891 FORMAT (10X,16HPAERO3 card no. ,i8,31H referenced by caero3 card no. ,&
              i8,15H does NOT EXIST)
 1001 CALL mesage (-61,0,nam)
 GO TO 1000
END SUBROUTINE apd3
