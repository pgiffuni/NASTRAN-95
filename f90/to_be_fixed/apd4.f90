SUBROUTINE apd4
     
 EXTERNAL        orf
 INTEGER :: nam(2),iz(1),back,pspa,ic(16),eid,pid,cidbx,silb,  &
     scr1,ecta,bgpa,gpla,useta,sila,acpt,buf10,ca4s,  &
     ca4e,pa4s,pa4e,auset(6,2),silc,orf,usa,uk,eidb,  &
     sildx(4),acsix(4),cid(5),necta(2)
 REAL :: vx1(3),vx2(3),acpl(3,3),rb1(3)
 DIMENSION       ai(6),head(9),ihead(9)
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
     (necta(1),eidb),(necta(2),cid(1)),  &
     (acsix(2),vx2(1)),(z(1),iz(1)),(eid,ic(1)),  &
     (sildx(1),icid),(sildx(3),silc),  &
     (ai(1),dy),(ai(2),bloc),(ai(3),d),(ai(4),ca),  &
     (ai(5),gap),(ai(6),nsize),(head(1),ihead(1))
 DATA     nam  / 4HAPD4,4H    /
 
 lca  = 16
 nc4  = ((ca4e-ca4s)+1)/lca
 ncam = ncam + nc4
 
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
 auset(5,2) = uk
 ihead(1) = 4
 silc = silb
 
!     LOOP ON NC4 MOVING CAERO4 CARD TO COMMON
 
 DO  i = 1,nc4
   ntot = 0
   n = (i-1)*lca - 1
   DO  j = 1,lca
     ic(j) = iz(ca4s+n+j)
   END DO
   mcstm = mcstm + 1
   iz(ca4s+n+2) = mcstm
   iz(ca4s+n+8) = 4
   acsix(1) = mcstm
   
!     FIND PAERO4 CARD
   
   CALL apdoe (pid,iz,pa4s,pa4e,ipid,npc)
   IF (ipid == 0) GO TO 999
   
!     FIND NUMBER OF STRIPS
   
   ispan = nspan
   iast  = 0
   IF (nspan /= 0) GO TO 20
   CALL apdoe (nchord,iz,naef1,naef2,iast,nspan)
   IF (iast == 0) GO TO 998
   nspan = nspan - 1
   iast  = iast + 1
   20 iz(ca4s+n+4) = nspan
   iz(ca4s+n+5) = 1
   ipp = ipid + 5
   npc = npc - 4
   npc = npc/3
   IF (npc < nspan ) GO TO 997
   ihead(8) = nspan
   
!     GET POINTS IN PROPER COORD SYSTEM
   
   CALL apdcs
   head(9) = 1.0/SQRT(1.0+((xp4+.25*(x43-x12))/yp4)**2)
   IF (next+6*nspan > left) GO TO 996
   ioc = next
   
!     GENERATE DATA FOR BOXES
   
   ncrdp = ncrd
   fsj1  = apdf(z(iast),1,ispan)
   yj1   = fsj1*yp4
   dj1   = fsj1*xp4
   cj1   = x12 + fsj1*(x43-x12)
   xij1  = dj1
   xi1j1 = dj1 + cj1
   eidb  = eid - 1
   DO  j = 1,nspan
     yj    = yj1
     fsj1  = apdf(z( iast),j+1,ispan)
     yj1   = fsj1*yp4
     dj1   = fsj1*xp4
     cj1   = x12 + fsj1*(x43-x12)
     dy    = (yj1 - yj)
     ya    = .5*dy + yj
     ysp   = ya
     cloc  = x12 - (x12-x43)*ya/yp4
     bloc  = cloc*.5
     doc   = z(ipp)
     caoc  = z(ipp+1)
     gapoc = z(ipp+2)
     ipp   = ipp + 3
     d     = doc*cloc
     ca    = caoc*cloc
     gap   = gapoc*cloc
     nsize = 2
     IF (caoc /= 0.0) nsize = 3
     nj    = nj + nsize
     nk    = nk + nsize
     ntot  = ntot + nsize
     kk    = 0
     DO  k = 1,6
       z(ioc+j+kk) =  ai(k)
       kk = kk + nspan
     END DO
     
!     EXTERNAL ID S
     
     eidb   = eidb + 1
     cid(1) = cidbx + 1 + 2*(j-1)
     cid(2) = cid(1) + 1
     cid(3) = cid(1) + 2
     cid(4) = cid(3) + 1
     cid(5) = eidb
     ncid   = cid(4)
     
!     BGPDT , SPL, AND USET
     
     xij    = xij1
     xi1j   = xi1j1
     xij1   = dj1
     xi1j1  = dj1 + cj1
     xic    = (xij+xij1+bloc)*.5
     vx1(3) = 0
     IF (j /= 1) GO TO 310
     ASSIGN 300 TO back
     icid   = cid(1)
     vx1(1) = xij
     vx1(2) = yj
     kk     = 1
     GO TO 340
     300 ASSIGN 310 TO back
     icid   = cid(2)
     vx1(1) = xi1j
     vx1(2) = yj
     kk     = 1
     GO TO 340
     310  ASSIGN 320 TO back
     icid   = cid(3)
     vx1(1) = xij1
     vx1(2) = yj1
     kk     = 1
     GO TO 340
     320 ASSIGN 330 TO back
     icid   = cid(4)
     vx1(1) = xi1j1
     vx1(2) = yj1
     kk     = 1
     GO TO 340
     330 ASSIGN 360 TO back
     icid   = cid(5)
     vx1(1) = xic
     IF (nsize == 3) auset(6,2) = uk
     vx1(2) = ysp
     kk     = 2
     340 CALL gmmats (acpl,3,3,0,vx1,3,1,0,vx2)
     DO  k = 1,3
       vx2(k) = vx2(k) + rb1(k)
     END DO
     CALL WRITE (bgpa,acsix,4,0)
     CALL WRITE (gpla,icid,1,0)
     CALL WRITE (useta,auset(1,kk),6,0)
     
!     SIL AND EQEXIN
     
     ncrd   = ncrd + 1
     silc   = silc + 6
     isiln  = isiln +6
     sildx(4) = isiln
     luseta = silc
     sildx(2) = 10*silc + 1
     CALL WRITE (sila,silc,1,0)
     CALL WRITE (scr2,isiln,1,0)
     CALL WRITE (scr2,silc,1,0)
     CALL WRITE (scr1,icid,2,0)
     GO TO back, (300,310,320,330,360)
     
!     ECT
     
     360 cid(1) = iapd(1,j  ,1,ncrdp)
     cid(2) = iapd(2,j  ,1,ncrdp)
     cid(4) = iapd(1,j+1,1,ncrdp)
     cid(3) = iapd(2,j+1,1,ncrdp)
     cid(5) = cid(3) + 1
     CALL WRITE (ecta,necta(1),6,0)
     auset(6,2) = usa
   END DO
   cidbx = ncid
   
!     PUT OUT ACPT REC
   
   ihead(2) = ntot
   ihead(3) = iz(ipid+1)
   lcla     = ihead(3)
   ihead(4) = iz(ipid+2)
   ihead(5) = iz(ipid+3)
   icirc    = ihead(5)
   ihead(6) = iz(ipid+4)
   ihead(7) = 0
   il       = 0
   in       = nspan + 1
   
!     PROPERTY DATA
   
   IF (lcla == 0 .AND. icirc == 0) GO TO 70
   IF (lcla == 0) GO TO 50
   CALL apdoe (ihead(4),iz,naef1,naef2,il,nw)
   IF (il == 0) GO TO 994
   IF (MOD(nw,in) /= 0) GO TO 994
   ihead(7) = nw/in
   GO TO 70
   50 IF (icirc == 0) GO TO 70
   CALL apdoe (ihead(6),iz,naef1,naef2,il,nw)
   IF (il == 0) GO TO 995
   in = 2 + 2*icirc
   IF (MOD(nw,in) /= 0) GO TO 995
   ihead(7) = nw/in
   70 CALL WRITE (acpt,ihead,9,0)
   CALL WRITE (acpt,z(ioc+1),nspan*6,0)
   IF (il /= 0) CALL WRITE (acpt,z(il+1),nw,0)
   CALL WRITE (acpt,0,0,1)
 END DO
 silb = silc
 1001 RETURN
 
!     ERROR MESSAGES
 
 994 i = ihead(4)
 9941 WRITE  (NOT,9942) ufm,i,eid
 9942 FORMAT (a23,' 2429, WRONG NUMBER OF WORDS OR CARD NOT FOUND FOR ',  &
     'CARD ID',i9, /29X,'ASSOCIATED WITH CAERO4 ID',i9)
 GO TO 1000
 995 i = ihead(6)
 GO TO 9941
 996 CALL mesage (-8,0,nam)
 997 i = pid
 GO TO 9941
 998 i = nchord
 GO TO 9941
 999 CALL emsg (0,2323,1,2,0)
 WRITE  (NOT,891) pid,eid
 891 FORMAT (10X,16HPAERO4 card no. ,i8,31H referenced by caero4 card n  &
     o. ,i8,15H does NOT EXIST)
 GO TO 1000
 1000 CALL mesage (-61,0,nam)
 GO TO 1001
END SUBROUTINE apd4
