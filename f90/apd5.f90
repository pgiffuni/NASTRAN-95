SUBROUTINE apd5
     
 EXTERNAL        orf
 INTEGER :: nam(2),iz(1),back,pspa,ic(16),eid,pid,cidbx,silb,  &
     scr1,ecta,bgpa,gpla,useta,sila,acpt,buf10,ca5s,  &
     ca5e,pa5s,pa5e,auset(6,2),silc,orf,usa,uk,eidb,  &
     sildx(4),acsix(4),cid(5),necta(2)
 REAL :: vx1(3),vx2(3),acpl(3,3),rb1(3)
 DIMENSION       ai(3),head(10),ihead(10)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nk,nj,luseta
 COMMON /apd1c / eid,pid,cp,nspan,nchord,nthry,nthick,igid,  &
     x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm,  &
     ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd,scr1,  &
     scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,sila,  &
     cstma,acpt,buf10,buf11,buf12,next,left,isiln, ncam,naef1,naef2,  &
     nca1,nca2,ca2s,ca2e,ca3s,ca3e,ca4s,ca4e,  &
     npa1,npa2,pa2s,pa2e,pa3s,pa3e,pa4s,pa4e,ca5s,ca5e, pa5s,pa5e
 COMMON /apd1d / icpl(14),yp4,s1,c1,xp2,xp3,xp4,ra1(3)
 COMMON /zzzzzz/ z(1)
 COMMON /bitpos/ ibit(64)
 COMMON /two   / itwo(32)
 COMMON /system/ sysbuf,NOT
 EQUIVALENCE     (icpl(3),rb1(1)),(icpl(6),acpl(1,1)),  &
     (necta(1),eidb),(necta(2),cid(1)),  &
     (acsix(2),vx2(1)),(z(1),iz(1)),(eid,ic(1)),  &
     (sildx(1),icid),(sildx(3),silc), (ai(1),dy),(ai(2),bloc),(ai(3),ca),  &
     (head(1),ihead(1))
 DATA    nam   / 4HAPD5,4H    /
 
 lca  = 16
 nc5  = ((ca5e-ca5s)+1)/lca
 ncam = ncam + nc5
 
!     INITIAL SETUP
 
 i17  = ibit(17)
 i18  = ibit(18)
 i19  = ibit(19)
 i20  = ibit(20)
 pspa = orf(itwo(i17),itwo(i20))
 usa  = orf(pspa,itwo(i18))
 uk   = orf(itwo(i19),itwo(i20))
 DO  j = 1,2
   DO  i = 1,6
     auset(i,j) = usa
   END DO
 END DO
 auset(3,2) = uk
 auset(5,2) = uk
 ihead(1)   = 5
 silc = silb
 
!     LOOP ON NC5 MOVING CAERO5 CARD TO COMMON
 
 DO  i = 1,nc5
   ntot = 0
   n = (i-1)*lca - 1
   DO  j = 1,lca
     ic(j) = iz(ca5s+n+j)
   END DO
   mcstm = mcstm + 1
   iz(ca5s+n+2) = mcstm
   iz(ca5s+n+8) = 5
   acsix(1) = mcstm
   
!     FIND PAERO5 CARD
   
   CALL apdoe (pid,iz,pa5s,pa5e,ipid,npc)
   IF (ipid == 0) GO TO 999
   
!     FIND NUMBER OF STRIPS
   
   ispan = nspan
   iast  = 0
   IF (nspan /= 0) GO TO 20
   CALL apdoe (nchord,iz,naef1,naef2,iast,nspan)
   IF (iast == 0) GO TO 998
   nspan = nspan - 1
   iast  = iast + 1
   20 iz(ca5s+n+4) = nspan
   iz(ca5s+n+5) = 1
   ipp = ipid + 7
   npc = npc - 6
   IF (npc < nspan) GO TO 997
   ihead(9) = nspan
   
!     GET POINTS IN PROPER COORD SYSTEM
   
   CALL apdcs
   head(10) = SQRT(1.0+(xp4/yp4)**2)
   IF (next+3*nspan > left) GO TO 996
   ioc = next
   
!     GENERATE DATA FOR BOXES
   
   ncrdp= ncrd
   fsj1 = apdf(z(iast),1,ispan)
   yj1  = fsj1*yp4
   dj1  = fsj1*xp4
   cj1  = x12 + fsj1*(x43-x12)
   xij1 = dj1
   xi1j1= dj1 + cj1
   eidb = eid - 1
   DO  j = 1,nspan
     yj   = yj1
     fsj1 = apdf(z(iast),j+1,ispan)
     yj1  = fsj1*yp4
     dj1  = fsj1*xp4
     cj1  = x12 + fsj1*(x43-x12)
     dy   = yj1 - yj
     ya   = .5*dy + yj
     ysp  = ya
     cloc = x12 - (x12-x43)*ya/yp4
     bloc = cloc*.5
     caoc = z(ipp)
     ipp  = ipp + 1
     ca   = caoc*cloc
     nsize= 2
     IF (caoc /= 0.0) nsize = 3
     nj   = nj + nsize
     nk   = nk + nsize
     ntot = ntot + nsize
     kk   = 0
     DO  k = 1,3
       z(ioc+j+kk) = ai(k)
       kk = kk + nspan
     END DO
     
!     EXTERNAL ID S
     
     eidb   = eidb   + 1
     cid(1) = cidbx  + 1 + 2*(j-1)
     cid(2) = cid(1) + 1
     cid(3) = cid(1) + 2
     cid(4) = cid(3) + 1
     cid(5) = eidb
     ncid   = cid(4)
     
!     BGPDT, SPL, AND USET
     
     xij  = xij1
     xi1j = xi1j1
     xij1 = dj1
     xi1j1= dj1 + cj1
     xic  = (xij+xij1+bloc)*.5
     vx1(3) = 0
     IF (j /= 1) GO TO 310
     ASSIGN 300 TO back
     icid   = cid(1)
     vx1(1) = xij
     vx1(2) = yj
     kk = 1
     GO TO 340
     300 ASSIGN 310 TO back
     icid   = cid(2)
     vx1(1) = xi1j
     vx1(2) = yj
     kk     = 1
     GO TO 340
     310 ASSIGN 320 TO back
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
     
     ncrd  = ncrd + 1
     silc  = silc + 6
     isiln = isiln +6
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
   ihead(4) = nthry
   IF (nthry == 1) head(10) = 0.0
   ihead(5) = nthick
   ihead(6) = iz(ipid+1)
   ihead(7) = iz(ipid+3)
   ihead(8) = iz(ipid+5)
   
!     PROPERTY DATA
   
   il = 0
   in = nspan + 1
   
!     ALPHAS
   
   CALL apdoe (iz(ipid+2),iz,naef1,naef2,il,nw)
   IF (il == 0) GO TO 994
   IF (ihead(6) == 1) GO TO 50
   ihead(3) = nw/in
   IF (MOD(nw,in) /= 0) GO TO 994
   ihead(6) = nspan
   GO TO 60
   50 ihead(3) = nw/2
   IF (MOD(nw,2) /= 0) GO TO 994
   
!     INTEGRALS
   
   60 INT = 0
   itn = 0
   IF (nthick == 0) GO TO 70
   CALL apdoe (nthick,iz,naef1,naef2,INT,nwi)
   IF (INT == 0) GO TO 995
   IF (ihead(7) == 0 .AND. nwi < 6) GO TO 995
   IF (ihead(7) /= 0 .AND. nwi /= 12) GO TO 995
   GO TO 90
   
!     TAUS
   
   70 CALL apdoe (z(ipid+6),iz,naef1,naef2,itn,nwt)
   IF (INT == 0) GO TO 993
   IF (ihead(8) == 0) GO TO 993
   IF (ihead(8) == 1) GO TO 80
   IF (nwt /= 3*nspan) GO TO 993
   ihead(8) = nspan
   GO TO 90
   80 IF (nwt /= 3) GO TO 993
   
!     THICKNESSES
   
   90 itk = 0
   IF (nthick /= 0 .AND. ihead(7) == 0) GO TO 99
   CALL apdoe (z(ipid+4),iz,naef1,naef2,itk,nwtk)
   IF (itk == 0) GO TO 992
   IF (INT == 0) GO TO 92
   IF (ihead(7) /= 1) ihead(7) = nspan
   IF (ihead(7) == nspan .AND. nwtk < nspan) GO TO 992
   GO TO 99
   92 IF (nwtk < 2) GO TO 992
   IF (ihead(8) == nspan .AND. nwtk < 2*nspan) GO TO 992
   99 CALL WRITE (acpt,ihead,10,0)
   CALL WRITE (acpt,z(ioc+1),nspan*3,0)
   CALL WRITE (acpt,z(il+1),nw,0)
   IF (INT /= 0) CALL WRITE (acpt,z(INT+1),nwi,0)
   IF (INT /= 0 .AND. itk /= 0) CALL WRITE (acpt,z(itk+1),nwtk,0)
   IF (itn /= 0) CALL WRITE (acpt,z(itn+1),nwt,0)
   IF (itn /= 0) CALL WRITE (acpt,z(itk+1),nwtk,0)
   CALL WRITE (acpt,0,0,1)
 END DO
 silb = silc
 1001 RETURN
 
!     ERROR MESSAGES
 
 992 i = iz(ipid+4)
 GO TO 9941
 993 i = iz(ipid+6)
 GO TO 9941
 994 i = ihead(6)
 9941 WRITE  (NOT,9942) ufm,i,eid
 9942 FORMAT (a23,' 2429, WRONG NUMBER OF WORDS OR CARD NOT FOUND FOR ',  &
     'CARD ID',i9, /29X,'ASSOCIATED WITH CAERO5 ID',i9)
 GO TO 1000
 995 i = nthick
 GO TO 9941
 996 CALL mesage (-8,0,nam)
 997 i = pid
 GO TO 9941
 998 i = nchord
 GO TO 9941
 999 CALL emsg (0,2323,1,2,0)
 WRITE  (NOT,891) pid,eid
 891 FORMAT (10X,16HPAERO5 card no. ,i8,31H referenced by caero5 card no. ,&
              i8,15H does NOT EXIST)
 GO TO 1000
 1000 CALL mesage (-61,0,nam)
 GO TO 1001
END SUBROUTINE apd5
