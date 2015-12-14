SUBROUTINE apd2 (iopt,cao1,cao2,ncore,id)
     
 
 INTEGER, INTENT(IN OUT)                  :: iopt
 INTEGER, INTENT(IN OUT)                  :: cao1(1)
 INTEGER, INTENT(IN OUT)                  :: cao2(1)
 INTEGER, INTENT(IN)                      :: ncore
 INTEGER, INTENT(IN)                      :: id
 INTEGER :: pa2s,pa2e,ca2s,ca2e,iz(1),nam(2),iax(1),  pc,ppc,bet,TYPE(3),  &
     cp,acsid,eid,eidb,cid(5),cidbx,auset(6,2),silb,  &
     uk,usa,necta(6),key(5),sildx(2),acsix(4),back,  &
     scr1,scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,  &
     sila,cstma,acpt,buf10,buf11,buf12,acsib,pid
 REAL :: rb1(3),acpl(3,3),vx1(3),vx2(3)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / nk,nj,luseta
 COMMON /system/ sysbuf,NOT
 COMMON /apd1c / eid,pid,cp,nspan,nchord,lspan,lchord,igid,  &
     x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm,  &
     ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd,  &
     scr1,scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,  &
     sila,cstma,acpt,buf10,buf11,buf12,next,left,isiln, ncam,naef1,naef2,  &
     nca1,nca2,ca2s,ca2e,ca3s,ca3e,ca4s,ca4e,  &
     npa1,npa2,pa2s,pa2e,pa3s,pa3e,pa4s,pa4e
 COMMON /apd1d / icpl(14),yp4,s1,c1,xp2,xp3,xp4,ra1(3)
 COMMON /apd12c/ key,auset,usa,uk,ncam2,nasb
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1)),(icpl(3),rb1(1)),(icpl(6),acpl(1,1)),  &
     (necta(1),eidb),(necta(2),cid(1)),(key(2),np),  &
     (key(3),nstrip),(key(4),ntp),(eid,iax(1)),  &
     (sildx(1),icid),(acsix(1),acsib),(acsix(2),vx2(1))
 DATA    TYPE  / 1HZ,2HZY,1HY /
 DATA    nam   / 4HAPD2,4H    /
 DATA    pio180/ .0174532925  /
 
 ibc   = 0
 nb    = 0
 IF (iopt == 1) GO TO 200
 
!     PROCESS CAERO2 WITHOUT CAERO1 ATTACHED
 
 np    = 0
 nstrip= 0
 ntp   = 0
 nas   = 0
 ipc   = 1
 ids   = cao2(1)
 10 IF (cao2(ipc) < 0) GO TO 191
 11 pc    = cao2(ipc+1) - 1
 IF (ibc /= nb) GO TO 31
 
!     LOOP OVER ALL CAERO2 WITH CURRENT ID TO SET UP POINTERS
 
 nb    = 0
 ibc   = 0
 nbz   = 0
 nby   = 0
 ntz   = 0
 nty   = 0
 ntzs  = 0
 ntys  = 0
 nbea1 = 0
 nsbea = 0
 nfl   = 0
 nt121 = 0
 nt122 = 0
 k     = ipc
 12 IF (cao2(k) /= ids) GO TO 19
 nb    = nb + 1
 l     = cao2(k+1) - 1
 DO  m = 1,7
   iax(m) = iz(l+m)
 END DO
 ASSIGN 14 TO iret
 GO TO 29
 14 CONTINUE
 SELECT CASE ( bet )
   CASE (    1)
     GO TO 15
   CASE (    2)
     GO TO 15
   CASE (    3)
     GO TO 16
 END SELECT
 15 nbz   = nbz  + 1
 ntz   = ntz  + nint
 ntzs  = ntzs + nsb
 IF (bet == 1) GO TO 17
 16 nby   = nby  + 1
 nty   = nty  + nint
 ntys  = ntys + nsb
 17 CONTINUE
 nbea1 = nbea1 + nint
 nsbea = nsbea + nsb
 nt121 = nt121 + nth1
 nt122 = nt122 + nth2
 nfl   = nfl + kt1
 k     = k + 2
 IF (k > ncam2*2) GO TO 19
 GO TO 12
 19 ids   = cao2(k)
 nto   = ntp + ntz + nty
 nas   = nasb
 
!     NOW SET UP POINTERS TO BUILD ACPT IN CORE
 
 i = ncore
 iz(i   ) = 2
 iz(i+ 1) = ntp
 iz(i+ 2) = ntp*2
 iz(i+ 3) = np
 iz(i+ 4) = nb
 iz(i+ 5) = ntp
 iz(i+ 6) = nbz
 iz(i+ 7) = nby
 iz(i+ 8) = ntz
 iz(i+ 9) = nty
 iz(i+10) = nto
 iz(i+11) = ntzs
 iz(i+12) = ntys
 iz(i+13) = nstrip
 inc    = i + 14
 inb    = inc  + np
 inas   = inb  + np
 inbea1 = inas + np
 inbea2 = inbea1 + nb
 insbea = inbea2 + nb
 izb    = insbea + nb
 iyb    = izb  + nb
 iavr   = iyb  + nb
 iarb   = iavr + nb
 infl   = iarb + nb
 ixle   = infl + nb
 ixte   = ixle + nb
 int121 = ixte + nb
 int122 = int121 + nb
 izs    = int122 + nb
 iys    = izs + nb + nstrip
 iee    = iys + nb + nstrip
 isg    = iee + nstrip
 icg    = isg + nstrip
 ix     = icg + nstrip
 idelx  = ix  + ntp   + nbea1
 ixic   = idelx + ntp + nbea1
 ixlam  = ixic  + ntp
 iao    = ixlam + ntp
 ixis1  = iao   + nsbea
 ixis2  = ixis1 + nsbea
 iaop   = ixis2 + nsbea
 iria   = iaop  + nsbea
 inasb  = iria  + nbea1
 ifla1  = inasb + nas
 ifla2  = ifla1 + nfl
 ith1a  = ifla2 + nfl
 ith2a  = ith1a + nt121
 nwr    = ith2a + nt122 - ncore
 na     = ith2a + nt122 - 1
 i      = na + np*6 +1
 IF (i > left) CALL mesage (-8,0,nam)
 
!     IF PANELS EXIST INSERT DATA FROM SCRATCH FILES
 
 IF(np == 0) GO TO 31
 nass = na
 CALL WRITE (scr3,0,0,1)
 CALL WRITE (scr4,0,0,1)
 CALL WRITE (scr5,0,0,1)
 CALL CLOSE (scr3,1)
 CALL CLOSE (scr4,1)
 CALL CLOSE (scr5,1)
 CALL gopen (scr3,z(buf10),0)
 CALL gopen (scr4,z(buf11),0)
 CALL gopen (scr5,z(buf12),0)
 DO  i = 1,np
   CALL fread (scr5,iz(inc),1,0)
   CALL fread (scr5,iz(inb),1,0)
   CALL fread (scr5,k,1,0)
   DO  j = 1,6
     iz(na+j) = iz(k+j)
   END DO
   inc = inc + 1
   inb = inb + 1
   na  = na  + 6
 END DO
 DO  i = 1,nstrip
   CALL fread (scr3,iz(iys),1,0)
   CALL fread (scr3,iz(izs),1,0)
   CALL fread (scr3,iz(iee),1,0)
   CALL fread (scr3,iz(isg),1,0)
   CALL fread (scr3,iz(icg),1,0)
   iys = iys + 1
   izs = izs + 1
   iee = iee + 1
   isg = isg + 1
   icg = icg + 1
 END DO
 DO  i = 1,ntp
   CALL fread (scr4,iz(ixic),1,0)
   CALL fread (scr4,iz(idelx),1,0)
   CALL fread (scr4,iz(ixlam),1,0)
   z(ix) = z(ixic) + .5*z(idelx)
   ixic  = ixic  + 1
   idelx = idelx + 1
   ixlam = ixlam + 1
   ix    = ix    + 1
 END DO
 CALL CLOSE (scr3,1)
 CALL CLOSE (scr4,1)
 CALL CLOSE (scr5,1)
 
!     FILL IN ASSOCIATED BODIES
 
 na  = nass
 DO  i = 1,np
   l   = 0
   loop25:  DO  j = 1,6
     IF (iz(na+j) == 0) CYCLE loop25
     l   = l + 1
     ibt = ipc
     DO  k = 1,nb
       m   = cao2(ibt+1)
       IF (iz(m) /= iz(na+j)) GO TO 28
       iz(inasb) = k
       inasb = inasb +1
       CYCLE loop25
       28 ibt = ibt + 2
     END DO
     GO TO 880
   END DO loop25
   iz(inas) = l
   inas = inas + 1
   na   = na + 6
 END DO
 31 CONTINUE
 ibc  = ibc + 1
 
!     MOVE TO COMMON
 
 DO  j = 1,16
   iax(j) = iz(j+pc)
 END DO
 iz(pc+2) = acsid
 acsib = acsid
 x4    = x1
 y4    = y1 + 1.0
 z4    = z1
 x43   = x12
 igid  =-igid
 CALL apdcs
 igid  =-igid
 
!     MOVE AERO CORD SYS TO ICPL
 
 IF (acsid == 0) GO TO 35
 DO  i = 1,14
   icpl(i) = iz(iacs+i-1)
 END DO
 35 CONTINUE
 ASSIGN 85 TO iret
 GO TO 29
 
!     FIND PAERO2 CARD
 
 29 CONTINUE
 IF (pa2s == 0) GO TO 990
 DO  j = pa2s,pa2e,15
   IF (pid == iz(j)) GO TO 40
 END DO
 GO TO 990
 40 ppc = j
 
!     GET BODY TYPE AND NUMBER OF ELEMENTS
 
 nsb  = nspan
 nint = nchord
 bet  = iz(ppc+1)
 DO  j = 1,3
   IF (bet == TYPE(j)) EXIT
 END DO
 60 bet  = j
 lth1 = iz(ppc+7)
 lth2 = iz(ppc+8)
 nth1 = 0
 nth2 = 0
 kt1  = 0
 IF (lspan == 0) GO TO 70
 CALL apdoe (lspan,iz,naef1,naef2,ispan,jspan)
 IF (ispan == 0) GO TO 950
 nsb  = jspan - 1
 70 IF (lchord == 0) GO TO 79
 CALL apdoe (lchord,iz,naef1,naef2,ichord,jchord)
 IF (ichord == 0) GO TO 960
 nint = jchord - 1
 79 IF (nint == 0) GO TO 80
 kt1  = kt1 + 1
 IF (iz(ppc+ 9) == 0) GO TO 920
 IF (iz(ppc+11) == 0) GO TO 75
 kt1  = kt1 + 1
 IF (iz(ppc+13) == 0) GO TO 75
 kt1  = kt1 + 1
 75 IF (lth1 == 0) GO TO 940
 CALL apdoe (lth1,iz,naef1,naef2,ith1,nth1)
 IF (ith1 == 0) GO TO 940
 IF (lth2 == 0) GO TO 80
 CALL apdoe (lth2,iz,naef1,naef2,ith2,nth2)
 IF (ith2 == 0) GO TO 930
 80 IF (nsb < 2) GO TO 970
 GO TO iret, (14,85)
 
!     PUT IN TERMS FOR SOME BODY ARRAYS
 
 85 iz(inbea1) = nint
 IF(ibc > 1 .AND. bet < iz(inbea2-1)) GO TO 870
 iz(inbea2) = bet
 iz(insbea) = nsb
 z(izb)  = ra1(3)
 z(iyb)  = ra1(2)
 z(izs)  = ra1(3)
 z(iys)  = ra1(2)
 z(iavr) = z(ppc+3)
 z(iarb) = z(ppc+4)
 iz(infl)= kt1
 iz(int121) = nth1
 iz(int122) = nth2
 inbea1  = inbea1 + 1
 inbea2  = inbea2 + 1
 insbea  = insbea+1
 izb  = izb + 1
 iyb  = iyb + 1
 izs  = izs + 1
 iys  = iys + 1
 iavr = iavr+ 1
 iarb = iarb+ 1
 infl = infl+ 1
 int121 = int121 + 1
 int122 = int122 + 1
 
!     ADD SOME MISC ARRAYS
 
 IF (nth1 == 0) GO TO 89
 DO  i = 1,nth1
   z(ith1a) = z(ith1+i)*pio180
   ith1a = ith1a + 1
 END DO
 IF (nth2 == 0) GO TO 88
 DO  i = 1,nth2
   z(ith2a) = z(ith2+i)*pio180
   ith2a = ith2a + 1
 END DO
 88 k = ppc + 9
 IF (iz(k) /= 1 .AND. iz(k+1) /= nint .AND. nth2 == 0) GO TO 910
 DO  i = 1,kt1
   iz(ifla1) = iz(k)
   iz(ifla2) = iz(k+1)
   k = k + 2
   IF (iz(ifla1) > iz(ifla2)) GO TO 910
   IF (iz(ifla2) > nint) GO TO 910
   IF (i == 1) GO TO 82
   IF (iz(ifla1) <= iz(ifla2-1)) GO TO 910
   82 ifla1 = ifla1 + 1
   ifla2 = ifla2 + 1
 END DO
 89 lrsb = iz(ppc+5)
 lrib = iz(ppc+6)
 IF (lrsb == 0) GO TO 91
 CALL apdoe (lrsb,iz,naef1,naef2,irsb,nrsb)
 IF (irsb ==     0) GO TO 900
 IF (nrsb /= nsb+1) GO TO 900
 91 IF (lrib ==     0) GO TO 92
 CALL apdoe (lrib,iz,naef1,naef2,irib,nrib)
 IF (irib ==      0) GO TO 890
 IF (nrib /= nint+1) GO TO 890
 92 CONTINUE
 width = z(ppc+3)
 
!     GENERATE ELEMENTS
 
 eidb   = eid - 1
 cidbx  = cidbx + 1
 vx1(2) = ra1(2)
 vx1(3) = ra1(3)
 
!     PUT IN PROPER MASKS FOR USET
 
 IF (bet == 1) GO TO 90
 auset(2,2) = uk
 auset(6,2) = uk
 IF (bet == 2) GO TO 90
 auset(3,2) = usa
 auset(5,2) = usa
 90 CONTINUE
 
!     BUMP NJ AND NK
 
 nja = nsb + nint
 nka = nsb*2
 nj  = nj + nja
 nk  = nk + nka
 iz(ncore+1) = iz(ncore+1) + nja
 iz(ncore+2) = iz(ncore+2) + nka
 IF (bet /= 2) GO TO 94
 nj  = nj + nja
 nk  = nk + nka
 iz(ncore+1) = iz(ncore+1) + nja
 iz(ncore+2) = iz(ncore+2) + nka
 94 i = 1
 95 eidb = eidb + 1
 cid(1) = cidbx
 cidbx  = cidbx + 1
 cid(2) = cidbx
 cid(5) = eidb
 
!     GRID POINTS IN AERO SYSTEM
 
 IF (i /= 1) GO TO 110
 ASSIGN 110 TO back
 icid = cid(1)
 IF (lspan == 0) vx1(1) = ra1(1) + (x12/nsb)*(i-1)
 IF (lspan /= 0) vx1(1) = ra1(1) + z(ispan+i)*x12
 oldx = vx1(1)
 z(ixle)  = oldx
 z(ixis1) = oldx
 ixis1 = ixis1 + 1
 kk = 1
 GO TO 130
 110 ASSIGN 120 TO back
 icid = cid(2)
 IF (lspan == 0) vx1(1) = ra1(1) + (x12/nsb)*i
 IF (lspan /= 0) vx1(1) = ra1(1) + z(ispan+i+1)*x12
 z(ixte ) = vx1(1)
 z(ixis2) = vx1(1)
 ixis2 = ixis2 + 1
 IF (i /= 1) z(ixis1) = oldx
 IF (i /= 1) ixis1 = ixis1 + 1
 kk = 1
 GO TO 130
 120 ASSIGN 160 TO back
 
!     A0 AND AOP
 
 z(iao ) = width
 z(iaop) = 0.0
 IF (lrsb == 0) GO TO 125
 z(iao ) = (z(irsb+i  ) + z(irsb+i+1))*.5
 z(iaop) = (z(irsb+i+1) - z(irsb+i))/(vx1(1)-oldx)
 125 iao  = iao  + 1
 iaop = iaop + 1
 temp = (vx1(1)+oldx)/2.0
 oldx = vx1(1)
 vx1(1) = temp
 icid = cid(5)
 kk   = 2
 
!     CONVERT TO BASIC
 
 130 IF (acsid == 0) GO TO 140
 CALL gmmats (acpl,3,3,0,vx1,3,1,0,vx2)
 DO  k = 1,3
   vx2(k) = vx2(k) + rb1(k)
 END DO
 GO TO 150
 140 DO  k = 1,3
   vx2(k) = vx1(k)
 END DO
 
!     PUT OUT BGPDT GPL USET
 
 150 CALL WRITE (bgpa,acsix,4,0)
 CALL WRITE (gpla,icid,1,0)
 CALL WRITE (useta,auset(1,kk),6,0)
 
!     BUMP POINTERS
!     PUT OUT SIL EQEXIN SILGA
 
 ncrd   = ncrd + 1
 silb   = silb + 6
 isiln  = isiln+ 6
 luseta = silb
 sildx(2) = 10*silb + 1
 CALL WRITE (sila,silb,1,0)
 CALL WRITE (scr2,isiln,1,0)
 CALL WRITE (scr2,silb,1,0)
 CALL WRITE (scr1,icid,2,0)
 GO TO back, (110,120,160)
 
!     PUT OUT ECT
 
 160 cid(1) = ncrd - 3
 IF (i == 1) cid(1) = cid(1) + 1
 cid(2) = ncrd - 1
 cid(3) = cid(1)
 cid(4) = cid(2)
 cid(5) = ncrd
 CALL WRITE (ecta,necta,6,0)
 i  = i + 1
 IF (i <= nsb) GO TO 95
 
!     INTEFERENCE CALCULATIONS AND ARRAYS
 
 IF (nint == 0) GO TO 170
 p1 = 1.0/nint
 DO  j = 1,nint
   z(iria) = width
   IF (lrib /= 0) z(iria) = .5*(z(irib+j)+z(irib+j+1))
   iria = iria + 1
   d1   = p1*(j-1)
   d2   = p1*j
   IF (lchord /= 0) d1 = z(ichord+j  )
   IF (lchord /= 0) d2 = z(ichord+j+1)
   z(idelx) = x12*(d2-d1)
   z(ix) = ra1(1) + x12*(d1+d2)/2.0
   IF (j ==    1) z(ixle) = ra1(1) + d1*x12
   IF (j == nint) z(ixte) = ra1(1) + d2*x12
   idelx = idelx + 1
   ix = ix + 1
 END DO
 170 CONTINUE
 ixle = ixle + 1
 ixte = ixte + 1
 iz(pc+ 4) = nsb
 iz(pc+ 5) = 1
 iz(pc+ 8) = 2
 iz(pc+16) = bet
 IF (bet == 1) GO TO 190
 auset(2,2) = usa
 auset(6,2) = usa
 auset(3,2) = uk
 auset(5,2) = uk
 190 IF (ibc == nb) CALL WRITE (acpt,iz(ncore),nwr,1)
 191 IF (iopt == 1) GO TO 230
 ipc = ipc + 2
 IF (ipc < ncam2*2) GO TO 10
 GO TO 1000
 
!     CAERO2 WITH CAERO1 ATTACHED
 
 200 ipc = 1
 ids = id
 210 IF (cao2(ipc) == id) GO TO 11
 220 ipc = ipc + 2
 IF (ipc < ncam2*2) GO TO 210
 GO TO 1000
 230 cao2(ipc) = -cao2(ipc)
 GO TO 220
 1000 RETURN
 
!     ERROR MESSAGES
 
 912 CALL mesage (-61,0,nam)
 870 WRITE  (NOT,8777) ufm,eid
 8777 FORMAT (a23,' 2273, CAERO2',i9,' NOT INPUT IN Z, ZY, Y SEQUENCE.')
 GO TO 912
 880 WRITE  (NOT,8888) ufm,iz(na+j),cao2(ibt)
 8888 FORMAT (a23,' 2274, ASSOCIATED BODY',i9,' WAS NOT FOUND WITH ',  &
     'CAERO2 GROUP',i9,1H.)
 GO TO 912
 890 j = lrib
 GO TO 941
 900 j = lrsb
 GO TO 941
 910 WRITE  (NOT,9111) ufm,eid
 9111 FORMAT (a23,' 2275, CAERO2',i9,' HAS INCONSISTENT USE FOR THI OR',  &
     ' THN, OR LTH2 IS REQUIRED.')
 GO TO 912
 920 WRITE  (NOT,9222) ufm,eid
 9222 FORMAT (a23,' 2276, THI1 AND THN1 REQUIRED FOR CAERO2',i9,1H.)
 GO TO 912
 930 j = lth2
 GO TO 941
 940 j = lth1
 941 WRITE  (NOT,9999) ufm,j,eid
 9999 FORMAT (a23,' 2429, WRONG NUMBER OF WORDS OR CARD NOT FOUND FOR',  &
     ' CARD ID',i9, /28X,'ASSOCIATED WITH CAERO2 ID',i9)
 GO TO 912
 950 CALL emsg (0,2326,1,2,0)
 WRITE  (NOT,951) eid,lspan
 951 FORMAT (10X,'CAERO2 ELEMENT NO.',i9,' REFERENCES AEFACT CARD NO.',  &
     i9,' WHICH DOES NOT EXIST.')
 GO TO 912
 960 CALL emsg (0,2327,1,2,0)
 WRITE (NOT,951) eid,lchord
 GO TO 912
 970 WRITE  (NOT,971) ufm,eid
 971 FORMAT (a23,' 2277, CAERO2 BODY',i9,' DOES NOT HAVE ENOUGH ',  &
     'SLENDER ELEMENTS.')
 GO TO 912
 990 CALL emsg (0,2323,1,2,0)
 WRITE  (NOT,991) pid,eid
 991 FORMAT (10X,'PAERO2 CARD NO.',i9,' REFERENCED BY CAERO2 CARD NO.',  &
     i9,' BUT DOES NOT EXIST.')
 GO TO 912
END SUBROUTINE apd2
