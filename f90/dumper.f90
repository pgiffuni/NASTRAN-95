SUBROUTINE dumper
     
!     THIS SUBROUTINE DUMPS THE OSCAR
 
 EXTERNAL         lshift,rshift,andf
 INTEGER :: ixtra(3),con1,con2
 INTEGER :: recno,dmapno,op,oscar(1),os(5),rshift,iname(2),  &
     tp,ap,andf,vps,ptype,bl,el,ml,cl,ceitbl
 DIMENSION        ra(4),roscar(1),loco(300),avps(1),ihd(96)
 DOUBLE PRECISION :: dprec,dprec1
 COMMON /output/ ititle(96),ihead(96)
 COMMON /zzzzzz/ core(1)
 COMMON /xgpi2 / lmpl,mplpnt,mpl(1)
 COMMON /xvps  / vps(1)
 COMMON /xceitb/ ceitbl(1)
 COMMON /system/ sysbuf,op,junk5(6),nlpp,junk6(2),nlines
 COMMON / xgpic/ junk22(28),nosgn
 EQUIVALENCE     (vps(1),avps(1)), (dprec,ra(1)), (dprec1,ra(3)),  &
     (oscar(1),roscar(1),os(5)), (core(1),os(1)), (iosbot,os(3))
 DATA   mask1,   mask2,      mask3,        mask4,      mask5     /  &
     32767,   32768, 1073676288,   1073741824,     983040     /
 DATA   con1,    con2 /     4HCONS,4HTANT                        /
 DATA   ihd/2*4H    ,4H COS,4HMIC ,4H/ na,4HSTRA,4HN dm,4HAP c,  &
     4HOMPI,4HLER ,4H- os,4HCAR ,4HLIST,4HING ,82*4H    /
 DATA   ixtra/4H(con,4HTINU,4HED)  /
 DATA   ion, ioff /  4HON  ,4HOFF  /
 
 10 FORMAT (20X,2A4,5H(i  ),2X,i10)
 20 FORMAT (20X,2A4,5H(r  ),2X,e15.6)
 30 FORMAT (20X,2A4,5H(bcd),5X,2A4)
 40 FORMAT (20X,2A4,5H(rdp),2X,d24.15)
 
!     INITIALIZE LOCO ARRAY - POINTS TO FIRST WORD IN MPL FOR MOD I
 
 j = 1
 i = 1
 50 loco(i) = j
 j = j + mpl(j)
 IF (j > lmpl) GO TO 60
 i = i + 1
 GO TO 50
 60 CONTINUE
 
 i = 1
 DO  k=1,96
   ihead(k) = ihd(k)
 END DO
 CALL page
 DO  k=1,3
   ihead(k+14) = ixtra(k)
 END DO
 
!     PROCESS ENTRY HEADER
 
 90 IF (mi == 11) GO TO 910
 nwe   = oscar(i  )
 recno = oscar(i+1)
 mi    = rshift(oscar(i+2),16)
 msave = loco(mi)
 itype = oscar(i+2) - lshift(rshift(oscar(i+2),16),16)
 iexflg= ioff
 IF (oscar(i+5) < 0) iexflg = ion
 dmapno = andf(nosgn,oscar(i+5))
 nlines = nlines + 4
 IF (nlines < nlpp) GO TO 100
 CALL page
 nlines = nlines + 4
 100 CONTINUE
 WRITE  (op,110)
 110 FORMAT (/1X,18(4H****))
 WRITE  (op,120) recno,itype,iexflg,oscar(i+3),oscar(i+4),dmapno
 120 FORMAT (2X,20HOSCAR record NUMBER ,i3,5X,14HMODULE TYPE = ,i2,  &
     5X,16HEXECUTE flag -- , a4, /2X,  &
     15HMODULE NAME -  ,2A4,5X,21HDMAP instruction no. ,i3)
 i   = i + 6
 nwe = nwe - 6
 SELECT CASE ( itype )
   CASE (    1)
     GO TO 130
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 800
   CASE (    4)
     GO TO 540
 END SELECT
 130 io  = 1
 nip = oscar(i)
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 140
 CALL page
 nlines = nlines + 2
 140 CONTINUE
 WRITE  (op,150) nip
 150 FORMAT (/10X,29HSUMMARY of INPUT DATA blocks(,i2,2H ) )
 j = 1
 160 iname(1) = oscar(i+1)
 iname(2) = oscar(i+2)
 ntu = andf(oscar(i+3),mask1)
 tp  = rshift(andf(oscar(i+3),mask2),15)
 ltu = rshift(andf(oscar(i+3),mask3),16)
 ap  = rshift(andf(oscar(i+3),mask4),30)
 IF (iname(1) == 0 .AND. io == 1) GO TO 190
 IF (iname(1) == 0 .AND. io == 0) GO TO 220
 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 170
 CALL page
 nlines = nlines + 1
 170 CONTINUE
 WRITE  (op,180) iname(1),iname(2),ap,ltu,tp,ntu
 180 FORMAT (20X,2A4,3X,i1,1H/,i5,1H/,i1,1H/,i5)
 GO TO 250
 190 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 200
 CALL page
 nlines = nlines + 1
 200 CONTINUE
 WRITE  (op,210) j
 210 FORMAT (20X,24H********INPUT DATA BLOCK,i3,8H is null)
 GO TO 250
 220 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 230
 CALL page
 nlines = nlines + 1
 230 CONTINUE
 WRITE  (op,240) j
 240 FORMAT (20X,25H********output DATA BLOCK,i3,8H is null)
 250 i = i + 3
 j = j + 1
 IF (j <= nip) GO TO 160
 IF (itype == 2 ) io = 0
 
!     PROCESS OUTPUT DATA BLOCKS
 
 IF (io == 0) GO TO 280
 io  = 0
 i   = i + 1
 nip = oscar(i)
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 260
 CALL page
 nlines = nlines + 2
 260 CONTINUE
 WRITE  (op,270) nip
 270 FORMAT (/10X,30HSUMMARY of output DATA blocks(,i2,2H ))
 j = 1
 GO TO 160
 
!     PROCESS PARAMETER SECTION
 
 280 i = i + 2
 nparm = oscar(i)
 IF (nparm == 0) GO TO 530
 j = 1
 mplp = msave + 7
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 290
 CALL page
 nlines = nlines + 2
 290 CONTINUE
 WRITE  (op,300) nparm
 300 FORMAT (/10X,22HSUMMARY of parameters(,i2,2H ))
 310 IF (oscar(i+1) > 0.0) THEN
   GO TO   320
 ELSE
   GO TO   440
 END IF
 
!     SEARCH MPL FOR TYPE OF VARIABLE
 
 320 iname(1) = con1
 iname(2) = con2
 kk = IABS(mpl(mplp))
 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 330
 CALL page
 nlines = nlines + 1
 330 CONTINUE
 SELECT CASE ( kk )
   CASE (    1)
     GO TO 340
   CASE (    2)
     GO TO 360
   CASE (    3)
     GO TO 370
   CASE (    4)
     GO TO 390
   CASE (    5)
     GO TO 400
   CASE (    6)
     GO TO 420
 END SELECT
 340 WRITE (op,10) iname(1),iname(2),oscar(i+2)
 350 i = i + 2
 IF (mpl(mplp) > 0 ) mplp = mplp+1
 mplp = mplp+1
 j = j+1
 IF (j > nparm) GO TO 530
 GO TO 310
 360 WRITE (op,20) iname(1),iname(2),roscar(i+2)
 GO TO 350
 370 WRITE (op,30) iname(1),iname(2),oscar(i+2),oscar(i+3)
 380 i = i + 3
 IF (mpl(mplp) > 0) mplp = mplp+2
 mplp = mplp + 1
 j = j + 1
 IF (j > nparm) GO TO 530
 GO TO 310
 390 ra(1) = roscar(i+2)
 ra(2) = roscar(i+3)
 WRITE (op,40) iname(1),iname(2),dprec
 GO TO 380
 400 WRITE  (op,410) iname(1),iname(2),roscar(i+2),roscar(i+3)
 410 FORMAT (20X,2A4,5H(csp),2X,2E15.6)
 GO TO 380
 420 ra(1) = roscar(i+2)
 ra(2) = roscar(i+3)
 ra(3) = roscar(i+4)
 ra(4) = roscar(i+5)
 WRITE  (op,430) iname(1),iname(2),dprec,dprec1
 430 FORMAT (20X,2A4,5H(cdp),2X,2D24.15)
 i = i + 5
 IF (mpl(mplp) > 0) mplp = mplp+4
 mplp = mplp+1
 j = j + 1
 IF (j > nparm) GO TO 530
 GO TO 310
 440 ivps = andf(nosgn,oscar(i+1))
 iname(1) = vps(ivps-3)
 iname(2) = vps(ivps-2)
 ptype = rshift(andf(vps(ivps-1),mask5),16)
 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 450
 CALL page
 nlines = nlines + 1
 450 CONTINUE
 SELECT CASE ( ptype )
   CASE (    1)
     GO TO 460
   CASE (    2)
     GO TO 470
   CASE (    3)
     GO TO 480
   CASE (    4)
     GO TO 490
   CASE (    5)
     GO TO 500
   CASE (    6)
     GO TO 510
 END SELECT
 460 WRITE (op,10) iname(1),iname(2),vps(ivps)
 GO TO 520
 470 WRITE (op,20) iname(1),iname(2),avps(ivps)
 GO TO 520
 480 WRITE (op,30) iname(1),iname(2),vps(ivps),vps(ivps+1)
 GO TO 520
 490 ra(1) = avps(ivps  )
 ra(2) = avps(ivps+1)
 WRITE (op,40) iname(1),iname(2),dprec
 GO TO 520
 500 WRITE (op,410) iname(1),iname(2),avps(ivps),avps(ivps+1)
 GO TO 520
 510 ra(1) = avps(ivps  )
 ra(2) = avps(ivps+1)
 ra(3) = avps(ivps+2)
 ra(4) = avps(ivps+3)
 WRITE (op,430) iname(1),iname(2),dprec,dprec1
 520 i = i + 1
 j = j + 1
 IF (mpl(mplp) > 0) mplp = mplp + ptype/3 + 1
 IF (ptype == 6) mplp = mplp + 1
 mplp = mplp + 1
 IF (j > nparm) GO TO 530
 GO TO 310
 
!     HAVE COMPLETED FUNCTIONAL MODULE
 
 530 i = i + 2
 IF (itype == 2) i = i - 1
 GO TO 90
 
!     PROCESS EXECUTIVE MODULES
 
 540 IF (mi - 3 > 0) THEN
   GO TO   600
 END IF
 
!     PROCESS CHKPNT
 
 550 ndb = oscar(i)
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 560
 CALL page
 nlines = nlines + 2
 560 CONTINUE
 WRITE  (op,570) ndb
 570 FORMAT (/10X,31HDATA blocks TO be checkpointed(,i2,2H ))
 ist  = i + 1
 ifin = ist + 2 * ndb - 1
 npage = (10+ndb)/10+1
 nlines = nlines + npage
 IF (nlines < nlpp) GO TO 580
 CALL page
 nlines = nlines + npage
 580 CONTINUE
 IF (ndb /= 0) WRITE(op,590) (oscar(k),k=ist,ifin)
 i = i + 2*ndb+1
 590 FORMAT ((20X,10(2A4,2X)),/)
 GO TO 90
 600 IF (mi - 8 > 0) THEN
   GO TO   670
 END IF
 
!     PROCESS SAVE
 
 610 nparm  = oscar(i)
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 620
 CALL page
 nlines = nlines + 2
 620 CONTINUE
 WRITE  (op,630) nparm
 630 FORMAT (/10X,23HPARAMETERS TO be saved(,i2,2H ))
 640 FORMAT (20X,2A4,2X,i5)
 j = 1
 650 ivps = oscar(i+1)
 iname(1) = vps(ivps-3)
 iname(2) = vps(ivps-2)
 nlines   = nlines + 1
 IF (nlines < nlpp) GO TO 660
 CALL page
 nlines = nlines + 1
 660 CONTINUE
 WRITE (op,640) iname(1),iname(2),oscar(i+2)
 j = j + 1
 i = i + 2
 IF (j <= nparm) GO TO 650
 i = i + 1
 GO TO 90
 670 ndb = oscar(i)
 nwe = nwe - 1
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 680
 CALL page
 nlines = nlines + 2
 680 CONTINUE
 IF (mi ==  9) WRITE (op,690) ndb
 IF (mi == 10) WRITE (op,700) ndb
 690 FORMAT (/10X,25HDATA blocks TO be purged( ,i2,2H ))
 700 FORMAT (/10X,26HDATA blocks TO be equived(,i2,2H ))
 ist  = i + 1
 ifin = ist + 2*ndb - 1
 IF (mi /= 10) GO TO 730
 ntu = rshift(oscar(ist+2),16)
 ltu = oscar(ist+2) - lshift(ntu,16)
 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 710
 CALL page
 nlines = nlines + 1
 710 CONTINUE
 WRITE  (op,720) oscar(ist),oscar(ist+1),ntu,ltu
 720 FORMAT (20X,19HPRIMARY DATA BLOCK ,2A4,3X,i5,1H/,i5)
 ist  = ist  + 3
 ifin = ifin + 1
 nwe  = nwe  - 3
 730 CONTINUE
 npage  = (10+ndb)/10+1
 nlines = nlines + npage
 IF (nlines < nlpp) GO TO 740
 CALL page
 nlines = nlines + npage
 740 CONTINUE
 WRITE  (op,750)(oscar(k),k=ist,ifin)
 750 FORMAT ((20X,10(2A4,2X)),/)
 nwe = nwe - 2*ndb + 2
 IF (mi == 9) nwe = nwe - 2
 ivps = oscar(ifin+1)
 nlines = nlines + 1
 IF (nlines < nlpp) GO TO 760
 CALL page
 nlines = nlines + 1
 760 CONTINUE
 IF (ivps < 0) WRITE (op,770)
 770 FORMAT (20X,35HDEFAULT PARAMETER - always negative)
 IF (ivps < 0) GO TO 790
 WRITE  (op,780) vps(ivps-3),vps(ivps-2)
 780 FORMAT (20X,21HCONTROL PARAMETER is ,2A4)
 790 CONTINUE
 i = i + 2*ndb + 2
 IF (mi == 10 ) i = i + 1
 nwe = nwe - 1
 IF (nwe > 0) GO TO 670
 GO TO 90
 
!     PROCESS CONTROL INSTRUCTIONS
 
 800 irn = rshift(oscar(i),16)
 IF (mi == 11 .OR. mi == 12) GO TO 810
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 810
 CALL page
 nlines = nlines + 2
 810 CONTINUE
 IF (mi /= 11 .AND. mi /= 12) WRITE (op,820) irn
 820 FORMAT (/10X,25HRE-ENTRY record NUMBER = ,i4)
 IF (mi == 6) GO TO 900
 iw = oscar(i) - lshift(irn,16)
 IF (mi /= 7) GO TO 860
 
!     CONDITIONAL INSTRUCTION
 
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 840
 CALL page
 nlines = nlines + 2
 840 CONTINUE
 WRITE  (op,850) vps(iw-3),vps(iw-2)
 850 FORMAT (/10X,21HPARAMETER for cond = ,2A4)
 GO TO 900
 860 bl = rshift(ceitbl(iw-1),16)
 el = ceitbl(iw-1) - lshift(bl,16)
 ml = rshift(ceitbl(iw),16)
 cl = ceitbl(iw  ) - lshift(ml,16)
 nlines = nlines + 2
 IF (nlines < nlpp) GO TO 870
 CALL page
 nlines = nlines + 2
 870 CONTINUE
 IF (mi == 5) WRITE (op,880) bl,el,ml,cl,ceitbl(iw+1), ceitbl(iw+2)
 IF (mi == 11 .OR. mi == 12) WRITE (op,890) el,ml,cl
 880 FORMAT (/20X,i5,1H/,i5,5X,i5,1H/,i5,5X,2A4)
 890 FORMAT (/20X,5X,1H/,i5,5X,i5,1H/,i5)
 900 i = i + 1
 GO TO 90
 910 RETURN
END SUBROUTINE dumper
