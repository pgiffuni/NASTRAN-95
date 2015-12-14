SUBROUTINE pltopr
     
 INTEGER :: prnt,ploter,pltype,eof,itype(4),list(20),skip,  &
     plttyp(3,3),a1(9),a2(12),a3(25),a4(13),a5(17),  &
     a6(6),b1(16),b2(16),b5(10),c1(16),c2(17),c3(24),  &
     film(2),paper(2),ilay(20),d1(11),d2(6),d3(29),  &
     prj(3,3),plus,minus,xyz(3),symm(4),BLANK,e1(12),  &
     e2(20),e3(14),e4(18),f1(22),f2(25),f3(12),f4(23),  &
     stress(34),dist(12),way(4),g1(13),g2(11),g3(16),  &
     g4(6),g5(11),g6(8),g7(4),displa(10),SPACE,layer,  &
     camera,bframs,pltmdl,tapden,paptyp,pensiz,penclr,  &
     axis,daxis,prject,origin,strain(4),p6(23),org, h1(11)
 REAL :: alst(20),cscale
 COMMON /BLANK / skpcom(20),prnt
 COMMON /xxparm/ pbufsz,camera,bframs,pltmdl(2),tapden,npens,  &
     papsiz(2),paptyp(2),pensiz(8),penclr(8,2),penpap,  &
     scale(2),skpscl(3),axis(3),daxis(3),vangle(9),  &
     vantx1,r0,s0l,s0r,t0,d0,vantx2(2),prject,s0s,  &
     for,org,norg,origin(11),origx3(11,4),  &
     xy(11,3),ncntr,cntr(50),icntvl,iwhere,idirec, skip23(23),layer
 COMMON /pltdat/ model,ploter,skpplt(17),cscale,skpa(2),cntsin,  &
     skpb(3),nopens,skpc(2),pltype,skpd(2),eof,cntin3
 EQUIVALENCE     (list(1),alst(1))
 
 DATA    nskip , skip    /1        ,4H(1X)                     /
 DATA    film  , paper   /4HFILM   ,1H       ,4HPAPE   ,1HR    /
 
!     PLOTTER TYPE FORMATS.
 
 DATA    np6   / 23 /
 DATA    p6    / 4H(10X  ,4H,38H   ,4HTHE    ,4HFOLL   ,4HOWIN  ,  &
     4HG pl  ,4HOTS    ,4HARE    ,4HFOR    ,4HA na  ,  &
     4HSTPL  ,4HT ,2   ,4HA4,a   ,4H2,8H   ,4HPLOT  ,  &
     4HTER   ,4H,2A4   ,4H,17H   ,4HTYPI   ,4HNG c  ,  &
     4HAPAB  ,4HILIT   ,4HY,/)   /
 DATA    plttyp/ 4HMICR  ,4HOFIL   ,1HM      ,  &
     4H  ta  ,4HBLE    ,1H       , 4H  dr  ,4HUM     ,1H       /
 DATA    itype / 4HWITH  ,4H       , 4HWITH  ,4HOUT    /
 
!     GENERAL PLOTTER FORMATS.
 
 DATA    na1   / 9 /
 DATA    a1    / 4H(//,  ,4H25H    ,4HP l    ,4HO t    ,4HT e   ,  &
     4HR     ,4H d a   ,4H t a   ,4H,/)    /
 DATA    na2   / 12 /
 DATA    a2    / 4H(10X  ,4H,27H   ,4HTHE    ,4HPLOT   ,4H tap  ,  &
     4HE is  ,4H wri   ,4HTTEN   ,4H at,   ,4HI4,4  , 4HH bp  ,4HI,/)   /
 DATA    na3   / 25 /
 DATA    a3    / 4H(10X  ,4H,89H   ,4HTHE    ,4HPLOT   ,4HS ar  ,  &
     4HE se  ,4HPARA   ,4HTED    ,4HBY e   ,4HND-o  ,  &
     4HF-fi  ,4HLE m   ,4HARKS   ,4H...t   ,4HWO e  ,  &
     4HND-o  ,4HF-fi   ,4HLE m   ,4HARKS   ,4H fol  ,  &
     4HLOW   ,4HTHE    ,4HLAST   ,4H plo   ,4HT,/)  /
 DATA    na4   / 13 /
 DATA    a4    / 4H(10X  ,4H,41H   ,4HAN e   ,4HND-o   ,4HF-fi  ,  &
     4HLE m  ,4HARK    ,4HFOLL   ,4HOWS    ,4HTHE   ,  &
     4HLAST  ,4H plo   ,4HT,/)   /
 DATA    na5   / 17 /
 DATA    a5    / 4H(10X  ,4H,56H   ,4HTHE    ,4HFIRS   ,4HT co  ,  &
     4HMMAN  ,4HD fo   ,4HR ea   ,4HCH p   ,4HLOT   ,  &
     4HCONT  ,4HAINS   ,4H the   ,4H plo   ,4HT nu  , 4HMBER  ,4H,/)    /
 DATA    na6   / 6  /
 DATA    a6    / 4H(10X  ,4H,9HC   ,4HSCAL   ,4HE =    ,4H,f5.  , 4H2,/)  /
 
!     TABLE PLOTTER FORMATS.
 
 DATA    nb1   / 16 /
 DATA    b1    / 4H(10X  ,4H,30H   ,4HSET    ,4HTHE    ,4HX +   ,  &
     4HY sc  ,4HALE    ,4HFACT   ,4HORS    ,4HAT,f  ,  &
     4H6.1,  ,4H12H    ,4HCOUN   ,4HTS/i   ,4HNCH,  , 4H/)    /
 DATA    nb2   / 16 /
 DATA    b2    / 4H(10X  ,4H,12H   ,4HPAPE   ,4HR si   ,4HZE =  ,  &
     4H,f5.  ,4H1,2H   ,4H x,f   ,4H5.1,   ,4H16H,  ,  &
     4H  pa  ,4HPER    ,4HTYPE   ,4H = ,   ,4H2A4,  , 4H/)    /
 DATA    nb5   / 10 /
 DATA    b5    / 4H(10X  ,4H,3HP   ,4HEN,i   ,4H2,7H   ,4H - s  ,  &
     4HIZE,  ,4HI2,2   ,4HH, ,   ,4H2A4,   ,4H/)    /
 
!     ELECTRONIC PLOTTER FORMATS.
 
 DATA    nc1   / 16 /
 DATA    c1    / 4H(10X  ,4H,37H   ,4HTHE    ,4HFOLL   ,4HOWIN   ,  &
     4HG pl  ,4HOTS    ,4HARE    ,4HREQU   ,4HESTE   ,  &
     4HD on  ,4H ,a4   ,4H,a1,   ,4H5H o   ,4HNLY,   , 4H/)    /
 DATA    nc2   / 17 /
 DATA    c2    / 4H(10X  ,4H,54H   ,4HTHE    ,4HFOLL   ,4HOWIN   ,  &
     4HG pl  ,4HOTS    ,4HARE    ,4HREQU   ,4HESTE   ,  &
     4HD on  ,4H bot   ,4HH fi   ,4HLM +   ,4H pap   , 4HER,/  ,4H)      /
 DATA    nc3   / 24 /
 DATA    c3    / 4H(10X  ,4H,i1,   ,4H79H    ,4HBLAN   ,4HK fr   ,  &
     4HAMES  ,4H wil   ,4HL be   ,4H ins   ,4HERTE   ,  &
     4HD on  ,4H fil   ,4HM on   ,4HLY b   ,4HETWE   ,  &
     4HEN e  ,4HACH    ,4HOF t   ,4HHE f   ,4HOLLO   ,  &
     4HWING  ,4H plo   ,4HTS,/   ,4H)      /
 
!     ENGINEERING DATA FORMATS.
 
 DATA    nd1   / 11 /
 DATA    d1    / 4H(//3  ,4H3H e   ,4H n g   ,4H i n   ,4H e e   ,  &
     4H r i  ,4H n g   ,4H       ,4HD a    ,4HT a,   , 4H/)    /
 DATA    nd2   / 6  /
 DATA    d2    / 4H(10X  ,4H,3A4   ,4H,11H   ,4H pro   ,4HJECT   , 4HION)  /
 DATA    nd3   / 29 /
 DATA    d3    / 4H(10X  ,4H,29H   ,4HROTA   ,4HTION   ,4HS (d   ,  &
     4HEGRE  ,4HES)    ,4H- ga   ,4HMMA    ,4H=,f7   ,  &
     4H.2,8  ,4HH, b   ,4HETA    ,4H=,f7   ,4H.2,9   ,  &
     4HH, a  ,4HLPHA   ,4H =,f   ,4H7.2,   ,4H10H,   ,  &
     4H  ax  ,4HES =   ,4H ,2A   ,4H1,2(   ,4H1H,,   ,  &
     4H2A1)  ,4H,2H,   ,4H ,4A   ,4H4)     /
 DATA    prj   , plus    ,minus    ,xyz      ,symm     ,BLANK    /  &
     4HORTH  ,4HOGRA   ,4HPHIC   ,4HPERS   ,4HPECT   ,  &
     4HIVE   ,4HSTER   ,4HEOSC   ,4HOPIC   ,1H+,1H-  ,  &
     1HX     ,1HY      ,1HZ      ,4HANTI   ,4HSYMM   ,  &
     4HETRI  ,1HC      ,1H       /
 
!     ORTHOGRAPHIC + PERSPECTIVE ENGINEERING DATA FORMATS.
 
 DATA    ne1   / 12 /
 DATA    e1    / 4H(10X  ,4H,29H   ,4HSCAL   ,4HE (o   ,4HBJEC   ,  &
     4HT-TO  ,4H-plo   ,4HT si   ,4HZE)    ,4H=,1P   , 4H,e13  ,4H.6)    /
 DATA    ne2   / 20 /
 DATA    e2    / 4H(10X  ,4H,29H   ,4HVANT   ,4HAGE    ,4HPOIN   ,  &
     4HT (i  ,4HNCHE   ,4HS) -   ,4H ro    ,4H=,1P   ,  &
     4H,e13  ,4H.6,6   ,4HH, s   ,4H0 =,   ,4HE13.   ,  &
     4H6,6H  ,4H, t0   ,4H =,e   ,4H13.6   ,4H)      /
 DATA    ne3   / 14 /
 DATA    e3    / 4H(10X  ,4H,38H   ,4HPROJ   ,4HECTI   ,4HON p   ,  &
     4HLANE  ,4H sep   ,4HARAT   ,4HION    ,4H(inc   ,  &
     4HHES)  ,4H =,1   ,4HP,e1   ,4H3.6)   /
 DATA    ne4   / 18 /
 DATA    e4    / 4H(10X  ,4H,6HO   ,4HRIGI   ,4HN,i8   ,4H,11H   ,  &
     4H   -  ,4H   x   ,4H0 =,   ,4H1P,e   ,4H14.6   ,  &
     4H,6H,  ,4H y0    ,4H=,e1   ,4H4.6,   ,4H5X,8   ,  &
     4HH(in  ,4HCHES   ,4H))      /
 
!     STEREO ENGINEERING DATA FORMATS.
 
 DATA    nf1   / 22 /
 DATA    f1    / 4H(10X  ,4H,30H   ,4HSCAL   ,4HES -   ,4H (mo   ,  &
     4HDEL-  ,4HTO-p   ,4HLOT    ,4HSIZE   ,4H =,1   ,  &
     4HP,e1  ,4H3.6,   ,4H25H,   ,4H  ob   ,4HJECT   ,  &
     4H-TO-  ,4HMODE   ,4HL si   ,4HZE =   ,4H,e13   , 4H.6,1  ,4HH))    /
 DATA    nf2   / 25 /
 DATA    f2    / 4H(10X  ,4H,29H   ,4HVANT   ,4HAGE    ,4HPOIN   ,  &
     4HT (i  ,4HNCHE   ,4HS) -   ,4H r0    ,4H=,1P   ,  &
     4H,e13  ,4H.6,9   ,4HH, s   ,4H0(l)   ,4H =,e   ,  &
     4H13.6  ,4H,9H,   ,4H s0(   ,4HR) =   ,4H,e13   ,  &
     4H.6,6  ,4HH, t   ,4H0 =,   ,4HE13.   ,4H6)     /
 DATA    nf3   / 12 /
 DATA    f3    / 4H(10X  ,4H,28H   ,4HOCUL   ,4HAR s   ,4HEPAR   ,  &
     4HATIO  ,4HN (i   ,4HNCHE   ,4HS) =   ,4H,1P,   , 4HE13.  ,4H6)     /
 DATA    nf4   / 23 /
 DATA    f4    / 4H(10X  ,4H,6HO   ,4HRIGI   ,4HN,i8   ,4H,14H   ,  &
     4H   -  ,4H   x   ,4H0(l)   ,4H =,1   ,4HP,e1   ,  &
     4H4.6,  ,4H9H,    ,4HX0(r   ,4H) =,   ,4HE14.   ,  &
     4H6,6H  ,4H, y0   ,4H =,e   ,4H14.6   ,4H,5X,   ,  &
     4H8H(i  ,4HNCHE   ,4HS))    /
 
!     CONTOUR PLOTTING DATA FORMATS
 
 DATA    ng1   / 13 /  &
     g1    / 4H(//4  ,4H2H c   ,4H o n   ,4H t o   ,4H u r   ,  &
     4H   p  ,4H l o   ,4H t t   ,4H i n   ,4H g     ,  &
     4H d a  ,4H t a   ,4H,/)    /
 DATA    ng2   / 11 /  &
     g2    / 4H(9X,  ,4H32HA   ,4HBOVE   ,4H plo   ,4HT is   ,  &
     4H a c  ,4HONTO   ,4HUR p   ,4HLOT    ,4HOF ,   , 4H4A4)  /
 DATA    ng3   / 16 /  &
     g3    / 4H(9X,  ,4H52HT   ,4HHE c   ,4HONTO   ,4HUR v   ,  &
     4HALUE  ,4HS ar   ,4HE ca   ,4HLCUL   ,4HATED   ,  &
     4H at   ,4HFIBR   ,4HE di   ,4HSTAN   ,4HCE ,   , 4H3A4)  /
 DATA    ng4   / 6  /  &
     g4    / 4H(9X,  ,4H4HIN   ,4H a,2   ,4HA4,6   ,4HHSYS   , 4HTEM)  /
 DATA    ng5   / 11 /  &
     g5    / 4H(//,  ,4H51X,   ,4H28HT   ,4HABLE   ,4H  of   ,  &
     4H  pl  ,4HOTTI   ,4HNG     ,4HSYMB   ,4HOLS,   , 4H/)    /
 DATA    ng6   / 8  /  &
     g6    / 4H(5(5  ,4HX,13   ,4HHSYM   ,4HBOL    ,4H val   ,  &
     4HUE,6  ,4HX),/   ,4H)      /
 DATA    ng7   / 4  / g7    / 4H(5(i  ,4H9,1P   ,4H,e15   ,4H.6))   /
 
 DATA    nh1   / 11 /
 DATA    h1    / 4H(//5  ,4H0X,2   ,4H9HPL   ,4HOT m   ,4HODUL   ,  &
     4HE me  ,4HSSAG   ,4HES c   ,4HONTI   ,4HNUE    , 4H,/)   /
 
 DATA   strain / 4HSTRA  ,4HIN e   ,4HNERG   ,4HIES    /,  &
     dist   / 4H z2   ,2*1H     ,4H z1    ,2*1H     ,4HMAX    ,  &
     4H- z1  ,4H,z2    ,4HAVER   ,4H-z1,   ,4HZ2     /,  &
     way    / 4H loc  ,4HAL     ,4H com   ,4HMON    /
 
 DATA   SPACE  / 4H      /  &
     displa / 4HDEFO  ,4HRMAT   ,4HION    ,1HX,1HY,1HZ, 3HMAG , 3*0     /
 
!                              1              3
 DATA   stress /       4HSTRE,4HSS,  ,4HSHEA,4HR -  ,
!               5 (1)          7 (2)          9 (3)         11 (4)  &
 4HMAJO,4HR-pr ,4HMINO,4HR-pr ,4HMAXI,4HMUM  ,4HNORM,4HAL x,
!              13 (5)         15 (6)         17     18      19  &
 4HNORM,4HAL y ,4HNORM,4HAL z ,4HXY  ,4HXZ   ,4HYZ   ,
!              20 (14)        22 (15)        24(16)         26 (17)  &
 4HNORM,4HAL 1 ,4HNORM,4HAL 2 ,4HSHEA,4HR 12 ,4HSHEA,4HR 1Z,
!              28 (18)        30 (19)        32             34  &
 4HSHEA,4HR 1Z ,4HBOND,4HSH12 ,4HLAYE,4HR nu ,4HMBER /
 
 DATA   ilay   / 4H  1 ,4H  2 ,4H  3 ,4H  4 ,4H  5 ,4H  6 ,  &
     4H  7 ,4H  8 ,4H  9 ,4H 10 ,4H 11 ,4H 12 ,  &
     4H 13 ,4H 14 ,4H 15 ,4H 16 ,4H 17 ,4H 18 , 4H 19 ,4H 20 /
 
 IF (ncntr > 0) GO TO 201
 
!     PRINT THE PLOTTER ID.
 
 list(1) = 0
 CALL WRITE  (prnt,list,1,0)
 CALL wrtprt (prnt,list,a1,na1)
 
!     NASTRAN GENERAL PURPOSE PLOTTER.
 
 list(1) = 5
 j = IABS(pltype)
 DO  i = 1,3
   list(i+1) = plttyp(i,j)
 END DO
 mm = 1
 IF (pltype < 0) mm = 3
 list(5) = itype(mm  )
 list(6) = itype(mm+1)
 CALL wrtprt (prnt,list,p6,np6)
 
!     GENERAL PLOTTER INFORMATION.
 
 IF (tapden <= 0) GO TO 151
 list(1) = 1
 list(2) = tapden
 CALL wrtprt (prnt,list,a2,na2)
 151 IF (eof /= 0) GO TO 152
 CALL wrtprt (prnt,0,a3,na3)
 GO TO 154
 152 CALL wrtprt (prnt,0,a4,na4)
 154 CALL wrtprt (prnt,0,a5,na5)
 list(1) = 1
 alst(2) = cscale
 CALL wrtprt (prnt,list,a6,na6)
 IF (IABS(pltype)-2 < 0.0) THEN
   GO TO   170
 ELSE IF (IABS(pltype)-2 == 0.0) THEN
   GO TO   160
 ELSE
   GO TO   163
 END IF
 
!     TABLE PLOTTER INFORMATION.
 
 160 list(1) = 1
 alst(2) = cntsin
 CALL wrtprt (prnt,list,b1,nb1)
 163 list(1) = 4
 alst(2) = papsiz(1)
 alst(3) = papsiz(2)
 list(4) = paptyp(1)
 list(5) = paptyp(2)
 CALL wrtprt (prnt,list,b2,nb2)
 
 list(1) = 4
 n = MIN0(npens,nopens)
 DO  i = 1,n
   list(2) = i
   list(3) = pensiz(i)
   IF (list(3) < 0) CYCLE
   list(4) = penclr(i,1)
   list(5) = penclr(i,2)
   IF (list(4) == BLANK .AND. list(5) == BLANK) CYCLE
   CALL wrtprt (prnt,list,b5,nb5)
 END DO
 CALL wrtprt (prnt,0,skip,nskip)
 GO TO 180
 
!     ELECTRONIC PLOTTER INFORMATION.
 
 170 IF (camera-2 < 0.0) THEN
   GO TO   171
 ELSE IF (camera-2 == 0.0) THEN
   GO TO   172
 ELSE
   GO TO   174
 END IF
 171 list(2) = film(1)
 list(3) = film(2)
 GO TO 173
 172 list(2) = paper(1)
 list(3) = paper(2)
 173 list(1) = 2
 CALL wrtprt (prnt,list,c1,nc1)
 GO TO 175
 174 CALL wrtprt (prnt,0,c2,nc2)
 175 IF (camera == 2 .OR. bframs == 0) GO TO 180
 list(1) = 1
 list(2) = bframs
 CALL wrtprt (prnt,list,c3,nc3)
 
!     ENGINEERING DATA.
 
 180 CALL wrtprt (prnt,0,d1,nd1)
 list(1) = 3
 DO  i = 1,3
   list(i+1) = prj(i,prject)
 END DO
 CALL wrtprt (prnt,list,d2,nd2)
 
 list(1) = 13
 alst(2) = vangle(3)
 IF (vangle(2) > -1.e10) GO TO 1815
 IF (prject /= 2) vangle(2) = vangle(4)
 IF (prject == 2) vangle(2) = vangle(5)
 1815 alst(3) = vangle(2)
 alst(4) = vangle(1)
 DO  i = 1,3
   j = 2*i + 3
   k = IABS(axis(i))
   list(j) = plus
   IF (axis(i) < 0) list(j) = minus
   list(j+1) = xyz(k)
 END DO
 n = 1
 IF (axis(1) == daxis(1)) n = 2
 list(14) = BLANK
 j = 1
 DO  i = n,4
   list(j+10) = symm(i)
   j = j + 1
 END DO
 CALL wrtprt (prnt,list,d3,nd3)
 IF (prject == 3) GO TO 195
 
!     ORTHOGRAPHIC + PERSPECTIVE ENGINEERING DATA.
 
 list(1) = 1
 alst(2) = scale(1)/cntsin
 CALL wrtprt (prnt,list,e1,ne1)
 IF (prject == 1) GO TO 191
 list(1) = 3
 alst(2) = r0
 alst(3) = s0l
 alst(4) = t0
 CALL wrtprt (prnt,list,e2,ne2)
 list(1) = 1
 alst(2) = d0
 CALL wrtprt (prnt,list,e3,ne3)
 
 191 CALL wrtprt (prnt,0,skip,nskip)
 list(1) = 3
 DO  i = 1,org
   list(2) = origin(i)
   alst(3) = xy(i,1)/cntsin
   alst(4) = xy(i,3)/cntsin
   CALL wrtprt (prnt,list,e4,ne4)
 END DO
 GO TO 260
 
!     STEREO ENGINEERING DATA.
 
 195 list(1) = 2
 alst(2) = scale(1)/cntin3
 alst(3) = scale(2)
 CALL wrtprt (prnt,list,f1,nf1)
 list(1) = 4
 alst(2) = r0
 alst(3) = s0l
 alst(4) = s0r
 alst(5) = t0
 CALL wrtprt (prnt,list,f2,nf2)
 list(1) = 1
 alst(2) = d0
 CALL wrtprt (prnt,list,e3,ne3)
 alst(2) = s0s
 CALL wrtprt (prnt,list,f3,nf3)
 
 CALL wrtprt (prnt,0,skip,nskip)
 list(1) = 4
 DO  i = 1,org
   list(2) = origin(i)
   alst(3) = xy(i,1)/cntsin
   alst(4) = xy(i,2)/cntsin
   alst(5) = xy(i,3)/cntsin
   CALL wrtprt (prnt,list,f4,nf4)
 END DO
 GO TO 260
 
!     CONTOUR PLOTTING DATA
 
 201 list(1) = 0
 CALL wrtprt (prnt,list,g1,ng1)
 list(1) = 4
 IF (icntvl > 9 .AND. icntvl < 14) GO TO 210
 
!     STRESS CONTOURS
 
 i = 1
 IF (icntvl > 6 .OR. icntvl == 3) i = 3
 IF (icntvl >= 14 .AND.icntvl <= 19) i = 1
 IF (icntvl /= 20) GO TO 203
 
!     STRAIN CONTOURS
 
 list(1) = 4
 list(2) = strain(1)
 list(3) = strain(2)
 list(4) = strain(3)
 list(5) = strain(4)
 CALL wrtprt (prnt,list,g2,ng2)
 GO TO 205
 
 203 list(2) = stress(i)
 list(3) = stress(i+1)
 i = icntvl*2 + 3
 IF (icntvl > 13 .AND. icntvl < 20) i = (icntvl-14)*2 + 20
 IF (icntvl > 6 .AND. icntvl <= 9) i = icntvl + 10
 list(4) = stress(i)
 list(5) = SPACE
 IF (icntvl < 7 .OR. icntvl >= 14) list(5) = stress(i+1)
 CALL wrtprt (prnt,list,g2,ng2)
 
!     ADDING LAYER NUMBER TO OUTPUT WHEN REQUESTED
 
 IF (icntvl < 14 .OR. icntvl == 20) GO TO 204
 list(1) = 4
 list(2) = stress(32)
 list(3) = stress(33)
 list(4) = stress(34)
 list(5) = ilay(layer)
 CALL wrtprt (prnt,list,g2,ng2)
 GO TO 205
 
 204 list(1) = 3
 i = iwhere
 IF (iwhere <= 0) i = 0
 i = i*3 + 1
 list(2) = dist(i)
 list(3) = dist(i+1)
 list(4) = dist(i+2)
 CALL wrtprt (prnt,list,g3,ng3)
 
 205 j = 2*(idirec-1) + 1
 GO TO 220
 
!     DISPLACEMENT CONTOURS
 
 210 i = 1
 list(2) = displa(i  )
 list(3) = displa(i+1)
 list(4) = displa(i+2)
 list(5) = displa(icntvl-6)
 CALL wrtprt (prnt,list,g2,ng2)
 j = 3
 220 IF (icntvl < 4 .OR. icntvl == 13) GO TO 235
 list(1) = 2
 list(2) = way(j  )
 list(3) = way(j+1)
 CALL wrtprt (prnt,list,g4,ng4)
 235 list(1) = 0
 CALL wrtprt (prnt,list,g5,ng5)
 CALL wrtprt (prnt,list,g6,ng6)
 l = (ncntr-1)/10 + 1
 list(1) = 2*l
 k = MIN0(ncntr,10)
 DO  j = 1,k
   n = j + (l-1)*10
   m = 2
   DO  i = j,n,10
     IF (i > ncntr) GO TO 245
     list(m  ) = i
     alst(m+1) = cntr(i)
     m = m + 2
   END DO
   GO TO 247
   245 list(1) = list(1) - 2
   l = l - 1
   247 CALL wrtprt (prnt,list,g7,ng7)
 END DO
 
 260 list(1) = 0
 CALL wrtprt (prnt,list,h1,nh1)
 RETURN
END SUBROUTINE pltopr
