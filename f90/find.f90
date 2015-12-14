SUBROUTINE find (mode,buf1,buf4,setid,x)
     
 
 INTEGER, INTENT(IN OUT)                  :: mode
 INTEGER, INTENT(IN OUT)                  :: buf1
 INTEGER, INTENT(IN)                      :: buf4
 INTEGER, INTENT(IN)                      :: setid(1)
 INTEGER, INTENT(IN OUT)                  :: x(1)
 INTEGER :: awrd(2), bufsiz,ERR(3),for,fscale,fvp,  &
     gpset,origin,org,parm,prject,prnt,region,set,  &
     setd, tra,word, hset,orig,poin,regi,  &
     scal,vant,msg1(20),msg3(21),msg6(20),NAME(2)
 REAL :: imsep,MAX,maxdef,MIN,mm17p5
 DOUBLE PRECISION :: dwrd
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ bufsiz, nout
 COMMON /BLANK / ngp,skp11,nsets,prnt,skp12,ngpset,skp13(4),  &
     parm,gpset,skp2(8),merr,setd
 COMMON /xxparm/ pltbuf,pltter(5),nopens,papsiz(2),penpap(27),  &
     scale,objmod,fscale,maxdef,defmax,axis(6),view(9),  &
     fvp,skpvp1(4),d0,skpvp2(2),prject,s0s,for,org,  &
     norg,origin(11),edge(11,4),xy(11,3)
 COMMON /pltdat/ skpplt(2),reg(4),axymax(14),skpa(3),cntchr(2)
 COMMON /rstxxx/ cstm(3,3),MIN(3),MAX(3),d(3),aver(3)
 EQUIVALENCE     (word,awrd(1),iwrd,fwrd,dwrd)
 DATA    NAME  / 4H  fi, 4HND  /
 DATA    mm17p5, rdist,  sqrt3 / .688975, 29., 1.732051/,  &
     orig  / 4HORIG/, regi / 4HREGI/, scal / 4HSCAL/,  &
     hset  / 3HSET /, vant / 4HVANT/, poin / 4HPOIN/
 DATA    nmsg1 , msg1  / 20,  &
     4H(34X, 4H,45H, 4HAN a, 4HTTEM, 4HPT h, 4HAS b,  &
     4HEEN , 4HMADE, 4H TO , 4HDEFI, 4HNE m, 4HORE ,  &
     4HTHAN, 4H ,i2, 4H,17H, 4H dis, 4HTINC, 4HT OR, 4HIGIN, 4HS)     /
 DATA    nmsg3 , msg3  / 21,  &
     4H(25X, 4H,27H, 4HAN u, 4HNREC, 4HOGNI, 4HZABL,  &
     4HE re, 4HQUES, 4HT (,, 4H2A4,, 4H37H), 4H has,  &
     4H bee, 4HN sp, 4HECIF, 4HIED , 4HON a, 4H -fi, 4HND- , 4HCARD, 4H)   /
 DATA    nmsg6 , msg6  / 20,  &
     4H(33X, 4H,71H, 4HMAXI, 4HMUM , 4HDEFO, 4HRMAT,  &
     4HION , 4HCARD, 4H nee, 4HDED , 4H- 5 , 4HPER ,  &
     4HCENT, 4H of , 4HMAXI, 4HMUM , 4HDIME, 4HNSIO, 4HN us, 4HED.)  /
 
 CALL rdmodx (parm,mode,word)
 set    = setd
 region = 0
 reg(1) = 0.
 reg(2) = 0.
 reg(3) = 1.
 reg(4) = 1.
 ratio  = 0.
 nogo   = 0
 IF (mode < 0) GO TO 480
 
!     INTERPRET THE REQUESTS ON THE -FIND- CARD.
 
 10 IF (mode <= 0) CALL rdmode (*10,*20,*480,mode,word)
 20 CALL rdword (mode,word)
 
!     IS AN ORIGIN TO BE FOUND
 
 30 IF (word /= orig) GO TO 90
 IF (mode /=    0) GO TO 10
 ASSIGN 40 TO tra
 GO TO 400
 40 IF (org == 0) GO TO 70
 DO  j = 1,org
   IF (origin(j) == iwrd) GO TO 80
 END DO
 IF (org  < norg) GO TO 70
 IF (prnt <    0) GO TO 60
 ERR(1) = 1
 ERR(2) = norg
 CALL wrtprt (merr,ERR,msg1,nmsg1)
 60 org = norg
 i   = org + 1
 edge(i,1) = 0.0
 edge(i,2) = 0.0
 edge(i,3) = 1.0
 edge(i,4) = 1.0
 70 org = org + 1
 origin(org) = iwrd
 j   = org
 80 for = j
 GO TO 10
 
!     IS A REGION SPECIFIED
 
 90 IF (word /= regi) GO TO 200
 IF (mode /=    0) GO TO 10
 region = 1
 ASSIGN 110 TO tra
 j = 0
 100 j = j + 1
 GO TO 440
 110 reg(j) = AMIN1(1.,ABS(fwrd))
 IF (j-4 < 0) THEN
   GO TO   100
 ELSE
   GO TO    10
 END IF
 
!     IS THE SCALE TO BE FOUND
 
 200 IF (word /= scal) GO TO 220
 fscale = 1
 IF (mode /= 0) GO TO 10
 ASSIGN 210 TO tra
 GO TO 440
 210 ratio = fwrd
 GO TO 10
 
!     IS THERE A SET ON THE FIND CARD
 
 220 IF (word /= hset) GO TO 300
 IF (mode /=    0) GO TO 10
 ASSIGN 230 TO tra
 GO TO 400
 230 DO  j = 1,nsets
   IF (iwrd == setid(j)) GO TO 260
 END DO
 WRITE  (nout,250) uwm,iwrd
 250 FORMAT (a25,' 700, SET',i9,' REQUESTED ON FIND CARD HAS NOT BEEN',  &
     ' DEFINED. DEFAULT SET',i9,' USED')
 nogo = 1
 GO TO 10
 260 set  = j
 GO TO 10
 
!     IS THE VANTAGE POINT TO BE FOUND
 
 300 IF (word /= vant) GO TO 320
 IF (mode == 0) CALL rdmode (*10,*310,*480,mode,word)
 310 CALL rdword (mode,word)
 IF (word /= poin) GO TO 30
 fvp = 1
 GO TO 10
 
!     UNRECOGNIZABLE OPTION ON THE FIND CARD
 
 320 IF (prnt < 0) GO TO 10
 ERR(1) = 2
 ERR(2) = awrd(1)
 ERR(3) = awrd(2)
 CALL wrtprt (merr,ERR,msg3,nmsg3)
 GO TO 10
 
!     READ AN INTEGER FROM THE FIND CARD
 
 400 CALL rdmode (*410,*10,*480,mode,word)
 410 IF (mode == -1) GO TO 430
 IF (mode == -4) GO TO 420
 iwrd = fwrd
 GO TO 430
 420 iwrd = dwrd
 430 GO TO tra, (40,230)
 
!     READ A REAL NUMBER FROM THE FIND CARD
 
 440 CALL rdmode (*450,*10,*480,mode,word)
 450 IF (mode == -4) GO TO 460
 IF (mode /= -1) GO TO 470
 fwrd = iwrd
 GO TO 470
 460 fwrd = dwrd
 470 GO TO tra, (110,210)
 
!     END OF THE FIND CARD
 
 480 IF (org > 0) GO TO 485
 
!     ALLOW NO ORIGIN REQUEST ON FIRST FIND CARD
!     ORIGIN ID IS ZERO
 
 org = 1
 origin(1) = 0
 region = 1
 485 IF (for    == 0) GO TO 500
 IF (region == 0) GO TO 490
 edge(for,1) = reg(1)
 edge(for,2) = reg(2)
 edge(for,3) = reg(3)
 edge(for,4) = reg(4)
 GO TO 500
 490 reg(1) = edge(for,1)
 reg(2) = edge(for,2)
 reg(3) = edge(for,3)
 reg(4) = edge(for,4)
 500 reg(1) = reg(1)*axymax(1)
 IF (reg(2) /= 0.) GO TO 510
 reg(2) = 4.*cntchr(2)
 GO TO 520
 510 reg(2) = reg(2)*axymax(2)
 520 reg(3) = reg(3)*axymax(1) - cntchr(1)*8.
 reg(4) = reg(4)*axymax(2) - cntchr(2)
 
!     CALCULATE THE ROTATION MATRIX + ROTATE THE CO-ORDINATES OF THE SET
 
 CALL gopen (gpset,x(buf4),0)
 i = 1
 CALL fwdrec (*810,gpset)
 IF (set == 1) GO TO 540
 DO  i = 2,set
   CALL fwdrec (*810,gpset)
 END DO
 
!     READ NGPSET
 
 540 CALL fread (gpset,ngpset,1,0)
 
!     CHECK CORE
 
 icrq = 3*ngpset + ngp - buf4 - bufsiz - 1
 IF (icrq > 0) GO TO 800
 CALL fread (gpset,x,ngp,0)
 CALL CLOSE (gpset,1)
 CALL fndset (x,x(ngp+1),buf1,0)
 DO  i = 1,3
   MIN(i) = +1.e+20
   MAX(i) = -1.e+20
 END DO
 CALL proces (x(ngp+1))
 IF (maxdef /= 0.0 .OR. prnt >= 0) GO TO 560
 
!     DEFORMED PLOTS AND MAXDEF WAS NOT SPECIFIED
 
 ERR(1) = 0
 CALL wrtprt (merr,ERR,msg6,nmsg6)
 maxdef = AMAX1(d(2),d(3))
 IF (maxdef <= 0.0) maxdef = 1.0
 maxdef = 0.05*maxdef
 560 CONTINUE
 SELECT CASE ( prject )
   CASE (    1)
     GO TO 600
   CASE (    2)
     GO TO 570
   CASE (    3)
     GO TO 700
 END SELECT
 
!     PERSPECTIVE PROJECTION (FIND VANTAGE POINT IF REQUESTED)
 
 570 DO  i = 1,3
   MIN(i) = +1.e+20
   MAX(i) = -1.e+20
 END DO
 CALL perpec (x(ngp+1),0)
 fvp = 0
 
!     ORTHOGRAPHIC OR PERSPECTIVE PROJECTION
 
!     FIND SCALE FACTOR (IF REQUESTED).
 
 600 IF (fscale == 0) GO TO 630
 a = d(2) + 2.*maxdef*sqrt3
 IF (a == 0.0) GO TO 610
 a = (reg(3)-reg(1))/a
 610 b = d(3) + 2.*maxdef*sqrt3
 IF (b == 0.0) GO TO 620
 b = (reg(4)-reg(2))/b
 620 scale = AMIN1(a,b)
 IF (scale <= 0.) scale = AMAX1(a,b)
 IF (scale <= 0.) scale = 1.
 IF (ratio /= 0.) scale = ratio*scale
 
!     FIND ORIGIN -FOR- IF REQUESTED
 
 630 IF (for == 0) GO TO 830
 xy(for,1) = aver(2)*scale - (reg(1)+reg(3))/2.
 xy(for,3) = aver(3)*scale - (reg(2)+reg(4))/2.
 GO TO 830
 
!     STEREO PROJECTION
 
!     FIND SCALE FACTORS (IF REQUESTED).
 
 700 IF (fscale == 0) GO TO 710
 diam = SQRT(d(1)**2 + d(2)**2 + d(3)**2)
 a = sqrt3*maxdef
 IF (d(2)+a >= diam .OR. d(3)+a >= diam) diam = diam + maxdef
 IF (diam == 0.0) diam = 1.e-5
 objmod = 10./diam
 scale  = AMIN1(reg(3)-reg(1),reg(4)-reg(2))/mm17p5
 IF (ratio /= 0.) scale=ratio*scale
 
!     FIND VANTAGE POINT (IF REQUESTED)
 
 710 CALL perpec (x(ngp+1),0)
 fvp = 0
 
!     FIND ORIGIN -FOR- IF REQUESTED
 
 IF (for == 0) GO TO 830
 imsep     = s0s*(rdist-d0)/(2.*rdist)
 xy(for,1) = scale*(aver(2)*objmod-imsep) - (reg(1)+reg(3))/2.
 xy(for,2) = scale*(aver(2)*objmod+imsep) - (reg(1)+reg(3))/2.
 xy(for,3) = scale*(aver(3)*objmod)       - (reg(2)+reg(4))/2.
 GO TO 830
 
 800 CALL mesage (-8,icrq,NAME)
 
 810 WRITE  (nout,820) ufm,setid(set)
 820 FORMAT (a23,' 703, SET',i9,' REQUESTED ON FIND CARD NOT IN ',  &
     'GPSETS FILE.')
 nogo = 1
 CALL CLOSE (gpset,1)
 GO TO 840
 
 830 fscale = 0
 for    = 0
 840 IF (nogo /= 0) CALL mesage (-37,0,NAME)
 RETURN
END SUBROUTINE find
