SUBROUTINE xdph
     
!     DATA POOL HOUSEKEEPER (XDPH)
 
!     THIS SUBROUTINE SCANS THE DATA POOL DICT AND TO DETERMINE THE
!     NUMBER AND SIZE OF ANY FILES NO LONGER NEEDED.  IF A SUFFICIENT
!     QUANTITY IS NOT NEEDED, THE FILE IS RECOPIED WITH THE DEAD FILES
!     DELETED.
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        rshift,andf,orf
 DIMENSION       ndpd(1),ndph(2),fequ(1),fntu(1),fon(1),ford(1),  &
     minp(1),mlsn(1),mout(1),mscr(1),sal(1),sdbn(1), sntu(1),sord(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ ibufsz,outtap,dum(36),nbpc,nbpw,ncpw
 COMMON /xfiat / fiat(1),fmxlg,fculg,FILE(1),fdbn(2),fmat(1)
 COMMON /xfist / fist(2)
 COMMON /xpfist/ npfist
 COMMON /xxfiat/ exfiat
 COMMON /xdpl  / dpd(1),dmxlg,dculg,ddbn(2),dfnu(1)
 COMMON /zzzzzz/ endsfa(1)
 COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(1),fpun(1),  &
     fcum(1),fcus(1),fknd(1)
 EQUIVALENCE     (dpd(1),dnaf),(fiat(1),funlg),(FILE(1),fequ(1)),  &
     (FILE(1),ford(1)),(endsfa(1),ndpd(1))
 EQUIVALENCE     (md(2),mlsn(1)),(md(3),minp(1)),(md(4),mout(1)),  &
     (md(5),mscr(1)), (sos(1),slgn)  ,(sos(2),sdbn(1)),(sos(4),sal(1)),  &
     (sos(4),sntu(1)),(sos(4),sord(1)),  &
     (comm(1),almsk),(comm(2),apndmk),(comm(3),cursno),  &
     (comm(4),entn1),(comm(5),entn2 ),(comm(6),entn3 ),  &
     (comm(7),entn4),(comm(8),flag  ),(comm(9),fnx   ),  &
     (comm(10),lmsk),(comm(11),lxmsk),  &
     (comm(13),rmsk),(comm(14),rxmsk),(comm(15),s    ),  &
     (comm(16),scornt),(comm(17),tapmsk), (comm(18),thcrmk),(comm(19),zap),  &
     (xf1at(1),fntu(1)),(xf1at(1),fon(1))
 
 DATA    nconst/ 100    /
 DATA    scrn1 / 4HSCRA /, scrn2 /4HTCH* /
 DATA    pool  , npol   /  4HPOOL,4HNPOL /,  ndph / 4HXDPH,4H    /
 
 
 flag = 0
 100 lmt3 = dculg*entn4
 lmt  = (dculg-1)*entn4 + 1
 ncnt = 0
 ngcnt= 0
 trial= dnaf - 1
 
!     COUNT DEAD FILE SIZE, PUT SIZE IN NCNT
 
 DO  i = 1,lmt3,entn4
   IF (ddbn(i) /= 0 .OR. ddbn(i+1) /= 0) GO TO 159
   IF (dfnu(i) >= 0) GO TO 130
   
!     DEAD FILE IS EQUIV
   
   flag = -1
   kk = andf(rmsk,dfnu(i))
   DO  j = 1,lmt3,entn4
     IF (dfnu(j) >= 0 .OR. i == j) CYCLE
     IF (kk /= andf(rmsk,dfnu(j))) CYCLE
     IF (ddbn(j) /= 0 .OR. ddbn(j+1) /= 0) GO TO 145
     dfnu(j) = 0
   END DO
   130 IF (kk  == trial) GO TO 140
   IF (dfnu(i) == 0) GO TO 150
   ncnt = ncnt + rshift(andf(lmsk,dfnu(i)),16)
   GO TO 150
   140 dnaf = trial
   145 dfnu(i) = 0
   150 IF (i /= lmt) CYCLE
   dculg = dculg - 1
   flag  = -1
   GO TO 100
   
!     COUNT GOOD STUFF ALSO
   
   159 ngcnt = ngcnt + rshift(andf(lmsk,dfnu(i)),16)
 END DO
 
!     CHECK FOR BREAKING OF EQUIV
 
 IF (flag == 0) GO TO 200
 loop180:  DO  i = 1,lmt3,entn4
   IF (dfnu(i) >= 0) CYCLE loop180
   kk = andf(rmsk,dfnu(i))
   DO  j = 1,lmt3,entn4
     IF (dfnu(j) >= 0  .OR.  i == j) CYCLE
     IF (kk == andf(rmsk,dfnu(j))) CYCLE loop180
   END DO
   dfnu(i) = andf(almsk,dfnu(i))
 END DO loop180
 
!     IS NCNT OF SUFFICIENT SIZE TO WARRANT RECOPYING POOL
 
 200 CALL sswtch (3,ix)
 IF (ix /= 1) GO TO 211
 CALL page1
 WRITE  (outtap,201) ncnt
 201 FORMAT (21H0DPH dead FILE count=,i6)
 WRITE  (outtap,202)(dpd(ix),ix=1,3)
 202 FORMAT (16H0DPD before dph ,3I4)
 ii = dculg*3 + 3
 DO  ix = 4,ii,3
   iprt1 = rshift(dpd(ix+2),nbpw-1)
   iprt2 = rshift(andf(lxmsk,dpd(ix+2)),16)
   iprt3 = andf(rxmsk,dpd(ix+2))
   203 FORMAT (1H ,2A4,3I6)
   WRITE  (outtap,203) dpd(ix),dpd(ix+1),iprt1,iprt2,iprt3
 END DO
 
!     RECOPY POOL IF THERE ARE MORE THAN 500,000 WORD DEAD AND
!     THE GOOD STUFF IS TWICE AS BIG AS THE DEAD STUFF
 
 211 IF (ncnt > nconst .AND. ncnt > 2*ngcnt) GO TO 230
 IF (ncnt > 0 .AND. dculg+5 >= dmxlg) GO TO 230
 RETURN
 
!     RECOPY POOL, SWITCH POOL FILE POINTERS
 
 230 lmt2 = funlg*entn1
 kk   = andf(thcrmk,scrn2)
 DO  i = 1,lmt2,entn1
   IF (fdbn(i) == 0 .AND. fdbn(i+1) == 0) GO TO 270
   IF (fdbn(i) == scrn1 .AND. andf(thcrmk,fdbn(i+1)) == kk) GO TO 270
 END DO
 
!     NO FILE AVAILABLE TO COPY ONTO, FORGET IT
 
 RETURN
 
!     SET-UP FOR A RECOPY
 
 270 isav = i
 CALL OPEN (*900,pool,endsfa,0)
 fnx = 1
 fist(2*npfist+4) = isav + 2
 fist(2) = npfist + 1
 fist(2*npfist+3) = npol
 CALL OPEN (*900,npol,endsfa(ibufsz+1),1)
 m = 2*ibufsz
 i = m + 1
 istart = i
 m = m + dculg*3 + 3
 iwkbuf = korsz(endsfa) - m
 IF (iwkbuf < 100) CALL mesage (-8,0,ndph)
 m = m + 1
 nfile = 1
 nculg = 0
 DO  j = 1,lmt3,entn4
   IF (ddbn(j) == 0 .AND. ddbn(j+1) == 0) CYCLE
   IF (ddbn(j) == 63 .AND. ddbn(j+1) == 63) CYCLE
   
!     RECOPY DICTIONARY
   
   ndpd(i  ) = ddbn(j  )
   ndpd(i+1) = ddbn(j+1)
   ndpd(i+2) = orf(andf(lxmsk,dfnu(j)),nfile)
   IF (dfnu(j) >= 0) GO TO 290
   ndpd(i+2) = orf(s,ndpd(i+2))
   kk = andf(rmsk,dfnu(j))
   DO  k  = 1,lmt3,entn4
     IF (dfnu(k) >= 0 .OR. j == k) CYCLE
     IF (kk /= andf(rmsk,dfnu(k))) CYCLE
     i = i + 3
     nculg   = nculg + 1
     ndpd(i) = ddbn(k)
     ddbn(k) = 63
     ndpd(i+1) = ddbn(k+1)
     ddbn(k+1) = 63
     ndpd(i+2) = ndpd(i-1)
   END DO
   290 i = i + 3
   nculg = nculg + 1
   
!     RECOPY NECESSARY FILE
   
   fn = andf(rmsk,dfnu(j))
   CALL xfilps (fn)
   CALL cpyfil (pool,npol,endsfa(m),iwkbuf,flag)
   CALL eof (npol)
   nfile = nfile + 1
   fnx   = fn + 1
 END DO
 
!     COPY TEMPORARY DPD INTO ACTUAL DPD
 
 i  = i - 1
 ix = 0
 DO  j = istart,i
   ix = ix + 1
   ddbn(ix) = ndpd(j)
 END DO
 dnaf = nfile
 dculg= nculg
 CALL CLOSE (pool,1)
 CALL CLOSE (npol,1)
 fnx = 1
 
!     COPY POOL BACK TO POOL UNIT
 
 CALL OPEN (*900,npol,endsfa,0)
 CALL OPEN (*900,pool,endsfa(ibufsz+1),1)
 nfile = nfile - 1
 DO  ix = 1,nfile
   CALL cpyfil (npol,pool,endsfa(m),iwkbuf,flag)
   CALL eof (pool)
 END DO
 CALL CLOSE (pool,1)
 CALL CLOSE (npol,1)
 
!     THE FOLLOWING 3 LINES OF CODE WILL FREE DISK AREA ON SOME CONFIG.
 
 CALL OPEN  (*900,npol,endsfa,1)
 CALL WRITE (npol,ndph,2,1)
 CALL CLOSE (npol,1)
 CALL sswtch (3,ix)
 IF (ix /= 1) RETURN
 
 WRITE  (outtap,500) (dpd(ix),ix=1,3)
 500 FORMAT (15H0DPD after dph ,3I4)
 ii = dculg*3 + 3
 DO  ix = 4,ii,3
   iprt1 = rshift(dpd(ix+2),nbpw-1)
   iprt2 = rshift(andf(lxmsk,dpd(ix+2)),16)
   iprt3 = andf(rxmsk,dpd(ix+2))
   WRITE (outtap,203) dpd(ix),dpd(ix+1),iprt1,iprt2,iprt3
 END DO
 RETURN
 
 900 WRITE  (outtap,901) sfm
 901 FORMAT (a25,' 1041, OLD/NEW POOL COULD NOT BE OPENED.')
 CALL mesage (-37,0,ndph)
 RETURN
END SUBROUTINE xdph
