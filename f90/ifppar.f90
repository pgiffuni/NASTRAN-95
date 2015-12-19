SUBROUTINE ifppar
     
!     SUBROUTINE TO TEST FOR PARAM CARD PARAMETERS REQUIRED BY VARIOUS
!     RIGID FORMATS.
 
 LOGICAL :: abort,hfreq,lfreq,lmode,nodje,p1,p2,p3,ptot,  &
     ctype,kindx,nsegs,ltest,queue
 INTEGER :: rf,app,hfre,ctyp,que,appr(4)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /system/ n1,nout,abort,n2(17),iapp,n3(3),rf
 COMMON /ifpdta/ idta(509),nparam
 COMMON /zzzzzz/ ibuff(1)
 DATA    appr  / 4HDMAP,   4HDISP, 4HHEAT,   4HAERO        /
 DATA    hfre  / 4HHFRE /, lfre /  4HLFRE /, lmod  /4HLMOD /
 DATA    nodj  / 4HNODJ /, ip1  /  4HP1   /, ip2   /4HP2   /
 DATA    ip3   / 4HP3   /, que  /  4HQ                     /
 DATA    ctyp  / 4HCTYP /, kind /  4HKIND /, nseg  /4HNSEG /
 DATA    hfreq / .false./, lfreq/  .false./, lmode /.false./
 DATA    nodje / .false./, ctype/  .false./, kindx /.false./
 DATA    nsegs / .false./, p1   /  .false./, p2    /.false./
 DATA    queue / .false./, p3   /  .false.                 /
 
 app = IABS(iapp)
 
!     NO PARAMS REQD FOR HEAT APPROACH,
!     DMAPS DISP 1 THRU 9, DISP 13, AND DISP 16 THRU 19, AND
!     AERO RF 9
 
 IF (app == 1 .OR.  app == 3) GO TO 9999
 IF (app == 2 .AND. (rf <= 9 .OR. rf == 13 .OR. rf >= 16)) GO TO 9999
 IF (app == 4 .AND.  rf == 9) GO TO 9999
 
!     FATAL ERROR IF NO PARAMS ENTERED AS REQUIRED
 
 IF (nparam == 0) GO TO 9800
 
!     LOOP TO TEST PARAMS IN PVT FOR PRESENCE OF REQUIRED ONES.
 
 ipm = 1
 500 ipn = 2*n1 + ipm
 
 IF (rf >= 14) GO TO 1000
 IF (ibuff(ipn) == hfre) hfreq = .true.
 IF (ibuff(ipn) == lfre) lfreq = .true.
 IF (ibuff(ipn) == lmod .AND. ibuff(ipn+2) /= 0) lmode = .true.
 
 IF (app /= 4) GO TO 2000
 IF (ibuff(ipn) == nodj .AND. ibuff(ipn+2) /= 0) nodje = .true.
 IF (ibuff(ipn) == ip1) p1 = .true.
 IF (ibuff(ipn) == ip2) p2 = .true.
 IF (ibuff(ipn) == ip3) p3 = .true.
 IF (ibuff(ipn) == que  .AND. ibuff(ipn+2) /= 0) queue = .true.
 GO TO 2000
 
 1000 IF (ibuff(ipn) == ctyp .AND. ibuff(ipn+2) /= 0) ctype = .true.
 IF (ibuff(ipn) == nseg .AND. ibuff(ipn+2) /= 0) nsegs = .true.
 IF (ibuff(ipn) == kind .AND. ibuff(ipn+2) /= 0) kindx = .true.
 
 2000 ipm = ipm + 4
 IF (ibuff(ipn+2) >= 3 .AND. ibuff(ipn+2) <= 5) ipm = ipm + 1
 IF (ibuff(ipn+2) >= 6) ipm = ipm + 3
 IF (ipm < nparam) GO TO 500
 
!     TEST TO VERIFY THAT ALL REQUIRED PARAMS ARE PRESENT
 
 IF (rf == 14 .OR. rf == 15) GO TO 4000
 IF (lmode .AND.  .NOT.(hfreq.OR.lfreq)) GO TO 3000
 IF (hfreq .AND. lfreq .AND. .NOT.lmode) GO TO 3000
 
!     SOMETING AMISS - - IS AN LMODES, HFREQ, OR LFREQ MISSING
 
 IF (.NOT.(lmode .OR. (hfreq .AND. lfreq))) GO TO 9810
 
!     IS LMODES PRESENT WITH HFREQ AND/OR LFREQ
 
 IF (lmode .AND. (hfreq .OR. lfreq)) GO TO 9820
 
 3000 IF (app /= 4) GO TO 9999
 
!     TEST FOR CORRECT NODJE SETUP FOR AERO RF 10 AND 11
 
 ptot = p1 .AND. p2 .AND. p3
 IF (nodje .AND. ptot) GO TO 3500
 IF (nodje .AND. .NOT.ptot) GO TO 9830
 IF ((p1.OR.p2.OR.p3) .AND. .NOT.nodje) GO TO 9840
 
!     TEST FOR Q REQUIRED BY AERO RF 11
 
 3500 IF (rf == 10) GO TO 9999
 IF (queue) GO TO 9999
 GO TO 9870
 
!     TEST FOR CTYPE, NSEGS, OR KINDEX REQD BY DISP RF 14 AND 15.
 
 4000 ltest = ctype .AND. nsegs
 IF (.NOT.ltest) GO TO 9850
 4100 IF (rf == 14) GO TO 9999
 
 IF (kindx) GO TO 9999
 GO TO 9860
 
!     SET UP ERROR MESSAGE
 
 9800 ASSIGN 9900 TO ierr
 msgno = 340
 GO TO  9890
 9810 ASSIGN 9910 TO ierr
 msgno = 341
 GO TO  9890
 9820 ASSIGN 9920 TO ierr
 msgno = 342
 GO TO  9895
 9830 ASSIGN 9930 TO ierr
 msgno = 343
 GO TO  9890
 9840 ASSIGN 9940 TO ierr
 msgno = 344
 GO TO  9895
 9850 ASSIGN 9950 TO ierr
 msgno = 345
 GO TO  9890
 9860 ASSIGN 9960 TO ierr
 msgno = 346
 GO TO  9890
 9870 ASSIGN 9970 TO ierr
 msgno = 347
 
 9890 CALL page2 (3)
 WRITE  (nout,9891) ufm,msgno
 9891 FORMAT (a23,i4)
 abort = .true.
 GO TO 9898
 9895 CALL page2 (3)
 WRITE  (nout,9896) uwm,msgno
 9896 FORMAT (a25,i4)
 9898 GO TO ierr, (9900,9910,9920,9930,9940,9950,9960,9970)
 
 9900 WRITE  (nout,9905) appr(app),rf
 9905 FORMAT (' PARAM CARDS REQUIRED BY ',a4,' RIGID FORMAT',i3,  &
     ' NOT FOUND IN BULK DATA.')
 GO TO 9999
 
 9910 WRITE  (nout,9915) appr(app),rf
 9915 FORMAT (' LMODES OR HFREQ/LFREQ PARAM REQUIRED BY ',a4,  &
     ' RIGID FORMAT',i3,' NOT IN BULK DATA OR TURNED OFF.')
 GO TO 3000
 
 9920 WRITE  (nout,9925)
 9925 FORMAT (' LMODES PARAM FOUND IN BULK DATA WITH HFREQ OR LFREQ.',  &
     '  LMODES TAKES PRECEDENCE.')
 GO TO 3000
 
 9930 WRITE  (nout,9935) rf
 9935 FORMAT (' NODJE PARAM SPECIFIED FOR AERO RIGID FORMAT',i3,  &
     ' BUT P1, P2, OR P3 OMITTED.')
 GO TO 3500
 
 9940 WRITE  (nout,9945)
 9945 FORMAT (' P1, P2, OR P3 PARAM FOUND IN BULK DATA BUT NODJE ',  &
     'MISSING OR TURNED OFF.')
 GO TO 3500
 
 9950 WRITE  (nout,9955) rf
 9955 FORMAT (' CTYPE OR NSEGS PARAM REQUIRED BY DISPLACEMENT RIGID ',  &
     'FORMAT',i3,' MISSING OR INCORRECT.')
 GO TO 4100
 
 9960 WRITE  (nout,9965)
 9965 FORMAT (' KINDEX PARAM REQUIRED BY DISPLACEMENT RIGID FORMAT 15',  &
     ' MISSING OR TURNED OFF.')
 GO TO 9999
 
 9970 WRITE  (nout,9975)
 9975 FORMAT (' DYNAMIC PRESSURE (Q) PARAM REQUIRED BY AERO RIGID FORM',  &
     'AT 11 NOT IN BULK DATA.')
 
 9999 RETURN
END SUBROUTINE ifppar
