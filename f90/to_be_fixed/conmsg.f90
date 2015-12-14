SUBROUTINE conmsg (mesage, nwords, idummy)
     
 
 INTEGER, INTENT(IN OUT)                  :: mesage(1)
 INTEGER, INTENT(IN OUT)                  :: nwords
 INTEGER, INTENT(IN)                      :: idummy
 INTEGER :: fchar
 REAL :: inctim, modtim
 CHARACTER (LEN=7) :: machos
 CHARACTER (LEN=11) :: mchnam
 CHARACTER (LEN=41) :: ahead
 CHARACTER (LEN=8) :: ctime
 
 CHARACTER (LEN=1) :: (LEN=12) real_clock(3), time
 INTEGER :: values(8)
 CHARACTER (LEN=1) :: 4) cvalues(8)
 EQUIVALENCE(values,cvalues)
 
 DIMENSION icrdat(3)
 DIMENSION idate(3), itime(3)
 
 COMMON /chmach/ mchnam, machos
 COMMON /logout/ lout
 COMMON /system/ isystm(175)
 
 EQUIVALENCE (isystm( 15), idate(1)), (isystm( 18), cpustr),  &
     (isystm( 42), icrdat), (isystm( 75), cputim),  &
     (isystm(151), nllog ), (isystm(152), loglin),  &
     (isystm(159), logpag), (isystm(160), oldcpu)
 
 DATA idsms,iwrtt,iaudt,impya /4HDSMS, 4HWRTT, 4HAUDT, 4HMPYA/
 DATA modtim /0.0/
 DATA idash /4H----/
 
!   ASSEMBLE PAGE HEADING
 
 ahead = ' '
 ncmnam = INDEX(mchnam,' ') - 1
 IF (ncmnam <= -1) ncmnam = 11
 ncmos  = INDEX(machos,' ') - 1
 IF (ncmos <= -1) ncmos = 7
 fchar = (18 - ncmnam - ncmos) + 1
 ahead(fchar:fchar+6)='LOG OF '
 fchar = fchar + 7
 WRITE (ahead(fchar:fchar+1),15) icrdat(3)
 15 FORMAT (a2)
 fchar = fchar + 3
 ahead(fchar:41) = mchnam(1:ncmnam) // ' ' // machos(1:ncmos) //  &
     ' NASTRAN JOB'
 
 imodtm = 0
 IF (idummy == 111111 .OR. idummy == 222222) imodtm = idummy/111111
 IF (loglin < nllog.AND.loglin > 0) GO TO 300
 IF (loglin == 0) WRITE (lout, 2000) idate, ahead
 IF (loglin == 0) WRITE (lout, 2055)
 IF (mesage(1) == idsms.AND.nwords == 1) RETURN
 IF (mesage(1) == iwrtt.AND.nwords == 1) RETURN
 IF (mesage(1) == iaudt.AND.nwords == 1) RETURN
 IF (mesage(1) == impya.AND.nwords == 1) RETURN
 300 CALL nastim (itime(1), itime(2), itime(3), cputmm)
 WRITE (ctime,2056) itime
 CALL date_and_time(real_clock(1), time, real_clock(3), values)
 
 2056 FORMAT (2( i2,':'),i2)
 IF (ctime(4:4) == ' ') ctime(4:4) = '0'
 IF (ctime(7:7) == ' ') ctime(7:7) = '0'
 cputmm = cputmm + oldcpu - cpustr
 inctim = cputmm - cputim
 IF (cputim == 0.0) inctim = 0.0
 IF (imodtm == 1) modtim = 0.0
 IF (imodtm == 2) modtim = cputmm - modtim
 mwords = MIN0 (nwords, 15)
 IF (imodtm /= 2) WRITE (lout,2100) time(1:2),time(3:4),time(5:10),  &
     cputmm, inctim,(mesage(i), i = 1, mwords)
 IF (imodtm == 2) WRITE (lout,2110) time(1:2),time(3:4),time(5:10),  &
     cputmm, inctim, modtim, (mesage(i), i = 1, mwords)
 CALL flush(4)
 
 loglin = loglin + 1
 cputim = cputmm
 IF (imodtm == 1) modtim = cputmm
 RETURN
 
 2000 FORMAT (1H1,  77(1H*)/ 1X , 1H*,  75X, 1H*/  &
     1X , 1H*, 7X, 'DATE ', 2(i2, '/'), i2, 7X, a41, 7X, 1H*/  &
     1X , 1H*,  75X, 1H*/ 1X ,  77(1H*)/)
 2055 FORMAT (1X, 2X, 'WALL', 15X, 'TOTAL', 7X,  &
     'INCREMENTAL', 6X, 'MODULE', 14X,  &
     'MODULE/'/ 1X, 2X, 'CLOCK', 15X,  &
     'CPU', 12X, 'CPU', 12X,  &
     'CPU', 13X, 'SUBROUTINE'/  &
     1X, 2X, 'TIME', 14X, 'SECONDS', 8X,  &
     'SECONDS', 8X, 'SECONDS', 13X,  &
     'STATUS'// 1X,  78(1H-)/)
 2100 FORMAT (1X, a2,':',a2,':',a6, 4X, f10.3,  5X, f10.3, 15X,  &
     5X, a4, 2X, 2A4, 2X, 12A4)
 2110 FORMAT (1X, a2,':',a2,':',a6, 4X, f10.3,  5X, f10.3,  5X, f10.3,  &
     5X, a4, 2X, 2A4, 2X, 12A4)
END SUBROUTINE conmsg
