SUBROUTINE xyout (iopt,buf,rbuf)
     
!     THIS SUBROUTINE IS CALLED BY XYTRAN AND OUTPUTS TO PRINTER AND
!     PUNCH
 
 EXTERNAL        lshift,rshift
 LOGICAL :: PRINT,punch
 INTEGER :: buf(300),names(44),TYPE(6),plt(2),imtd(6), itype(4),rshift
 REAL :: rbuf(300)
 COMMON /machin/ mach,ihalf
 COMMON /BLANK / icom1,dum(4),icard
 COMMON /system/ sysbuf,l,d1(6),maxlns,d2(2),line,d3(78),lpch
 COMMON /output/ ihead(96)
 DATA    names / 4HDISP ,4HLACE ,4HMENT ,4H     ,  &
     4HVELO ,4HCITY ,4H     ,4H     , 4HACCE ,4HLERA ,4HTION ,4H     ,  &
     4HS p  ,4HC f  ,4H     ,4H     , 4HLOAD ,4H     ,4H     ,4H     ,  &
     4HELEM ,4HENT- ,4HSTRE ,4HSS   , 4HELEM ,4HENT- ,4HFORC ,4HE    ,  &
     4HS-di ,4HSPLA ,4HCEME ,4HNT   , 4HS-ve ,4HLOCI ,4HTY   ,4H     ,  &
     4HS-ac ,4HCELE ,4HRATI ,4HON   ,  &
     4HNONL ,4HINEA ,4HR-fo ,4HRCE                  /
 DATA    TYPE  / 4HWHOL ,4HE    ,4HUPPE ,4HR    ,4HLOWE ,4HR    /
 DATA    irand / 4HRAND /
 DATA    ivg   / 4HVG   /
 DATA    plt   / 4HNAST ,4HPLT  /
 DATA    imtd  / 4HFILM ,1H     ,4HTABL ,1HE    ,4HDRUM ,1H     /
 DATA    itype / 4HWITH ,4H     , 4HWITH ,4HOUT  /
 
 IF (icom1 == ivg) GO TO 86
 
!     BRANCH ON OPTION
 
 IF (iopt < 0) THEN
   GO TO    10
 ELSE
   GO TO    90
 END IF
 
!     PRINT XY-OUTPUT SUMMARY
 
 
!     FILL OUT HEADING
 
 10 DO  i = 1,96
   ihead(i) = buf(i+50)
 END DO
 CALL page1
 WRITE (l,150)
 IF (icom1 == irand) GO TO 30
 WRITE (l,170)buf(1)
 GO TO 40
 30 WRITE (l,160) rbuf(1)
 WRITE (l,161) rbuf(42)
 40 itempv = 4*buf(6) - 3
 
!     PRINT TYPE OF PLOT
 
 IF (buf(245)-2 < 0.0) THEN
   GO TO    41
 ELSE IF (buf(245)-2 == 0.0) THEN
   GO TO    42
 ELSE
   GO TO    43
 END IF
 41 WRITE (l,460)
 GO TO 45
 42 WRITE (l,470)
 GO TO 45
 43 WRITE (l,480)
 
!     PRINT DATA TYPE AND CURVE
 
 45 icomp = buf(5)
 IF (buf(6) /= 6 .AND. buf(6) /= 7) icomp = buf(5) - 2
 IF (buf(7) < 0.0) THEN
   GO TO    70
 ELSE IF (buf(7) == 0.0) THEN
   GO TO    60
 END IF
 50 WRITE (l,200) names(itempv),names(itempv+1),names(itempv+2),  &
     names(itempv+3),buf(4),icomp
 itemp = 3
 GO TO 72
 60 WRITE (l,180) names(itempv),names(itempv+1),names(itempv+2),  &
     names(itempv+3),buf(4),icomp
 itemp = 1
 GO TO 72
 70 WRITE (l,190) names(itempv),names(itempv+1),names(itempv+2),  &
     names(itempv+3),buf(4),icomp
 itemp  = 5
 72 icount = icard + 1
 WRITE (l,210)
 IF (buf(288) > 0) WRITE (l,230)
 IF (buf(290) > 0) WRITE (l,240) icount
 
!     PLOTTER INFORMATION
 
 IF (buf(289) <= 0) GO TO 84
 WRITE (l,220)
 j = rshift(buf(284),ihalf)
 model = buf(284) - lshift(j,ihalf) - 100
 m = 1
 IF (model < 0) m = 3
 
!   . NASPLOT...
 
 k = 2*IABS(model) - 1
 WRITE (l,380) plt(1),plt(2),imtd(k),imtd(k+1),itype(m),itype(m+1)
 IF (buf(283) <= 0) buf(283) = 1
 
!     WRITE CSCALE DATA OUT
 
 WRITE (l,490) rbuf(282)
 IF (IABS(model)-2 < 0) THEN
   GO TO    81
 ELSE
   GO TO    82
 END IF
 
!   . CAMERA, DENSITY...
 
 81 IF (buf(287) >= 3) WRITE (l,410)
 IF (buf(287) == 2) WRITE (l,430)
 IF (buf(287) <= 1) WRITE (l,420)
 WRITE (l,450) buf(283)
 GO TO 83
 
!   . PAPER SIZE
!     (THE LOGIC HERE IS SIMILAR TO THAT IN SUBROUTINE PLTSET)
 
 82 IF (IABS(model) == 2) GO TO 822
 
!   . DRUM PLOTTERS
 
 IF (rbuf(285) <= 0.0) rbuf(285) = 30.0
 IF (rbuf(286) <= 0.0) rbuf(286) = 30.0
 GO TO 824
 
!   . TABLE PLOTTERS
 
 822 IF (rbuf(285) <=  0.0) rbuf(285) = 11.0
 IF (rbuf(285) > 30.0) rbuf(285) = 30.0
 IF (rbuf(286) <=  0.0) rbuf(286) = 8.5
 824 IF (rbuf(286) > 30.0) rbuf(286) = 30.0
 WRITE (l,390) rbuf(285),rbuf(286)
 
!   . PEN SIZE
 
 WRITE (l,440) buf(283)
 83 WRITE (l,250) buf(3),TYPE(itemp),TYPE(itemp+1),buf(2)
 
!  .  PAPER PLOT
 
 84 IF (buf(289) > 0 .AND. buf(289) /= 2) GO TO 85
 WRITE (l,400) buf(281)
 
 85 CONTINUE
 WRITE (l,260) (buf(j),j=147,174),(buf(j),j=179,206), (buf(j),j=211,238)
 WRITE (l,270)
 WRITE (l,290) rbuf( 11),rbuf( 12)
 WRITE (l,300) rbuf(293),rbuf(294)
 WRITE (l,310) rbuf(295),rbuf(296)
 WRITE (l,280) rbuf(291),rbuf(292)
 WRITE (l,300) rbuf(297),rbuf(298)
 WRITE (l,310) rbuf(299),rbuf(300)
 WRITE (l,320)
 IF (buf(288) > 0) WRITE (l,330)
 86 itempv = 4*buf(6) - 3
 IF (buf(7) < 0.0) THEN
   GO TO    89
 ELSE IF (buf(7) == 0.0) THEN
   GO TO    88
 END IF
 87 itemp = 3
 GO TO 891
 88 itemp = 1
 GO TO 891
 89 itemp = 5
 891 iprint= 0
 id    = buf(4)
 icomp = buf(5)
 IF (buf(6) /= 6 .AND. buf(6) /= 7) icomp = buf(5) - 2
 PRINT = .false.
 punch = .false.
 IF (buf(290) > 0) punch = .true.
 IF (buf(288) > 0) PRINT = .true.
 IF (.NOT.PRINT) RETURN
 line = maxlns + 1
 RETURN
 
!     PRINT AND OR PUNCH OUTPUT
 
 90 iprint = iprint + 1
 IF (.NOT.punch) GO TO 100
 icard = icard + 1
 WRITE (lpch,370) iprint,rbuf(1),rbuf(2),icard
 100 IF (.NOT.PRINT) RETURN
 IF (line < maxlns) GO TO 110
 CALL page1
 WRITE (l,340) names(itempv),names(itempv+1),names(itempv+2),  &
     names(itempv+3),id,icomp,TYPE(itemp),TYPE(itemp+1)
 110 line = line + 1
 IF (.NOT.punch) GO TO 120
 WRITE (l,350) iprint,rbuf(1),rbuf(2),icard
 RETURN
 120 WRITE (l,350) iprint,rbuf(1),rbuf(2)
 RETURN
 
 150 FORMAT (///44X,33HX y - o u t p u t   s u m m a r y)
 160 FORMAT (//5X,24HROOT mean square value =,1P,e15.6)
 161 FORMAT (6X,38HFREQUENCY of zero crossings (n zero) =,1P,e15.6)
 170 FORMAT (//5X,7HSUBCASE,i10)
 180 FORMAT (6X,4A4,5HCURVE,i9,1H(,i2,1H))
 190 FORMAT (6X,4A4,5HCURVE,i9,4H(--,,i2,1H))
 200 FORMAT (6X,4A4,5HCURVE,i9,1H(,i2,4H,--))
 210 FORMAT (1H )
 220 FORMAT (6X,44HXY-pairs within frame limits will be plotted)
 230 FORMAT (6X,46HXY-pairs between xmin AND xmax will be printed)
 240 FORMAT (6X,64HXY-pairs between xmin AND xmax will be punched bigin&
     &ning on card,i8)
 250 FORMAT (//5X,13HTHIS is curve,i4,4H of ,a4,a2,5HFRAME,i5)
 260 FORMAT (//5X,14HCURVE  title =,28A4,/6X,14HX-axis title =,28A4,  &
     /6X,14HY-axis title =,28A4)
 270 FORMAT (/////5X,62HTHE following information is for the above defi&
     &ned curve only. )
 280 FORMAT (//5X,36HWITHIN the x-limits of all DATA (x =,1P,e14.6,  &
     8H TO  x =,1P,e14.6,1H))
 290 FORMAT (///6X,36HWITHIN the frame x-limits       (x =,1P,e14.6,  &
     8H TO  x =,1P,e14.6,1H))
 300 FORMAT (//30X,22HTHE smallest y-value =,1P,e14.6,7H at x =,e15.6)
 310 FORMAT (//30X,22HTHE largest  y-value =,1P,e14.6,7H at x =,e15.6, //)
 320 FORMAT (//45X,27HE n d   o f   s u m m a r y)
 330 FORMAT (//25X,69HP r i n t e d   d a t a   f o r   t h i s   c u r  &
     &v e   f o l l o w s)
 340 FORMAT (//5X,4A4,12HCURVE   id =,i9,5X,11HCOMPONENT =,i3,5X,a4,a2,  &
     5HFRAME,///27X,12HPRINT NUMBER,10X,7HX-value,14X,  &
     7HY-value,14X,11HCARD NUMBER  )
 350 FORMAT (28X,i7,1P,e25.6,e21.6,10X,i8)
 370 FORMAT (i10,10X,1P,2E20.6,12X,i8)
 380 FORMAT (6X,21HPLOTTER specified is ,3A4,a1,9H plotter ,2A4,  &
     18HTYPING capability.)
 390 FORMAT (6X,11HPAPER size ,f5.2,3H x ,f5.2,18H inches specified.)
 400 FORMAT (6X,38HTHIS curve will be paper-plotted frame,i5)
 410 FORMAT (6X,36HCAMERA 3 used. (paper AND 35MM film))
 420 FORMAT (6X,26HCAMERA 2 used. (35MM film))
 430 FORMAT (6X,22HCAMERA 1 used. (paper))
 440 FORMAT (6X,9HPENSIZE =,i3)
 450 FORMAT (6X,9HDENSITY =,i3)
 460 FORMAT (6X,8HRESPONSE)
 470 FORMAT (6X,38HPOWER-spectral-density-FUNCTION (psdf))
 480 FORMAT (6X,15HAUTOCORRELATION)
 490 FORMAT (6X,9HCSCALE = ,f5.2)
END
