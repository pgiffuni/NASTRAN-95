SUBROUTINE ttlpge (topt)
     
 
 
 INTEGER, INTENT(OUT)                     :: topt
 INTEGER :: idate(3),card(20), fchar
 CHARACTER (LEN=7) :: machos
 CHARACTER (LEN=11) :: mchnam
 CHARACTER (LEN=15) :: vn
 CHARACTER (LEN=28) :: mchttl
 COMMON /chmach/ mchnam, machos
 COMMON /machin/ machx
 COMMON /system/ ksystm(100)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (ksystm( 2),nout), (ksystm(42),idate(1)),  &
     (ksystm( 9),nlpp), (ksystm(11),   ipage), (ksystm(91),lpch)
 
!     ASSEMBLE MCHTTL AND VN LINE
 
 mchttl = ' '
 vn = ' '
 ncmnam = INDEX(mchnam,' ') - 1
 IF (ncmnam <= -1) ncmnam = 11
 ncmos  = INDEX(machos,' ') - 1
 IF (ncmos <= -1) ncmos = 7
 fchar = (11 - ncmnam)/2 + 1
 mchttl(fchar:fchar+ncmnam+16) = mchnam(1:ncmnam) // ' COMPUTER SYSTEMS'
 fchar = (7 - ncmos)/2 + 1
 vn(fchar:fchar+ncmos+7) = machos(1:ncmos) // ' VERSION'
 
!     SET TOPT DEFAULT TO +2 FOR THE MAIN FRAMES, OR TO -1 FOR UNIX
!     BASE WORKSTATION
 
 IF (topt /= -9) GO TO 1
 topt = +2
 IF (machx >= 6 .AND. machx <= 20) topt = -1
 1  CONTINUE
 
!     BRANCH ON OPTION
 
!     TOPT = 1, PRINT ONE NASTRAN LOGO TITLE PAGE
!          = 2, PRINT TWO NASTRAN LOGO TITLE PAGES
!          = 3, PRINT DUMMY MESSAGE AND ONE SHORT TITLE PAGE
!          = 4, READ AND PRINT ONE LINE USER INPUT CARD AND PRINT ONE
!               NASTRAN SHORT TITLE PAGE
!          = 0, OR .GE.5, NO TITLE PAGE PRINTED
!          = NEGATIVE INTEGER, PRINT ONE NASTRAN SHORT TITLE PAGE
 
 5 IF (topt /= 2  .AND. topt /= 1) GO TO 110
 
!     TOPT = 1, OR 2
 
 DO  i = 1,topt
   IF (ipage <= 0 .OR. i == 2) WRITE (nout,10)
   IF (nlpp  > 48) WRITE (nout,20)
   WRITE (nout,30) mchttl,vn
   WRITE (nout,50)
   WRITE (nout,60) idate(2),idate(3)
   WRITE (nout,70)
   WRITE (nout,75)
   WRITE (nout,80)
   WRITE (nout,85)
   WRITE (nout,90)
   WRITE (nout,95)
   10   FORMAT (1H1)
   20   FORMAT (///)
   30   FORMAT (34X,17(1HM), /28X,29(1HM),  &
       /25X,35(1HM), /22X,20(1HM),1X,20(1HM),22X,1H/,6X,a28,  &
       /20X,45(1HM),18X,2H//,9X,a20)
   40   FORMAT (1H+,93X,a4,10H version -,i5,1HK)
   50   FORMAT (18X,16(1HM),2X,31(1HM),14X,3H///, /16X,53(1HM),10X,4(1H/),  &
       /14X,13(1HM),9X,35(1HM),6X,5(1H/))
   60   FORMAT (13X,12(1HM),2X, 9(1HM),2X,34(1HM),3X, 6(1H/),9X,  &
       3X,18HSYSTEM release  - , a3,a2,4H ed.)
   70   FORMAT (12X,12(1HM),1X,13(1HM),3X,15(1HM),2X,15(1HM),6(1H/),  &
       /11X,12(1HM),1X,17(1HM),2X,28(1HM),6(1H/),  &
       /10X,13(1HM),1X,19(1HM),2X,24(1HM),6(1H/),  &
       /9X,5(1HM),2X,7(1HM),1X,13(1HM),1X,7(1HM),2X,19(1HM),8(1H/) ,      2HMM,  &
       /9X,14(1HM),1X,23(1HM),2X,14(1HM),8(1H/),1H-,4(1HM),  &
       43X,1H*,1X,1H*,1X,1H*,  &
       /8X,16(1HM),1X,24(1HM),1X,9(1HM),9(1H/),2H--,7(1HM), 41X,1H*,5X,1H*)
   75   FORMAT (8X,16(1HM),1X,25(1HM),2X,4(1HM),10(1H/),2H--,9(1HM),  &
       41X,1H*,2X,1HR,2X,1H*,  &
       /8X,16(1HM),1X,27(1HM),1X,1HM,8(1H/),4HMM--,11(1HM), 41X,1H*,5X,1H*,  &
       /7X,8(1HM),4X,6(1HM),4X,5(1HM),5X,10(1HM),8X,4H//mm,11X,  &
       2HMM,3X,6(1HM),7X,5(1HM),8X,4(1HM),6X,4(1HM),2X,1H*,1X,1H*, 1X,1H*,  &
       /7X,9(1HM),4X,6(1HM),2X,7(1HM),4X,6(1HM),14H///   /// m  m,  &
       25HM- mmm   mmm mmm  m   mmm,7X,4(1HM),9X,4(1HM),6X,2HMM)
   80   FORMAT (7X,9(1HM),5X,5(1HM),2X,6(1HM),3H  m,3X,8(1H/),3X,4(1HM),  &
       5H mm--,5(1HM),3X,7(1HM),21H  m    mmm     mm mmm,8X, 5(1HM),5X,2HMM,  &
       /7X,9(1HM),2X,1HM,4X,5HMMM  ,4(1HM),6H// ///,3X,5(1H/),  &
       13HMMM   mmmm-- ,6(1HM),3X,7(1HM),41H  m    mmm     m   mmm  &
       mm mmmm   mm, /7X,9(1HM),2X,2HMM,4X,2HMM,3X,4(1H/),2X,3H///,4X,8(1HM),  &
       4X,4H--m ,7(1HM),3X,7(1HM),2X,1HM,3X,3HMMM,5X,2HMM,3X,  &
       4(1HM),6X,2HMM,2X,4(1HM),2X,2HMM)
   85   FORMAT (7X,9(1HM),2X,4(1HM),6X,11H/ /// ///mm,4X,8(1HM),4H---m,  &
       4X,6(1HM),3X,7(1HM),2X,6(1HM),6X,1HM,5X,4(1HM),5X,2HMM, 4X,6(1HM),  &
       /7X,9(1HM),2X,5(1H/),5X,4H// m,11X,6(1HM),3H---,4(1HM),4X,  &
       5(1HM),3X,7(1HM),7H  m mmm,6X,11(1HM),5X,2HMM,5X,5(1HM),  &
       /7X,2HMM,7(1H/),2X,6(1HM),4X,3HMMM,2X,7(1HM),4X,7HMMM----,  &
       4HMMMM,4X,6HM mmmm,3X,7(1HM),8H  m  mmm,4X,2HMM,7X,  &
       4(1HM),4X,2HMM,6X,4(1HM),  &
       /5X,4(1H/),6(1HM),4X,7(1HM),2X,2HMM,4X,5(1HM),5X,5H----m,  &
       9X,6HMM mmm,5X,5(1HM),3X,2HMM,3X,2HMM,2X,4(1HM),5X,  &
       6(1HM),2X,4(1HM),7X,2HMM)
   90   FORMAT (3X,2H//,3X,26(1HM),1X,6(1HM),4(1H-),16(1HM),1X,15(1HM),  &
       6X,3HMMM, /8X, 27(1HM),7H mm----,19(1HM),1X,15(1HM),  &
       /8X, 27(1HM),3H---,23(1HM),1X,15(1HM),  &
       /9X, 24(1HM),7H---mm  ,22(1HM),1X,13(1HM),  &
       /9X, 22(1HM),2H--,6(1HM),4X,19(1HM),1X,5(1HM),2X,6(1HM),  &
       /10X,19(1HM),3H---,7(1HM),4X,19(1HM),1X,12(1HM),  &
       /11X, 9(1HM),1X,6(1HM),2H--,  33(1HM), 1X, 11(1HM),  &
       /12X,13(1HM),3H---,33(1HM),1X,11(1HM))
   95   FORMAT (13X,11(1HM),2H--, 22(1HM),2X, 9(1HM), 2X, 11(1HM),  &
       /14X, 8(1HM),2H--,26(1HM), 9X,12(1HM),  &
       /16X, 5(1HM),2H--,46(1HM),24X,14HDISTRIBUTED by,  &
       /18X, 4HMM--,13(1HM),2X,30(1HM), /19X, 1H-,   45(1HM),5X,  &
       51HCOMPUTER software management AND information center, 9H (cosmic),  &
       /18X,1H-,3X,41(1HM),26X,22HUNIVERSITY of  georgia,  &
       /17X,1H-,7X,35(1HM),29X,22HATHENS, georgia  30602, /28X,29(1HM),  &
       /1X,14X, 19X,17(1HM),28X,40HPHONE: (706)542-3265   fax: (706)542-480  &
       ,1H7)
 END DO
 GO TO 240
 
 110  IF (topt   < 0.0) THEN
   GO TO   160
 ELSE IF (topt   == 0.0) THEN
   GO TO   240
 END IF
 120  IF (topt-4 < 0.0) THEN
   GO TO   130
 ELSE IF (topt-4 == 0.0) THEN
   GO TO   210
 ELSE
   GO TO   240
 END IF
 
!     TOPT = 3
 
 130  WRITE  (nout, 10)
 WRITE  (nout,140)
 140  FORMAT (' THIS COMMENT CAN BE USED TO IDENTIFY LOCAL FIXES - ',  &
     'TO CHANGE, UPDATE DECK TTLPGE.')
 GO TO 160
 
!     TOPT = NEGATIVE (AND 3, AND 4)
 
 160  IF (ipage <= 0) CALL page1
 WRITE  (nout,170) mchttl
 WRITE  (nout,180) vn,idate(2),idate(3)
 WRITE  (nout,190)
 170  FORMAT (//////34X,4H****, /32X,1H*,6X,1H*, /31X,1H*,8X,1H*,  &
     /31X,16H*  n a s t r a n, /31X,1H*,8X,1H*, /32X,1H*,6X,1H*, /34X,4H****,  &
     ///25X,a28)
 180  FORMAT(27X,a20,//26X,17HSYSTEM release - ,a3,a2, 4H ed.)
 190  FORMAT (/32X,'DISTRIBUTED BY', //9X,'COMPUTER SOFTWARE MANAGE',  &
     'MENT AND INFORMATION CENTER (COSMIC)', /17X,'UNIVERSITY ',  &
     'OF GEORGIA, ATHENS, GEORGIA 30602', /17X,  &
     'PHONE: (706)542-3265', 6X, 'FAX: (706)542-4807')
 GO TO 240
 
!     TOPT = 4
 
 210  WRITE  (nout,10)
 CALL xread (*240,card)
 WRITE  (nout,220) card
 220  FORMAT (1X,20A4)
 GO TO 160
 
!     CALL NSINFO TO PRINTOUT INSTALLATION-CENTER-TO-USER MESSAGES,
!     FROM THE THIRD SECTION OF THE NASINFO FILE
 
 240  CALL nsinfo (3)
 
 RETURN
END SUBROUTINE ttlpge
