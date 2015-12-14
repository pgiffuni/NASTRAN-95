SUBROUTINE pltset
     
!     COMMENTS FROM G.C. -
!     THE DRIVER FOR DMAP MODULE PLTSET IS DPLTST
!     THIS ROUTINE HAS NOTHING TO DO WITH DPLTST.  IT IS CALLED ONLY
!     BY PARAM (IN MODULE PLOT), XYPLOT, AND SEEMAT
 
 
 LOGICAL :: tapbit
 INTEGER :: chrwrd,pbfsiz,pbufsz,pdata,pltdat,pltype,  &
     ploter,plt1,plt2,pltnum,offscl
 REAL :: xymax(2),cntchr(2),xysize(2)
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim
 COMMON /BLANK / skp4(4),pltnum
 COMMON /xxparm/ pbufsz,skparm(6),papsiz(2),skp235(226),offscl
 COMMON /pltdat/ model,ploter,reg(2,2),axymax(2),xyedge(11),chrscl,  &
     pdata(20),pltdat(20,1)
 COMMON /machin/ mach
 COMMON /system/ ksystm(65)
 EQUIVALENCE     (pdata(1),xymax(1)) ,  (pdata(3),cntsin)  ,  &
     (pdata(4),cntchr(1)),  (pdata(10),pltype) ,  &
     (pdata(12),pbfsiz)  ,  (nout  ,ksystm( 2)),  &
     (chrwrd,ksystm(41)) ,  (itrack,ksystm(59))
 DATA    xysize/ 11.0, 8.5 /,  plt1,plt2 / 4HPLT1, 4HPLT2  /
 
!     INITIALIZE -PDATA-
 
 DO  i = 1,20
   pdata(i) = pltdat(i,ploter)
 END DO
 
!     PLT2 FILE WAS HARD CODED INTO THE 11TH WORD OF PLTDAT(11,PLOTER)
!     BY PLOTBD. IF USER REQUESTS PLT1 FILE, WE MUST MAKE A SWITCH HERE
 
 IF (.NOT.tapbit(plt2) .AND. tapbit(plt1)) pdata(11) = plt1
 IF (pltnum == 0 .AND. offscl == 0) WRITE (nout,110) uim,pdata(11)
 110 FORMAT (a29,', PLOT FILE GOES TO ',a4)
 
 IF (offscl == 0) offscl = 1
 
!     SCALE THE CHARACTURE SIZE BEFORE SETTING BORDERS
 
 cntchr(1) = chrscl*cntchr(1)
 cntchr(2) = chrscl*cntchr(2)
 pbufsz    = pbfsiz/chrwrd
 
!     FOR UNIVAC 9 TRACK CALCOMP PLOT TAPES QUARTER WORD MODE WILL
!     BE USED LIMITING THE NUMBER OF CHARACTERS PER WORD TO 4
!     ITRACK = 2 FOR 9 TRACK TAPES - OTHERWISE 1 FOR 7 TRACK TAPES
!     THE DEFAULT IS FOR 7 TRACK TAPES
 
!     IF (MACH.EQ.3 .AND. ITRACK.EQ.2) PBUFSZ = PBFSIZ/4
 
!     SINCE GENERAL PLOTTER IS THE ONLY ONE SUPPORTED BY NASTRAN, THE
!     PLOT BUFFER FOR UNIVAC MUST BE 500 WORDS FOR BOTH FORTRAN V AND
!     ASCII FORTRAN. (SEE PROG. MANUAL PAGE 6.10-15)
 
 IF (mach == 3) pbufsz = pbfsiz/6
 
 pltype = model
 
!     INITIALIZE PAPER SIZE AND BORDERS
 
 DO  i = 1,2
   IF (IABS(pltype)-2 < 0.0) THEN
     GO TO   121
   ELSE IF (IABS(pltype)-2 == 0.0) THEN
     GO TO   125
   ELSE
     GO TO   124
   END IF
   121 IF (pltype > 0.0) THEN
     GO TO   123
   END IF
   
!     CRT PLOTTERS
   
   122 axymax(i) = xymax(i) - cntchr(i)
   xyedge(i) = cntchr(i)*.5
   GO TO 129
   123 axymax(i) = xymax(i)
   xyedge(i) = 0.
   GO TO 129
   
!     DRUM PLOTTERS
   
   124 IF (papsiz(i) <= 0.0) papsiz(i) = xymax(i)/cntsin
   SELECT CASE ( i )
     CASE (    1)
       GO TO 127
     CASE (    2)
       GO TO 126
   END SELECT
   
!     TABLE PLOTTERS
   
   125 IF (papsiz(i) <= 0.0) papsiz(i) = xysize(i)
   
   126 IF (cntsin*papsiz(i) > xymax(i)) papsiz(i) = xymax(i)/cntsin
   127 axymax(i) = cntsin*papsiz(i) - cntsin
   xyedge(i) = cntsin*.5
   129 reg(i,1)  = 0.
   reg(i,2)  = axymax(i)
 END DO
 
 RETURN
END SUBROUTINE pltset
