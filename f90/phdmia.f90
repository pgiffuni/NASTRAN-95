SUBROUTINE phdmia
     
!     PUNCH SINGLE- OR DOUBLE-FIELD DMI CARDS FOR REAL, SINGLE-
!     PRECISION MATRICES.
 
!  $MIXED_FORMAT
 
 LOGICAL :: first
 INTEGER :: FMT(8),iqx(8),kfmt(22),ret,NAME(2),ret1,lfmt(23),  &
     h1,h2,h3,c1,c2,c3,h1a,c1a,erno
 REAL :: qx(8)
 COMMON /phdmix/  NAME,nam,ifo,itin,itout,ir,ic,noutpt,kpp,nlpp,  &
     erno,icol,iro,x,icard1
 COMMON /system/  dum90(90),np
 COMMON /mahcin/  mach
 EQUIVALENCE      (qx(1),iqx(1)), (lfmt(2),kfmt(1))
 DATA    dmi   ,  iz, p  , dmis  , s       /  &
     3HDMI ,  0 , 1H+, 4HDMI*, 1H*     /
 DATA    kfmt  /  3*4H$$$$,16*4H**** , 4HA1,a,4H2,i5,4H)   /
 DATA    kfmti ,  kfmtr1,kfmtr2/4HI8 ,,4HF8.1,4H,          /
 DATA    kdmti ,  kdmtr1,kdmtr2/4HI16,,4H1PE1,4H6.8,       /
 DATA    kfmtb ,  kfmt8 ,kdmtr0/4H    ,4H4X, ,4H  e1       /
 DATA    kdmtb ,  kdmt8 /4H    ,4H8X,      /
 DATA    h1    ,  h2    ,h3    ,c1    ,c2    ,c3    ,h1a   ,c1a   /  &
     4H(a4,,  4H4X, ,4H2A4,,4H(a1,,4HA2, ,4HI5, ,4H a4,,4H a1,/
 DATA   lfmt(1)/  4H(1X,/
 
!     IBM/AIX (MACH=8) DOES NOT LIKE THE NON-ANSI STANDARD FORMAT
!     1PE16.8 (THE STANDARD IS 1P,E16.8).
 
 IF (mach == 8) kdmtr1 = kdmtr0
 
!     CALLED INITIALLY FOR EACH MATRIX.
 
!     SET PUNCH UNIT TO 7 FOR IBM AND CDC AND TO 1 FOR UNIVAC
 
 erno  = 0
 nout  = noutpt
 kp    = kpp
 nkp   = 8/kp
 icard =-1
 icard1= 0
 SELECT CASE ( kp )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 20
 END SELECT
 10 dmips = dmi
 ps = p
 GO TO 30
 20 dmips = dmis
 ps = s
 DO  i = 12,19
   kfmt(i) = kdmtb
 END DO
 30 WRITE  (np,1) dmi,NAME,iz,ifo,itin,itout,ir,ic,p,nam,icard1
 1 FORMAT (a3,5X,2A4,4I8,8X,2I8,a1,a2,i5)
 IF (nout <= 0) GO TO 40
 WRITE  (nout,2) dmi,NAME,iz,ifo,itin,itout,ir,ic,p,nam,icard1
 2 FORMAT (1H1,/1X,a3,5X,2A4,4I8,8X,2I8,a1,a2,i5)
 l = 1
 40 RETURN
 
 
 ENTRY phdmib
!     ============
 
!     CALLED FOR FIRST NON-ZERO ELEMENT OF EACH COLUMN.
 
 iq   = 0
 irow = iro
 first=.true.
 iq   = iq + 1
 FMT(iq) = 0
 iqx(iq) = 0
 iq   = iq + 1
 FMT(iq) = 1
 iqx(iq) = icol
 iq   = iq + 1
 FMT(iq) = 1
 iqx(iq) = irow
 iq   = iq + 1
 FMT(iq) = 2
 qx(iq)  = x
 RETURN
 
 
 ENTRY phdmic
!     ============
 
!     CALLED FOR EACH NON-ZERO ELEMENT OF COLUMN EXCEPT FIRST ONE.
 
!     LOOK FOR FULL CARD
 
 IF (iq < nkp) GO TO 100
 ASSIGN 100 TO ret
 GO TO 700
 
!     DETERMINE IF NEW ENTRY IS CONSECUTIVE OR NON-CONSECUTIVE.
 
 100 IF (iro /= irow+1) GO TO 200
 irow = iro
 iq   = iq + 1
 FMT(iq) = 2
 qx(iq)  = x
 RETURN
 
 200 iq   = iq + 1
 irow = iro
 FMT(iq) = 1
 iqx(iq) = iro
 IF (iq < nkp) GO TO 300
 ASSIGN 300 TO ret
 GO TO 700
 300 iq = iq + 1
 FMT(iq) = 2
 qx(iq)  = x
 RETURN
 
 
 ENTRY phdmid
!     ============
 
!     ENTRY POINT FOR COLUMN TERMINATION CALL
 
 IF (iq <= 0) RETURN
 ASSIGN 500 TO ret
 GO TO 700
 500 RETURN
 
!     PUNCH CARD
 
 700 n = iq
 ASSIGN 800 TO ret1
 GO TO 1000
 800 IF (first) GO TO 900
 WRITE (np,kfmt,ERR=810) ps,nam,icard,(qx(l),l=1,iq),ps,nam,icard1
 810 lfmt(2) = c1a
 IF (nout <= 0) GO TO 850
 IF (l < nlpp) GO TO 830
 WRITE  (nout,3)
 3 FORMAT (1H1)
 l = 0
 830 WRITE (nout,lfmt,ERR=840) ps,nam,icard,(qx(l),l=1,iq), ps,nam,icard1
 840 l  = l + 1
 850 iq = 0
 GO TO 950
 900 WRITE (np,kfmt,ERR=910) dmips,NAME,(qx(l),l=2,iq),ps,nam,icard1
 910 lfmt(2) = h1a
 IF (nout <= 0) GO TO 940
 IF (l < nlpp) GO TO 920
 WRITE (nout,3)
 l  = 0
 920 WRITE (nout,lfmt,ERR=930) dmips,NAME,(qx(l),l=2,iq),ps,nam,icard1
 930 l  = l + 1
 940 first = .false.
 iq = 0
 950 GO TO ret, (100,300,500)
 
!     BUILD FORMAT FOR CARD IMAGE.
 
 1000 icard  = icard + 1
 icard1 = icard + 1
 IF (icard1 > 99999) GO TO 9901
 SELECT CASE ( kp )
   CASE (    1)
     GO TO 1001
   CASE (    2)
     GO TO 1101
 END SELECT
 1001 IF (first) GO TO 1005
 i1 = 1
 kfmt(1) = c1
 kfmt(2) = c2
 kfmt(3) = c3
 GO TO 1009
 1005 i1 = 2
 kfmt(1) = h1
 kfmt(2) = h2
 kfmt(3) = h3
 kfmt(4) = kfmtb
 kfmt(5) = kfmtb
 1009 DO  i = i1,n
   k  = FMT(i)
   IF (k == 2) GO TO 1020
   1010 kfmt(2*i+2) = kfmti
   kfmt(2*i+3) = kfmtb
   CYCLE
   1020 kfmt(2*i+2) = kfmtr1
   kfmt(2*i+3) = kfmtr2
 END DO
 IF (n >= nkp) GO TO 1999
 n1 = n + 1
 DO  i = n1,nkp
   kfmt(2*i+2) = kfmt8
   kfmt(2*i+3) = kfmt8
 END DO
 GO TO 1999
 1101 IF (first) GO TO 1105
 i1 = 1
 kfmt(1) = c1
 kfmt(2) = c2
 kfmt(3) = c3
 GO TO 1109
 1105 i1 = 2
 kfmt(1) = h1
 kfmt(2) = h2
 kfmt(3) = h3
 kfmt(4) = kdmt8
 kfmt(5) = kdmtb
 1109 DO  i = i1,n
   k  = FMT(i)
   IF (k == 2) GO TO 1120
   1110 kfmt(2*i+2) = kdmti
   kfmt(2*i+3) = kdmtb
   CYCLE
   1120 kfmt(2*i+2) = kdmtr1
   kfmt(2*i+3) = kdmtr2
 END DO
 IF (n >= nkp) GO TO 1999
 n1  = n + 1
 DO  i = n1,nkp
   kfmt(2*i+2) = kdmt8
   kfmt(2*i+3) = kdmt8
 END DO
 GO TO 1999
 1999 GO TO ret1, (800)
 
 
!     ERROR MESSAGES
 
 9901 erno = 1
 GO TO 9999
 
 
 9999 RETURN
 
END SUBROUTINE phdmia
