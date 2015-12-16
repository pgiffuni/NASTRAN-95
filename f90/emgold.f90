SUBROUTINE emgold
     
!     THIS IS A DRIVING ROUTINE OF THE -EMG- MODULE WHICH ALLOWS PIVOT-
!     POINT-LOGIC ELEMENT SUBROUTINES TO BE USED IN CONJUNCTION WITH THE
!     NON-PIVOT-POINT PROCESS.
 
 LOGICAL :: error,last,heat,kheat,lheat,hydro
 INTEGER :: outpt,sil,posvec,eltype,elid,elem,dict,estwds,  &
     estid,estbuf,filtyp,precis,mdict(15),kdict(15), bdict(15),flags,qp
 DOUBLE PRECISION :: dummy
 CHARACTER (LEN=31) :: sim
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm,swm,sim
 COMMON /machin/  mach,dum2(2),lqro
 COMMON /system/  ksystm(65)
 COMMON /gpta1 /  nelems,nlast,incr,elem(1)
 COMMON /emgdic/  eltype,ldict,nlocs,elid,estid
 COMMON /emgest/  estbuf(100)
 COMMON /sma1et/  kecpt(100)
 COMMON /sma2et/  mecpt(100)
 COMMON /sma1io/  smaio(36)
 COMMON /sma1bk/  s1dum(10)
 COMMON /sma2bk/  s2dum(10)
 COMMON /sma1dp/  kwork(700)
 COMMON /sma2dp/  mwork(700)
 COMMON /sma1cl/  iopt4,k4ggsw,knpvt,skip19(19),knogo,ksafe(200)
 COMMON /sma2cl/  ioptb,bggind,mnpvt,skip17(17),mnogo,skip2(2), msafe(200)
 COMMON /emg1bx/  nsils,posvec(10),ibloc,nbloc,irows,dict(15),  &
     filtyp,sil(10),last
 COMMON /emgprm/  skipxx(15),flags(3),precis,error,heat,icmbar,  &
     lcstm,lmat,lhmat,kflags(3),l38
 COMMON /hydroe/  hydro
 COMMON /sma1ht/  kheat
 COMMON /sma2ht/  lheat
 COMMON /iemgod/  dummy,ltypes
 EQUIVALENCE      (ksystm(2),outpt), (ksystm(40),nbpw),  &
     (smaio(11),ifkgg), (smaio(13),if4gg)
 
 ktemp = knogo
 mtemp = mnogo
 qp    = MOD(lqro/100,10)
 jltype= 2*(eltype-1) + precis
 kheat = heat
 lheat = heat
 izero = incr*(eltype - 1)
 IF (eltype == ltypes) GO TO 20
 CALL page2 (3)
 INDEX = izero
 IF (.NOT.heat) GO TO 3
 IF (eltype == 62 .OR. eltype == 63) INDEX = 15*incr
 3 IF (l38 == 1) WRITE (outpt,4) sim,elem(INDEX+1),elem(INDEX+2)
 4 FORMAT (a31,' 3107',/5X,'EMGOLD CALLED BY EMGPRO TO PROCESS ',2A4,  &
     ' ELEMENTS.')
 ltypes = eltype
 GO TO 20
 5 WRITE  (outpt,10) swm,elid,elem(izero+1),elem(izero+2)
 10 FORMAT (a27,' 3121, EMGOLD HAS RECEIVED A CALL FOR ELEMENT ID',i9,  &
     ' (ELEMENT TYPE ',2A4,2H)., /5X,'ELEMENT IGNORED AS THIS ',  &
     'ELEMENT TYPE IS NOT HANDLED BY EMGOLD.')
 GO TO 1220
 
 20 nsils = elem(izero+10)
 isil  = elem(izero+13)
 IF (elem(izero+9) /= 0) isil = isil - 1
 estwds = elem(izero+12)
 i1 = isil
 i2 = isil + nsils - 1
 l  = nsils
 
!     MOVE SILS TO SEPARATE ARRAY
 
!     SORT ARRAY OF SILS
 
!     POSITION VECTOR
 
 DO  i = i1,i2
   IF (estbuf(i) == 0) GO TO 72
   k = 1
   DO  j = i1,i2
     IF (estbuf(j) - estbuf(i) < 0.0) THEN
       GO TO    60
     ELSE IF (estbuf(j) - estbuf(i) == 0.0) THEN
       GO TO    50
     ELSE
       GO TO    70
     END IF
     50 IF (j >= i) CYCLE
     60 IF (estbuf(j) /= 0) k = k + 1
     70 CONTINUE
   END DO
   GO TO 74
   72 k = l
   l = l - 1
   74 posvec(k) = i - i1 + 1
   sil(k) = estbuf(i)
 END DO
 
!     ELIMINATE DUP SILS THAT MAY OCCUR,E.G. CHBDY WITH AMB.PTS.
 
 k = 1
 icount = 1
 DO  i = 2,nsils
   82 k = k + 1
   IF (k <= nsils) GO TO 84
   sil(i) = 0
   posvec(i) = 0
   CYCLE
   84 IF (sil(k) == sil(k-1)) GO TO 82
   sil(i) = sil(k)
   IF (sil(k) /= 0) icount = icount + 1
   posvec(i) = posvec(k)
 END DO
 nsils = icount
 
!     SETUP VALUES AND DICTIONARY IN /EMG1BX/ FOR EMG1B USE
 
 dict(1) = estid
 dict(2) = 1
 
!     PSUEDO SMA1-SMA2 FILE NUMBERS
 
 ifkgg = 201
 if4gg = 202
 
!     DICT(4) WILL BE RESET TO EITHER 1 OR 63 BY EMG1B
!     BASED ON INCOMING DATA TO EMG1B
 
 DO  i = 5,15
   dict(i) = 0
 END DO
 
!     CALL ELEMENT FOR EACH PIVOT ROW
 
 last  = .false.
 knogo = 0
 mnogo = 0
 DO  i = 1,nsils
   IF (i == nsils) last = .true.
   
!     STIFFNESS MATRIX
   
   IF (flags(1) == 0) GO TO 550
   
!     RESTORE K-DICTIONARY IF NECESSARY
   
   IF (i == 1) GO TO 110
   DO  l = 1,15
     dict(l) = kdict(l)
   END DO
   110 CONTINUE
   filtyp = 1
   
!     IOPT4 IS TURNED ON SO THAT DAMPING CONSTANTS ARE SENT TO EMG1B
!     IN ALL AVAILABLE CASES BY ELEMENT ROUTINES.  MATRIX DATA WILL BE
!     IGNORED BY EMG1B ON EMG1B CALLS SENDING DAMPING CONSTANTS.
!     DAMPING CONSTANTS WILL BE PLACED IN 5TH WORD OF ELEMENT DICTIONARY
!     ENTRY.
   
   iopt4  = 1
   k4ggsw = 0
   knpvt  = sil(i)
   
!     FULL 6X6 MATRIX FORCED FOR STIFFNESS WITH OLD ELEMENT ROUTINES
   
   dict(2) = 1
   IF (sil(i) /= 0) GO TO 115
   CALL emg1b (dummy,0,1,1,0)
   GO TO 520
   115 CONTINUE
   DO  l = 1,estwds
     kecpt(l) = estbuf(l)
   END DO
   hydro = .false.
   IF (eltype >= 76 .AND. eltype <= 79) hydro = .true.
   
!     CALL THE PROPER ELEMENT STIFFNESS ROUTINE
   
   local = jltype - 100
   IF (local > 0) THEN
     GO TO   140
   END IF
   
!     PAIRED -GO TO- ENTRIES PER ELEMENT SINGLE/DOUBLE PRECISION
   
!             1 CROD      2 C.....    3 CTUBE     4 CSHEAR    5 CTWIST
   130 GO TO (  5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!             6 CTRIA1    7 CTRBSC    8 CTRPLT    9 CTRMEM   10 CONROD  &
   ,      190,  190,    5,    5,  210,  210,    5,    5,    5,    5     &
   
!            11 ELAS1    12 ELAS2    13 ELAS3    14 ELAS4    15 CQDPLT  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,  280,  280     &
   
!            16 CQDMEM   17 CTRIA2   18 CQUAD2   19 CQUAD1   20 CDAMP1  &
   ,      290,  290,  300,  300,  310,  310,  320,  320,    5,    5     &
   
!            21 CDAMP2   22 CDAMP3   23 CDAMP4   24 CVISC    25 CMASS1  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            26 CMASS2   27 CMASS3   28 CMASS4   29 CONM1    30 CONM2  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5  &
   
!            31 PLOTEL   32 C.....   33 C.....   34 CBAR     35 CCONEAX  &
   ,      520,  520,    5,    5,    5,    5,    5,    5,  340,  345  &
   
!            36 CTRIARG  37 CTRAPRG  38 CTORDRG  39 CTETRA   40 CWEDGE  &
   ,      350,  350,  370,  370,    5,    5,  400,  400,  410,  410  &
   
!            41 CHEXA1   42 CHEXA2   43 CFLUID2  44 CFLUID3  45 CFLUID4  &
   ,      420,  420,  430,  430,  440,  440,  450,  450,  460,  460  &
   
!            46 CFLMASS  47 CAXIF2   48 CAXIF3   49 CAXIF4   50 CSLOT3  &
   ,      520,  520,  440,  440,  450,  450,  460,  460,  470,  470  &
    ), jltype
   
!            51 CSLOT4   52 CHBDY    53 CDUM1    54 CDUM2    55 CDUM3
   140 GO TO (480,  480,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            56 CDUM4    57 CDUM5    58 CDUM6    59 CDUM7    60 CDUM8  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            61 CDUM9    62 CQDMEM1  63 CQDMEM2  64 CQDMEM3  65 CIHEX1  &
   ,        5,    5,  292,  292,  295,  294,    5,    5,    5,    5 &
   
!            66 CIHEX2   67 CIHEX3   68 CQUADTS  69 CTRIATS  70 CTRIAAX  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            71 CTRAPAX  72 CAERO1   73 CTRIM6   74 CTRPLT1  75 CTRSHL  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            76 CFHEX1   77 CFHEX2   78 CFTETRA  79 CFWEDGE  80 CIS2D8  &
   ,      420,  420,  430,  430,  400,  400,  410,  410,    5,    5 &
   
!            81 CELBOW   82 FTUBE    83 CTRIA3   84 CPSE2    85 CPSE3  &
   ,      390,  390,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            86 CPSE4  &
   ,        5,    5 &
    ), local
   
   
!     IN -HEAT- FORMULATIONS SOME ELEMENTS ARE IGNORED (OPTION(1)=HEAT)
!     IN STRUCTURE PROBLEMS SOME ELEMENTS ARE IGNORED (OPTION(1)=STRUCT)
   
   190 CALL ktriqd (1)
   GO TO 520
   210 CALL ktrplt
   GO TO 520
   280 CALL kqdplt
   GO TO 520
   290 CALL kqdmem
   GO TO 520
   
!     REPLACE ELEMENT TYPE CQDMEM1 BY ELEMENT TYPE CQDMEM
!     IN -HEAT- FORMULATION
   
   292 IF (heat) GO TO 290
   GO TO 5
   
!     REPLACE ELEMENT TYPE CQDMEM2 BY ELEMENT TYPE CQDMEM
!     IN -HEAT- FORMULATION
   
   294 IF (heat) GO TO 290
   GO TO 5
   
!     REPLACE ELEMENT TYPE CQDMEM2 BY ELEMENT TYPE CQDMEM
!     IN -HEAT- FORMULATION
   
   295 IF (heat) GO TO 290
   GO TO 5
   300 CALL ktriqd (2)
   GO TO 520
   310 CALL ktriqd (4)
   GO TO 520
   320 CALL ktriqd (3)
   GO TO 520
   340 CALL kcones
   GO TO 520
   343 CALL kcone2
   GO TO 520
   345 IF (mach ==  3) GO TO 340
   IF (nbpw >= 60) GO TO 343
   IF (qp == 0) CALL kconed
   IF (qp /= 0) CALL kconeq
   GO TO 520
   350 IF (heat) GO TO 360
   IF (knogo == 2) CYCLE
   CALL ktrirg
   IF (knogo == 2) CYCLE
   GO TO 520
   360 CALL hring (3)
   GO TO 520
   370 IF (heat) GO TO 380
   CALL ktrapr
   GO TO 520
   390 CALL kelbow
   GO TO 520
   380 CALL hring (4)
   GO TO 520
   400 CALL ktetra (0,0)
   GO TO 520
   410 CALL ksolid (1)
   GO TO 520
   420 CALL ksolid (2)
   GO TO 520
   430 CALL ksolid (3)
   GO TO 520
   440 CALL kflud2
   GO TO 520
   450 CALL kflud3
   GO TO 520
   460 CALL kflud4
   GO TO 520
   470 CALL kslot (0)
   GO TO 520
   480 CALL kslot (1)
   GO TO 520
   
!     OUTPUT THE PIVOT ROW PARTITION NOW COMPLETED BY -EMG1B-
   
   520 CALL emg1b (0.0D0,-1111111,0,0,0.0D0)
   
!     SAVE K-DICTIONARY
   
   DO  l = 1,15
     kdict(l) = dict(l)
   END DO
   
!     MASS MATRIX M
   
   550 IF (flags(2) == 0) GO TO 1090
   IF (heat) GO TO 1090
   
!     RESTORE M-DICTIONARY IF NECESSARY
   
   IF (i == 1) GO TO 570
   DO  l = 1,15
     dict(l) = mdict(l)
   END DO
   570 CONTINUE
   filtyp = 2
   ioptb  = 0
   bggind =-1
   mnpvt  = sil(i)
   dict(2)= 1
   IF (sil(i) /= 0)  GO TO 575
   CALL emg1b (dummy,0,1,2,0)
   GO TO 1060
   575 CONTINUE
   DO  l = 1,estwds
     mecpt(l) = estbuf(l)
   END DO
   
!     CALL THE PROPER ELEMENT MASS ROUTINE.
   
   590 local = jltype - 100
   IF (local > 0) THEN
     GO TO   610
   END IF
   
!     PAIRED -GO TO- ENTRIES PER ELEMENT SINGLE/DOUBLE PRECISION
   
!             1 CROD      2 C.....    3 CTUBE     4 CSHEAR    5 CTWIST
   600 GO TO (  5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!             6 CTRIA1    7 CTRBSC    8 CTRPLT    9 CTRMEM   10 CONROD  &
   ,      670,  670,    5,    5,  710,  710,    5,    5,    5,    5 &
   
!            11 ELAS1    12 ELAS2    13 ELAS3    14 ELAS4    15 CQDPLT  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,  740,  740 &
   
!            16 CQDMEM   17 CTRIA2   18 CQUAD2   19 CQUAD1   20 CDAMP1  &
   ,      760,  760,  770,  770,  790,  790,  810,  810,    5,    5 &
   
!            21 CDAMP2   22 CDAMP3   23 CDAMP4   24 CVISC    25 CMASS1  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            26 CMASS2   27 CMASS3   28 CMASS4   29 CONM1    30 CONM2  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            31 PLOTEL   32 C.....   33 C.....   34 CBAR     35 CCONEAX  &
   ,     1060, 1060,    5,    5,    5,    5,    5,    5,  910,  910 &
   
!            36 CTRIARG  37 CTRAPRG  38 CTORDRG  39 CTETRA   40 CWEDGE  &
   ,      930,  930,  940,  940,    5,    5,  960,  960,  970,  970 &
   
!            41 CHEXA1   42 CHEXA2   43 CFLUID2  44 CFLUID3  45 CFLUID4  &
   ,      980,  980,  990,  990, 1000, 1000, 1010, 1010, 1020, 1020 &
   
!            46 CFLMASS  47 CAXIF2   48 CAXIF3   49 CAXIF4   50 CSLOT3  &
   ,     1030, 1030, 1000, 1000, 1010, 1010, 1020, 1020, 1040, 1040 &
    ), jltype
   
   
!            51 CSLOT4   52 CHBDY    53 CDUM1    54 CDUM2    55 CDUM3
   610 GO TO (1050,1050,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            56 CDUM4    57 CDUM5    58 CDUM6    59 CDUM7    60 CDUM8  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            61 CDUM9    62 CQDMEM1  63 CQDMEM2  64 CQDMEM3  65 CIHEX1  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            66 CIHEX2   67 CIHEX3   68 CQUADTS  69 CTRIATS  70 CTRIAAX  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            71 CTRAPAX  72 CAERO1   73 CTRIM6   74 CTRPLT1  75 CTRSHL  &
   ,        5,    5,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            76 CFHEX1   77 CFHEX2   78 CFTETRA  79 CFWEDGE  80 CIS2D8  &
   ,     1060, 1060, 1060, 1060, 1060, 1060, 1060, 1060,    5,    5 &
   
!            81 CELBOW   82 FTUBE    83 CTRIA3   84 CPSE2    85 CPSE3  &
   ,      950,  950,    5,    5,    5,    5,    5,    5,    5,    5 &
   
!            86 CPSE4  &
   ,        5,    5 &
    ), local
   
   
!     CONVENTIONAL MASS MATRIX GENERATION ROUTINE CALLED WHEN
!     ICMBAR .LT. 0
!     OTHERWISE CONSISTENT MASS MATRIX GENERATION ROUTINE CALLED
   
   670 IF (icmbar < 0) GO TO 680
   CALL mtriqd (1)
   GO TO 1060
   680 CALL masstq (5)
   GO TO 1060
   
   710 IF (icmbar < 0) GO TO 720
   CALL mtrplt
   GO TO 1060
   720 CALL masstq (3)
   GO TO 1060
   
   740 IF (icmbar < 0) GO TO 750
   CALL mqdplt
   GO TO 1060
   750 CALL masstq (7)
   GO TO 1060
   760 CALL masstq (1)
   GO TO 1060
   
   770 IF (icmbar < 0) GO TO 780
   CALL mtriqd (2)
   GO TO 1060
   780 CALL masstq (4)
   GO TO 1060
   
   790 IF (icmbar < 0) GO TO 800
   CALL mtriqd (4)
   GO TO 1060
   800 CALL masstq (1)
   GO TO 1060
   
   810 IF (icmbar < 0) GO TO 820
   CALL mtriqd (3)
   GO TO 1060
   820 CALL masstq (2)
   GO TO 1060
   910 CALL mcone
   GO TO 1060
   930 IF (mnogo == 2) CYCLE
   IF (heat) GO TO 935
   CALL mtrirg
   IF (mnogo == 2) CYCLE
   GO TO 1060
   935 CALL mring (3)
   GO TO 1060
   940 IF (heat) GO TO 945
   CALL mtrapr
   GO TO 1060
   945 CALL mring (4)
   GO TO 1060
   950 CALL melbow
   GO TO 1060
   960 CALL msolid (1)
   GO TO 1060
   970 CALL msolid (2)
   GO TO 1060
   980 CALL msolid (3)
   GO TO 1060
   990 CALL msolid (4)
   GO TO 1060
   1000 CALL mflud2
   GO TO 1060
   1010 CALL mflud3
   GO TO 1060
   1020 CALL mflud4
   GO TO 1060
   1030 CALL mfree
   GO TO 1060
   1040 CALL mslot (0)
   GO TO 1060
   1050 CALL mslot (1)
   GO TO 1060
   
!     OUTPUT THE PIVOT ROW PARTITION NOW COMPLETED BY -EMG1B-
   
   1060 CALL emg1b (0.0D0,-1111111,0,0,0.0D0)
   IF (heat) GO TO 1185
   
!     SAVE M-DICTIONARY
   
   DO  l = 1,15
     mdict(l) = dict(l)
   END DO
   
!     DAMPING MATRIX B
   
   1090 IF (flags(3) == 0) CYCLE
   IF (.NOT.heat) CYCLE
   
!     RESTORE B-DICTIONARY IF NECESSARY
   
   IF (i == 1) GO TO 1110
   DO  l = 1,15
     dict(l) = bdict(l)
   END DO
   1110 filtyp = 3
   ioptb  =-1
   bggind =-1
   mnpvt  = sil(i)
   dict(2)= 1
   IF (sil(i) /= 0) GO TO 1115
   CALL emg1b (dummy,0,1,3,0)
   GO TO 1180
   1115 DO  l = 1,estwds
     mecpt(l) = estbuf(l)
   END DO
   GO TO 590
   
!     OUTPUT THE PIVOT ROW PARTITION NOW COMPLETED BY -EMG1B-
   
   1180 CALL emg1b (0.0D0,-1111111,0,0,0.0D0)
   
!     SAVE DICTIONARY
   
   1185 DO  l = 1,15
     bdict(l) = dict(l)
   END DO
   
 END DO
 IF (knogo == 0) knogo = ktemp
 IF (mnogo == 0) mnogo = mtemp
 
 1220 RETURN
END SUBROUTINE emgold
