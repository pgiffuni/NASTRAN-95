SUBROUTINE emgpro (ibuf)
     
!     THIS ROUTINE OF THE -EMG- MODULE IS THE MAIN PROCESSOR.  IT WILL
!     PASS THE -EST- DATA BLOCK ONCE, ELEMENT TYPE BY ELEMENT TYPE.
 
!     ELEMENT TYPES CONTRIBUTING TO STIFFNESS, MASS, OR DAMPING MATRICES
!     WILL BE PROCESSED.
 
 
 INTEGER, INTENT(IN)                      :: ibuf(7)
 LOGICAL :: anycon, error, heat
 INTEGER :: z, est, cstm, dit, geom2, dictn, savjcr, elid,  &
     outpt, eor, subr(2), eltype, precis, estbuf,  &
     elem, estwds, estid, savncr, dosi(2), flags, sil(32), sysbuf, scr3, scr4, ret
 DOUBLE PRECISION :: dummy
 DIMENSION        iz(1), ipos(32), trim6(2), trpl1(2), trshl(2), estx(12)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm, uwm, uim, sfm, swm
 COMMON /BLANK /  nok, nom, nob, nok4gg, nokdgg, nocmas, ncpbar,  &
     ncprod, ncpqd1, ncpqd2, ncptr1, ncptr2, ncptub,  &
     ncpqdp, ncptrp, ncptrb, volume, surfac
 COMMON /system/  ksystm(65)
 COMMON /gpta1 /  nelem, last, incr, elem(1)
 COMMON /emgfil/  est, cstm, mpt, dit, geom2, mats(3), dictn(3)
 COMMON /emgprm/  icore, jcore, ncore, icstm, ncstm, imat, nmat,  &
     ihmat, nhmat, idit, ndit, icong, ncong, lcong,  &
     anycon, flags(3), precis, error, heat,  &
     icmbar, lcstm, lmat, lhmat, kflags(3), l38
 COMMON /emgest/  estbuf(200)
 COMMON /emgdic/  eltype, ldict, nlocs, elid, estid
 COMMON /zzzzzz/  z(1)
 COMMON /iemgot/  nval(3)
 COMMON /matout/  egnu(6), rho
 COMMON /iemgod/  dummy, ktypes
 COMMON /iemg1b/  icall, ilast
 COMMON /sma1cl/  kdummy(22), knogo
 COMMON /sma2cl/  mdummy(20), mnogo
 EQUIVALENCE      (ksystm( 2),outpt), (ksystm( 1),sysbuf ),  &
     (ksystm(55),iprec), (estbuf( 1),estx(1)), (iz    ( 1),z(1) )
 DATA    trim6,         trpl1,         trshl               /  &
     4HCTRI,  4HM6  , 4HCTRP,4HLT1 , 4HCTRS,4HHL       /
 DATA    scr3  ,  scr4  / 303, 304 /
 DATA    eor   ,  noeor / 1, 0     /, subr / 4HEMGP,4HRO   /
 DATA    dosi  /  4HDOUB, 4HSING   /
 
 iqdmm1 = 0
 iqdmm2 = 0
 nval(1)= 0
 nval(2)= 0
 nval(3)= 0
 ltypes = 0
 ktypes = 0
 dummy  = 0.0D0
 icall  = 0
 ilast  = 0
 
!     INITIALIZE /SMA1CL/ AND /SMA2CL/
 
 knogo = 0
 mnogo = 0
 kdummy(10) = 10
 mdummy(10) = 10
 
!     FOLLOWING CALL PREPS /GPTA1/ FOR DUMMY ELEMENTS
 
 CALL delset
 
!     DEFINE WORKING CORE BLOCK FOR RESET PURPOSES.
 
 ipr = precis
 IF (ipr /= 1) ipr = 0
 savjcr = jcore
 savncr = ncore
 estid  = 0
 lnum   = lcong/2
 
!     READ THE ELEMENT TYPE FROM THE EST.
 
 10 CALL READ (*1340,*1360,est,eltype,1,noeor,iwords)
 izero = incr*(eltype-1)
 
!     CHECK FOR ALLOWABLE ELEMENT TYPES
 
 IF (eltype == 2 .OR. eltype == 32 .OR. eltype == 33 .OR.  &
     eltype == 68 .OR. eltype == 69 .OR. eltype == 72) GO TO 15
 IF (eltype >= 1 .AND. eltype <= nelem) GO TO 40
 15 WRITE  (outpt,20) sfm,elem(izero+1),elem(izero+2),eltype
 20 FORMAT (a25,' 3105, EMGPRO FINDS ',2A4,' ELEMENTS (ELEM. TYPE ',  &
     i3,') UNDEFINED IN EST DATA BLOCK AND/OR ELEMENT ROUTINE.')
 30 CALL fwdrec (*1350,est)
 error = .true.
 GO TO 10
 
!     RESTORE CORE POINTERS
 
 40 jcore = savjcr
 ncore = savncr
 
!     CLEAR ESTBUF
 
 DO  i = 1,200
   estbuf(i) = 0
 END DO
 
!     SET VARIOUS PARAMETERS = FUNCTION OF THIS ELEMENT TYPE
 
!     TURN ON COUPLED MASS FLAG IF EITHER OF ALL-COUPLED-MASS-FLAG
!     OR SPECIFIC-TYPE-COUPLED-MASS-FLAG IS ON.
 
 IF (flags(2) == 0.0) THEN
   GO TO    53
 END IF
 51 IF (nocmas   < 0) THEN
   GO TO    53
 ELSE IF (nocmas   == 0) THEN
   GO TO    52
 ELSE
   GO TO    54
 END IF
 52 IF (eltype == 34) IF (ncpbar) 53,53,54
 IF (eltype ==  1) IF (ncprod) 53,53,54
 IF (eltype == 19) IF (ncpqd1) 53,53,54
 IF (eltype == 18) IF (ncpqd2) 53,53,54
 IF (eltype ==  6) IF (ncptr1) 53,53,54
 IF (eltype == 17) IF (ncptr2) 53,53,54
 IF (eltype ==  3) IF (ncptub) 53,53,54
 IF (eltype == 15) IF (ncpqdp) 53,53,54
 IF (eltype ==  8) IF (ncptrp) 53,53,54
 IF (eltype ==  7) IF (ncptrb) 53,53,54
 53 icmbar = -1
 GO TO 56
 
 54 icmbar = 1
 
 56 jltype = 2*eltype - ipr
 estwds = elem(izero+12)
 nsils  = elem(izero+10)
 isil   = elem(izero+13)
 IF (elem(izero+9) /= 0) isil = isil - 1
 i1 = isil
 i2 = isil + nsils - 1
 isave2 = 0
 IF (estwds. LE. 200) GO TO 70
 WRITE  (outpt,60) sfm,eltype
 60 FORMAT (a25,' 3106, EMGPRO FINDS THAT ELEMENT TYPE ',i3,  &
     ' HAS EST ENTRIES TOO LARGE TO HANDLE CURRENTLY.')
 GO TO 30
 
!     CHECK TO SEE IF ILLEGAL ELEMENTS ARE USED IN -HEAT- FORMULATION
 
 70 IF (.NOT.heat) GO TO 80
 IF (eltype == 1 .OR.  eltype == 3 .OR. eltype == 6) GO TO 80
 IF (eltype >= 9 .AND. eltype <= 14) GO TO 80
 IF (eltype >= 16 .AND. eltype <= 24) GO TO 80
 IF (eltype == 34 .OR.  eltype == 36 .OR. eltype == 37) GO TO 80
 IF (eltype >= 39 .AND. eltype <= 42) GO TO 80
 IF (eltype == 52 .OR.  eltype == 62 .OR. eltype == 63) GO TO 80
 IF (eltype >= 64 .AND. eltype <= 67) GO TO 80
 IF (eltype == 80 .OR.  eltype == 81 .OR. eltype == 83) GO TO 80
 
 WRITE  (outpt,75) ufm,elem(izero+1),elem(izero+2),eltype
 75 FORMAT (a23,' 3115, EMGPRO FINDS ',2A4,' ELEMENTS (ELEMENT TYPE ',  &
     i3,') PRESENT IN A HEAT FORMULATION.')
 GO TO 30
 
!     SET UP VARIABLES TO BE WRITTEN AS DICTIONARY 3-WORD HEADER
 
 80 nlocs = nsils
 ldict = nlocs + 5
 
!     READ AN ELEMENT EST ENTRY
 
 90 CALL READ (*1350,*1200,est,estbuf,estwds,noeor,iwords)
 elid  = estbuf(1)
 estid = estid + 1
 
!     CHECK TO SEE IF THIS ELEMENT IS CONGRUENT TO ANOTHER ALREADY
!     POSSESSING A DICTIONARY IN CORE.
 
 IF (.NOT.anycon) GO TO 150
 CALL bisloc (*150,elid,z(icong),2,lnum,j)
 
!     MATCH FOUND.  CHECK FOR DICTIONARY-TABLE ON PRIMARY.
 
 iprime = z(icong+j  )
 idprim = z(icong+j-1)
 100 IF (iprime < 0) THEN
   GO TO   120
 ELSE IF (iprime == 0) THEN
   GO TO   104
 ELSE
   GO TO   110
 END IF
 
!     SET UP ELEMENT MATRIX MAPPING ARRAY FOR LATER USE BY OTHER
!     ELEMENTS IN THIS CONGRUENT SET
 
 104 icg    = jcore
 jjcore = jcore  + 2*nsils + 5
 icrq   = jjcore - ncore
 IF (jjcore .GE .ncore) GO TO 1800
 jcore   = jjcore
 iz(icg) = idprim
 iz(icg+1) = nsils
 iz(icg+2) = 0
 iz(icg+3) = 0
 iz(icg+4) = 0
 igoto = 0
 GO TO 1380
 
!     IPRIME POINTS TO PRIMARY ID
 
 110 idprim = z(iprime  )
 iprime = z(iprime+1)
 GO TO 100
 
!     IPRIME IS NEGATIVE TABLE ADDRESS IMPLYING DICTIONARY EXISTS.
 
 120 IF (error) GO TO 150
 iprime =-iprime
 imatch = 0
 ibfind = 1
 j = 0
 125 j = j + 1
 iadd = z(iprime+j)
 IF (iadd > 0) THEN
   GO TO   130
 ELSE
   GO TO   140
 END IF
 
!     COPY DICTIONARY FROM CORE TO DICTIONARY FILE.
 
 130 z(iadd)  = estid
 flags(j) = flags(j) + 1
 CALL WRITE (dictn(j),z(iadd),5,noeor)
 iaddd  = iadd + 5
 IF (imatch == 1) GO TO 135
 IF (imatch == 2) GO TO 1600
 indcng = savjcr
 igoto  = 1
 131 IF (iz(indcng) == idprim) GO TO 1380
 jjcore = indcng + 2*iz(indcng+1) + 5
 IF (jjcore >= ncore) GO TO 1820
 indcng = jjcore
 GO TO 131
 133 DO  l = 1,nsils
   IF (ipos(l) /= iz(indcng+nsils+l+4)) GO TO 137
 END DO
 imatch = 1
 135 CALL WRITE (dictn(j),z(iaddd),nsils,noeor)
 GO TO 140
 137 imatch = 2
 GO TO 1600
 140 ibfind = ibfind + 2
 IF (j < 3) GO TO 125
 GO TO 90
 
!     BRANCH ON ELEMENT TYPE.  INDIVIDUAL ROUTINES WILL COMPUTE AND
!     OUTPUT ALL MATRIX TYPES DESIRED BASED ON FLAGS AVAILABLE TO THEM.
 
 150 IF (eltype == ltypes) GO TO 152
 ltypes = eltype
 IF (ltypes > nelem) GO TO 15
 CALL page2 (3)
 WRITE  (outpt,151) uim,dosi(ipr+1),elem(izero+1),elem(izero+2), eltype,elid
 151 FORMAT (a29,' 3113,', /5X,'EMG MODULE PROCESSING ',a4,  &
     'LE PRECISION ',2A4,' ELEMENTS (ELEMENT TYPE ',i3, ') STARTING WITH ID ',i8)
 IF (eltype >= 84 .AND. eltype <= 86) WRITE (outpt,2300)
 152 IF (l38 /= 1) GO TO 154
 CALL page2 (1)
 WRITE  (outpt,153) elid
 153 FORMAT (5X,'ELEMENT ',i8,' IS BEING PROCESSED')
 154 local = jltype - 100
 IF (local > 0) THEN
   GO TO   156
 END IF
 
!     PAIRED -GO TO- ENTRIES PER ELEMENT SINGLE/DOUBLE PRECISION
 
!             1 CROD      2 C.....    3 CTUBE     4 CSHEAR    5 CTWIST
 155 GO TO (210,  215,   15,   15,  230,  235,  240,  245,  250,  255,
 
!             6 CTRIA1    7 CTRBSC    8 CTRPLT    9 CTRMEM   10 CONROD  &
 260,  265,  270,  275,  280,  285,  290,  295,  210,  215,
 
!            11 ELAS1    12 ELAS2    13 ELAS3    14 ELAS4    15 CQDPLT  &
 320,  320,  325,  325,  335,  335,  345,  345,  350,  355,
 
!            16 CQDMEM   17 CTRIA2   18 CQUAD2   19 CQUAD1   20 CDAMP1  &
 360,  365,  370,  375,  380,  385,  390,  395,  405,  405,
 
!            21 CDAMP2   22 CDAMP3   23 CDAMP4   24 CVISC    25 CMASS1  &
 415,  415,  425,  425,  435,  435,  440,  445,  455,  455,
 
!            26 CMASS2   27 CMASS3   28 CMASS4   29 CONM1    30 CONM2  &
 465,  465,  475,  475,  485,  485,  490,  495,  500,  505,
 
!            31 PLOTEL   32 C.....   33 C.....   34 CBAR     35 CCONEAX  &
 510,  515,   15,   15,   15,   15,  540,  545,  550,  555,
 
!            36 CTRIARG  37 CTRAPRG  38 CTORDRG  39 CTETRA   40 CWEDGE  &
 560,  565,  570,  575,  580,  585,  590,  595,  600,  605,
 
!            41 CHEXA1   42 CHEXA2   43 CFLUID2  44 CFLUID3  45 CFLUID4  &
 610,  615,  620,  625,  630,  635,  640,  645,  650,  655,
 
!            46 CFLMASS  47 CAXIF2   48 CAXIF3   49 CAXIF4   50 CSLOT3  &
 660,  665,  670,  675,  680,  685,  690,  695,  700,  705
  ), jltype
 
 
!            51 CSLOT4   52 CHBDY    53 CDUM1    54 CDUM2    55 CDUM3
 156 GO TO (710,  715,  720,  725,  730,  730,  740,  740,  750,  750,
 
!            56 CDUM4    57 CDUM5    58 CDUM6    59 CDUM7    60 CDUM8  &
 760,  760,  770,  770,  780,  780,  790,  790,  800,  800,
 
!            61 CDUM9    62 CQDMEM1  63 CQDMEM2  64 CQUAD4   65 CIHEX1  &
 810,  810,  820,  825,  830,  835,  950,  955,  850,  855,
 
!            66 CIHEX2   67 CIHEX3   68 CQUADTS  69 CTRIATS  70 CTRIAAX  &
 850,  855,  850,  855,   15,   15,   15,   15,  880,  885,
 
!            71 CTRAPAX  72 CAERO1   73 CTRIM6   74 CTRPLT1  75 CTRSHL  &
 890,  895,   15,   15,  900,  905,  910,  915,  920,  925,
 
!            76 CFHEX1   77 CFHEX2   78 CFTETRA  79 CFWEDGE  80 CIS2D8  &
 610,  615,  620,  625,  590,  595,  600,  605,  930,  935,
 
!            81 CELBOW   82 FTUBE    83 CTRIA3   84 CPSE2    85 CPSE3  &
 940,  945,  840,  840,  960,  965,   90,   90,   90,   90,
 
!            86 CPSE4  &
 90,   90
  ), local
 
!     ==================================================================
!     A WALKING TOUR OF EMG TO COMPUTE STIFFNESS (K-) AMD MASS (M-)
!     MATRICES FOR AN 'OLD' ELEMENT SUCH AS CTRIA2.
!     SEE HOW EASY IT IS.      G.CHAN/UNISYS, 7/87
 
!       EMG SUPPORTING ROUTINES -
!       EMGTAB,EMGCNG,EMGCOR,EMGFIN,
!       EMGSOC (WHICH COMPUTES OFFSET BETWWEN /ZZEMGX/ AND /ZZEMII/ AND
!           /   SETS ICORE,JCORE,NCORE IN /EMGPRM/ FOR OPEN CORE USAGE)
!          /
!         /                                --->EMG1B---->EMGOUT
!        /                                /   OUTPUT PIVOT ROW PARTITION
!     EMG---->EMGPRO---->CTRIA2          /    AFTER KTRIQD IS DONE, AND
!               /      AN ENTRY POINT   /     ALSO AFTER MTRIQD
!            /ZZEMGX/     IN           /                             (*)
!                        OLDEL3--->EMGOLD--->KTRIQD--->KTRMEM--->KTRPLT
!                                    /                              /
!                                    ------->MTRIQD             /ZZEM14/
!            (*)                                 (&)          UNIT 14 IS
!             KTRPLT---------------->KTRBSC                   ALLOCATED
!           TO COMPUTE BENDING     TO COMPUTE MEMBRANE        TO CTRIA2
!           FOR CTRIA2             FOR CTRIA2                 BY TA1ABD
!                  /                      /
!                 ------->SMA1B<----------
!                           \
!                            ---->EMG1B---->EMGOUT
!                                          OUTPUT A K-MATRIX
!           (&)                            PARTITION
!             MTRIQD BRANCH, FOR M-MATRIX
!             FOR CTRIA2 ELEMENT, IS SIMILARLY
!             STRUCTURED AS THAT OF THE KTRIQD BRANCH
 
!             REPEAT DAMPING B-MATRIX IF NECESSARY
!             IF ELEMENT HAS HEAT CAPBABILITY - WHAT DO I DO NOW?
 
!     THIS SYMBOL '>' IS RIGHT ARROW HEAD, AND '<' IS LEFT ARROW HEAD
!     ==================================================================
 
 210 CALL rods
 GO TO 90
 215 CALL rodd
 GO TO 90
 230 CALL tubes
 GO TO 90
 235 CALL tubed
 GO TO 90
 240 CALL shears
 GO TO 90
 245 CALL sheard
 GO TO 90
 250 CALL twists
 GO TO 90
 255 CALL twistd
 GO TO 90
 260 CALL tria1s
 GO TO 300
 265 CALL tria1d
 GO TO 300
 270 CALL trbscs
 GO TO 300
 275 CALL trbscd
 GO TO 300
 280 CALL trplts
 GO TO 300
 285 CALL trpltd
 GO TO 300
 290 CALL trmems
 GO TO 300
 295 CALL trmemd
 300 kht = 7
 l   = 9
 IF (volume == 0 .AND. surfac == 0) GO TO 90
 CALL WRITE (scr4,elem(izero+1),2,0)
 CALL WRITE (scr4,estbuf(1),1,0)
 estx(5)   = estx(kht)
 estx(6)   = rho
 estbuf(7) = 3
 CALL WRITE (scr4,estbuf(5), 3,0)
 CALL WRITE (scr4,estbuf(2), 3,0)
 CALL WRITE (scr4,estbuf(l),12,1)
 GO TO 90
 310 kht = 8
 l   = 10
 315 IF (volume == 0 .AND. surfac == 0) GO TO 90
 CALL WRITE (scr4,elem(izero+1),2,0)
 CALL WRITE (scr4,estbuf(1),1,0)
 estx(5)   = estx(kht)
 estx(6)   = rho
 estbuf(7) = 4
 CALL WRITE (scr4,estbuf(5), 3,0)
 CALL WRITE (scr4,estbuf(2), 4,0)
 CALL WRITE (scr4,estbuf(l),16,1)
 GO TO 90
 320 nscal1 = 1
 GO TO 346
 325 nscal1 = 2
 GO TO 346
 335 nscal1 = 3
 GO TO 346
 345 nscal1 = 4
 346 nscal2 = 1
 GO TO 487
 350 CALL qdplts
 352 kht = 10
 l   = 14
 GO TO 315
 355 CALL qdpltd
 GO TO 352
 360 CALL qdmems
 GO TO 310
 365 CALL qdmemd
 GO TO 310
 370 CALL tria2s
 GO TO 300
 375 CALL tria2d
 GO TO 300
 380 CALL quad2s
 GO TO 310
 385 CALL quad2d
 GO TO 310
 390 CALL quad1s
 392 IF (estx(12) <= 0.0) estx(12) = estx(8)
 kht = 12
 l   = 14
 GO TO 315
 395 CALL quad1d
 GO TO 392
 405 nscal1 = 1
 GO TO 436
 415 nscal1 = 2
 GO TO 436
 425 nscal1 = 3
 GO TO 436
 435 nscal1 = 4
 436 nscal2 = 3
 GO TO 487
 440 CALL viscs
 IF (flags(3) == 0) WRITE (outpt,442) uwm
 442 FORMAT (a25,' 2422, VISC DATA NOT PROCESSED BY EMGPRO.')
 GO TO 90
 445 CALL viscd
 IF (flags(3) == 0) WRITE (outpt,442) uwm
 GO TO 90
 455 nscal1 = 1
 GO TO 486
 465 nscal1 = 2
 GO TO 486
 475 nscal1 = 3
 GO TO 486
 485 nscal1 = 4
 486 nscal2 = 2
 487 CALL scaled (nscal1,nscal2)
 GO TO 90
 490 CALL conm1s
 GO TO 90
 495 CALL conm1d
 GO TO 90
 500 CALL conm2s
 GO TO 90
 505 CALL conm2d
 GO TO 90
 510 CALL plotls
 GO TO 90
 515 CALL plotld
 GO TO 90
 540 CALL bars
 GO TO 90
 545 CALL bard
 GO TO 90
 550 CALL cones
 GO TO 90
 555 CALL coned
 GO TO 90
 560 CALL triars
 GO TO 90
 565 CALL triard
 GO TO 90
 570 CALL traprs
 GO TO 90
 575 CALL traprd
 GO TO 90
 580 CALL tordrs
 GO TO 90
 585 CALL tordrd
 GO TO 90
 590 CALL tetras
 GO TO 90
 595 CALL tetrad
 GO TO 90
 600 CALL wedges
 GO TO 90
 605 CALL wedged
 GO TO 90
 610 CALL hexa1s
 GO TO 90
 615 CALL hexa1d
 GO TO 90
 620 CALL hexa2s
 GO TO 90
 625 CALL hexa2d
 GO TO 90
 630 CALL flud2s
 GO TO 90
 635 CALL flud2d
 GO TO 90
 640 CALL flud3s
 GO TO 90
 645 CALL flud3d
 GO TO 90
 650 CALL flud4s
 GO TO 90
 655 CALL flud4d
 GO TO 90
 660 CALL flmass
 GO TO 90
 665 CALL flmasd
 GO TO 90
 670 CALL axif2s
 GO TO 90
 675 CALL axif2d
 GO TO 90
 680 CALL axif3s
 GO TO 90
 685 CALL axif3d
 GO TO 90
 690 CALL axif4s
 GO TO 90
 695 CALL axif4d
 GO TO 90
 700 CALL slot3s
 GO TO 90
 705 CALL slot3d
 GO TO 90
 710 CALL slot4s
 GO TO 90
 715 CALL slot4d
 GO TO 90
 720 CALL hbdys
 GO TO 90
 725 CALL hbdyd
 GO TO 90
 730 CALL kdum1
 GO TO 90
 740 CALL kdum2
 GO TO 90
 750 CALL kdum3
 GO TO 90
 760 CALL kdum4
 GO TO 90
 770 CALL kdum5
 GO TO 90
 780 CALL kdum6
 GO TO 90
 790 CALL kdum7
 GO TO 90
 800 CALL kdum8
 GO TO 90
 810 CALL kdum9
 GO TO 90
 820 IF (.NOT.heat) GO TO 822
 IF (iqdmm1 /= 0) GO TO 360
 ASSIGN 360 TO ret
 iqdmm1 = 1
 GO TO 1000
 822 CALL qdmm1s
 GO TO 310
 825 IF (.NOT.heat) GO TO 827
 IF (iqdmm1 /= 0) GO TO 365
 ASSIGN 365 TO ret
 iqdmm1 = 1
 GO TO 1000
 827 CALL qdmm1d
 GO TO 310
 830 IF (.NOT.heat) GO TO 832
 IF (iqdmm2 /= 0) GO TO 360
 ASSIGN 360 TO ret
 iqdmm2 = 1
 GO TO 1000
 832 CALL qdmm2s
 GO TO 310
 835 IF (.NOT.heat) GO TO 837
 IF (iqdmm2 /= 0) GO TO 365
 ASSIGN 365 TO ret
 iqdmm2 = 1
 GO TO 1000
 837 CALL qdmm2d
 GO TO 310
 840 CALL ftube
 GO TO 90
 850 CALL ihexs (eltype-64)
 GO TO 90
 855 CALL ihexd (eltype-64)
 GO TO 90
 880 CALL triaax
 GO TO 90
 885 CALL triaad
 GO TO 90
 890 CALL trapax
 GO TO 90
 895 CALL trapad
 GO TO 90
 900 CALL ktrm6s
 l = 14
 GO TO 927
 905 CALL ktrm6d
 l = 14
 GO TO 927
 910 CALL ktrpls
 l = 24
 GO TO 927
 915 CALL ktrpld
 l = 24
 GO TO 927
 920 CALL ktshls
 l = 28
 GO TO 927
 925 CALL ktshld
 l = 28
 927 IF (volume == 0.0 .AND. surfac == 0.0) GO TO 90
 estx(8) = elem(izero+1)
 estx(9) = elem(izero+2)
 IF (estx(11) <= 0.0) estx(11) = estx(10)
 IF (estx(12) <= 0.0) estx(12) = estx(10)
 thk = (estx(10) + estx(11) + estx(12))/3.
 estbuf(10) = estbuf(1)
 estx  (11) = thk
 estx  (12) = rho
 estbuf(13) = 6
 CALL WRITE (scr4,estbuf(8), 6,0)
 CALL WRITE (scr4,estbuf(2), 6,0)
 CALL WRITE (scr4,estbuf(l),24,1)
 GO  TO  90
 930 CALL is2d8s
 GO TO 90
 935 CALL is2d8d
 GO TO 90
 940 CALL elbows
 GO TO 90
 945 CALL elbowd
 GO TO 90
 950 CALL quad4s
 GO TO 90
 955 CALL quad4d
 GO TO 90
 960 CALL tria3s
 GO TO 90
 965 CALL tria3d
 GO TO 90
 
!     PRINT WARNING MESSAGE TO INDICATE THAT QDMEM1 ELEMENTS
!     (ELEMENT TYPE 62) AND QDMEM2 ELEMENTS (ELEMENT TYPE 63)
!     ARE REPLACED BY QDMEM ELEMENTS (ELEMENT TYPE 16) IN
!     -HEAT- FORMULATION
 
 1000 INDEX  = 15*incr
 index1 = 16
 CALL page2 (3)
 WRITE  (outpt,1100) uwm,elem(izero+1),elem(izero+2),eltype,  &
     elem(INDEX+1),elem(INDEX+2),index1
 1100 FORMAT (a25,' 3144, EMGPRO FINDS ',2A4,' ELEMENTS (ELEMENT TYPE ',  &
     i3,') PRESENT IN A HEAT FORMULATION AND IS',/5X,'REPLACING'  &
     ,      ' THE SAME BY ',2A4,' ELEMENTS (ELEMENT TYPE ',i3,2H).)
 GO TO ret, (360,365)
 
!     ALL ELEMENTS OF THIS ELEMENT TYPE PROCESSED.
!     COMPLETE DICTIONARY RECORD FOR ELEMENT TYPE.
 
 1200 IF (error) GO TO 1310
 DO  i = 1,3
   IF (flags(i) <= 0) CYCLE
   flags(i) = -flags(i)
   CALL WRITE (dictn(i),0,0,eor)
 END DO
 
!     FOR SAFETY AND IF CONGRUENCY EXISTS CLEAR OFF ANY TABLE POINTERS
!     ON PRIMARY-IDS IN THE CONGRUENCY LIST
 
 1310 IF (.NOT.anycon) GO TO 10
 DO  i = icong,ncong,2
   IF (z(i+1) < 0) z(i+1) = 0
 END DO
 GO TO 10
 
!     ALL ELEMENT TYPES HAVE BEEN PROCESSED.
 
 1340 CONTINUE
 IF (knogo > 0 .OR. mnogo > 0) CALL mesage (-61,0,0)
 RETURN
 
!     IMPROPER ENCOUNTER OF AN -EOF-
 
 1350 jfile = est
 1355 CALL mesage (-2,jfile,subr)
 
!     IMPROPER ENCOUNTER OF AN -EOR-
 
 1360 jfile = est
 CALL mesage (-3,jfile,subr)
 
!     FILE NOT IN FIST
 
 1370 CALL mesage (-1,jfile,subr)
 
!     COMPUTE MAPPING DATA FOR CONGRUENT ELEMENTS
 
 1380 l1 = nsils
 DO  l = i1,i2
   IF (estbuf(l) == 0) GO TO 1420
   m = 1
   DO  n = i1,i2
     IF (estbuf(n)-estbuf(l) < 0.0) THEN
       GO TO  1400
     ELSE IF (estbuf(n)-estbuf(l) == 0.0) THEN
       GO TO  1390
     ELSE
       GO TO  1410
     END IF
     1390 IF (n >= l) CYCLE
     1400 IF (estbuf(n) /= 0) m = m + 1
   END DO
   GO TO 1425
   1420 m  = l1
   l1 = l1 - 1
   1425 ipos(m) = l - i1 + 1
   sil(m) = estbuf(l)
 END DO
 IF (igoto == 1) GO TO 133
 DO  l = i1,i2
   l1 = l - i1 + 1
   DO  n = 1,nsils
     IF (estbuf(l) /= sil(n)) CYCLE
     iz(icg+l1+4) = n
     EXIT
   END DO
   1450 iz(icg+nsils+l1+4) = ipos(l1)
 END DO
 GO TO 150
 
!     CHECK IF THE ELEMENT MATRIX IS DIAGONAL
 
 1600 IF (iz(iadd+1) /= 2) GO TO 1604
 
!     ELEMENT MATRIX IS DIAGONAL.
!     RE-WRITE ONLY THE ELEMENT DICTIONARY FOR A CONGRUENT ELEMENT.
 
 DO  l = 1,nsils
   m = ipos(l)
   n = iz(indcng+m+4) - 1
   CALL WRITE (dictn(j),z(iaddd+n),1,noeor)
 END DO
 GO TO 140
 
!     ELEMENT MATRIX IS SQUARE.
!     PICK UP ELEMENT MATRIX DATA FOR A CONGRUENT ELEMENT THAT HAS
!     ALREADY BEEN PROCESSED AND STORE IT ON SCR3.
 
 1604 ibuf1 = ncore - sysbuf - 2
 icrq  = jcore - ibuf1
 IF (icrq > 0) GO TO 1840
 ibuf3 = ibuf1 - 1
 icrq  = jcore - ibuf3 + 36*nsils*iprec
 IF (icrq > 0) GO TO 1840
 ifile = mats(j)
 IF (iz(indcng+j+1) /= 0) GO TO 1640
 CALL savpos (ifile,isave1)
 CALL CLOSE  (ifile,1)
 ibuf2 = ibuf(ibfind+1)
 CALL gopen  (ifile,z(ibuf2),0)
 CALL filpos (ifile,z(iaddd))
 IF (isave2 /= 0) GO TO 1605
 CALL gopen (scr3,z(ibuf1),1)
 GO TO 1607
 1605 jfile = scr3
 CALL OPEN (*1370,scr3,z(ibuf1),3)
 1607 jfile = ifile
 DO  l1 = 1,nsils
   CALL READ  (*1355,*1610,ifile,z(jcore),ibuf3,eor,n)
   1610 CALL WRITE (scr3,z(jcore),n,eor)
   IF (l1 == 1) CALL savpos (scr3,iz(indcng+j+1))
 END DO
 CALL filpos (ifile,isave1)
 CALL skprec (ifile,1)
 CALL CLOSE  (ifile,2)
 CALL OPEN   (*1370,ifile,z(ibuf2),3)
 CALL savpos (scr3,isave2)
 CALL CLOSE  (scr3,1)
 
!     ELEMENT MATRIX DATA IS AVAILABLE ON SCR3.  REARRANGE IT IN
!     THE REQUIRED ORDER AND WRITE IT ON THE OUTPUT DATA BLOCK.
 
 1640 CALL gopen (scr3,z(ibuf1),0)
 jfile = scr3
 DO  l = 1,nsils
   CALL filpos (scr3,iz(indcng+j+1))
   m = ipos(l)
   n = iz(indcng+m+4) - 1
   CALL skprec (scr3,n)
   CALL READ (*1355,*1650,scr3,z(jcore),ibuf3,eor,n)
   1650 nnwrds = n/(nsils*iprec)
   nnwrds = SQRT(nnwrds+0.5)
   nwords = nnwrds*iprec
   jjcore = jcore
   DO  l2 = 1,nnwrds
     DO  l1 = 1,nsils
       m = ipos(l1)
       n = iz(indcng+m+4) - 1
       CALL WRITE (ifile,z(jjcore+n*nwords),nwords,noeor)
     END DO
     jjcore = jjcore + nwords*nsils
   END DO
   CALL WRITE  (ifile,0,0,1)
   CALL savpos (ifile,isave1)
   CALL WRITE  (dictn(j),isave1,1,noeor)
 END DO
 CALL filpos (scr3,isave2)
 CALL skprec (scr3,1)
 CALL CLOSE  (scr3,2)
 GO TO 140
 1800 WRITE (outpt,2000) uim,idprim
 WRITE (outpt,2400) icrq
 GO TO 1850
 1820 WRITE (outpt,2100) swm,estid
 GO TO 1850
 1840 WRITE (outpt,2200) uim,estid
 WRITE (outpt,2400) icrq
 1850 CALL page2 (4)
 GO TO 150
 
 2000 FORMAT (a29,' 2382, ELEMENT MATRICES FOR ELEMENTS CONGRUENT TO ',  &
     'ELEMENT ID =',i10, /5X,'WILL BE RE-COMPUTED AS THERE IS',  &
     ' INSUFFICIENT CORE AT THIS TIME TO HOLD CONGRUENCY ', 'MAPPING DATA.')
 2100 FORMAT (a27,' 2383, UNABLE TO LOCATE CONGRUENCY MAPPING DATA FOR',  &
     ' ELEMENT ID =',i10,1H., /5X,'ELEMENT MATRICES FOR THIS ',  &
     'ELEMENT WILL, THEREFORE, BE RE-COMPUTED.')
 2200 FORMAT (a29,' 2384, CONGRUENCY OF ELEMENT ID =',i10,  &
     ' WILL BE IGNORED AND ITS ELEMENT MATRICES', /5X,  &
     'WILL BE RE-COMPUTED AS THERE IS INSUFFICIENT CORE AT ',  &
     'THIS TIME TO PERFORM CONGRUENCY MAPPING COMPUTATIONS.')
 2300 FORMAT (5X,'(STEPPING THRU ONLY. NO REAL COMPUTATION HERE FOR ',  &
     'THIS DIFFERENTIAL STIFFNESS ELEMENT)')
 2400 FORMAT (5X,'ADDITIONAL CORE NEEDED =',i9,' WORDS.')
 
END SUBROUTINE emgpro
