SUBROUTINE ofp1c (line)
     
!     THIS SUBROUTINE WAS FORMED ONLY TO REDUCE THE SIZE OF OFP1 FOR
!     COMPILATION PURPOSES.  IT IS CALLED ONLY BY OFP1.
!     THIS ROUTINE WAS PART OF OFP1B BEFORE.
 
 
 INTEGER, INTENT(IN)                      :: line
 COMMON /system/ ibuf,l
!ZZ   COMMON /ZZOFPX/ L123(1)
 COMMON /zzzzzz/ l123(20000)
 
 
!WKBR NCL93012 3/94      IF (LINE .GT. 467) GO TO 100
!WKBR SPR94001 7/94      IF (LINE .GT. 470) GO TO 100
 IF (line > 474) GO TO 100
 local = line - 380
 GO TO (381,382,383,384,385,386,387,388,389,390,  &
     391,392,393,394,395,396,397,398,399,400,  &
     401,402,403,404,405,406,407,408,409,410,  &
     411,412,413,414,415,416,417,418,419,420,  &
     421,422,423,424,425,426,427,428,429,430,  &
     431,432,433,434,435,436,437,438,439,440,  &
     441,442,443,444,445,446,447,448,449,450,  &
     451,452,453,454,455,456,457,458,459,460,
!WKBR NCL93012 3/94 8       461,462,463,464,465,466,467), LOCAL
!WKBD SPR94001 7/94 8       461,462,463,464,465,466,467,100,469,470), LOCAL
!WKBNB SPR94001 7/94  &
 461,462,463,464,465,466,467,100,469,470, 471,472,473,474), local
!WKBNE SPR94001 7/94
 
 100 WRITE  (l,110) line
 110 FORMAT ('0*** OFP ERROR/OFP1C,  LINE=',i9)
 CALL mesage (-61,0,0)
 
 381 WRITE (l,881)
 GO TO 1000
 382 WRITE (l,882)
 GO TO 1000
 383 WRITE (l,883)
 GO TO 1000
 384 WRITE (l,884)
 GO TO 1000
 385 WRITE (l,885)
 GO TO 1000
 386 WRITE (l,886)
 GO TO 1000
 387 WRITE (l,887)
 GO TO 1000
 388 WRITE (l,888)
 GO TO 1000
 389 WRITE (l,889)
 GO TO 1000
 390 WRITE (l,890)
 GO TO 1000
 391 WRITE (l,891)
 GO TO 1000
 392 WRITE (l,892)
 GO TO 1000
 393 WRITE (l,893)
 GO TO 1000
 394 WRITE (l,894)
 GO TO 1000
 395 WRITE (l,895)
 GO TO 1000
 396 WRITE (l,896)
 GO TO 1000
 397 WRITE (l,897)
 GO TO 1000
 398 WRITE (l,898)
 GO TO 1000
 399 WRITE (l,899)
 GO TO 1000
 400 WRITE (l,900)
 GO TO 1000
 401 WRITE (l,901)
 GO TO 1000
 402 WRITE (l,902)
 GO TO 1000
 403 WRITE (l,903)
 GO TO 1000
 404 WRITE (l,904)
 GO TO 1000
 405 WRITE (l,905)
 GO TO 1000
 406 WRITE (l,906)
 GO TO 1000
 407 WRITE (l,907)
 GO TO 1000
 408 WRITE (l,908)
 GO TO 1000
 409 WRITE (l,909)
 GO TO 1000
 410 WRITE (l,910)
 GO TO 1000
 411 WRITE (l,911)
 GO TO 1000
 412 WRITE (l,912)
 GO TO 1000
 413 WRITE (l,913)
 GO TO 1000
 414 WRITE (l,914)
 GO TO 1000
 415 WRITE (l,915)
 GO TO 1000
 416 WRITE (l,916)
 GO TO 1000
 417 WRITE (l,917)
 GO TO 1000
 418 WRITE (l,918)
 GO TO 1000
 419 WRITE (l,919)
 GO TO 1000
 420 WRITE (l,920)
 GO TO 1000
 421 WRITE (l,921)
 GO TO 1000
 422 WRITE (l,922)
 GO TO 1000
 423 WRITE (l,923)
 GO TO 1000
 424 WRITE (l,924)
 GO TO 1000
 425 WRITE (l,925)
 GO TO 1000
 426 WRITE (l,926)
 GO TO 1000
 427 WRITE (l,927)
 GO TO 1000
 428 WRITE (l,928)
 GO TO 1000
 429 WRITE (l,929)
 GO TO 1000
 430 WRITE (l,930)
 GO TO 1000
 431 WRITE (l,931)
 GO TO 1000
 432 WRITE (l,932)
 GO TO 1000
 433 WRITE (l,933)
 GO TO 1000
 434 WRITE (l,934)
 GO TO 1000
 435 WRITE (l,935)
 GO TO 1000
 436 WRITE (l,936)
 GO TO 1000
 437 WRITE (l,937)
 GO TO 1000
 438 WRITE (l,938)
 GO TO 1000
 439 WRITE (l,939)
 GO TO 1000
 440 WRITE (l,940)
 GO TO 1000
 441 WRITE (l,941)
 GO TO 1000
 442 WRITE (l,942)
 GO TO 1000
 443 WRITE (l,943)
 GO TO 1000
 444 WRITE (l,944)
 GO TO 1000
 445 WRITE (l,945)
 GO TO 1000
 446 WRITE (l,946)
 GO TO 1000
 447 WRITE (l,947)
 GO TO 1000
 448 WRITE (l,948)
 GO TO 1000
 449 WRITE (l,949)
 GO TO 1000
 450 WRITE (l,950)
 GO TO 1000
 451 WRITE (l,951)
 GO TO 1000
 452 WRITE (l,952)
 GO TO 1000
 453 WRITE (l,953)
 GO TO 1000
 454 WRITE (l,954)
 GO TO 1000
 455 WRITE (l,955)
 GO TO 1000
 456 WRITE (l,956)
 GO TO 1000
 457 WRITE (l,957)
 GO TO 1000
 458 WRITE (l,958)
 GO TO 1000
 459 WRITE (l,959)
 GO TO 1000
 460 WRITE (l,960)
 GO TO 1000
 461 WRITE (l,961)
 GO TO 1000
 462 WRITE (l,962)
 GO TO 1000
 463 WRITE (l,963)
 GO TO 1000
 464 WRITE (l,964)
 GO TO 1000
 465 WRITE (l,965)
 GO TO 1000
 466 WRITE (l,966)
 GO TO 1000
 467 WRITE (l,967)
 GO TO 1000
!WKBNB NCL93012 3/94
 469   WRITE ( l,969)
 GO TO 1000
 470   WRITE ( l,970)
 GO TO 1000
!WKBNE NCL93012 3/94
!WKBNB SPR94001 7/94
 471   WRITE (l,971)
 GO TO 1000
 472   WRITE (l,972)
 GO TO 1000
 473   WRITE (l,973)
 GO TO 1000
 474   WRITE (l,974)
 GO TO 1000
!WKBNE SPR94001 7/94
 1000 RETURN
 
!     ******************************************************************
 
 881 FORMAT (4X,'S T R A I N S / C U R V A T U R E S   I N   G E N E ',  &
     'R A L   Q U A D R I L A T E R A L   E L E M E N T S',6X, '( C Q U A D 2 )')
 882 FORMAT (4X,'S T R A I N S / C U R V A T U R E S   I N   G E N E ',  &
     'R A L   Q U A D R I L A T E R A L   E L E M E N T S',6X, '( C Q U A D 1 )')
!WKBRB NCL93012 3/94
!  883 FORMAT (2X,7HELEMENT,24X,37HSTRNS./CURVS. IN ELEMENT COORD SYSTEM,
!     1       6X,38HPRIN. STRNS./CURVS. (ZERO SHEAR/TWIST),7X,7HMAXIMUM)
!WKBRE NCL93012 3/94
 883 FORMAT (2X,7HELEMENT,8X,'STRAIN',8X  &
     ,      37HSTRNS./curvs. in element coord system  &
     ,      6X,38HPRIN. strns./curvs. (zero shear/twist),7X,7HMAXIMUM)
 884 FORMAT (4X,3HID.,6X,15HID./output code,5X,8HNORMAL-x,7X,  &
     8HNORMAL-y,6X,8HSHEAR-xy,7X,5HANGLE,9X,5HMAJOR,11X,5HMINOR,  &
     7X,11HSHEAR/twist)
 885 FORMAT (2X,7HELEMENT,4X,16HMAT. coord. sys.,4X,'STRNS./CURVS. ',  &
     ' IN MATERIAL COORD SYSTEM',5X,  &
     38HPRIN. strns./curvs. (zero shear/twist),7X,7HMAXIMUM)
 886 FORMAT (33X,'S T R A I N S / C U R V A T U R E S   A T   G R I D',  &
     '   P O I N T S')
 887 FORMAT (2X,7H point ,4X,16HMAT. coord. sys.,6X,  &
     33HSTRESSES inmaterial coord system , 12X,  &
     31HPRINCIPAL stresses (zero shear), 12X,3HMAX)
 888 FORMAT (2X,7H point ,4X,16HMAT. coord. sys.,4X,  &
     38HSTRNS./curvs. in material coord system, 5X,  &
     38HPRIN. strns./curvs. (zero shear/twist), 7X,7HMAXIMUM)
 889 FORMAT (50X,30H(in element coordinate system),/)
 890 FORMAT (50X,31H(in material coordinate system),/)
!WKBRB NCL93012 3/94
!  891 FORMAT (4X,3HID.,26X,8HNORMAL-X, 7X,8HNORMAL-Y, 6X,8HSHEAR-XY,
!     1       7X,5HANGLE, 9X,5HMAJOR, 11X,5HMINOR, 7X,11HSHEAR/TWIST)
!WKBRE NCL93012 3/94
 891 FORMAT (4X,3HID.,9X,'CURVATURE',7X  &
     ,      8HNORMAL-x, 7X,8HNORMAL-y, 6X,8HSHEAR-xy  &
     ,      7X,5HANGLE, 9X,5HMAJOR, 11X,5HMINOR, 7X,11HSHEAR/twist)
 892 FORMAT (4X,'C O M P L E X   F O R C E S   I N   A X I S - S Y M ',  &
     'M E T R I C   T R I A N G U L A R   R I N G   E L E M E ',  &
     'N T S   (CTRIAAX)',/)
 893 FORMAT (2X,'C O M P L E X   S T R E S S E S   I N   A X I S - S ',  &
     'Y M M E T R I C   T R I A N G U L A R   R I N G   E L E ',  &
     'M E N T S   (CTRIAAX)',/)
 894 FORMAT (3X,'C O M P L E X   F O R C E S   I N   A X I S - S Y M ',  &
     'M E T R I C   T R A P E Z O I D A L   R I N G   E L E M ',  &
     'E N T S   (CTRAPAX)',/)
 895 FORMAT (' C O M P L E X   S T R E S S E S   I N   A X I S - S Y ',  &
     ' M M E T R I C   T R A P E Z O I D A L   R I N G   E L E',  &
     ' M E N T S   (CTRAPAX)',/)
 896 FORMAT (3X,'SUBCASE   HARMONIC    POINT',12X,'RADIAL',12X,  &
     'CIRCUMFERENTIAL',12X,'AXIAL',16X,'CHARGE', /14X,  &
     'NUMBER     ANGLE',13X,'(R)',17X,'(THETA-T)',16X,'(Z)')
 897 FORMAT (' SUBCASE   HARMONIC    POINT    RADIAL      AXIAL     ',  &
     'CIRCUM.     SHEAR      SHEAR      SHEAR      F L U X   ',  &
     'D E N S I T I E S', /11X,'NUMBER      ANGLE     (R)',9X,  &
     '(Z)     (THETA-T)    (ZR)       (RT)       (ZT)',8X,  &
     '(R)        (Z)        (T)')
 898 FORMAT ('   FREQUENCY  HARMONIC    POINT            RADIAL',12X,  &
     'CIRCUMFERENTIAL',12X,'AXIAL',16X,'CHARGE', /14X,  &
     'NUMBER     ANGLE',13X,'(R)',17X,'(THETA-T)',16X,'(Z)')
 899 FORMAT (' FREQUENCY HARMONIC    POINT    RADIAL      AXIAL     ',  &
     'CIRCUM.     SHEAR      SHEAR      SHEAR      F L U X   ',  &
     'D E N S I T I E S', /10X,'NUMBER      ANGLE     (R)',9X,  &
     '(Z)     (THETA-T)    (ZR)       (RT)       (ZT)',8X,  &
     '(R)        (Z)        (T)')
 900 FORMAT (4X,'TIME     HARMONIC    POINT            RADIAL',12X,  &
     'CIRCUMFERENTIAL',12X,'AXIAL',16X,'CHARGE', /14X,  &
     'NUMBER     ANGLE',13X,'(R)',17X,'(THETA-T)',16X,'(Z)')
 901 FORMAT (2X,'TIME     HARMONIC    POINT    RADIAL      AXIAL     ',  &
     'CIRCUM.     SHEAR      SHEAR      SHEAR      F L U X   ',  &
     'D E N S I T I E S', /11X,'NUMBER      ANGLE     (R)',9X,  &
     '(Z)     (THETA-T)    (ZR)       (RT)       (ZT)',8X,  &
     '(R)        (Z)        (T)')
 902 FORMAT (5X,4HTIME,7X,8HHARMONIC,8X,2HT1,13X,2HT2,13X,2HT3,13X,  &
     2HR1,13X,2HR2,13X,2HR3)
 903 FORMAT (4X,7HSUBCASE,5X,8HHARMONIC,8X,2HT1,13X,2HT2,13X,2HT3,  &
     13X,2HR1,13X,2HR2,13X,2HR3)
 904 FORMAT (3X,9HFREQUENCY,4X,8HHARMONIC,8X,2HT1,13X,2HT2,13X,2HT3,  &
     13X,2HR1,13X,2HR2,13X,2HR3)
 905 FORMAT (19X,'F I N I T E   E L E M E N T   M A G N E T I C   F I',  &
     ' E L D   A N D   I N D U C T I O N',/)
 906 FORMAT (4X,'ELEMENT-ID   EL-TYPE         X-FIELD',10X,'Y-FIELD',  &
     10X,'Z-FIELD        X-INDUCTION      Y-INDUCTION',6X, 'Z-INDUCTION')
 907 FORMAT (28X,'G R I D   P O I N T   S T R E S S E S   F O R   I S',  &
     ' 2 D 8   E L E M E N T S',/)
 908 FORMAT (2X,7HELEMENT,3X,5HNO.of,4X,5HNO.of,7X,4HGRID,3X,6HCOORD.)
 909 FORMAT (4X,3HID.,4X,9HGRID pts.,1X,8HSTRESSES,3X,5HPOINT,2X,  &
     7HSYS id.,5X,5HSIG-x,8X,5HSIG-y,8X,6HTAU-xy)
 910 FORMAT (12X,5HNO.of,4X,5HNO.of,13X,6HCOORD.)
 911 FORMAT (4X,4HTIME,3X,9HGRID pts.,1X,8HSTRESSES,2X,7HGRID pt,1X,  &
     7HSYS id.,5X,5HSIG-x,8X,5HSIG-y,8X,6HTAU-xy)
 912 FORMAT (20X,'C O M P L E X   G R I D   P O I N T   S T R E S S E',  &
     ' S   F O R   I S 2 D 8   E L E M E N T S',/)
 913 FORMAT (2X,7HELEMENT,3X,5HNO.of,4X,5HNO.of,13X,6HCOORD.,/4X,3HID.,  &
     4X,9HGRID pts.,1X,8HSTRESSES,2X,7HGRID pt,1X,7HSYS id.,11X,  &
     5HSIG-x,22X,5HSIG-y,22X,6HTAU-xy)
 914 FORMAT (12X,5HNO.of,4X,5HNO.of,13X,6HCOORD., /1X,9HFREQUENCY,1X,  &
     9HGRID pts.,1X,8HSTRESSES,2X,7HGRID pt,1X,7HSYS id.,11X,  &
     5HSIG-x,22X,5HSIG-y,22X,6HTAU-xy)
 915 FORMAT (12X,5HNO.of,4X,5HNO.of,13X,6HCOORD., /2X,7HSUBCASE,2X,  &
     9HGRID pts.,1X,8HSTRESSES,2X,7HGRID pt,1X,7HSYS id.,5X,  &
     5HSIG-x,8X,5HSIG-y,8X,6HTAU-xy)
 916 FORMAT (26X,'F O R C E S    I N    C U R V E D    B E A M    E L',  &
     ' E M E N T S',8X,'( C E L B O W )',/)
 917 FORMAT (5X,7HELEMENT,11X,16H-bending moment-,21X,7H-shear-,18X,  &
     13H-axial force-,7X,8H-torque-)
 918 FORMAT (7X,3HID.,7X,13HPLANE-1 END-a,2X,13HPLANE-2 END-a,5X,  &
     13HPLANE-1 END-a,2X,13HPLANE-2 END-a,13X,5HEND-a,13X, 5HEND-a)
 919 FORMAT (25X,5HEND-b,10X,5HEND-b,13X,5HEND-b,28X,5HEND-b,13X, 5HEND-b)
 920 FORMAT (26X,'S T R E S S E S    I N    C U R V E D    B E A M   ',  &
     ' E L E M E N T S',8X,'( C E L B O W )',/)
 921 FORMAT (23X,16H-bending moment-,21X,7H-shear-,18X,  &
     13H-axial force-,7X,8H-torque-)
 922 FORMAT (4X,4HTIME,9X,13HPLANE-1 END-a,2X,13HPLANE-2 END-a,5X,  &
     13HPLANE-1 END-a,8X,7HPLANE-2,13X,5HEND-a,13X,5HEND-a,  &
     /25X,5HEND-b,10X,5HEND-b,13X,5HEND-b,28X,5HEND-b,13X, 5HEND-b)
 923 FORMAT (6X,7HSUBCASE,4X,13HPLANE-1 END-a,2X,13HPLANE-2 END-a,5X,  &
     13HPLANE-1 END-a,2X,13HPLANE-2 END-a,13X,5HEND-a,13X,  &
     5HEND-a, /25X,5HEND-b,10X,5HEND-b,13X,5HEND-b,28X,5HEND-b, 13X,5HEND-b)
 924 FORMAT (19X,'C O M P L E X   F O R C E S   I N   C U R V E D   ',  &
     'B E A M   E L E M E N T S   ( C E L B O W )',/)
 925 FORMAT (7X,9HFREQUENCY,21X,14HBENDING-moment,19X,11HSHEAR-force,  &
     22X,5HAXIAL,10X,6HTORQUE, /35X,7HPLANE 1,8X,7HPLANE 2,11X,  &
     7HPLANE 1,8X,7HPLANE 2,13X,5HFORCE)
 926 FORMAT (18X,'C O M P L E X   S T R E S S E S   I N   C U R V E D',  &
     '   B E A M   E L E M E N T S   ( C E L B O W )',/)
 927 FORMAT (7X,9HELEMENT  ,21X,14HBENDING-moment,19X,11HSHEAR-force,  &
     22X,5HAXIAL,10X,6HTORQUE, /35X,7HPLANE 1,8X,7HPLANE 2,11X,  &
     7HPLANE 1,8X,7HPLANE 2,13X,5HFORCE)
 928 FORMAT (23X,'F O R C E S   I N   F L U I D   H E X A H E D R A L',  &
     '   E L E M E N T S   ( C F H E X 2 )',/)
 929 FORMAT (23X,'F O R C E S   I N   F L U I D   H E X A H E D R A L',  &
     '   E L E M E N T S   ( C F H E X 1 )',/)
 930 FORMAT (19X,'F O R C E S   I N   F L U I D   T E T R A H E D R A',  &
     ' L   E L E M E N T S   ( C F T E T R A )',/)
 931 FORMAT (26X,'F O R C E S   I N   F L U I D   W E D G E   E L E M',  &
     ' E N T S    ( C F W E D G E )',/)
 932 FORMAT (24X,'P O W E R   C O N V E C T E D   B Y   F T U B E   ',  &
     'E L E M E N T S   ( C F T U B E )',/)
 933 FORMAT (47X,4HTIME,26X,5HPOWER)
 934 FORMAT (45X,10HELEMENT-id ,22X,5HPOWER)
 935 FORMAT (2X,7HELEMENT,3X,16HMAT. coord. sys.,30X,  &
     42H- stresses in material coordinate system -, /4X,  &
     3HID., 5X,16HID./output coded, 14X,8HNORMAL-x, 26X,  &
     8HNORAML-y, 25X,8HSHEAR-xy )
 936 FORMAT (16X,16HMAT. coord. sys., 30X,  &
     42H- stresses in material coordinate system -, /4X,  &
     9HFREQUENCY, 3X,15HID./output code,  &
     14X,8HNORAML-x, 26X,8HNORMAL-y, 25X,8HSHEAR-xy )
 937 FORMAT (50X, 29H(in stress coordinate system),/)
 938 FORMAT (2X,7HELEMENT,6X,5HFIBRE,15X,'STRESSES IN STRESS COORD. ',  &
     'SYSTEM',13X,31HPRINCIPAL stresses (zero shear),12X,3HMAX)
 939 FORMAT (4X,3HID.,7X,8HDISTANCE,11X,8HNORMAL-x,7X,8HNORMAL-y,6X,  &
     8HSHEAR-xy,7X,5HANGLE,9X,5HMAJOR,11X,5HMINOR,10X,5HSHEAR)
 940 FORMAT (20X,'F O R C E S   I N   G E N E R A L   Q U A D R I ',  &
     'L A T E R A L   E L E M E N T S     ( Q U A D 4 )',/)
 941 FORMAT (6X,'ELEMENT',12X,'- MEMBRANE  FORCES -',22X,'- BENDING',  &
     '   MOMENTS -',11X,'- TRANSVERSE SHEAR FORCES -')
 942 FORMAT (8X,'ID',10X,2HFX,12X,2HFY,12X,3HFXY,11X,2HMX,12X,2HMY,  &
     12X,3HMXY,11X,2HVX,12X,2HVY)
 943 FORMAT (19X,5HFIBRE,11X,32HSTRESSES in stress coord. system,13X,  &
     31HPRINCIPAL stresses (zero shear),10X,7HMAXIMUM, /7X,  &
     4HTIME,7X,8HDISTANCE,7X,8HNORMAL-x,7X,8HNORMAL-y,6X,  &
     8HSHEAR-xy,7X,5HANGLE,9X,5HMAJOR,11X,5HMINOR,10X,5HSHEAR)
 944 FORMAT (19X, 5HFIBRE, 11X, 32HSTRESSES in stress coord. system,  &
     13X, 31HPRINCIPAL stresses (zero shear), 10X, 7HMAXIMUM,  &
     /5X, 7HSUBCASE, 6X, 8HDISTANCE, 7X, 8HNORMAL-x, 7X,  &
     8HNORMAL-y, 6X, 8HSHEAR-xy, 7X, 5HANGLE, 9X, 5HMAJOR,  &
     11X, 5HMINOR, 10X, 5HSHEAR)
 945 FORMAT (6X,' TIME  ',18X,'- MEMBRANE  FORCES -',22X,'- BENDING',  &
     '   MOMENTS -',11X,'- TRANSVERSE SHEAR FORCES -')
 946 FORMAT (26X,2HFX,12X,2HFY,12X,3HFXY,11X,2HMX,12X,2HMY,12X,3HMXY,  &
     11X,2HVX,12X,2HVY)
 947 FORMAT (6X,'SUBCASE',18X,'- MEMBRANE  FORCES -',22X,'- BENDING',  &
     '   MOMENTS -',11X,'- TRANSVERSE SHEAR FORCES -')
 948 FORMAT (6X,'C O M P L E X   S T R E S S E S   I N   G E N E R A ',  &
     'L   Q U A D R I L I A T E R A L   E L E M E N T S   ', '( C Q U A D 4 )')
 949 FORMAT (9H  element,7X,5HFIBRE,38X,'- STRESSES IN STRESS COORDI',  &
     'NATE SYSTEM -', /4X,3HID.,8X,8HDISTANCE,18X,8HNORMAL-x,  &
     26X,8HNORMAL-y,25X,8HSHEAR-xy)
 950 FORMAT (20X,5HFIBRE,38X,'- STRESSES IN STRESS COORDINATE SYSTEM -'  &
     ,      /4X,9HFREQUENCY,6X,8HDISTANCE,18X,8HNORMAL-x,26X,  &
     8HNORMAL-y,25X,8HSHEAR-xy)
 951 FORMAT (6X,7HELEMENT,15X,6HCENTER,22X,7HEDGE  1,14X,7HEDGE  2,14X,  &
     7HEDGE  3,14X,7HEDGE  4, /8X,3HID.,9X,'R ------- PHI ----',  &
     '-- Z',4X,4(8X,13HS ------- phi))
 952 FORMAT (28X,6HCENTER,22X,7HEDGE  1,14X,7HEDGE  2,14X,7HEDGE  3,  &
     14X,7HEDGE 4, /7X,4HTIME,9X,22HR ------- phi ------ z,4X,  &
     4(8X,13HS ------- phi))
 953 FORMAT (29X,6HCENTER,21X,7HEDGE  1,14X,7HEDGE  2, 14X,7HEDGE  3,  &
     14X,7HEDGE  4,/4X,9HFREQUENCY,7X,22HR ------- phi ------ z,  &
     4X,4(8X,13HS ------- phi))
 954 FORMAT (9X,'C O M P L E X   S T R E S S E S   I N   T R I A N G ',  &
     'U L A R   M E M B R A N E   E L E M E N T S   ', '( C T R I M 6 )')
 955 FORMAT (11X,'C O M P L E X   F O R C E S   I N   T R I A N G U L',  &
     ' A R   M E M B R A N E   E L E M E N T S   ( C T R I M 6 )')
 956 FORMAT (9X,'C O M P L E X   S T R E S S E S   I N   T R I A N G ',  &
     'U L A R   B E N D I N G   E L E M E N T S   ', '( C T R P L T 1 )')
 957 FORMAT (11X,'C O M P L E X   F O R C E S   I N   T R I A N G U L',  &
     ' A R   B E N D I N G   E L E M E N T S   ( C T R P L T 1 )')
 958 FORMAT (12X,'C O M P L E X   S T R E S S E S   I N   T R I A N G',  &
     ' U L A R   S H E L L   E L E M E N T S   ( C T R S H L )')
 959 FORMAT (14X,'C O M P L E X   F O R C E S   I N   T R I A N G U L',  &
     ' A R   S H E L L   E L E M E N T S   ( C T R S H L )')
 960 FORMAT (9X,'C O M P L E X   F O R C E S   I N   G E N E R A L   ',  &
     'Q U A D R I L A T E R A L   E L E M E N T S   ', '( C Q U A D 4 )')
 961 FORMAT (3X,'FREQUENCY',14X,'- MEMBRANE  FORCES -',23X,'- BENDING',  &
     '   MOMENTS -',10X,'- TRANSVERSE SHEAR FORCES -',  &
     /22X,2HFX,12X,2HFY,11X,3HFXY,13X,2HMX,12X,2HMY,11X,3HMXY, 13X,2HVX,12X,2HVY)
 962 FORMAT (16X,4HGRID,11X,35HSTRESSES in basic coordinate system,13X,  &
     12HDIR. cosines, /3X,9HFREQUENCY,3X,5HPOINT,5X,8HNORMAL-x,  &
     9X,8HNORMAL-y,9X,8HNORMAL-z,9X,8HSHEAR-xy,9X,8HSHEAR-yz,9X, 8HSHEAR-zx)
 963 FORMAT (22X,'F O R C E S   I N   G E N E R A L   T R I A N G ',  &
     'U L A R   E L E M E N T S     ( C T R I A 3 )',/)
 964 FORMAT (12X,'C O M P L E X   F O R C E S   I N   G E N E R A L  ',  &
     ' T R I A N G U L A R   E L E M E N T S   ( C T R I A 3 )')
 965 FORMAT (21X,'S T R E S S E S   I N   G E N E R A L   T R I A N G',  &
     ' U L A R   E L E M E N T S',6X,'( C T R I A 3 )')
 966 FORMAT (9X,'C O M P L E X   S T R E S S E S   I N   G E N E R A ',  &
     'L   T R I A N G U L A R   E L E M E N T S   ', '( C T R I A 3 )')
 967 FORMAT (107X,22HOCTAHEDRAL    pressure, /6X,10H   subcase,8X,  &
     8HSIGMA-xx,6X,8HSIGMA-yy,6X,8HSIGMA-zz,7X,6HTAU-yz,8X,  &
     6HTAU-xz,8X,6HTAU-xy,8X,5HTAU-0,10X,1HP)
!WKBNB NCL93012 3/94
 969 FORMAT ( 4X, 'S T R A I N S / C U R V A T U R E S   I N   G E N '  &
     ,'E R A L   Q U A D R I L A T E R A L   E L E M E N T S ',6X  &
     ,'( Q U A D 4 )' )
 970 FORMAT ( 4X, 'S T R A I N S / C U R V A T U R E S   I N   G E N '  &
     ,'E R A L   T R I A N G U L A R   E L E M E N T S ',6X ,'( T R I A 3 )' )
!WKBNE NCL93012 3/94
!WKBNB SPR94001 7/94
 971 FORMAT (' SUBCASE',5X,'STRESS',15X,'RADIAL',16X  &
     ,       'CIRCUMFERENTIAL',16X,'AXIAL',21X,'SHEAR')
 972 FORMAT (5X,'NO ',6X,'POINT',17X,'(X)',21X,'(THETA)',21X,'(Z)',23X,  &
     '(ZX)')
 973 FORMAT (' SUBCASE   CORNER',18X,'RADIAL',26X,'CIRCUMFERENTIAL',  &
     26X,'AXIAL')
 974 FORMAT ('     NO     POINT',20X,'(X)',31X,'(THETA)',31X,'(Z)')
!WKBNE SPR94001 7/94
END SUBROUTINE ofp1c
