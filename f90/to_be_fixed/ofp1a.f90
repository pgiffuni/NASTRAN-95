SUBROUTINE ofp1a (line)
     
!     THIS ROUTINE WAS NAMED OFP1 BEFORE.
 
 
 INTEGER, INTENT(IN)                      :: line
 INTEGER :: l123(5),id(50),of(6)
 REAL :: fid(50),rt(8,15),sectn(2)
 COMMON /system/ ksys(80)
 COMMON /zzzzzz/ core(1)
 COMMON /ofp1id/ id22,m
 EQUIVALENCE     (ksys(2),l), (ksys(12),linet), (ksys(33),iflag),  &
     (ksys(3),nogo), (core(1),of(1),l123(1)), (fid(1) ,id(1), of(6))
 
 DATA   rt/ 4HSING, 4HULAR, 4HITIE, 4HS en, 4HCOUN, 4HTERE, 4HD.  , 4H    ,  &
     4H4 sh, 4HIFT , 4HPTS., 4HPER , 4HROOT, 4H exc, 4HEEDE, 4HD.  ,  &
     4HALL , 4HEIGE, 4HNVAL, 4HUES , 4HFOUN, 4HD in, 4H ran, 4HGE. ,  &
     4H3X e, 4HST.r, 4HOOTS, 4H in , 4HRANG, 4HE sp, 4HECIF, 4HIED.,  &
     4HNO m, 4HORE , 4HEIGE, 4HNVAL, 4HUES , 4HIN p, 4HROBL, 4HEM. ,  &
     4HNO. , 4HOF r, 4HOOTS, 4H des, 4HIRED, 4H wer, 4HE fo, 4HUND.,  &
     4H1 OR, 4H mor, 4HE ro, 4HOT o, 4HUTSI, 4HDE f, 4HR.ra, 4HNGE.,  &
     4HINSU, 4HFFIC, 4HIENT, 4H tim, 4HE fo, 4HR NE, 4HXT r, 4HOOT.,  &
     4HUNAB, 4HLE t, 4HO co, 4HNVER, 4HGE. , 4H    , 4H    , 4H    ,  &
     4HNORM, 4HAL t, 4HERMI, 4HNATI, 4HON  , 4H    , 4H    , 4H    ,  &
     4HEIGE, 4HNVAL, 4HUES , 4HOUTS, 4HIDE , 4HFREQ, 4H. ra, 4HNGE ,  &
     4HINSU, 4HFFIC, 4HIENT, 4H tim, 4HE re, 4HMAIN, 4HING , 4H    ,  &
     4HFEWE, 4HR th, 4HAN r, 4HEQUE, 4HSTED, 4H roo, 4HTS f, 4HOUND,  &
     4HROOT, 4HS fo, 4HUND , 4HWITH, 4H req, 4H. ac, 4HCURA, 4HCY  ,  &
     4HNO r, 4HOOTS ,4H fou, 4HND, , 4HNONE, 4H pas, 4HSED , 4HTEST/
 DATA     sectn / 4H.3.3, 4H.7.3  /     , twopi /6.283185307    /
 
 local = line - 100
 IF (local > 0) THEN
   GO TO  2004
 END IF
 2003 GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,  &
     23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  &
     44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,  &
     65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,  &
     86,87,88,89,90,91,92,93,94,95,96,97,98,99,100), line
 2004 GO TO (  101,102,103,104,105,106,107,108,109,110,111,112,113,114,  &
     115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,  &
     131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,  &
     147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,  &
     163,164,165,166,167,168,169,170,171,172,173,174), local
 
 1 WRITE (l,501)
 GO TO 1000
 2 WRITE (l,502)
 GO TO 1000
 3 WRITE (l,503)
 GO TO 1000
 4 WRITE (l,504)
 GO TO 1000
 5 WRITE (l,505)
 IF (id22 == -999) GO TO 1000
 IF (id22 < 0) THEN
   GO TO  2053
 ELSE IF (id22 == 0) THEN
   GO TO  2052
 END IF
 2051 WRITE (l,5050)
 WRITE (l,5051) id22
 GO TO 2054
 2052 IF (.NOT.((m >= 3 .AND. m <= 7) .OR. m == 10 .OR. m == 13)) GO TO 1000
 WRITE (l,5050)
 WRITE (l,5052)
 GO TO 2054
 2053 IF (m /= 8 .AND. m /= 12) GO TO 1000
 WRITE (l,5050)
 WRITE (l,5053)
 2054 WRITE (l,5054)
 linet = linet + 10
 id22  =-999
 GO TO 1000
 6 WRITE (l,506) id(5)
 GO TO 1000
 7 WRITE (l,507)
 GO TO 1000
 8 WRITE (l,508)
 GO TO 1000
 9 WRITE (l,509)
 GO TO 1000
 10 WRITE (l,510)
 GO TO 1000
 11 WRITE (l,511)
 GO TO 1000
 12 WRITE (l,512)
 GO TO 1000
 13 WRITE (l,513)
 GO TO 1000
 14 WRITE (l,514)
 GO TO 1000
 15 WRITE (l,515)
 GO TO 1000
 16 WRITE (l,516)
 GO TO 1000
 17 WRITE (l,517)
 GO TO 1000
 18 WRITE (l,518)
 GO TO 1000
 19 WRITE (l,519)
 GO TO 1000
 20 WRITE (l,520)
 GO TO 1000
 21 WRITE (l,521)
 GO TO 1000
 22 WRITE (l,522)
 GO TO 1000
 23 WRITE (l,523)
 GO TO 1000
 24 WRITE (l,524)
 GO TO 1000
 25 WRITE (l,525)
 GO TO 1000
 26 WRITE (l,526)
 GO TO 1000
 27 WRITE (l,527)
 GO TO 1000
 28 WRITE (l,528)
 GO TO 1000
 29 WRITE (l,529)
 GO TO 1000
 
!     PROCESS SPC AND MPC SET IDS PROPERLY TO ACCOUNT FOR AXISYMMETRIC
!     PROBLEMS
 
 30 DO  j = 3,4
   IF (id(j) <  100000000) CYCLE
   id(j) = id(j) - 100000000
   IF (id(j) <  100000000) CYCLE
   id(j) = id(j) - 100000000
 END DO
 WRITE (l,530) id(3),id(4)
 linet = linet + 1
 GO TO 1000
 31 WRITE (l,531)
 GO TO 1000
 32 WRITE (l,532)
 GO TO 1000
 33 WRITE (l,533)
 GO TO 1000
 34 WRITE (l,534)
 GO TO 1000
 35 WRITE (l,535)
 GO TO 1000
 36 WRITE (l,536)
 GO TO 1000
 37 WRITE (l,537)
 GO TO 1000
 38 WRITE (l,538)
 GO TO 1000
 39 WRITE (l,539)
 GO TO 1000
 40 WRITE (l,540)
 GO TO 1000
 41 WRITE (l,541)
 GO TO 1000
 42 WRITE (l,542)
 GO TO 1000
 43 WRITE (l,543)
 GO TO 1000
 44 WRITE (l,544)
 GO TO 1000
 45 WRITE (l,545)
 GO TO 1000
 46 WRITE (l,546)
 GO TO 1000
 47 WRITE (l,547)
 GO TO 1000
 48 WRITE (l,548)
 GO TO 1000
 49 WRITE (l,549)
 GO TO 1000
 50 WRITE (l,550)
 GO TO 1000
 51 WRITE (l,551)
 GO TO 1000
 52 WRITE (l,552)
 GO TO 1000
 53 WRITE (l,553)
 GO TO 1000
 54 WRITE (l,554)
 GO TO 1000
 55 WRITE (l,555)
 GO TO 1000
 56 WRITE (l,556)
 GO TO 1000
 57 WRITE (l,557)
 GO TO 1000
 58 WRITE (l,558)
 GO TO 1000
 59 WRITE (l,559)
 GO TO 1000
 60 WRITE (l,560)
 GO TO 1000
 61 WRITE (l,561)
 GO TO 1000
 62 WRITE (l,562)
 GO TO 1000
 63 WRITE (l,563)
 GO TO 1000
 64 WRITE (l,564)
 GO TO 1000
 65 WRITE (l,565)
 GO TO 1000
 66 WRITE (l,566)
 GO TO 1000
 67 WRITE (l,567)
 GO TO 1000
 68 WRITE (l,568)
 GO TO 1000
 69 WRITE (l,569)
 GO TO 1000
 70 WRITE (l,570)
 GO TO 1000
 71 WRITE (l,571)
 GO TO 1000
 72 WRITE (l,572)
 GO TO 1000
 73 WRITE (l,573)
 GO TO 1000
 74 WRITE (l,574)
 GO TO 1000
 75 WRITE (l,575)
 GO TO 1000
 76 WRITE (l,576)
 GO TO 1000
 77 WRITE (l,577)
 GO TO 1000
 78 WRITE (l,578)
 GO TO 1000
 79 WRITE (l,579)
 GO TO 1000
 80 WRITE (l,580)
 GO TO 1000
 81 WRITE (l,581)
 GO TO 1000
 82 WRITE (l,582)
 GO TO 1000
 83 WRITE (l,583)
 GO TO 1000
 84 WRITE (l,584)
 GO TO 1000
 85 WRITE (l,585)
 GO TO 1000
 86 WRITE (l,586)
 GO TO 1000
 87 WRITE (l,587)
 GO TO 1000
 88 WRITE (l,588)
 GO TO 1000
 89 WRITE (l,589)
 GO TO 1000
 90 IF (id(16) == 1) GO TO 905
 WRITE (l,590)
 GO TO 1000
 905 WRITE (l,5905)
 GO TO 1000
 91 WRITE (l,591) (id(k),k=11,15),id(17),fid(18),(id(k),k=19,21)
 m = id(17)
 IF (m >= 8) nogo = 14
 IF (id(16) /= 1) GO TO 911
 IF (m == 2 .OR. m > 3) nogo = 14
 IF (m == 0) m = 10
 IF (m == 1) m = 13
 IF (m == 3) m = 14
 911 IF (m > 0) WRITE (l,5911) (rt(k,m),k=1,8),sectn(1)
 id22 = id(22)
 GO TO 1000
 92 WRITE (l,592)
 GO TO 1000
 93 WRITE (l,593) (id(k),k=11,17),fid(18),(id(k),k=19,21)
 m = id(17)
 933 IF (m >= 3) nogo = 14
 IF (m == 1) m = 6
 IF (m == 2) m = 11
 IF (m == 3) m = 8
 IF (m == 4) m = 1
 GO TO 911
 94 WRITE (l,594)
 GO TO 1000
 95 WRITE (l,595)
 GO TO 1000
 96 WRITE (l,596)
 sectn(1) = sectn(2)
 GO TO 1000
 97 WRITE (l,597)
 GO TO 1000
 98 WRITE (l,598) (id(k),k=11,18)
 m = id(18)
 GO TO 933
 99 IF (id(3) == 4) GO TO 1007
 WRITE (l,599) (id(k),k=11,16)
 m = id(16)
 IF (id(17) /= 1) GO TO 911
 IF (m > 2) nogo = 14
 IF (m == 0) m = 10
 IF (m == 1) m = 13
 IF (m == 2) m = 15
 GO TO 911
 
!     ID(3)=2, ID(17)=0, METHOD IS COMPLEX INV
!     ID(3)=2, ID(17)=1, METHOD IS COMPLEX FEER
!     ID(3)=4, ID(17)=0, METHOD IS COMPLEX HESS
 
 100 sectn(1) = sectn(2)
 IF (id(17) == 1) GO TO 1005
 IF (id( 3) == 4) GO TO 1006
 WRITE (l,600)
 GO TO 1000
 1005 WRITE (l,6005)
 GO TO 1000
 1006 WRITE (l,6006)
 GO TO 1000
 1007 WRITE (l,6007) id(11),id(12),id(18)
 m = id(18)
 GO TO 933
 101 f = SQRT(ABS(fid(6)))/twopi
 WRITE (l,601) fid(6),f
 linet = linet + 1
 GO TO 1000
 102 f = SQRT(ABS(fid(6)))/twopi
 WRITE (l,602) fid(6),f
 linet = linet + 1
 GO TO 1000
 103 WRITE (l,603) fid(5)
 GO TO 1000
 104 CONTINUE
 WRITE (l,604) fid(5)
 GO TO 1000
 105 WRITE (l,605) fid(5)
 GO TO 1000
 106 CONTINUE
 WRITE (l,606) fid(5)
 GO TO 1000
 107 WRITE (l,607) id(5)
 GO TO 1000
 108 WRITE (l,608) id(5)
 GO TO 1000
 109 cycfrq = ABS(fid(7)) / twopi
 WRITE (l,609) fid(6),fid(7),cycfrq
 GO TO 1000
 110 CONTINUE
 cycfrq = ABS(fid(7)) / twopi
 WRITE (l,610) fid(6),fid(7),cycfrq
 GO TO 1000
 111 WRITE (l,611)
 GO TO 1000
 112 WRITE (l,612)
 GO TO 1000
 113 WRITE (l,613)
 GO TO 1000
 114 WRITE (l,614)
 GO TO 1000
 115 WRITE (l,615)
 GO TO 1000
 116 WRITE (l,616)
 GO TO 1000
 117 WRITE (l,617)
 GO TO 1000
 118 WRITE (l,618)
 GO TO 1000
 119 WRITE (l,619)
 GO TO 1000
 120 WRITE (l,620)
 GO TO 1000
 121 WRITE (l,621)
 GO TO 1000
 122 WRITE (l,622)
 GO TO 1000
 123 WRITE (l,623)
 GO TO 1000
 124 WRITE (l,624) id(5)
 GO TO 1000
 125 WRITE (l,625)
 GO TO 1000
 126 WRITE (l,626)
 GO TO 1000
 127 WRITE (l,627)
 GO TO 1000
 128 WRITE (l,628)
 GO TO 1000
 129 WRITE (l,629)
 GO TO 1000
 130 WRITE (l,630)
 GO TO 1000
 131 WRITE (l,631)
 GO TO 1000
 132 WRITE (l,632)
 GO TO 1000
 133 WRITE (l,633)
 GO TO 1000
 134 WRITE (l,634)
 GO TO 1000
 135 WRITE (l,635)
 GO TO 1000
 136 WRITE (l,636)
 GO TO 1000
 137 WRITE (l,637)
 GO TO 1000
 138 WRITE (l,638)
 GO TO 1000
 139 WRITE (l,639)
 GO TO 1000
 140 WRITE (l,640)
 GO TO 1000
 141 CONTINUE
 GO TO 1000
 142 WRITE (l,642)
 GO TO 1000
 143 WRITE (l,643)
 GO TO 1000
 144 WRITE (l,644)
 GO TO 1000
 145 WRITE (l,645)
 GO TO 1000
 146 WRITE (l,646)
 GO TO 1000
 147 CONTINUE
 GO TO 1000
 148 WRITE (l,648)
 GO TO 1000
 149 WRITE (l,649)
 GO TO 1000
 150 WRITE (l,650)
 GO TO 1000
 151 WRITE (l,651)
 GO TO 1000
 152 WRITE (l,652)
 GO TO 1000
 153 WRITE (l,653)
 GO TO 1000
 154 WRITE (l,654)
 GO TO 1000
 155 WRITE (l,655)
 GO TO 1000
 156 WRITE (l,656)
 GO TO 1000
 157 WRITE (l,657)
 GO TO 1000
 158 WRITE (l,658)
 GO TO 1000
 159 WRITE (l,659)
 GO TO 1000
 160 WRITE (l,660)
 GO TO 1000
 161 WRITE (l,661)
 GO TO 1000
 162 WRITE (l,662)
 GO TO 1000
 163 WRITE (l,663)
 GO TO 1000
 164 WRITE (l,664)
 GO TO 1000
 165 WRITE (l,665)
 GO TO 1000
 166 WRITE (l,666)
 GO TO 1000
 167 WRITE (l,667)
 GO TO 1000
 168 WRITE (l,668)
 GO TO 1000
 169 WRITE (l,669)
 GO TO 1000
 170 WRITE (l,670)
 GO TO 1000
 171 WRITE (l,671)
 GO TO 1000
 172 WRITE (l,672)
 GO TO 1000
 173 WRITE (l,673)
 GO TO 1000
 174 WRITE (l,674)
 1000 CONTINUE
 RETURN
 
 
 501 FORMAT (45X,37HD i s p l a c e m e n t   v e c t o r)
 502 FORMAT (6X,16HPOINT id.   TYPE,10X,2HT1,13X,2HT2,13X,2HT3,13X,  &
     2HR1,13X,2HR2,13X,2HR3)
 503 FORMAT (46X, 31HR e a l   e i g e n v a l u e s,/)
 504 FORMAT (3X, 4HMODE, 4X, 10HEXTRACTION, 7X, 10HEIGENVALUE, 12X,  &
     6HRADIAN, 14X, 6HCYCLIC, 2X, 2(9X, 11HGENERALIZED))
 505 FORMAT (4X,3HNO., 7X, 5HORDER, 30X, 9HFREQUENCY, 11X, 9HFREQUENCY,  &
     12X,4HMASS, 14X, 9HSTIFFNESS,/)
 5050 FORMAT (/37X,16(4H****), /37X,1H*,62X,1H*, /37X,1H*)
 5051 FORMAT (1H+,45X,'NASTRAN INFORMATION MESSAGE 3307, POTENTIALLY',  &
9X,1H*, /37X,1H*,i10,' EIGENVALUE(S) AT LOW FREQ. END NOT',  &
    ' FOUND',11X,1H*)
5052 FORMAT (1H+,39X,'NASTRAN INFORMATION MESSAGE 3308, LOWEST EIGEN',  &
    'VALUE FOUND',3X,1H*, /37X,1H*,2X,'AS INDICATED BY THE ',  &
    'STURM''S SEQUENCE OF THE DYNAMIC MATRIX',2X,1H*)
5053 FORMAT (1H+,42X,'NASTRAN WARNING MESSAGE 3309, ALL LOWER EIGEN',  &
    'VALUES',6X,1H*, /37X,1H*, 5X,'NOT NECESSARY FOUND.',37X, 1H*)
5054 FORMAT (37X,1H*,62X,1H*, /37X,1H*,8X,  &
    43H(this message can be suppressed by diag 37),11X,1H*, /37X,16(4H****),/)
506 FORMAT (41X,39HR e a l   e i g e n v e c t o r   n o .,i11)
507 FORMAT (7X,7HELEMENT,11X,5HAXIAL,37X,7HELEMENT,11X,5HAXIAL)
508 FORMAT (9X,3HID.,13X,5HFORCE,10X,6HTORQUE,23X,3HID.,13X,5HFORCE,  &
    10X,6HTORQUE)
509 FORMAT (12H0    element,9X,17HBEND-moment END-a,12X,17HBEND-moment  &
    END-b,16X,9H- shear -,15X,5HAXIAL)
510 FORMAT (4X,9H   id.   ,3(6X,7HPLANE 1,7X,7HPLANE 2,2X),7X,5HFORCE,  &
    9X,6HTORQUE)
511 FORMAT (7X,7HELEMENT,11X,5HFORCE,10X,5HFORCE,22X,7HELEMENT,11X,  &
    5HFORCE,10X,5HFORCE)
512 FORMAT (9X,3HID.,12X,7HPTS 1,3,8X,7HPTS 2,4,23X,3HID.,12X,  &
    7HPTS 1,3,8X,7HPTS 2,4)
513 FORMAT (7X,7HELEMENT,10X,6HMOMENT,9X,6HMOMENT,22X,7HELEMENT,10X,  &
    6HMOMENT,9X,6HMOMENT)
514 FORMAT (1H0,8X,7HELEMENT,2(11X,11HBEND-moment),10X,12HTWIST-moment  &
    ,      2(13X,5HSHEAR,4X))
515 FORMAT (11X,3HID.,17X,1HX,21X,1HY,43X,1HX,21X,1HY)
516 FORMAT (6X,3(7HELEMENT,9X,5HFORCE,12X),7HELEMENT,9X,5HFORCE)
517 FORMAT (8X,3(3HID.,30X),3HID.)
518 FORMAT (2(7X,7HELEMENT,7X,5HAXIAL,7X,6HSAFETY,6X,9HTORSIONAL,5X,  &
    6HSAFETY))
519 FORMAT (2(9X,3HID.,8X,6HSTRESS,7X,6HMARGIN,8X,6HSTRESS,6X, 6HMARGIN))
520 FORMAT (2X,7HELEMENT,8X,3HSA1,12X,3HSA2,12X,3HSA3,15X,1HS,14X,  &
    6HSA-MAX,9X,6HSA-MIN,11X,6HM.s.-t)
521 FORMAT (4X,3HID.,10X,3HSB1,12X,3HSB2,12X,3HSB3,30X,6HSB-MAX,9X,  &
    6HSB-MIN,11X,6HM.s.-c)
522 FORMAT (2(9X,7HELEMENT,12X,3HMAX,12X,3HAVG,8X,6HSAFETY))
523 FORMAT (2(11X,3HID.,13X,5HSHEAR,10X,5HSHEAR,7X,6HMARGIN))
524 FORMAT (2(11X,3HID.,40X,6HMARGIN))
525 FORMAT (2X,7HELEMENT,11X,32HSTRESSES in element coord system,12X,  &
    9HPRINCIPAL,11X,18HPRINCIPAL stresses,12X,3HMAX)
526 FORMAT (4X,3HID.,11X,8HNORMAL-x,7X,8HNORMAL-y,7X,8HSHEAR-xy,6X,  &
    12HSTRESS angle,9X,5HMAJOR,10X,5HMINOR,10X,5HSHEAR)
527 FORMAT (2X,7HELEMENT,6X,5HFIBRE,15X,'STRESSES IN ELEMENT COORD ',  &
    'SYSTEM',13X,'PRINCIPAL STRESSES (ZERO SHEAR)',12X,3HMAX)
528 FORMAT (4X,3HID.,7X,8HDISTANCE,11X,8HNORMAL-x,7X,8HNORMAL-y,6X,  &
    8HSHEAR-xy,7X,5HANGLE,9X,5HMAJOR,11X,5HMINOR,10X,5HSHEAR)
529 FORMAT (6X,3(7HELEMENT,9X,6HSTRESS,11X),7HELEMENT,9X,6HSTRESS)
530 FORMAT (30X,'G R I D   P O I N T   S I N G U L A R I T Y   ',  &
    'T A B L E',6X,3HSPC,i9,3X,3HMPC,i9)
531 FORMAT (8X,5HPOINT,10X,11HSINGULARITY,18X,'LIST OF COORDINATE ',  &
    'COMBINATIONS THAT WILL REMOVE SINGULARITY')
532 FORMAT (9X,3HID.,3X,4HTYPE,7X,5HORDER,7X,21HSTRONGEST combination,  &
    15X,18HWEAKER combination,17X,19HWEAKEST combination)
533 FORMAT (53X,21HL o a d   v e c t o r)
534 FORMAT (2X,7HELEMENT,8X,3HSA1,12X,3HSA2,12X,3HSA3,12X,3HSA4,11X,  &
    5HAXIAL,10X,6HSA-MAX,9X,17HSA-MIN     m.s.-t)
535 FORMAT (4X,3HID.,10X,3HSB1,12X,3HSB2,12X,3HSB3,12X,3HSB4,11X,  &
    6HSTRESS,9X,6HSB-MAX,9X,17HSB-MIN     m.s.-c)
536 FORMAT (43X,'F O R C E S   I N   R O D   E L E M E N T S',5X,  &
    '( C R O D )')
537 FORMAT (33X,'F O R C E S   I N   B E A M   E L E M E N T S',8X,  &
    '( C B E A M )')
538 FORMAT (27X,'F O R C E S   A C T I N G   O N   S H E A R   ',  &
    'P A N E L   E L E M E N T S   ( C S H E A R )')
539 FORMAT (37X,'F O R C E S   I N   T W I S T   P A N E L S',6X,  &
    '( C T W I S T )')
540 FORMAT (21X,'F O R C E S   I N   B A S I C   B E N D I N G   ',  &
    'T R I A N G L E S',7X,'( C T R B S C )')
541 FORMAT (30X,'F O R C E S   I N   S C A L A R   S P R I N G S',8X,  &
    '( C E L A S 1 )')
542 FORMAT (30X,'F O R C E S   I N   S C A L A R   S P R I N G S',8X,  &
    '( C E L A S 2 )')
543 FORMAT (30X,'F O R C E S   I N   S C A L A R   S P R I N G S',8X,  &
    '( C E L A S 3 )')
544 FORMAT (30X,'F O R C E S   I N   S C A L A R   S P R I N G S',8X,  &
    '( C E L A S 4 )')
545 FORMAT (31X,'F O R C E S   O F   S I N G L E - P O I N T   ',  &
    'C O N S T R A I N T')
546 FORMAT (43X,'F O R C E S   I N   R O D   E L E M E N T S',5X,  &
    '( C O N R O D )')
547 FORMAT (33X,'F O R C E S   I N   B A R   E L E M E N T S',9X,  &
    '( C B A R )')
548 FORMAT (17X,'F O R C E S   I N   B E N D I N G   Q U A D R I L A',  &
    ' T E R A L S',9X,'( C Q D P L T )')
549 FORMAT (17X,'F O R C E S   I N   G E N E R A L   Q U A D R I L A',  &
    ' T E R A L   E L E M E N T S     ( C Q U A D 1 )')
550 FORMAT (17X,'F O R C E S   I N   G E N E R A L   Q U A D R I L A',  &
    'T E R A L   E L E M E N T S     ( C Q U A D 2 )')
551 FORMAT (21X,'F O R C E S   I N   G E N E R A L   T R I A N G U L',  &
    ' A R   E L E M E N T S',8X,'( C T R I A 1 )')
552 FORMAT (21X,'F O R C E S   I N   G E N E R A L   T R I A N G U L',  &
    ' A R   E L E M E N T S',8X,'( C T R I A 2 )')
553 FORMAT (27X,'F O R C E S   I N   B E N D I N G   T R I A N G L E',  &
    ' S       ( C T R P L T )')
554 FORMAT (33X,'F O R C E S   I N   R O D   E L E M E N T S     ',  &
    '( C T U B E )')
555 FORMAT (37X,'S T R E S S E S   I N   R O D   E L E M E N T S',6X,  &
    '( C R O D )')
556 FORMAT (34X,'S T R E S S E S   I N   B E A M   E L E M E N T S',  &
    8X,'( C B E A M )')
557 FORMAT (40X,'S T R E S S E S   I N   S H E A R   P A N E L S',6X,  &
    '( C S H E A R )')
558 FORMAT (40X,'S T R E S S E S   I N   T W I S T   P A N E L S',7X,  &
    '( C T W I S T )')
559 FORMAT (22X,'S T R E S S E S   I N   T R I A N G U L A R   ',  &
    'M E M B R A N E S      ( C T R M E M )')
560 FORMAT (19X,'S T R E S S E S   I N   B A S I C   B E N D I N G  ',  &
    ' T R I A N G L E S',8X,'( C T R B S C )')
561 FORMAT (30X,'S T R E S S E S   I N   S C A L A R   S P R I N G S',  &
    8X,'( C E L A S 1 )')
562 FORMAT (30X,'S T R E S S E S   I N   S C A L A R   S P R I N G S',  &
    8X,'( C E L A S 2 )')
563 FORMAT (30X,'S T R E S S E S   I N   S C A L A R   S P R I N G S',  &
    8X,'( C E L A S 3 )')
564 FORMAT (33X,'S T R E S S E S   I N   B A R   E L E M E N T S',10X,  &
    '( C B A R )')
565 FORMAT (37X,'S T R E S S E S   I N   R O D   E L E M E N T S',6X,  &
    '( C O N R O D )')
566 FORMAT (21X,'S T R E S S E S   I N   Q U A D R I L A T E R A L  ',  &
    ' M E M B R A N E S      ( C Q D M E M )')
567 FORMAT (18X,'S T R E S S E S   I N   B E N D I N G   Q U A D R I',  &
    ' L A T E R A L S',13X,'( C Q D P L T )')
568 FORMAT (18X,'S T R E S S E S   I N   G E N E R A L   Q U A D R I',  &
    ' L A T E R A L   E L E M E N T S',6X,'( C Q U A D 1 )')
569 FORMAT (18X,'S T R E S S E S   I N   G E N E R A L   Q U A D R I',  &
    ' L A T E R A L   E L E M E N T S',6X,'( C Q U A D 2 )')
570 FORMAT (18X,'S T R E S S E S   I N   G E N E R A L   T R I A N G',  &
    ' U L A R   E L E M E N T S',7X,'( C T R I A 1 )')
571 FORMAT (18X,'S T R E S S E S   I N   G E N E R A L   T R I A N G',  &
    ' U L A R   E L E M E N T S',7X,'( C T R I A 2 )')
572 FORMAT (24X,'S T R E S S E S   I N   B E N D I N G   T R I A N G',  &
    ' L E S',8X,'( C T R P L T )')
573 FORMAT (36X,'S T R E S S E S   I N   R O D   E L E M E N T S',6X,  &
    '( C T U B E )')
574 FORMAT (20X,'S T R E S S E S   F O R   T H E   T R I A N G U L A',  &
    ' R   R I N G S',5X,'( C T R I A R G )')
575 FORMAT (5X,3HEL ,13X,6HRADIAL,20X,15HCIRCUMFERENTIAL,20X,5HAXIAL,  &
    25X,5HSHEAR)
576 FORMAT (5X,3HID ,15X,3H(x),25X,7H(theta),25X,3H(z),27X,4H(zx))
577 FORMAT (18X,'S T R E S S E S   F O R   T H E   T R A P E Z O I D',  &
    ' A L   R I N G S',5X,'( C T R A P R G )')
578 FORMAT (5X,3HEL ,5X,6HSTRESS,15X,6HRADIAL,16X,15HCIRCUMFERENTIAL,  &
    16X,5HAXIAL,21X,5HSHEAR)
579 FORMAT (5X,3HID ,6X,5HPOINT,17X,3H(x),21X,7H(theta),21X,3H(z),23X, 4H(zx))
580 FORMAT (11X,'S T R E S S   R E S U L T A N T S   F O R   T H E  ',  &
    ' T O R O I D A L   R I N G S     ( C T O R D R G )')
581 FORMAT (5X, 3HEL , 8H  stress, 15X, 17HMEMBRANE (forces), 26X,  &
    17HFLEXURE (moments), 23X, 5HSHEAR)
582 FORMAT (5X,2HID,9H    point,8X,10HTANGENTIAL,10X,'CIRCUMFERENTIAL'  &
    ,      8X,10HTANGENTIAL,11X,15HCIRCUMFERENTIAL,10X,7H(force))
583 FORMAT (22X,'F O R C E S   F O R   T H E   T R I A N G U L A R  ',  &
    ' R I N G S     ( C T R I A R G )')
584 FORMAT (5X,12HEL    corner,18X,6HRADIAL,26X,15HCIRCUMFERENTIAL,  &
    26X,5HAXIAL)
585 FORMAT (5X,12HID     point,20X,3H(x),31X,7H(theta),31X,3H(z))
586 FORMAT (21X,'F O R C E S   F O R   T H E   T R A P E Z O I D A L',  &
    '   R I N G S     ( C T R A P R G )')
587 FORMAT (23X,'F O R C E S   F O R   T H E   T O R O I D A L   ',  &
    'R I N G S     ( C T O R D R G )')
588 FORMAT (5X,12HEL    corner,9X,6HRADIAL,8X,15HCIRCUMFERENTIAL,7X,  &
    5HAXIAL,13X,6HMOMENT,9X,13HDIRECT strain,7X,9HCURVATURE)
589 FORMAT (5X,12HID     point,11X,3H(x),13X,7H(theta),12X,3H(z),15X,  &
    4H(zx),14X,4H(xi),13X,7H(xi,xi))
590 FORMAT (30X,'E I G E N V A L U E   A N A L Y S I S   S U M M A R',  &
    ' Y     (INVERSE POWER METHOD)')
5905 FORMAT (30X,'E I G E N V A L U E   A N A L Y S I S   S U M M A R',  &
    ' Y',9X,'(FEER METHOD)')
591 FORMAT (1H0, /1H0,39X,32HNUMBER of eigenvalues extracted ,6(2H .),  &
    i10,/1H0,39X,30HNUMBER of starting points used,7(2H .),i10,  &
    /1H0,39X,30HNUMBER of starting point moves,7(2H .),i10,  &
    /1H0,39X,36HNUMBER of triangular decompositions ,4(2H .),  &
    i10,/1H0,39X,34HTOTAL NUMBER of vector iterations ,5(2H .),  &
    i10,//1H0,39X,22HREASON for termination,11(2H .),i10,1H*,/,  &
    /1H0,39X,36HLARGEST off-diagonal modal mass term,4(2H .),  &
    e10.2, /1H0,77X,3(2H .),i10, /50X,9HMODE pair ,10(2H .),  &
    /78X,3(2H .),i10, /1H0,39X,'NUMBER OF OFF-DIAGONAL MODAL ',  &
    'MASS', /45X,23HTERMS failing criterion,8(2H .),i10)
5911 FORMAT (/1H0,39X,3H(* ,8A4, /41X,'SEE NASTRAN U.M. VOL II, ',  &
    'SECTION 2',a4,1H))
592 FORMAT (26X,'E I G E N V A L U E   A N A L Y S I S   S U M M A R',  &
    ' Y       (DETERMINANT METHOD)')
593 FORMAT (1H0, /1H0,39X,32HNUMBER of eigenvalues extracted ,6(2H .),  &
    i9,/1H0,39X,44HNUMBER of passes through starting points . .  &
    ,      i9,/1H0,39X,26HNUMBER of criteria changes,9(2H .),i9,  &
    /1H0,39X,30HNUMBER of starting point moves,7(2H .),i9,  &
    /1H0,39X,36HNUMBER of triangular decompositions ,4(2H .),  &
    i9,/1H0,39X,44HNUMBER of failures TO iterate TO a root  . .  &
    ,      i9, //1H0,39X,22HREASON for termination,11(2H .),i9,1H*,  &
    //,1H0,39X,36HLARGEST off-diagonal modal mass term,4(2H .),  &
    e9.2,/1H0,77X,3(2H .),i9,/50X, 9HMODE pair ,10(2H .), /78X,  &
    3(2H .),i9, /1H0,39X,33HNUMBER of off-diagonal modal mass,  &
    /45X,23HTERMS failing criterion,8(2H .),i9)
594 FORMAT (10X,14HSTARTING point,6X,6HLAMBDA,9X,'RADIAN FREQUENCY  ',  &
    '  CYCLIC FREQUENCY    DETERMINANT',9X,'SCALE FACTOR',/)
595 FORMAT (1H0,40X,'S W E P T   D E T E R M I N A N T   F U N C T I',  &
    ' O N',/)
596 FORMAT (20X,'C O M P L E X   E I G E N V A L U E   A N A L Y S I',  &
    ' S   S U M M A R Y     (DETERMINANT METHOD)')
597 FORMAT (42X,5H- p -,35X,10H- det(p) -, /10X,14HSTARTING point,10X,  &
    4HREAL,13X,4HIMAG,20X,9HMAGNITUDE,9X,5HPHASE,5X, 12HSCALE factor)
598 FORMAT (1H0, /1H0,39X,32HNUMBER of eigenvalues extracted ,6(2H .),  &
    i9,/1H0,39X,44HNUMBER of passes through starting points . .  &
    ,      i9,/1H0,39X,26HNUMBER of criteria changes,9(2H .),i9,  &
    /1H0,39X,30HNUMBER of starting point moves,7(2H .),i9,  &
    /1H0,39X,36HNUMBER of triangular decompositions ,4(2H .),  &
    i9,/1H0,39X,44HNUMBER of failures TO iterate TO a root  . .  &
    ,      i9,/1H0,39X,36HNUMBER of predictions outside region,4(2H .)  &
    ,      i9,/1H0,/1H0,39X,22HREASON for termination,11(2H .),i9,1H*)
599 FORMAT (1H0, /1H0, /1H0,35X,32HNUMBER of eigenvalues extracted ,  &
    9(2H .),i9, /1H0,35X,30HNUMBER of starting points used,  &
    10(2H .),i9, /1H0,35X,  &
    50HNUMBER of starting point OR shift point moves  . .,i9,  &
    /1H0,35X,42HTOTAL NUMBER of triangular decompositions ,  &
    4(2H .),i9, /1H0,35X,34HTOTAL NUMBER of vector iterations ,  &
    8(2H .),i9, /1H0, /1H0,35X,22HREASON for termination, 14(2H .),i9,1H*)
600 FORMAT (19X,'C O M P L E X   E I G E N V A L U E   A N A L Y S I',  &
    ' S   S U M M A R Y   (INVERSE POWER METHOD)')
6005 FORMAT (23X,'C O M P L E X   E I G E N V A L U E   A N A L Y S I',  &
    ' S   S U M M A R Y     (FEER METHOD)')
6006 FORMAT (20X,'C O M P L E X   E I G E N V A L U E   A N A L Y S I',  &
    ' S   S U M M A R Y     (HESSENBERG METHOD)')
6007 FORMAT (1H0, /1H0, /,1H0,35X,32HNUMBER of eigenvalues extracted ,  &
    9(2H .),i9, /,1H0,35X,30HNUMBER of eigenvalues desired ,  &
    10(2H .),i9, /,1H0,35X,22HREASON for termination,14(2H .), i9,1H*)
601 FORMAT (6X, 'EIGENVALUE =',e14.6,  &
    4X, '(CYCLIC FREQUENCY =', e14.6, ' HZ)'/)
602 FORMAT (6X, 'EIGENVALUE =',e14.6,  &
    4X, '(CYCLIC FREQUENCY =', 1P,e14.6, ' HZ)'/)
603 FORMAT (6X,11HFREQUENCY =,e14.6)
604 FORMAT (6X,11HFREQUENCY =,1P,e14.6)
605 FORMAT (6X,6HTIME =,e14.6)
606 FORMAT (6X,6HTIME =,1P,e14.6)
607 FORMAT (6X,10HPOINT-id =,i8)
608 FORMAT (6X,12HELEMENT-id =,i8,/)
609 FORMAT (6X,20HCOMPLEX eigenvalue =,e14.6,1H,,e14.6,  &
    4X, '(CYCLIC FREQUENCY =', e14.6, 'HZ)')
610 FORMAT (6X,20HCOMPLEX eigenvalue =,1P,e14.6,1H,,1P,e14.6,  &
    4X, '(CYCLIC FREQUENCY =', 1P,e14.6, 'HZ)')
611 FORMAT (6X,16HFREQUENCY   TYPE,10X,2HT1,13X,2HT2,13X,2HT3,13X,  &
    2HR1,13X,2HR2,13X,2HR3)
612 FORMAT (6X,16H time       TYPE,10X,2HT1,13X,2HT2,13X,2HT3,13X,  &
    2HR1,13X,2HR2,13X,2HR3)
613 FORMAT (48X,30HV e l o c i t y    v e c t o r )
614 FORMAT (44X,38HA c c e l e r a t i o n    v e c t o r )
615 FORMAT (41X,45HN o n - l i n e a r - f o r c e   v e c t o r )
616 FORMAT (40X,'C O M P L E X   E I G E N V A L U E   S U M M A R Y')
617 FORMAT (1H0,16X,19HROOT     extraction,18X,10HEIGENVALUE,21X,  &
    9HFREQUENCY,14X,7HDAMPING)
618 FORMAT (18X,3HNO.,8X,5HORDER,13X,6H(REAL),11X,6H(imag),16X,  &
    8H(cycles),12X,11HCOEFFICIENT)
619 FORMAT (39X,'C O M P L E X   D I S P L A C E M E N T   V E C T O R  &
    '      )
620 FORMAT (43X,'C O M P L E X   V E L O C I T Y   V E C T O R')
621 FORMAT (39X,'C O M P L E X   A C C E L E R A T I O N   V E C T O R  &
    '      )
622 FORMAT (25X,'C O M P L E X   F O R C E S   O F   S I N G L E   ',  &
    'P O I N T   C O N S T R A I N T')
623 FORMAT (47X,'C O M P L E X   L O A D   V E C T O R')
624 FORMAT (39X,'C O M P L E X   E I G E N V E C T O R   NO.',i11)
625 FORMAT (58X,'(REAL/IMAGINARY)')
626 FORMAT (57X,'(MAGNITUDE/PHASE)')
627 FORMAT (27X,'C O M P L E X   S T R E S S E S   I N   B A R   E L',  &
    ' E M E N T S   ( C B A R )')
628 FORMAT (23X,'C O M P L E X   S T R E S S E S   I N   S C A L A R',  &
    '   S P R I N G S   ( C E L A S 1 )')
629 FORMAT (23X,'C O M P L E X   S T R E S S E S   I N   S C A L A R',  &
    '   S P R I N G S   ( C E L A S 2 )')
630 FORMAT (23X,'C O M P L E X   S T R E S S E S   I N   S C A L A R',  &
    '   S P R I N G S   ( C E L A S 3 )')
631 FORMAT (25X,'C O M P L E X   S T R E S S E S   I N   R O D   E L',  &
    ' E M E N T S   ( C O N R O D )')
632 FORMAT (14X,'C O M P L E X   S T R E S S E S   I N   Q U A D R I',  &
    ' L A T E R A L   M E M B R A N E S   ( C Q D M E M )')
633 FORMAT (16X,'C O M P L E X   S T R E S S E S   I N   B E N D I N',  &
    ' G   Q U A D R I L A T E R A L S   ( C Q D P L T )')
634 FORMAT (6X,'C O M P L E X   S T R E S S E S   I N   G E N E R A L'  &
    ,      '   Q U A D R I L I A T E R A L   E L E M E N T S   ',  &
    '( C Q U A D 1)')
635 FORMAT (6X,'C O M P L E X   S T R E S S E S   I N   G E N E R A L'  &
    ,      '   Q U A D R I L I A T E R A L   E L E M E N T S   ',  &
    '( C Q U A D 2 )')
636 FORMAT (27X,'C O M P L E X   S T R E S S E S   I N   R O D   E L',  &
    ' E M E N T S   ( C R O D )')
637 FORMAT (25X,'C O M P L E X   S T R E S S E S   I N   S H E A R  ',  &
    ' P A N E L S   ( C S H E A R )')
638 FORMAT (14X,'C O M P L E X   S T R E S S E S   I N   B A S I C  ',  &
    ' B E N D I N G   T R I A N G L E S   ( C T R B S C )')
639 FORMAT (10X,'C O M P L E X   S T R E S S E S   I N   G E N E R A',  &
    ' L   T R I A N G U L A R   E L E M E N T S   ', '( C T R I A 1 )')
640 FORMAT (11X,'C O M P L E X   S T R E S S E S   I N   G E N E R A',  &
    ' L  T R I A N G U L A R   E L E M E N T S   ', '( C T R I A 2 )')
642 FORMAT (17X,'C O M P L E X   S T R E S S E S   I N   T R I A N G',  &
    ' U L A R   M E M B R A N E S   ( C T R M E M )')
643 FORMAT (20X,'C O M P L E X   S T R E S S E S   I N   B E N D I N',  &
    ' G   T R I A N G L E S   ( C T R P L T )')
644 FORMAT (26X,'C O M P L E X   S T R E S S E S   I N   R O D   ',  &
    'E L E M E N T S   ( C T U B E )')
645 FORMAT (25X,'C O M P L E X   S T R E S S E S   I N   T W I S T  ',  &
    ' P A N E L S   ( C T W I S T )')
646 FORMAT (29X,'C O M P L E X   F O R C E S   I N   B A R   E L E M',  &
    ' E N T S   ( C B A R )')
648 FORMAT (25X,'C O M P L E X   F O R C E S   I N   S C A L A R   ',  &
    'S P R I N G S   ( C E L A S 1 )')
649 FORMAT (25X,'C O M P L E X   F O R C E S   I N   S C A L A R   ',  &
    'S P R I N G S   ( C E L A S 2 )')
650 FORMAT (25X,'C O M P L E X   F O R C E S   I N   S C A L A R   ',  &
    'S P R I N G S   ( C E L A S 3 )')
651 FORMAT (25X,'C O M P L E X   F O R C E S   I N   S C A L A R   ',  &
    'S P R I N G S   ( C E L A S 4 )')
652 FORMAT (27X,'C O M P L E X   F O R C E S   I N   R O D   E L E M',  &
    ' E N T S   ( C O N R O D )')
653 FORMAT (17X,'C O M P L E X   F O R C E S   I N   B E N D I N G  ',  &
    ' Q U A D R I L A T E R A L S   ( C Q D P L T )')
654 FORMAT (9X,'C O M P L E X   F O R C E S   I N   G E N E R A L   ',  &
    'Q U A D R I L A T E R A L   E L E M E N T S   ', '( C Q U A D 1 )')
655 FORMAT (9X,'C O M P L E X   F O R C E S   I N   G E N E R A L   ',  &
    'Q U A D R I L A T E R A L   E L E M E N T S   ', '( C Q U A D 2 )')
656 FORMAT (29X,'C O M P L E X   F O R C E S   I N   R O D   E L E M',  &
    ' E N T S   ( C R O D )')
657 FORMAT (7X,'C O M P L E X   F O R C E S   A C T I N G   O N   ',  &
    'S H E A R   P A N E L   E L E M E N T S   (C S H E A R)')
658 FORMAT (16X,'C O M P L E X   F O R C E S   I N   B A S I C   B E',  &
    ' N D I N G   T R I A N G L E S   ( C T R B S C )')
659 FORMAT (12X,'C O M P L E X   F O R C E S   I N   G E N E R A L  ',  &
    ' T R I A N G U L A R   E L E M E N T S   ( C T R I A 1 )')
660 FORMAT (12X,'C O M P L E X   F O R C E S   I N   G E N E R A L  ',  &
    ' T R I A N G U L A R   E L E M E N T S   ( C T R I A 2 )')
661 FORMAT (22X, 'C O M P L E X   F O R C E S   I N   B E N D I N G ',  &
    '  T R I A N G L E S   ( C T R P L T )')
662 FORMAT (28X,'C O M P L E X   F O R C E S   I N   R O D   E L E M',  &
    ' E N T S   ( C T U B E )')
663 FORMAT (27X,'C O M P L E X   F O R C E S   I N   T W I S T   P A',  &
    ' N E L S   ( C T W I S T )')
664 FORMAT (12X,7HELEMENT,20X,4(8HLOCATION,7X ),6X,7HAVERAGE, /14X,  &
    3HID.,26X,1H1,14X,1H2,14X,1H3,14X,1H4,13X,12HAXIAL stress)
665 FORMAT (17X,7HELEMENT,29X,5HAXIAL,39X,9HTORSIONAL, /19X,3HID.,30X,  &
    6HSTRESS,41X,6HSTRESS)
666 FORMAT (17X,7HELEMENT,28X,7HMAXIMUM,39X,7HAVERAGE, /19X,3HID.,31X,  &
    5HSHEAR,41X,5HSHEAR)
667 FORMAT (17X,7HELEMENT,29X,5HAXIAL,41X,6HTORQUE, /19X,3HID.,31X, 5HFORCE)
668 FORMAT (9H  element,7X,5HFIBRE,37X,'- STRESSES IN ELEMENT COORDI',  &
    'NATE SYSTEM -', /4X,3HID.,8X,8HDISTANCE,18X,8HNORMAL-x,  &
    26X,8HNORMAL-y,25X,8HSHEAR-xy)
669 FORMAT (13X,7HELEMENT,33X,'- STRESSES IN ELEMENT COORDINATE SYST',  &
    'EM -', /15X,3HID.,18X,8HNORMAL-x,26X,8HNORMAL-y,26X, 8HSHEAR-xy)
670 FORMAT (2(16X,7HELEMENT,35X), /2(18X,3HID.,20X,5HFORCE,12X))
671 FORMAT (2(16X,7HELEMENT,35X), /2(18X,3HID.,19X,6HSTRESS,12X))
672 FORMAT (17X,7HELEMENT,29X,5HFORCE,42X,5HFORCE)
673 FORMAT (19X,3HID.,30X,7HPTS 1,3,40X,7HPTS 2,4)
674 FORMAT (17X,7HELEMENT,28X,6HMOMENT,41X,6HMOMENT)

END SUBROUTINE ofp1a