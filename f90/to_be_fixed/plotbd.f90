BLOCK DATA plotbd
!PLOTBD
 IMPLICIT INTEGER (a-z)
 INTEGER :: char1(60,3),char2(60,1),chr19(2,79),chram(2,88),  &
     chrnz(2,84),chlpqm(2,52),chrsym(2,19),  &
     npens(20,2),pltype(20,2),pbfsiz(20,2),eof(20,2)
 REAL :: chrscl,cntr,DATA,d02,d03,edge,g,maxdef,papsiz, scale,s0s,vangle
 COMMON /char94/ CHAR(60,4)
 COMMON /chrdrw/ lstchr,chrind(60),chr(2,350)
 COMMON /xxparm/ bufsiz,
!    1    ... PLOTTING DATA  &
 camera,bframs,pltmdl(2),tapden,
!    2    ... PEN + PAPER DATA  &
 nopens,papsiz(2),paptyp(2),pensiz(8),penclr(8,2),penpap,
!    3    ... SCALING DATA  &
 scale(2),fscale,maxdef,defmax,
!    4    ... VIEWING DATA  &
 axis(3),daxis(3),vangle(5),view(4),
!    5    ... VANTAGE POINT,               PROJECTION,OCULAR SEPARATION  &
 fvp,r0,s0l,s0r,t0,d0,d02,d03,prject,    s0s,
!    6    ... ORIGIN DATA  &
 forg,org,norg,origin(11),edge(11,4),xy(11,3),
!    7    ... CONTOUR PLOTTING DATA  &
 ncntr,cntr(50),icntvl,where,DIRECT,subcas,flag,value, lasset,
!    8    ... DATA FOR USER PLOT TITLE CARD  &
 fpltit,pltitl(17),color,layer,
!    9    ... OFFSET SCALE (WILL BE SET TO 1 BY PLTSET)  &
 offscl
 COMMON /pltdat/ model,ploter,xymin(2),xymax(2),axymax(2),  &
     xyedge(11),chrscl,pltdat(20),DATA(20,2)
 COMMON /symbls/ nsym,symbol(20,2)
 COMMON /pltscr/ ncor,pltsc(50)
 COMMON /drwaxs/ g(12)
 
! ... EQUIV FOR   /CHAR94/...
 EQUIVALENCE (CHAR(1,1),char1(1,1))  , (CHAR(1,4),char2(1,1))
 
! ... EQUIV FOR   /CHRDRW/...
 EQUIVALENCE (chr(1,  1),chr19(1,1)) , (chr(1, 80),chram(1,1)) ,  &
     (chr(1,168),chrnz(1,1)) , (chr(1,252),chlpqm(1,1)), (chr(1,304),chrsym(1,1))
 
! ... EQUIV FOR   /PLTDAT/...
 EQUIVALENCE (DATA( 7,1),npens(1,1)) , (DATA(10,1),pltype(1,1)),  &
     (DATA(12,1),pbfsiz(1,1)), (DATA(13,1),eof(1,1))
 
 DATA char1 / 1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1HA,1HB,1HC,1HD,1HE,  &
     1HF,1HG,1HH,1HI,1HJ,1HK,1HL,1HM,1HN,1HO,1HP,1HQ,1HR,1HS,1HT,  &
     1HU,1HV,1HW,1HX,1HY,1HZ,1H(,1H),1H+,1H-,1H*,1H/,1H=,1H.,1H,,  &
     1H$,1H-,1H ,12*0,
 
! ... THE FOLLOWING ARE NUMERIC EQUIVALENTS OF 7094 BINARY CHARACTERS.
 
!    5  ... 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F  &
 00,01,02,03,04,05,06,07,08,09,17,18,19,20,21,22,
!    6  ... G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V  &
 23,24,25,33,34,35,36,37,38,39,40,41,50,51,52,53,
!    7  ... W, X, Y, Z, (, ),  , -, *, /, =, ., ,, $, -,BLANK  &
 54,55,56,57,60,28,16,32,44,49,11,27,59,43,12,48,
!    8  . EOR,EOF, SPECIAL, FILLER  &
 58, 15, 63,42,26,   7*0,
 
! ... THE FOLLOWING ARE NUMBERIC EQUIVALENTS OF 7094 BCD CHARACTERS.
 
!    9  ... 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F  &
 10,01,02,03,04,05,06,07,08,09,49,50,51,52,53,54,
!    O  ... G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V  &
 55,56,57,33,34,35,36,37,38,39,40,41,18,19,20,21,
!   11  ... W, X, Y, Z, (, ),  , -, *, /, =, ., ,, $, -,BLANK  &
 22,23,24,25,28,60,48,32,44,17,11,59,27,43,12,16,
!   12  . EOR,EOF, SPECIAL, FILLER  &
 26, 15, 31,42,58,   7*0/
 
! ... THE FOLLOWING ARE NUMERIC VALUES ON CDC 6600 TO PRODUCE 7094 BCD
!     CHARACTERS.
 
 DATA char2 /
!    1  ... 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F  &
 27,28,29,30,31,32,33,34,35,36,01,02,03,04,05,06,
!    2  ... G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V  &
 07,08,09,10,11,12,13,14,15,16,17,18,19,20,21,22,
!    3  ... W, X, Y, Z, (, ),  , -, *, /, =, ., ,, $, -,BLANK  &
 23,24,25,26,41,42,37,38,39,40,44,47,46,43,52,45,
!    4  . EOR,EOF, SPECIAL, FILLER  &
 50, 49, 55,54,58,   7*0/
 
! ... DATA FOR DRAWING 6X6 CHARACTERS (8 UNITS WIDE - 16 UNITS HIGH).
 
 
!     THE FOLLOWING ARE INDICES USED TO DRAW CHARACTERS.
 
 DATA lstchr,chrind /  52,
!    1   0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F  &
 -25,001,006,014,027,031,041,052,055,071,080,086,098,106,113,120,
!    2   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V  &
 126,136,142,148,155,160,163,168,172,181,188,199,208,220,225,231,
!    3   W   X   Y   Z   (   )   +   -   *   /   =   .   ,   $   -  DOT  &
 234,239,243,248,252,256,260,264,266,274,276,280,285,287,302,304,
!    4  CIRCLE SQUARE DIAMOND TRIANGLE END FILLER  &
 -25,   309,    314,    319,   323, 7*0 /
 
! ... DATA FOR DRAWING CHARACTERS 1 TO 9.
 
 DATA chr19/ 2,5, 3,6, 3,0, 2,0, 4,0,  &
     0,5, 1,6, 4,6, 5,5, 5,4, 0,1, 0,0, 5,0,  &
     0,5, 1,6, 4,6, 5,5, 5,4, 4,3, 2,3, 4,3, 5,2, 5,1, 4,0, 1,0, 0,1,  &
     4,0, 4,6, 0,2, 5,2, 5,6, 0,6, 0,3, 1,4, 4,4, 5,3, 5,1, 4,0, 1,0, 0,1,  &
     4,6, 1,6, 0,5, 0,1, 1,0, 4,0, 5,1, 5,2, 4,3, 1,3, 0,2, 0,6, 5,6, 2,0,  &
     4,3, 5,4, 5,5, 4,6, 1,6, 0,5, 0,4, 1,3, 4,3, 5,2, 5,1, 4,0,  &
     1,0, 0,1, 0,2, 1,3, 5,0, 5,5, 4,6, 2,6, 1,5, 1,4, 2,3, 4,3, 5,4/
 
! ... DATA FOR DRAWING CHARACTERS A TO M.
 
 DATA chram / 0,0, 3,6, 5,2, 1,2, 5,2, 6,0,  &
     0,0, 0,6, 4,6, 5,5, 5,4, 4,3, 0,3, 4,3, 5,2, 5,1, 4,0, 0,0,  &
     5,5, 4,6, 1,6, 0,5, 0,1, 1,0, 4,0, 5,1, 5,4, 4,6, 0,6, 0,0, 4,0, 5,2, 5,4,  &
     5,6, 0,6, 0,3, 3,3, 0,3, 0,0, 5,0, 5,6, 0,6, 0,3, 3,3, 0,3, 0,0,  &
     5,5, 4,6, 1,6, 0,5, 0,1, 1,0, 4,0, 5,1, 5,3, 3,3,  &
     0,6, 0,0, 0,3, 5,3, 5,0, 5,6, 2,6, 4,6, 3,6, 3,0, 2,0, 4,0,  &
     3,6, 5,6, 4,6, 4,1, 3,0, 1,0, 0,1, 0,6, 0,0,-5,0, 0,3, 5,6,  &
     0,6, 0,0, 5,0, 0,0, 0,6, 3,0, 6,6, 6,0/
 
! ... DATA FOR DRAWING CHARACTERS N TO Z.
 
 DATA chrnz / 0,0, 0,6, 5,0, 5,6,  &
     6,5, 5,6, 1,6, 0,5, 0,1, 1,0, 5,0, 6,1, 6,5,  &
     0,0, 0,6, 4,6, 5,5, 5,4, 4,3, 0,3,  &
     6,5, 5,6, 1,6, 0,5, 0,1, 1,0, 5,0, 6,1, 6,5,-4,2, 6,0,  &
     0,0, 0,6, 4,6, 5,5, 5,4, 4,3, 0,3, 3,3, 5,0,  &
     5,5, 4,6, 1,6, 0,5, 0,4, 1,3, 4,3, 5,2, 5,1, 4,0, 1,0, 0,1,  &
     0,6, 3,6, 3,0, 3,6, 6,6, 0,6, 0,1, 1,0, 4,0, 5,1, 5,6,  &
     0,6, 3,0, 6,6, 0,6, 1,0, 3,4, 5,0, 6,6,  &
     0,6, 6,0,-6,6, 0,0, 0,6, 3,3, 3,0, 3,3, 6,6,  &
     0,6, 6,6, 0,0, 6,0/
 
! ... DATA FOR DRAWING CHARACTERS ( TO -.
 
 DATA chlpqm / 5,6, 3,4, 3,2, 5,0,  &
     1,6, 3,4, 3,2, 1,0, 3,5, 3,1,-1,3, 5,3,  &
     1,3, 5,3, 1,5, 5,1,-3,5, 3,1,-5,5, 1,1,-5,3, 1,3,  &
     0,0, 6,6, 1,4, 4,4,-1,2, 4,2,  &
     2,0, 2,1, 3,1, 3,0, 2,0, 1,0, 3,2,  &
     6,5, 5,6, 1,6, 0,5, 0,4, 1,3, 5,3, 6,2, 6,1, 5,0, 3,0, 3,6, 3,0, 1,0, 0,1,  &
     3,6, 3,4/
 
! ... DATA FOR DRAWING DOT, SQUARE, DIAMOND, TRIANGLE.
 
 DATA chrsym / 3,4, 2,3, 3,2, 4,3, 3,4,  &
     0,0, 0,6, 6,6, 6,0, 0,0, 3,6, 0,3, 3,0, 6,3, 3,6,  &
     0,0, 3,6, 6,0, 0,0/
 
 DATA bufsiz / 0 /,
!    1 ... CAMERA 2,  1 BLANK FRAME,  PLOTTER MODEL --M,1--  &
 camera,bframs,pltmdl,tapden / 2, 1, 1HM, 1, 0 /,
!    2 ... PAPER = DEFAULT,VELLUM...PEN SIZE = 1, COLOR = BLACK  &
 nopens,papsiz,paptyp,pensiz,penclr /  &
     8, 2*0., 4HVELL, 2HUM, 8*1, 8*4HBLAC, 8*1HK /,
!    3 ... FIND THE SCALES, MAX DEFORMATION = 0  &
 scale(2),fscale,maxdef / 1.,1,0. /,
!    4 ... AXES = +X,+Y,+Z, VIEW ANGLES  &
 axis,daxis,vangle / 1,2,3,1,2,3, 0.,-1.e10,34.27,23.17,0./,
!    5 ... FIND VANTAGE POINT, ORTHOGRAPIC PROJECTION, PLANE+OCULAR SEP.  &
 fvp,prject,d02,d03,s0s / 1,1,1.,2.,2.756  /,
!    6 ... LEFT=BOTTOM=0, RIGHT=TOP=1.  &
 norg,org,forg,edge / 10,0,1,22*0.,22*1.   /,
!    7 ... NCNTR=10=NO. CONTOURS, CNTR=LIST CONTOUR VALUES, ICNTVL=
!          MAJOR PRIN. STRESS, WHERE = Z1, DIRECT = COMMON  &
 ncntr,cntr,icntvl,where,DIRECT,flag,lasset/  &
     10 ,50*0.0,1,     1,    2,     0,   0     /,
!    8 ... DATA FOR USER PLOT TITLE CARD  &
 fpltit,pltitl / 0, 17*4H      /,
!    9 ... OFFSET SCALE (AND ALSO PLOT TAPE MESSAGE CONTROL)  &
 offscl / 0    /
 
! ... PLOTTER DATA.
 
 DATA  model,ploter,chrscl / -1, 1, 1.0 /
 
! ... 1  NASTRAN GENERAL PURPOSE MICROFILM PLOTTER.
 
 DATA  DATA( 1,1) /1023.0   /, DATA( 2,1) /1023.0   /,  &
     DATA( 3,1) / 146.1429/, DATA( 4,1) /   8.0   /,  &
     DATA( 5,1) /  16.0   /, DATA( 6,1) /1023.0   /,  &
     DATA( 8,1) /   0.0   /, DATA( 9,1) /   0.0   /,  &
     DATA(11,1) /4HPLT2   /, DATA(14,1) /1484.761 /,  &
     DATA(15,1) /   0.0   /, DATA(16,1) /   0.0   /,  &
     DATA(17,1) /   0.0   /, DATA(18,1) /   0.0   /,  &
     DATA(19,1) /   0.0   /, DATA(20,1) /   0.0   /
 
! ... 2  NASTRAN GENERAL PURPOSE TABLE OR DRUM PLOTTER
 
 DATA  DATA( 1,2) /3000.0   /, DATA( 2,2) /3000.0   /,  &
     DATA( 3,2) / 100.0   /, DATA( 4,2) /   8.0   /,  &
     DATA( 5,2) /  16.0   /, DATA( 6,2) /3000.0   /,  &
     DATA( 8,2) /   0.0   /, DATA( 9,2) /   0.0   /,  &
     DATA(11,2) /4HPLT2   /, DATA(14,2) / 100.0   /,  &
     DATA(15,2) /   0.0   /, DATA(16,2) /   0.0   /,  &
     DATA(17,2) /   0.0   /, DATA(18,2) /   0.0   /,  &
     DATA(19,2) /   0.0   /, DATA(20,2) /   0.0   /
 
 DATA npens(1,1),pltype(1,1),pbfsiz(1,1),eof(1,1)/ 64,-1,3000,1 /,  &
     npens(1,2),pltype(1,2),pbfsiz(1,2),eof(1,2)/ 64,-2,3000,1 /
 
! ... SYMBOL DATA.
 
 DATA nsym,symbol /  9,
!           X, *, +, -, DOT, CIRCLE, SQUARE, DIAMOND, TRIANGLE  &
 34,41,39,40,  48,     49,     50,      51,       52, 11*0,  &
     34,41,39,40,  48,     49,     50,      51,       52, 11*0/
 
! ... PLOTTER SCRATCH AREA
 
!          NCOR = ARRAY LENGTH
 DATA ncor,pltsc  / 50,50*0 /
 
! ... DATA FOR DRAWING A X-Y-Z COORDINATE TRIAD IN /DRWAXS/
!     G   - X,Y,Z COORD. POINT DATA AND SYMBOLS
 
 DATA    g     / 9*0.0, 1HX, 1HY, 1HZ /
 
END SUBROUTINE plot
