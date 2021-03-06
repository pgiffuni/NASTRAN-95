BLOCK DATA ofp1bd
!OFP1BD
 INTEGER :: d1   ,d201 ,d401 ,d601 ,d801 ,d1001,d1201,d1401,d1601,  &
     d1801,d2001,d2201,d2401,d2601,d2801,d3001,d3201,d3401,  &
     d3601,d3801,d4001,d4201
 COMMON /ofpbd1/       d1(200), d201(200), d401(200), d601(200),  &
     d801(200),d1001(200),d1201(200),d1401(200),d1601(200),  &
     d1801(200),d2001(200),d2201(200),d2401(200),d2601(200),  &
     d2801(200),d3001(200),d3201(200),d3401(200),d3601(200),  &
     d3801(200),d4001(200),d4201( 90)
!*****
!     DATA RECORD DEFINITION DATA IS IN THE D-ARRAY...
!*****
!     WHEN ADDING STRINGS REMEMBER THAT
!     FIRST WORD OF EACH STRING = NUMBER OF LINES OF OUTPUT THE FORMAT
!                                 STRING WILL PRODUCE
!*****
!     POINTERS TO THE OFP5BD FORMAT BLOCKS
!     NEGATIVE POINTERS REFER TO ESINGL ARRAYS FOR SPACING OR BCD WORDS
!     POSITIVE POINTERS REFER TO THE E ARRAYS FOR DATA PRINT FORMAT
!     SUMMARY OF FORMAT BLOCKS IN OFP5BD -
!     FORMAT INDEX     FORMAT INDEX     FORMAT INDEX     FORMAT INDEX
!     ------ -----     ------ -----     ------ -----     ------ -----
!     E-ARRAY
!      E15.6    1       F20.4   26         I10   51     '0',I20   76
!      E16.6    2       F16.4   27       I7,1X   52      I10,5X   77
!      E17.6    3       F22.4   28       3X,A4   53               78
!      E18.6    4       E27.6   29     '0',I13   54       I8,2X   79
!      E19.6    5       F12.5   30      1X,I20   55        F8.3   80
!      E20.6    6       E13.5   31    5X,A1,3X   56     '0',I27   81
!      E21.6    7       F13.3   32      1X,I22   57      '0',I5   82
!      E30.6    8       F18.4   33         I12   58      '0',I3   83
!      E26.6    9       F26.4   34      1X,I19   59          I4   84
!      E24.6   10       E14.5   35         I16   60       E11.4   85
!      F11.4   11       F14.3   36          I8   61          A4   86
!      F14.4   12        F5.2   37          I9   62      9E11.3   87
!      E28.6   13       E13.6   38         I11   63       F22.3   88
!      E37.6   14               39         I20   64      /E11.3   89
!      E22.6   15        E9.1   40         I19   65       F19.4   90
!      E14.6   16    6X,A1,3X   41      1X,I23   66        F8.2   91
!      F15.4   17         I15   42         I23   67       E12.5   92
!       F9.4   18       I9,1X   43         I28   68     '0',I12   93
!      F15.3   19      '0',I8   44    /1H ,I18   69       4X,I8   94
!      E23.6   20      1X,I13   45     '0',I15   70               95
!      E35.6   21       1X,I8   46     '0',I14   71               96
!      E25.5   22      '0',I7   47       F22.4   72               97
!      E50.6   23       6X,I8   48       F16.4   73               98
!      F46.4   24      1X,I15   49       F10.4   74               99
!              25      1X,I12   50     '0',I19   75              100
 
!     E-SINGL ARRAY
!          /   -1        /14X  -17        /28X  -33        'ZX'  -49
!        15X   -2         11X  -18        /15X  -34         'C'  -50
!        10X   -3        /24X  -19        /19X  -35        'LZ'  -51
!         5X   -4         '0'  -20        /21X  -36        'CP'  -52
!         1X   -5        ' /'  -21        /11X  -37        'MP'  -53
!       /10X   -6        'EN'  -22        /17X  -38        'C '  -54
!        16X   -7        'DA'  -23          2X  -39          3X  -55
!       '1 '   -8        'DB'  -24         'X'  -40        /30X  -56
!       '2 '   -9        /'0'  -25        'XY'  -41          9X  -57
!       '3 '  -10         23X  -26         'A'  -42        /23X  -58
!       '4 '  -11        /26X  -27        'LX'  -43          6X  -59
!       '5 '  -12         /9X  -28         'Y'  -44         39X  -60
!         7X  -13        /12X  -29        'YZ'  -45         24X  -61
!       /16X  -14        /' '  -30         'B'  -46              -62
!       /13X  -15        /20X  -31        'LY'  -47              -63
!         4X  -16        /32X  -32         'Z'  -48              -64
!*****
 
 DATA d1/  1,   45,   41,    1,    1,    1,    1,    1,    1,    0  &
     ,    1,   50,    3,   16,    1,   16,    1,   16,    1,    1  &
     ,    0,    0,    3,   47,    3,    1,    1,    1,    1,    1  &
     ,    1,   40,   -6,    1,    1,    1,    1,   -2,    1,    1  &
     ,   40,    0,  201,   45,   -4,    1,    1,  -18,   45,   -4  &
     ,    1,    1,    0,  401,   50,    4,   42,    4,   42,    4  &
     ,   42,    4,    0,    0,    1,   49,   15,   15,   15,   15  &
     ,   15,    0,    0,    0,  201,   48,    2,   -5,   40,    2  &
     ,   -5,   40,  -16,   46,    2,   -5,   40,    2,   -5,   40  &
     ,    0,    0,    3,   44,    2,    1,    1,    3,    3,    1  &
     ,   -4,   40,   -1,   -3,    1,    1,    1,   -7,    4,    1  &
     ,   -4,   40,    0,    0,  201,   49,   -4,    1,    1,   -5  &
     ,   40,   42,   -4,    1,    1,   -5,   40,    0,    0,    3  &
     ,   44,    2,    4,    1,    1,   11,    2,    2,    1,  -28  &
     ,    2,    4,    1,    1,   11,    2,    2,    1,    0,    0  &
     ,    0,    1,   46,    6,    1,    1,   12,    5,    1,    1  &
     ,    0,    0,    1,   46,    7,    8,    8,    8,    0,    0  &
     ,    6,   44,  -13,   -8,    9,    9,    9,    9,  -14,   -9  &
     ,    9,    9,    9,    9,  -14,  -10,    9,    9,    9,    9  &
     ,  -14,  -11,    9,    9,    9,    9,  -14,  -12,    9,    9/
 DATA d201/9,    9,    0,    0,    0,    0,    0,    0,    0,    0  &
     ,    4,   44,  -16,   -8,    5,   10,    6,   10,    6,  -15  &
     ,   -9,    5,   10,    6,   10,    6,  -15,  -10,    5,   10  &
     ,    6,   10,    6,    0,    0,    4,   44,   -4,   -8,   13  &
     ,   14,   14,  -17,   -9,   13,   14,   14,  -17,  -10,   13  &
     ,   14,   14,    0,    0,    5,   44,   -4,   -8,   13,   14  &
     ,   14,  -17,   -9,   13,   14,   14,  -17,  -10,   13,   14  &
     ,   14,  -17,  -11,   13,   14,   14,    0,    0,    3,   44  &
     ,   -4,   -8,    5,    4,    4,    4,    4,    4,  -17,   -9  &
     ,    5,    4,    4,    4,    4,    4,    0,    1,   46,   51  &
     ,    6,    6,    6,    6,    6,    0,    0,    0,    0,    0  &
     ,    1,   77,   56,    1,    1,    1,    1,    1,    1,    0  &
     ,    1,   45,   53,   61,   63,   62,   62,   65,   62,   62  &
     ,   65,   62,   62,    0,    0,    1,   55,    5,    5,    5  &
     ,    5,   42,    0,    0,    1,   55,   15,    2,    9,   17  &
     ,   63,    0,    0,    1,    1,   56,    1,    1,    1,    1  &
     ,    1,    1,    0,    0,    1,   57,   58,   15,    3,   15  &
     ,   15,    0,    0,    3,   54,   41,    1,    1,    1,    1  &
     ,    1,    1,  -19,    1,    1,    1,    1,    1,    1,    0  &
     ,    0,    3,  -20,   16,   56,    1,    1,    1,    1,    1/
 DATA d401/1,  -19,    1,    1,    1,    1,    1,    1,    0,    0  &
     ,    3,   54,   41,    1,    1,    1,    1,    1,    1,  -31  &
     ,   17,   17,   17,   17,   17,   17,    0,    0,    3,  -20  &
     ,   16,   56,    1,    1,    1,    1,    1,    1,  -31,   17  &
     ,   17,   17,   17,   17,   17,    0,    0,    1,   59,    4  &
     ,  -21,   16,    4,  -21,   16,    4,  -21,   16,    0,    0  &
     ,    1,   59,    4,  -21,   18,   20,  -21,   18,   20,  -21  &
     ,   18,    0,    1,   66,    8,  -21,   16,    8,  -21,   16  &
     ,    0,    0,    1,   66,    8,  -21,   18,   21,  -21,   18  &
     ,    0,    0,  201,   57,    5,  -21,   16,   67,    5,  -21  &
     ,   16,    0,    0,  201,   57,    5,  -21,   18,   68,    5  &
     ,  -21,   18,    0,    0,    3,   44,    2,    4,  -21,   16  &
     ,    4,  -21,   16,    4,  -21,   16,   -1,   22,    4,  -21  &
     ,   16,    4,  -21,   16,    4,  -21,   16,    0,    0,    3  &
     ,   44,    2,    4,  -21,   18,   20,  -21,   18,   20,  -21  &
     ,   18,   -1,   22,    4,  -21,   18,   20,  -21,   18,   20  &
     ,  -21,   18,    0,    0,    6,   69,   -4,  -22,  -23,   15  &
     ,    1,    1,    1,    6,   -1,   23,    1,    1,    1,    6  &
     ,  -25,  -26,  -22,  -24,   15,    1,    1,    1,   -1,   23  &
     ,    1,    1,    1,    0,    6,   69,   -4,  -22,  -23,   15/
 DATA d601/1,    1,    1,    6,   -1,   24,   17,   17,   17,   26  &
     ,  -25,  -26,  -22,  -24,   15,    1,    1,    1,   -1,   24  &
     ,   17,   17,   17,    0,    0,    3,   71,    1,   16,    1  &
     ,   16,    1,   16,    1,    1,  -17,    2,   16,    1,   16  &
     ,    1,   16,    1,    1,    0,    0,    3,   71,    1,   16  &
     ,    1,   16,    1,   16,    1,    1,  -17,   -5,   11,   12  &
     ,   17,   12,   17,   12,   17,   17,    0,    3,   70,   15  &
     ,   15,   15,   15,   15,  -14,   15,   15,   15,   15,   15  &
     ,    0,    0,    3,   70,   15,   15,   15,   15,   15,  -29  &
     ,   28,   28,   28,   28,   28,    0,    0,  201,    1,    5  &
     ,    1,    9,    5,    1,    0,    0,    1,    1,    1,   16  &
     ,    1,   16,    1,   16,    1,    1,    0,    0,    1,    2  &
     ,   15,   15,   15,   15,   15,    0,    0,  401,    2,    1  &
     ,    4,    1,    4,    1,    4,    1,    0,    0,  201,    1  &
     ,    1,   -5,   40,    2,   -5,   40,   16,    1,   -5,   40  &
     ,    2,   -5,   40,    0,    0,    3,  -20,   16,    2,    1  &
     ,    1,    3,    4,    1,   -4,   40,  -17,    3,    1,    1  &
     ,    3,    4,    1,   -4,   40,    0,    0,  201,    2,    6  &
     ,    1,   -5,   40,    1,    6,    1,   -5,   40,    0,    0  &
     ,    3,  -20,   16,   16,   16,    1,    1,   11,    2,    2/
 DATA d801/1,  -17,    1,   16,    1,    1,   11,    2,    2,    1  &
     ,    0,    0,    1,    1,   16,   16,    1,   11,    6,    1  &
     ,    1,    0,    0,    0,    3,   16,    4,    1,    1,    1  &
     ,    1,    2,    1,   40,  -14,    2,    1,    1,    1,   -2  &
     ,    2,    1,   40,  -30,    0,    1,   10,    8,  -21,   16  &
     ,    8,  -21,   16,    0,    0,    1,   10,    8,  -21,   18  &
     ,   21,  -21,   18,    0,    3,  -20,    1,   15,   15,   15  &
     ,   15,   15,  -14,   15,   15,   15,   15,   15,    0,    0  &
     ,    3,  -20,    1,   15,   15,   15,   15,   15,  -29,   28  &
     ,   28,   28,   28,   28,    0,    0,  201,   20,    5,  -21  &
     ,   16,   20,    5,  -21,   16,    0,    0,  201,   20,    5  &
     ,  -21,   18,   13,    5,  -21,   18,    0,    0,    3,  -20  &
     ,   16,    1,   16,    1,   16,    1,   16,    1,    1,  -17  &
     ,    2,   16,    1,   16,    1,   16,    1,    1,    0,    0  &
     ,    3,  -20,   16,    1,   16,    1,   16,    1,   16,    1  &
     ,    1,  -17,   -5,   11,   12,   17,   12,   17,   12,   17  &
     ,   17,    0,    0,    3,  -20,   16,   16,    4,  -21,   16  &
     ,    4,  -21,   16,    4,  -21,   16,   -1,  -16,   22,    4  &
     ,  -21,   16,    4,  -21,   16,    4,  -21,   16,    0,    3  &
     ,  -20,   16,   16,    4,  -21,   18,   20,  -21,   18,   20/
 DATA d1001/-21,18,   -1,  -16,   22,    4,  -21,   18,   20,  -21  &
     ,   18,   20,  -21,   18,    0,    1,    6,    4,  -21,   16  &
     ,    4,  -21,   16,    4,  -21,   16,    0,    0,    1,    6  &
     ,    4,  -21,   18,   20,  -21,   18,   20,  -21,   18,    0  &
     ,    0,    6,   -1,    5,   -4,  -22,  -23,   15,    1,    1  &
     ,    1,    6,   -1,   23,    1,    1,    1,    6,  -25,  -26  &
     ,  -22,  -24,   15,    1,    1,    1,   -1,   23,    1,    1  &
     ,    1,    0,    0,    6,   -1,    5,   -4,  -22,  -23,   15  &
     ,    1,    1,    1,    6,   -1,   24,   17,   17,   17,   26  &
     ,  -25,  -26,  -22,  -24,   15,    1,    1,    1,   -1,   24  &
     ,   17,   17,   17,    0,    0,    3,   47,   73,   16,   16  &
     ,   16,   16,   74,   16,   16,   16,  -19,   16,   16,   16  &
     ,   16,   74,   16,   16,   16,    0,    1,   46,   72,   15  &
     ,    4,    4,    4,    4,    0,    0,    1,   60,    4,   16  &
     ,   16,   16,   16,   16,   16,   16,    0,    0,    1,    2  &
     ,    4,   16,   16,   16,   16,   16,   16,   16,    0,    0  &
     ,    0,    3,   81,    4,   16,   16,   16,   16,   16,  -32  &
     ,   16,   16,   16,   16,   16,   16,    0,    0,    3,  -20  &
     ,   29,    4,   16,   16,   16,   16,   16,  -32,   16,   16  &
     ,   16,   16,   16,   16,    0,    0,    3,   81,    4,   16/
 DATA d1201/16, 16,   16,   16,  -33,   12,   12,   12,   12,   12  &
     ,   12,    0,    0,    3,  -20,   29,    4,   16,   16,   16  &
     ,   16,   16,  -33,   12,   12,   12,   12,   12,   12,    0  &
     ,    1,   77,    1,   16,    1,   16,    1,   16,    1,    1  &
     ,    0,    1,   42,   31,   31,   31,   31,   31,   31,   31  &
     ,   31,   31,    0,    3,   71,   31,   31,   31,   31,   31  &
     ,   31,   31,   31,   31,  -34,   31,   31,   31,   31,   31  &
     ,   31,   31,   31,   31,    0,    3,   71,   31,   31,   31  &
     ,   31,   31,   31,   31,   31,   31,  -37,   32,   32,   32  &
     ,   32,   32,   32,   32,   32,   32,    0,    3,  -20,   16  &
     ,   31,   31,   31,   31,   31,   31,   31,   31,   31,  -34  &
     ,   31,   31,   31,   31,   31,   31,   31,   31,   31,    0  &
     ,    3,  -20,   16,   31,   31,   31,   31,   31,   31,   31  &
     ,   31,   31,  -37,   32,   32,   32,   32,   32,   32,   32  &
     ,   32,   32,    0,    1,    6,    4,    4,    4,    4,    4  &
     ,    4,    0,    3,  -20,    5,    4,    4,    4,    4,    4  &
     ,    4,  -31,    4,    4,    4,    4,    4,    4,    0,    3  &
     ,  -20,    5,    4,    4,    4,    4,    4,    4,  -14,   33  &
     ,   33,   33,   33,   33,   33,    0,    1,   59,    4,    4  &
     ,    4,    4,    4,    4,    0,    3,   75,    4,    4,    4/
 DATA d1401/4,   4,    4,  -31,    4,    4,    4,    4,    4,    4  &
     ,    0,    3,   75,    4,    4,    4,    4,    4,    4,  -14  &
     ,   33,   33,   33,   33,   33,   33,    0,    1,    7,    9  &
     ,    9,    9,    9,    0,    3,  -20,    6,    9,    9,    9  &
     ,    9,  -36,    9,    9,    9,    9,    0,    3,  -20,    6  &
     ,    9,    9,    9,    9,  -38,   34,   34,   34,   34,    0  &
     ,    1,   55,    9,    9,    9,    9,    0,    3,   76,    9  &
     ,    9,    9,    9,  -36,    9,    9,    9,    9,    0,    3  &
     ,   76,    9,    9,    9,    9,  -38,   34,   34,   34,   34  &
     ,    0,    1,    1,   16,   31,   31,   31,   31,   31,   31  &
     ,   31,   31,   31,    0,    0,  201,   77,    6,    1,   -5  &
     ,   40,   77,    6,    1,   -5,   40,    0,    0,    0,    0  &
     ,    3,   54,   31,   35,   35,   35,   35,   35,   35,   35  &
     ,  -31,   35,   35,   35,   35,   35,   35,   35,   35,    0  &
     ,    0,    3,  -20,   31,   31,   35,   35,   35,   35,   35  &
     ,   35,   35,  -31,   35,   35,   35,   35,   35,   35,   35  &
     ,   35,    0,    0,    5,   45,   31,   35,   35,   35,   35  &
     ,   35,   35,   35,  -15,   35,   35,   35,   35,   35,   35  &
     ,   35,   35,  -31,   35,   35,   35,   35,   35,   35,   35  &
     ,   35,  -31,   35,   35,   35,   35,   35,   35,   35,   35/
 DATA d1601/0,   0,    5,  -20,   31,   31,   35,   35,   35,   35  &
     ,   35,   35,   35,  -15,   35,   35,   35,   35,   35,   35  &
     ,   35,   35,  -31,   35,   35,   35,   35,   35,   35,   35  &
     ,   35,  -31,   35,   35,   35,   35,   35,   35,   35,   35  &
     ,    0,    0,    5,   45,   31,   35,   35,   35,   35,   35  &
     ,   35,   35,  -28,   36,   36,   36,   36,   36,   36,   36  &
     ,   36,  -31,   35,   35,   35,   35,   35,   35,   35,   35  &
     ,  -14,   36,   36,   36,   36,   36,   36,   36,   36,    0  &
     ,    0,    5,  -20,   31,   31,   35,   35,   35,   35,   35  &
     ,   35,   35,  -28,   36,   36,   36,   36,   36,   36,   36  &
     ,   36,  -31,   35,   35,   35,   35,   35,   35,   35,   35  &
     ,  -14,   36,   36,   36,   36,   36,   36,   36,   36,    0  &
     ,    1,   77,   16,   16,    1,   12,    6,    1,    1,    0  &
     ,    3,   77,    3,    1,    1,    1,    1,    2,    1,   40  &
     ,  -14,    2,    1,    1,    1,   -2,    2,    1,   40,  -30  &
     ,    0,    3,  -20,   61,   -4,   31,   35,   35,   35,   35  &
     ,   35,   35,   35,  -31,   35,   35,   35,   35,   35,   35  &
     ,   35,   35,    0,  201,   77,    1,   -5,   40,    2,   -5  &
     ,   40,   62,   -4,    1,   -5,   40,    2,   -5,   40,    0  &
     ,    0,  201,   77,    5,    1,   55,   -4,    5,    1,    0/
 DATA d1801/4,  44,   46,  -39,  -40,   16,  -39,  -41,   16,  -39  &
     ,  -42,   16,  -39,  -43,   -5,   37,   37,   37,    1,    1  &
     ,  -31,  -44,   16,  -39,  -45,   16,  -39,  -46,   16,  -39  &
     ,  -47,   -5,   37,   37,   37,  -31,  -48,   16,  -39,  -49  &
     ,   16,  -39,  -50,   16,  -39,  -51,   -5,   37,   37,   37  &
     ,    0,    3,   44,   46,    1,    3,    3,    3,    3,    3  &
     ,  -35,   16,    3,    3,    3,    3,    3,    0,    4,   44  &
     ,   46,  -39,  -40,   16,  -39,  -41,   16,  -39,  -42,   16  &
     ,  -39,  -43,   -5,   37,   37,   37,    1,    1,   -6,   79  &
     ,  -44,   16,  -39,  -45,   16,  -39,  -46,   16,  -39,  -47  &
     ,   -5,   37,   37,   37,  -31,  -48,   16,  -39,  -49,   16  &
     ,  -39,  -50,   16,  -39,  -51,   -5,   37,   37,   37,    0  &
     ,    3,   44,   46,    1,    3,    3,    3,    3,    3,  -28  &
     ,   79,   16,    3,    3,    3,    3,    3,    0,    1,   50  &
     ,   80,    1,    1,    4,    1,   15,    5,    0,    1,   50  &
     ,   80,    1,   16,   16,   16,    1,    1,   16,   40,    0  &
     ,    9,   44,   -4,  -52,    5,    4,    4,    4,    4,    4  &
     ,  -17,  -52,    5,    4,    4,    4,    4,    4,  -17,  -52  &
     ,    5,    4,    4,    4,    4,    4,  -17,  -52,    5,    4  &
     ,    4,    4,    4,    4,  -17,  -53,    5,    4,    4,    4/
 DATA d2001/4,   4,  -17,  -53,    5,    4,    4,    4,    4,    4  &
     ,  -17,  -53,    5,    4,    4,    4,    4,    4,  -17,  -53  &
     ,    5,    4,    4,    4,    4,    4,    0,    7,   44,   -4  &
     ,  -52,    5,    4,    4,    4,    4,    4,  -17,  -52,    5  &
     ,    4,    4,    4,    4,    4,  -17,  -52,    5,    4,    4  &
     ,    4,    4,    4,  -17,  -53,    5,    4,    4,    4,    4  &
     ,    4,  -17,  -53,    5,    4,    4,    4,    4,    4,  -17  &
     ,  -53,    5,    4,    4,    4,    4,    4,    0,    6,   82  &
     ,  -55,   -8,    1,    1,    1,    1,    1,    1,    1,    1  &
     ,  -28,   -9,    1,    1,    1,    1,    1,    1,    1,    1  &
     ,  -28,  -10,    1,    1,    1,    1,    1,    1,    1,    1  &
     ,  -28,  -11,    1,    1,    1,    1,    1,    1,    1,    1  &
     ,  -28,  -54,    1,    1,    1,    1,    1,    1,    1,    1  &
     ,    0,    5,   82,  -55,  -54,    1,    1,    1,    1,    1  &
     ,    1,    1,    1,  -28,   -8,    1,    1,    1,    1,    1  &
     ,    1,    1,    1,  -28,   -9,    1,    1,    1,    1,    1  &
     ,    1,    1,    1,  -28,  -10,    1,    1,    1,    1,    1  &
     ,    1,    1,    1,    0,    0,    2,   47,   72,   87,    0  &
     ,    4,   47,   72,   15,   15,   15,   15,  -56,   15,   15  &
     ,   15,   15,  -56,   15,   15,   15,   15,    0,    0,    0/
 DATA d2201/6,  47,   72,   87,  -56,   87,  -56,   87,  -56,   87  &
     ,  -56,   87,    0,    5,   47,   72,   15,   15,   15,   15  &
     ,  -56,   15,   15,   15,   15,  -56,   15,   15,   15,   15  &
     ,  -56,   15,   15,   15,   15,    0,    0,    0,    5,   47  &
     ,   72,    8,    8,    8,  -56,    8,    8,    8,  -56,    8  &
     ,    8,    8,  -56,    8,    8,    8,    0,    1,  -26,   57  &
     ,  -16,    4,   -4,   33,    0,    1,   -5,   51,   -5,   50  &
     ,  -16,   86,   86,  -16,    1,    1,    1,    1,    1,    1  &
     ,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0  &
     ,    9,   44,  -16,   -8,    1,    1,   16,   16,    1,   16  &
     ,   16,    1,  -34,    1,    1,   16,   16,    1,   16,   16  &
     ,    1,  -15,  -10,    1,    1,   16,   16,    1,   16,   16  &
     ,    1,  -34,    1,    1,   16,   16,    1,   16,   16,    1  &
     ,  -15,  -12,    1,    1,   16,   16,    1,   16,   16,    1  &
     ,  -34,    1,    1,   16,   16,    1,   16,   16,    1,  -15  &
     ,  -54,    1,    1,   16,   16,    1,   16,   16,    1,  -34  &
     ,    1,    1,   16,   16,    1,   16,   16,    1,    0,    0  &
     ,    4,   54,  -57,   -8,    6,    6,    6,    6,    6,  -58  &
     ,  -10,    6,    6,    6,    6,    6,  -58,  -12,    6,    6  &
     ,    6,    6,    6,    0,    0,    0,    0,    0,    0,    0/
 DATA d2401/5,  44,  -13,   -8,    1,    1,    1,    1,    1,    1  &
     ,    1,  -14,  -10,    1,    1,    1,    1,    1,    1,    1  &
     ,  -14,  -12,    1,    1,    1,    1,    1,    1,    1,  -14  &
     ,  -54,    1,    1,    1,    1,    1,    1,    1,    0,    0  &
     ,    1,   77,   -5,   15,   15,   15,   15,   15,    0,    0  &
     ,  401,   51,  -59,    1,   58,  -59,    1,   58,  -59,    1  &
     ,   58,  -59,    1,    0,    0,    3,  -20,   61,  -55,  -55  &
     ,    2,    1,    1,    3,    4,    1,   -4,   40,  -17,    3  &
     ,    1,    1,    3,    4,    1,   -4,   40,    0,    0,  201  &
     ,   62,  -13,    6,    1,   -5,   40,   61,  -13,    6,    1  &
     ,   -5,   40,    0,    0,    3,  -20,   61,  -55,  -55,   16  &
     ,   16,    1,    1,   11,    2,    2,    1,  -17,    1,   16  &
     ,    1,    1,   11,    2,    2,    1,    0,    0,    1,   62  &
     ,  -13,    4,   16,   16,   16,   16,   16,   16,   16,    0  &
     ,    0,    3,   44,   42,    5,    1,    1,   11,    2,    2  &
     ,    1,  -28,   42,    5,    1,    1,   11,    2,    2,    1  &
     ,    0,    0,    3,   44,   86,  -59,  -59,    4,    1,    1  &
     ,   11,    2,    2,    1,  -28,   86,  -59,  -59,    4,    1  &
     ,    1,   11,    2,    2,    1,    0,    9,   47,   72,   15  &
     ,   15,   15,   15,  -56,   15,   15,   15,   15,   -1,  -56/
 DATA d2601/15, 15,   15,   15,  -56,   15,   15,   15,   15,   -1  &
     ,  -56,   15,   15,   15,   15,  -56,   15,   15,   15,   15  &
     ,    0,    9,   47,   72,   15,   15,   15,   15,  -58,   88  &
     ,   88,   88,   88,   -1,  -56,   15,   15,   15,   15,  -58  &
     ,   88,   88,   88,   88,   -1,  -56,   15,   15,   15,   15  &
     ,  -58,   88,   88,   88,   88,    0,    3,   47,   72,   87  &
     ,  -56,   87,    0,   12,   47,   72,   15,   15,   15,   15  &
     ,  -56,   15,   15,   15,   15,   -1,  -56,   15,   15,   15  &
     ,   15,  -56,   15,   15,   15,   15,   -1,  -56,   15,   15  &
     ,   15,   15,  -56,   15,   15,   15,   15,   -1,  -56,   15  &
     ,   15,   15,   15,  -56,   15,   15,   15,   15,    0,   12  &
     ,   47,   72,   15,   15,   15,   15,  -58,   88,   88,   88  &
     ,   88,   -1,  -56,   15,   15,   15,   15,  -58,   88,   88  &
     ,   88,   88,   -1,  -56,   15,   15,   15,   15,  -58,   88  &
     ,   88,   88,   88,   -1,  -56,   15,   15,   15,   15,  -58  &
     ,   88,   88,   88,   88,    0,   15,   47,   72,   87,  -56  &
     ,   87,   -1,  -56,   87,  -56,   87,   -1,  -56,   87,  -56  &
     ,   87,   -1,  -56,   87,  -56,   87,   -1,  -56,   87,  -56  &
     ,   87,    0,    0,    0,    0,    0,    9,   89,   90,   15  &
     ,   15,   15,   15,  -56,   15,   15,   15,   15,   -1,  -56/
 DATA d2801/15, 15,   15,   15,  -56,   15,   15,   15,   15,   -1  &
     ,  -56,   15,   15,   15,   15,  -56,   15,   15,   15,   15  &
     ,    0,    9,   89,   90,   15,   15,   15,   15,  -58,   88  &
     ,   88,   88,   88,   -1,  -56,   15,   15,   15,   15,  -58  &
     ,   88,   88,   88,   88,   -1,  -56,   15,   15,   15,   15  &
     ,  -58,   88,   88,   88,   88,    0,    3,   89,   90,   87  &
     ,  -56,   87,    0,   12,   89,   90,   15,   15,   15,   15  &
     ,  -56,   15,   15,   15,   15,   -1,  -56,   15,   15,   15  &
     ,   15,  -56,   15,   15,   15,   15,   -1,  -56,   15,   15  &
     ,   15,   15,  -56,   15,   15,   15,   15,   -1,  -56,   15  &
     ,   15,   15,   15,  -56,   15,   15,   15,   15,    0,   12  &
     ,   89,   90,   15,   15,   15,   15,  -58,   88,   88,   88  &
     ,   88,   -1,  -56,   15,   15,   15,   15,  -58,   88,   88  &
     ,   88,   88,   -1,  -56,   15,   15,   15,   15,  -58,   88  &
     ,   88,   88,   88,   -1,  -56,   15,   15,   15,   15,  -58  &
     ,   88,   88,   88,   88,    0,   15,   89,   90,   87,  -56  &
     ,   87,   -1,  -56,   87,  -56,   87,   -1,  -56,   87,  -56  &
     ,   87,   -1,  -56,   87,  -56,   87,   -1,  -56,   87,  -56  &
     ,   87,    0,    4,   89,   90,   15,   15,   15,   15,  -56  &
     ,   15,   15,   15,   15,  -56,   15,   15,   15,   15,    0/
 DATA d3001/5,  89,   90,   15,   15,   15,   15,  -56,   15,   15  &
     ,   15,   15,  -56,   15,   15,   15,   15,  -56,   15,   15  &
     ,   15,   15,    0,    2,   89,   90,   87,    0,    6,   89  &
     ,   90,   87,  -56,   87,  -56,   87,  -56,   87,  -56,   87  &
     ,    0,    3,  -20,   92,   31,   91,   31,   31,   91,   31  &
     ,   91,   31,   91,   31,   91,  -15,   31,   91,   31,   31  &
     ,   91,   31,   91,   31,   91,   31,   91,    0,    3,  -20  &
     ,   92,   31,   91,   31,   31,   91,   31,   91,   31,   91  &
     ,   31,   91,  -28,   32,   91,   32,   32,   91,   32,   91  &
     ,   32,   91,   32,   91,    0,    0,    3,  -20,   92,   16  &
     ,   16,   16,    1,   16,   16,    1,   16,   -1,  -16,  -16  &
     ,   36,   36,   36,   19,   36,   36,   19,   36,    0,    3  &
     ,  -20,   92,   16,   16,   16,    1,   16,   16,    1,   16  &
     ,  -15,   16,   16,   16,    1,   16,   16,    1,   16,    0  &
     ,    8,   46,   -5,   84,  -13,   84,  -55,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,  -15,   -2,   61,   46,  -39/
 DATA d3201/85,-39,   85,  -39,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -39,   85,  -39,   85,    0,    8,   85,   -5,   84  &
     ,   -4,   84,  -55,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -39,   85,  -39  &
     ,   85,    0,    8,   46,   -5,   84,  -13,   84,  -55,   61  &
     ,   46,  -39,   85,  -21,   85,  -55,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -15,   -2,   61,   46,  -39,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,  -55,   85,  -21,   85,  -15  &
     ,   -2,   61,   46,  -39,   85,  -21,   85,  -55,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -21,   85,  -55,   85,  -21,   85,  -55,   85,  -21  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -55,   85,  -21,   85,  -15,   -2,   61/
 DATA d3401/46,-39,   85,  -21,   85,  -55,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -15,   -2,   61,   46,  -39,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,  -55,   85,  -21,   85,  -15  &
     ,   -2,   61,   46,  -39,   85,  -21,   85,  -55,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,    0,    8,   85,   -5,   84  &
     ,   -4,   84,  -55,   61,   46,  -39,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -55,   85,  -21,   85,  -15,   -2,   61  &
     ,   46,  -39,   85,  -21,   85,  -55,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -15,   -2,   61,   46,  -39,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,  -55,   85,  -21,   85,  -15  &
     ,   -2,   61,   46,  -39,   85,  -21,   85,  -55,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,  -15,   -2,   61,   46,  -39  &
     ,   85,  -21,   85,  -55,   85,  -21,   85,  -55,   85,  -21  &
     ,   85,  -15,   -2,   61,   46,  -39,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -55,   85,  -21,   85,  -15,   -2,   61  &
     ,   46,  -39,   85,  -21,   85,  -55,   85,  -21,   85,  -55  &
     ,   85,  -21,   85,  -15,   -2,   61,   46,  -39,   85,  -21  &
     ,   85,  -55,   85,  -21,   85,  -55,   85,  -21,   85,    0  &
     ,    3,   -1,   50,    4,    1,    4,    1,    4,    4,  -15  &
     ,    4,    1,    4,   -2,    4,    4,    0,    0,    0,    0/
 DATA d3601/6,  -1,   48,   -3,  -22,  -23,   15,    1,    1,    1  &
     ,    6,   -1,   24,   17,   17,   17,   26,   -1,  -61,  -22  &
     ,  -24,   15,    1,    1,    1,    6,   -1,   24,   17,   17  &
     ,   17,   26,    0,    3,   47,    3,    1,    1,    1,    1  &
     ,    1,    1,   40,   -6,    1,    1,    1,    1,    1,    1  &
     ,    1,   40,    0,    0,    3,    1,    1,    1,    4,    1  &
     ,    4,    4,  -15,    4,    1,    4,   -2,    4,    4,  -30  &
     ,    0,    3,   77,    3,    1,    1,    1,    1,    1,    1  &
     ,   40,  -14,    2,    1,    1,    1,    1,    1,    1,   40  &
     ,    0,    6,   -1,    5,   -4,  -22,  -23,    5,    1,    4  &
     ,    1,    4,    1,  -33,    5,    1,    4,    1,    4,    1  &
     ,  -25,  -26,  -22,  -24,    5,    1,    4,   -2,    4,    1  &
     ,  -33,    5,    1,    4,   -2,    4,    1,    0,    6,   -1  &
     ,    5,   -4,  -22,  -23,    5,    1,    4,    1,    4,    1  &
     ,  -33,   17,   17,   33,   17,   33,   17,  -25,  -26,  -22  &
     ,  -24,    5,    1,    4,   -2,    4,    1,  -33,   17,   17  &
     ,   33,   -2,   33,   17,    0,    6,   -1,    5,   -4,  -22  &
     ,  -23,   15,    1,    1,    1,    6,   -1,   23,    1,    1  &
     ,    1,    6,  -25,  -26,  -22,  -24,   15,    1,    1,    1  &
     ,    6,   -1,   23,    1,    1,    1,    6,    0,    6,   -1/
 DATA d3801/5,  -4,  -22,  -23,   15,    1,    1,    1,    6,   -1  &
     ,   24,   17,   17,   17,   26,  -25,  -26,  -22,  -24,   15  &
     ,    1,    1,    1,    6,   -1,   24,   17,   17,   17,   26  &
     ,    0,    6,   -1,   48,   -3,  -22,  -23,    5,    1,    4  &
     ,    1,    4,    1,  -33,    5,    1,    4,    1,    4,    1  &
     ,  -25,  -26,  -22,  -24,    5,    1,    4,   -2,    4,    1  &
     ,  -33,    4,    1,    4,   -2,    4,    1,    0,    6,   -1  &
     ,   48,   -3,  -22,  -23,    5,    1,    4,    1,    4,    1  &
     ,  -33,   17,   17,   33,   17,   33,   17,  -25,  -26,  -22  &
     ,  -24,    5,    1,    4,   -2,    4,    1,  -33,   17,   17  &
     ,   33,   -2,   33,   17,    0,    6,   -1,   48,   -3,  -22  &
     ,  -23,   15,    1,    1,    1,    6,   -1,   23,    1,    1  &
     ,    1,    6,   -1,  -61,  -22,  -24,   15,    1,    1,    1  &
     ,    6,   -1,   23,    1,    1,    1,    6,    0,    1,  -60  &
     ,   42,   -7,    2,    0,    1,  -60,    2,   -7,   16,    0  &
     ,    3,    4,   63,   -4,    4,  -21,   16,    4,  -21,   16  &
     ,    4,  -21,   16,   -1,   64,   -4,    4,  -21,   16,    4  &
     ,  -21,   16,    4,  -21,   16,    0,    3,   44,   63,   -4  &
     ,    4,  -21,   18,   20,  -21,   18,   20,  -21,   18,   -1  &
     ,   64,   -4,    4,  -21,   18,   20,  -21,   18,   20,  -21/
 DATA d4001/18,  0,    3,  -20,   16,   62,   -4,    4,  -21,   16  &
     ,    4,  -21,   16,    4,  -21,   16,   -1,   -2,   62,   -4  &
     ,    4,  -21,   16,    4,  -21,   16,    4,  -21,   16,    0  &
     ,    3,  -20,   16,   62,   -4,    4,  -21,   18,   20,  -21  &
     ,   18,   20,  -21,   18,   -1,   -2,   62,   -4,    4,  -21  &
     ,   18,   20,  -21,   18,   20,  -21,   18,    0,    0,    0  &
     ,    2,  -20,  -16,   61,   -5,   -5,   31,   -5,   31,   -5  &
     ,   31,   -5,   31,   -5,   31,   -5,   31,   -5,   31,   -5  &
     ,   31,    0,    2,  -20,   -5,   61,  -57,   -5,   31,   -5  &
     ,   31,   -5,   31,   -5,   31,   -5,   31,   -5,   31,   -5  &
     ,   31,   -5,   31,    0,    1,   -5,    1,  -55,   -5,   31  &
     ,   -5,   31,   -5,   31,   -5,   31,   -5,   31,   -5,   31  &
     ,   -5,   31,   -5,   31,    0,    1,   50,   31,   91,   31  &
     ,   31,   91,   31,   91,   31,   91,   31,   91,    0,    1  &
     ,   92,   31,   91,   31,   31,   91,   31,   91,   31,   91  &
     ,   31,   91,    0,    3,   93,   31,   91,   31,   31,   91  &
     ,   31,   91,   31,   91,   31,   91,  -15,   31,   91,   31  &
     ,   31,   91,   31,   91,   31,   91,   31,   91,    0,    3  &
     ,   93,   31,   91,   31,   31,   91,   31,   91,   31,   91  &
     ,   31,   91,  -28,   32,   91,   32,   32,   91,   32,   91/
 DATA d4201/32, 91,   32,   91,    0,    3,  -20,   94,   16,   16  &
     ,   16,    1,   16,   16,    1,   16,   -1,  -16,  -16,   36  &
     ,   36,   36,   19,   36,   36,   19,   36,    0,    3,  -20  &
     ,   94,   16,   16,   16,    1,   16,   16,    1,   16,  -15  &
     ,   16,   16,   16,    1,   16,   16,    1,   16,    0,    3  &
     ,  -20,   85,   61,    1,    3,    3,    3,    3,    3,  -36  &
     ,   16,    3,    3,    3,    3,    3,    0,    3,  -20,   85  &
     ,   61,    1,    3,    3,    3,    3,    3,  -37,   79,   16  &
     ,    3,    3,    3,    3,    3,    0,    0,    0,    0,    0/

END BLOCK DATA
