      SUBROUTINE XSEM00
C **********************************************************************
C THE PURPOSE OF THIS ROUTINE IS TO EXECUTE THE PREFACE AND THEN TO
C EXECUTE MODULES ACCORDING TO THE DMAP.  THE DMAP IS READ FROM THE
C OSCAR.  FOR EACH MODULE TO BE EXECUTED, THE FIST AND XVPS ARE SETUP.
C
CWKBD 5/95
C     INTEGER ORF
      INTEGER ANDF ,DATABF,ERRFLG,FIST  ,FISTNM,FSTRST,OPNTR
     1       ,PARML,PARAM ,PARMN ,POOL  ,RSHIFT,SCRTCH
     2       ,VPS  ,VPARML,TYPECD,VPSX  ,WORDB ,WORDE
     3       ,PLOTF,EXIT  ,SYSBUF,SUBNAM(2)
      INTEGER EQUIV(2), PURGE(2), XEQU, XPUR, XSAV, YCHK
CWKBI 5/95
      CHARACTER*4   WORDC
C
      LOGICAL LVAX
C
      DIMENSION SCRTCH(3),WORDB(4),WORDE(2),NUMBR(10)
C
      COMMON/MACHIN/MACH
      COMMON/SEM   /MASK  ,MASK2 ,MASK3 ,LINKNM(15)
     1
     F      /SYSTEM/SYSBUF,XX(20),LINKNO,XXX(16),NBPC,NBPW,NCPW,XXXX(53)
     F             ,ISPERLNK
     G
     H      /XLINK /LXLINK,MAXLNK,MXLINK(1)
     1
     2      /XFIST /FIST(2)
     3
     4      /XPFIST/FSTRST
     5
     6      /OSCENT/INOSCR(200)
     7
     8      /ZZZZZZ/DATABF(1)
     9
     A      /BLANK /PARAM(100)
     B
     C      /XVPS  /VPS(1)
     D
     E      /MSGX  /NMSG
C
      EQUIVALENCE (XX(1),NOUT),
     1            (xx(3),nin) ! input file number
      EQUIVALENCE (XX(19),PLOTF)
      equivalence (xx(17),itmbgn)
CWKBI 5/95
      EQUIVALENCE ( WORDC, WORDB )
C
      DATA POOL /4HPOOL/
     3,    SCRTCH  /4HSCRA,4HTCH0,4HTCH0/
     4,    NUMBR   /1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H0 /
     5,    WORDB   /4HSEM1,4HBEGN,4H    ,4H    /
     5,    WORDE   /4HBEGN,4HEND /
     6,    IBLNK   /4H    /
     6,    MODX    /   215/
     7,    EXIT    /4HEXIT/
      DATA SUBNAM  /4HXSEM,2H00/
      DATA EQUIV, PURGE /4HEQUI, 4HV   , 4HPURG, 4HE   /
      DATA XEQU , XPUR  /4HXEQU, 4HXPUR/
      DATA XSAV , YCHK  /4HXSAV, 4HXCHK/
!     Set varables
      integer error_id
	  integer nin, nout
      character(80) proj,ft05,ft06,output,infile
C*****
C INITIALIZE MACHINE DEPENDENT CONSTANTS
      CALL BTSTRP
!----------------------------------------------------------------------
! hgs 12/06/2104 - The NASA delivery uses redirected input and out put
!                  so these files are not explicitly opened. I need
!                  to change the stdin and stdout to allow the use of
!                  GDB in the script file. Therefore the following mod-
!                  ifications explicitly open the stdin and stdout files
!                  using the FTN5 and FTN6 ENV set by the script.
!
c
c     open inout and output files
c
      call getenv('PROJ',proj)
	  call getenv('FT05',ft05)
	  call getenv('FT06',ft06)
	  output = trim(proj)//'/'//trim(ft06)
	  error_id = 0
  101 continue
      ifile = nout
      open(nout,file=output,form='formatted',
     1status='unknown',iostat=ierr,err=102)
      go to 103
  102 continue
      error_id = -2
      go to 104
  103 continue
      ifile = nin
	  infile = trim(proj)//'/'//trim(ft05)
      open(nin,file=infile,form='formatted'
     1    ,status='unknown', iostat =ierr,err=106)
      go to 105
  106 continue
      error_id = -1
  104 continue
c
c     open error
c
      write(nout,*) 'Error in opening file =',ifile,' IOSTAT = ',ierr
	  select case(error_id)
	  case(-1)
	    write(nout,'(a)') 'File name: ',infile
	  case(-2)
	    write(nout,'(a)') 'File name: ',output
	  end select
      call pexit ! close down and exit with ierror
	  return
c
 105  continue
      LVAX = MACH.EQ.5
C*****
C EXECUTE PREFACE
C*****
      KSCR= LSHIFT(1,NBPW-4*NBPC)
      CALL TDATE(XX(14))
      CALL CONMSG(WORDB,2,1)
      CALL SEMINT ( 0 )
      ISPERLNK = 1
      WORDB(2) = WORDE(2)
      CALL CONMSG ( WORDB,2,1)
      IPLOT = PLOTF
      IF (PLOTF .LT. 0) PLOTF=1
      IBUF1 = KORSZ(DATABF)-SYSBUF
      GO TO 20
C*****
C RETURN HERE AFTER MODULE HAS EXECUTED
C*****
   10 IF (INOSCR(4).EQ.XSAV.OR.INOSCR(4).EQ.YCHK) GO TO 20
      WORDB(4) = WORDE(2)
C      CALL CONMSG(WORDB,4,0)
      CALL CONMSG(WORDB,4,222222)
   20 IF(NMSG .GT. 0) CALL MSGWRT
      CALL OPEN(*270,POOL,DATABF(IBUF1),2)
C*****
C READ THE OSCAR ENTRY
C*****
   30 CALL READ(*280,*40,POOL,INOSCR,200,1,ERRFLG)
      GO TO 290
   40 IF (INOSCR(6))50,30,30
C*****
C TRY AGAIN IF EXECUTE FLAG IS OFF
C*****
   50 CALL CLOSE(POOL,2)
      TYPECD= ANDF(INOSCR(3),MASK)
C*****
C NOW DETERMINE TYPE OF OSCAR FORMAT
C*****
      IF(TYPECD .GT. 2)  GO TO 200
C*****
C*****
C NOW PROCESSING TYPE O AND F
C*****
   60 MODNO= INOSCR(2)
      FIST(2)= FSTRST
      OPNTR = 7
      ASSIGN 110 TO MM
      FISTNM=101
C*****
C PROCESS FILES IN OSCAR ENTRY.
C*****
   70 J=INOSCR(OPNTR)
      OPNTR=OPNTR+1
      IF(J.EQ.0) GO TO 100
      DO 90 I=1,J
      CALL GNFIST(INOSCR(OPNTR),FISTNM,MODNO)
      IF(MODNO) 60,260,80
   80 OPNTR= OPNTR+ 3
   90 FISTNM=FISTNM+1
  100 GO TO MM,(110,120)
C*****
C SETUP TO PROCESS OUTPUT FILES
C*****
  110 IF(TYPECD.EQ.2) GO TO 120
      ASSIGN 120 TO MM
      FISTNM=201
      GO TO 70
C*****
C PROCESS SCRATCH FILES
C*****
  120 J1= INOSCR(OPNTR)
      IF(J1.EQ.0) GO TO 140
      FISTNM= 301
      SCRTCH(2) = SCRTCH(3)
      LL = 1
      L  = 0
      DO 130 J=1,J1
      L = L + 1
      IF ( L .EQ. 10 ) SCRTCH(2) = KHRFN1(SCRTCH(2),3,NUMBR(LL),1)
      SCRTCH(2) = KHRFN1(SCRTCH(2),4,NUMBR(L),1)
      CALL GNFIST(SCRTCH,FISTNM,MODNO)
      IF ( L .NE. 10 ) GO TO 125
      L  = 0
      LL = LL + 1
  125 IF(MODNO) 60,260,130
  130 FISTNM=FISTNM+1
  140 OPNTR=OPNTR+1
C*****
C NOW PROCESS PARAMETER LIST IN OSCAR
C  PARMN = NO. OF PARAMETERS TO PROCESS
C*****
      PARMN=INOSCR(OPNTR)
      IF(PARMN .EQ. 0)  GO TO 200
      II=1
      OPNTR= OPNTR+ 1
      DO 190 J2=1,PARMN
      IF(INOSCR(OPNTR))170,150,150
C*****
C NOW PROCESS CONSTANT PARAMETER
C*****
  150 PARML=INOSCR(OPNTR)
      OPNTR=OPNTR+1
      DO 160 J3=1,PARML
      PARAM(II)=INOSCR(OPNTR)
      II=II+1
  160 OPNTR=OPNTR+1
      GO TO 190
C*****
C MOVE VARIABLE INTO COMMON VIA VPS TABLE
C*****
  170 VPSX= ANDF(INOSCR(OPNTR),MASK3)
      OPNTR=OPNTR+1
      VPARML=VPS(VPSX-1)
      DO 180 J5=1,VPARML
      PARAM(II)=VPS(VPSX)
      II=II+1
  180 VPSX=VPSX+1
  190 CONTINUE
  200 MODX = RSHIFT(INOSCR(3),16)
C*****
C MODULE IS IN THIS LINK
C PRINT TIME MODULE BEGAN EXECUTION IF FUNCTIONAL MODULE
C*****
  245 WORDB(2) = INOSCR(4)
      WORDB(3) = INOSCR(5)
      IF (INOSCR(4).NE.XEQU.AND.INOSCR(4).NE.XPUR) GO TO 250
      IF (INOSCR(4).NE.XEQU) GO TO 248
      WORDB(2) = EQUIV(1)
      WORDB(3) = EQUIV(2)
      GO TO 250
  248 WORDB(2) = PURGE(1)
      WORDB(3) = PURGE(2)
  250 CALL TMTOGO (KTIME)
      IF (KTIME.LE.0.AND.WORDB(2).NE.EXIT)
     *   CALL MESAGE (-50, 0, WORDB(2))
      IF (INOSCR(4).EQ.XSAV.OR.INOSCR(4).EQ.YCHK) GO TO 1000
      WORDB(1) = IBLNK
      WORDB(4) = WORDE(1)
C
C     EXTRACT DMAP SEQUENCE NUMBER
C
      IDIN  = ANDF(INOSCR(6),MASK)
CWKBIB 5/95
      WRITE( WORDC, 251 ) IDIN
251   FORMAT( I4 )
CWKBIE 5/95
CWKBDB 5/95
C      DO 251  I =1,4
C      ICHR  = IDIN -(IDIN/10)*10 +1
C      L = NBPW-NBPC
C      IF (.NOT.LVAX)  WORDB(1) =
C     *    ORF(RSHIFT(WORDB(1),NBPC),LSHIFT(RSHIFT(NUMBR(ICHR),L),L))
C      IF (LVAX)  WORDB(1)=KHRFN1(WORDB(1),5-I,NUMBR(ICHR),1)
C      IDIN = IDIN/10
C      IF(IDIN .EQ. 0)  GO TO 252
C  251 CONTINUE
CWKBDE 5/95
  252 CONTINUE
C      CALL CONMSG(WORDB,4,0)
      CALL CONMSG(WORDB,4,111111)
      GO TO 1000
C*****
C                   E R R O R   M E S S A G E S
C*****
C MODULE REQUIREMENTS EXCEED AVAILABLE FILES
  260 INOSCR(6) = ANDF(INOSCR(6),MASK)
      CALL MESAGE(-18,INOSCR(6),INOSCR(4))
C
C UNEXPECTED ALTERNATE RETURN TAKEN WHILE ATTEMPTING TO OPEN POOL TAPE.
  270 CONTINUE
      KODE = 270
      GO TO 990
C
C OSCAR FILE POSITIONED INCORRECTLY - HIT EOF.
  280 CONTINUE
      KODE = 280
      GO TO 990
C
C OSCAR RECORD TOO LARGE FOR /OSCENT/
  290 CONTINUE
      KODE = 290
      GO TO 990
C
C LINK SPECIFICATIONS INCORRECT FOR THIS MODULE.
  940 CONTINUE
      WRITE (NOUT,945) WORDB,MODX
  945 FORMAT (/1X,4A4,I9)
      KODE = 940
      GO TO 990
C
C
  990 CONTINUE
      WRITE(NOUT,991) KODE
  991 FORMAT(64H0*** SYSTEM FATAL MESSAGE 1006, LINK DRIVER LOGIC ERROR
     *- CODE =,I4)
      CALL MESAGE(-37,0,SUBNAM)
C**********************************************************************
C EXECUTE MODULE
 1000 CALL SSWTCH ( 2, LDIAG )
C     IF ( LDIAG .NE. 0 .AND. MODX .GT. 14 ) CALL DBMDIA
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1( 940,  940, 2003,  940, 2005, 2006, 2007, 2008, 2009, 2010),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 2060),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109, 2110),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2111, 2112, 2113, 2114, 2115, 2116, 2117, 2118, 2119, 2120),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128, 2129, 2130),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2131, 2132, 2133, 2134, 2135, 2136, 2137, 2138, 2139, 2140),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2141, 2142, 2143, 2144, 2145, 2146, 2147, 2148, 2149, 2150),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2151, 2152, 2153, 2154, 2155, 2156, 2157, 2158, 2159, 2160),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2161, 2162, 2163, 2164, 2165, 2166, 2167, 2168, 2169, 2170),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2171, 2172, 2173, 2174, 2175, 2176, 2177, 2178, 2179, 2180),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2181, 2182, 2183, 2184, 2185, 2186, 2187, 2188, 2189, 2190),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2191, 2192, 2193, 2194, 2195, 2196, 2197, 2198, 2199, 2200),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE. 10 ) GO TO
     1(2201, 2202, 2203, 2204, 2205, 2206, 2207, 2208, 2209, 2210),MODX
      MODX = MODX - 10
      IF ( MODX .GE.  1  .AND. MODX .LE.  7 ) GO TO
     1(2211, 2212, 2213, 2214, 2215, 2216, 2217), MODX
      GO TO 940
 2003 CALL XCHK
      GO TO 10
 2005 CALL XCEI
      GO TO 10
 2006 CALL XCEI
      GO TO 10
 2007 CALL XCEI
      GO TO 10
 2008 CALL XSAVE
      GO TO 10
 2009 CALL XPURGE
      GO TO 10
 2010 CALL XEQUIV
      GO TO 10
 2011 CALL XCEI
      GO TO 10
 2012 CALL XCEI
      GO TO 10
 2013 CALL XCEI
      GO TO 10
 2014 CALL DADD
      GO TO 10
 2015 CALL DADD5
      GO TO 10
 2016 CALL AMG
      GO TO 10
 2017 CALL AMP
      GO TO 10
 2018 CALL APD
      GO TO 10
 2019 CALL BMG
      GO TO 10
 2020 CALL CASE
      GO TO 10
 2021 CALL CYCT1
      GO TO 10
 2022 CALL CYCT2
      GO TO 10
 2023 CALL CEAD
      GO TO 10
 2024 CALL CURV
      GO TO 10
 2025 CONTINUE
      GO TO 10
 2026 CALL DDR
      GO TO 10
 2027 CALL DDR1
      GO TO 10
 2028 CALL DDR2
      GO TO 10
 2029 CALL DDRMM
      GO TO 10
 2030 CALL DDCOMP
      GO TO 10
 2031 CALL DIAGON
      GO TO 10
 2032 CALL DPD
      GO TO 10
 2033 CALL DSCHK
      GO TO 10
 2034 CALL DSMG1
      GO TO 10
 2035 CALL DSMG2
      GO TO 10
 2036 CONTINUE
      GO TO 10
 2037 CALL DUMOD1
      GO TO 10
 2038 CALL DUMOD2
      GO TO 10
 2039 CALL DUMOD3
      GO TO 10
 2040 CALL DUMOD4
      GO TO 10
 2041 CONTINUE
      GO TO 10
 2042 CALL EMA1
      GO TO 10
C         SET LINKNO TO FLAG SUBROUTINE SMA1B TO CALL EMG1B
 2043 LINKNO = LINKNM(8)
      CALL EMG
      LINKNO = LINKNM(1)
      GO TO 10
 2044 CALL FA1
      GO TO 10
 2045 CALL FA2
      GO TO 10
 2046 CALL DFBS
      GO TO 10
 2047 CALL FRLG
      GO TO 10
 2048 CALL FRRD
      GO TO 10
 2049 CONTINUE
      GO TO 10
 2050 CALL GI
      GO TO 10
 2051 CALL GKAD
      GO TO 10
 2052 CALL GKAM
      GO TO 10
 2053 CALL GP1
      GO TO 10
 2054 CALL GP2
      GO TO 10
 2055 CALL GP3
      GO TO 10
 2056 CALL GP4
      GO TO 10
 2057 CALL GPCYC
      GO TO 10
 2058 CALL GPFDR
      GO TO 10
 2059 CALL DUMOD5
      GO TO 10
 2060 CALL GPWG
      GO TO 10
 2061 CONTINUE
      GO TO 10
 2062 CALL INPUT
      GO TO 10
 2063 CALL INPTT1
      GO TO 10
 2064 CALL INPTT2
      GO TO 10
 2065 CALL INPTT3
      GO TO 10
 2066 CALL INPTT4
      GO TO 10
 2067 CALL MATGEN
      GO TO 10
 2068 CALL MATGPR
      GO TO 10
 2069 CALL MATPRN
      GO TO 10
 2070 CALL PRTINT
      GO TO 10
 2071 CALL MCE1
      GO TO 10
 2072 CALL MCE2
      GO TO 10
 2073 CALL MERGE1
      GO TO 10
 2074 CONTINUE
      GO TO 10
 2075 CALL MODA
      GO TO 10
 2076 CALL MODACC
      GO TO 10
 2077 CALL MODB
      GO TO 10
 2078 CALL MODC
      GO TO 10
 2079 CALL DMPYAD
      GO TO 10
 2080 CALL MTRXIN
      GO TO 10
 2081 CALL OFP
      GO TO 10
 2082 CALL OPTPR1
      GO TO 10
 2083 CALL OPTPR2
      GO TO 10
 2084 CONTINUE
      GO TO 10
 2085 CALL OUTPT
      GO TO 10
 2086 CALL OUTPT1
      GO TO 10
 2087 CALL OUTPT2
      GO TO 10
 2088 CALL OUTPT3
      GO TO 10
 2089 CALL OUTPT4
      GO TO 10
 2090 CALL QPARAM
      GO TO 10
 2091 CALL PARAML
      GO TO 10
 2092 CALL QPARMR
      GO TO 10
 2093 CALL PARTN1
      GO TO 10
 2094 CONTINUE
      GO TO 10
 2095 CALL MRED1
      GO TO 10
 2096 CALL MRED2
      GO TO 10
 2097 CALL CMRD2
      GO TO 10
 2098 CALL PLA1
      GO TO 10
 2099 CALL PLA2
      GO TO 10
 2100 CALL PLA3
      GO TO 10
 2101 CALL PLA4
      GO TO 10
 2102 CONTINUE
      GO TO 10
 2103 CALL DPLOT
      GO TO 10
 2104 CALL DPLTST
      GO TO 10
 2105 CALL PLTTRA
      GO TO 10
 2106 CALL PRTMSG
      GO TO 10
 2107 CALL PRTPRM
      GO TO 10
 2108 CALL RANDOM
      GO TO 10
 2109 CALL RBMG1
      GO TO 10
 2110 CALL RBMG2
      GO TO 10
 2111 CALL RBMG3
      GO TO 10
 2112 CALL RBMG4
      GO TO 10
 2113 CONTINUE
      GO TO 10
 2114 CALL REIG
      GO TO 10
 2115 CALL RMG
      GO TO 10
 2116 CALL SCALAR
      GO TO 10
 2117 CALL SCE1
      GO TO 10
 2118 CALL SDR1
      GO TO 10
 2119 CALL SDR2
      GO TO 10
 2120 CALL SDR3
      GO TO 10
 2121 CALL SDRHT
      GO TO 10
 2122 CALL SEEMAT
      GO TO 10
 2123 CONTINUE
      GO TO 10
 2124 CALL SETVAL
      GO TO 10
 2125 CALL SMA1
      GO TO 10
 2126 CALL SMA2
      GO TO 10
 2127 CALL SMA3
      GO TO 10
 2128 CALL SMP1
      GO TO 10
 2129 CALL SMP2
      GO TO 10
 2130 CALL SMPYAD
      GO TO 10
 2131 CALL SOLVE
      GO TO 10
 2132 CONTINUE
      GO TO 10
 2133 CALL SSG1
      GO TO 10
 2134 CALL SSG2
      GO TO 10
 2135 CALL SSG3
      GO TO 10
 2136 CALL SSG4
      GO TO 10
 2137 CALL SSGHT
      GO TO 10
 2138 CALL TA1
      GO TO 10
 2139 CALL TABPCH
      GO TO 10
 2140 CONTINUE
      GO TO 10
 2141 CALL TABFMT
      GO TO 10
 2142 CALL TABPT
      GO TO 10
 2143 CONTINUE
      GO TO 10
 2144 CALL TIMTST
      GO TO 10
 2145 CALL TRD
      GO TO 10
 2146 CALL TRHT
      GO TO 10
 2147 CALL TRLG
      GO TO 10
 2148 CALL DTRANP
      GO TO 10
 2149 CALL DUMERG
      GO TO 10
 2150 CALL DUPART
      GO TO 10
 2151 CALL VDR
      GO TO 10
 2152 CALL VEC
      GO TO 10
 2153 CONTINUE
      GO TO 10
 2154 CALL XYPLOT
      GO TO 10
 2155 CALL XYPRPT
      GO TO 10
 2156 CALL XYTRAN
      GO TO 10
 2157 CONTINUE
      GO TO 10
 2158 CALL COMB1
      GO TO 10
 2159 CALL COMB2
      GO TO 10
 2160 CALL EXIO
      GO TO 10
 2161 CALL RCOVR
      GO TO 10
 2162 CALL EMFLD
      GO TO 10
 2163 CONTINUE
      GO TO 10
 2164 CALL RCOVR3
      GO TO 10
 2165 CALL REDUCE
      GO TO 10
 2166 CALL SGEN
      GO TO 10
 2167 CALL SOFI
      GO TO 10
 2168 CALL SOFO
      GO TO 10
 2169 CALL SOFUT
      GO TO 10
 2170 CALL SUBPH1
      GO TO 10
 2171 CALL PLTMRG
      GO TO 10
 2172 CONTINUE
      GO TO 10
 2173 CALL COPY
      GO TO 10
 2174 CALL SWITCH
      GO TO 10
 2175 CALL MPY3
      GO TO 10
 2176 CALL DDCMPS
      GO TO 10
 2177 CALL LODAPP
      GO TO 10
 2178 CALL GPSTGN
      GO TO 10
 2179 CALL EQMCK
      GO TO 10
 2180 CALL ADR
      GO TO 10
 2181 CALL FRRD2
      GO TO 10
 2182 CALL GUST
      GO TO 10
 2183 CALL IFT
      GO TO 10
 2184 CALL LAMX
      GO TO 10
 2185 CALL EMA
      GO TO 10
 2186 CALL ANISOP
      GO TO 10
 2187 CONTINUE
      GO TO 10
 2188 CALL GENCOS
      GO TO 10
 2189 CALL DDAMAT
      GO TO 10
 2190 CALL DDAMPG
      GO TO 10
 2191 CALL NRLSUM
      GO TO 10
 2192 CALL GENPAR
      GO TO 10
 2193 CALL CASEGE
      GO TO 10
 2194 CALL DESVEL
      GO TO 10
 2195 CALL PROLAT
      GO TO 10
 2196 CALL MAGBDY
      GO TO 10
 2197 CALL COMUGV
      GO TO 10
 2198 CALL FLBMG
      GO TO 10
 2199 CALL GFSMA
      GO TO 10
 2200 CALL TRAIL
      GO TO 10
 2201 CALL SCAN
      GO TO 10
 2202 CONTINUE
      GO TO 10
 2203 CALL PTHBDY
      GO TO 10
 2204 CALL VARIAN
      GO TO 10
 2205 CALL FVRST1
      GO TO 10
 2206 CALL FVRST2
      GO TO 10
 2207 CALL ALG
      GO TO 10
 2208 CALL APDB
      GO TO 10
 2209 CALL PROMPT
      GO TO 10
 2210 CALL OLPLOT
      GO TO 10
 2211 CALL INPTT5
      GO TO 10
 2212 CALL OUTPT5
      GO TO 10
 2213 CONTINUE
      GO TO 10
 2214 CALL QPARMD
      GO TO 10
 2215 CALL GINOFL
      GO TO 10
 2216 CALL DBASE
      GO TO 10
 2217 CALL NORMAL
      GO TO 10
      END
