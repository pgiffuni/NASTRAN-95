SUBROUTINE xsem00
! **********************************************************************
! THE PURPOSE OF THIS ROUTINE IS TO EXECUTE THE PREFACE AND THEN TO
! EXECUTE MODULES ACCORDING TO THE DMAP.  THE DMAP IS READ FROM THE
! OSCAR.  FOR EACH MODULE TO BE EXECUTED, THE FIST AND XVPS ARE SETUP.
 
!WKBD 5/95
!     INTEGER ORF
 INTEGER :: andf ,databf,errflg,fist  ,fistnm,fstrst,opntr  &
     ,parml,param ,parmn ,pool  ,rshift,scrtch  &
     ,vps  ,vparml,typecd,vpsx  ,wordb ,worde ,plotf,EXIT  ,sysbuf,subnam(2)
 INTEGER :: equiv(2), purge(2), xequ, xpur, xsav, ychk
!WKBI 5/95
 CHARACTER (LEN=4) :: wordc
 
 LOGICAL :: lvax
 
 DIMENSION scrtch(3),wordb(4),worde(2),numbr(10)
 
 COMMON/machin/mach
 COMMON/sem   /mask  ,mask2 ,mask3 ,linknm(15)   &
     /system/sysbuf,xx(20),linkno,xxx(16),nbpc,nbpw,ncpw,xxxx(53) ,isperlnk  &
      /xlink /lxlink,maxlnk,mxlink(1)  &
      /xfist /fist(2)  &
      /xpfist/fstrst  &
      /oscent/inoscr(200)  &
      /zzzzzz/databf(1)  &
      /BLANK /param(100)  &
      /xvps  /vps(1)  &
      /msgx  /nmsg
 
 EQUIVALENCE (xx(1),nout), (xx(3),nin) ! input file number
 EQUIVALENCE (xx(19),plotf)
 EQUIVALENCE (xx(17),itmbgn)
!WKBI 5/95
 EQUIVALENCE ( wordc, wordb )
 
 DATA pool /4HPOOL/ ,    scrtch  /4HSCRA,4HTCH0,4HTCH0/  &
     ,    numbr   /1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,1H0 /  &
     ,    wordb   /4HSEM1,4HBEGN,4H    ,4H    / ,    worde   /4HBEGN,4HEND /  &
     ,    iblnk   /4H    / ,    modx    /   215/  &
     ,    EXIT    /4HEXIT/
 DATA subnam  /4HXSEM,2H00/
 DATA equiv, purge /4HEQUI, 4HV   , 4HPURG, 4HE   /
 DATA xequ , xpur  /4HXEQU, 4HXPUR/
 DATA xsav , ychk  /4HXSAV, 4HXCHK/
!     Set varables
 INTEGER :: error_id
 INTEGER :: nin, nout
 CHARACTER (LEN=1) :: 80) proj,ft05,ft06,output,infile
!*****
! INITIALIZE MACHINE DEPENDENT CONSTANTS
 CALL btstrp
!----------------------------------------------------------------------
! hgs 12/06/2104 - The NASA delivery uses redirected input and out put
!                  so these files are not explicitly opened. I need
!                  to change the stdin and stdout to allow the use of
!                  GDB in the script file. Therefore the following mod-
!                  ifications explicitly open the stdin and stdout files
!                  using the FTN5 and FTN6 ENV set by the script.
!
 
!     open inout and output files
 
 CALL getenv('PROJ',proj)
 CALL getenv('FT05',ft05)
 CALL getenv('FT06',ft06)
 output = trim(proj)//'/'//trim(ft06)
 error_id = 0
 101 CONTINUE
 ifile = nout
 OPEN(nout,FILE=output,FORM='formatted', STATUS='unknown',IOSTAT=ierr,ERR=102)
 GO TO 103
 102 CONTINUE
 error_id = -2
 GO TO 104
 103 CONTINUE
 ifile = nin
 infile = trim(proj)//'/'//trim(ft05)
 OPEN(nin,FILE=infile,FORM='formatted'  &
     ,STATUS='unknown', IOSTAT =ierr,ERR=106)
 GO TO 105
 106 CONTINUE
 error_id = -1
 104 CONTINUE
 
!     open error
 
 WRITE(nout,*) 'Error in opening file =',ifile,' IOSTAT = ',ierr
 select case(error_id)
 case(-1)
 WRITE(nout,'(a)') 'File name: ',infile
 case(-2)
 WRITE(nout,'(a)') 'File name: ',output
END select
CALL pexit ! close down and exit with ierror
RETURN

105  CONTINUE
lvax = mach == 5
!*****
! EXECUTE PREFACE
!*****
kscr= lshift(1,nbpw-4*nbpc)
CALL tdate(xx(14))
CALL conmsg(wordb,2,1)
CALL semint ( 0 )
isperlnk = 1
wordb(2) = worde(2)
CALL conmsg ( wordb,2,1)
iplot = plotf
IF (plotf < 0) plotf=1
ibuf1 = korsz(databf)-sysbuf
GO TO 20
!*****
! RETURN HERE AFTER MODULE HAS EXECUTED
!*****
10 IF (inoscr(4) == xsav.OR.inoscr(4) == ychk) GO TO 20
wordb(4) = worde(2)
!      CALL CONMSG(WORDB,4,0)
CALL conmsg(wordb,4,222222)
20 IF(nmsg > 0) CALL msgwrt
CALL OPEN(*270,pool,databf(ibuf1),2)
!*****
! READ THE OSCAR ENTRY
!*****
30 CALL READ(*280,*40,pool,inoscr,200,1,errflg)
GO TO 290
40 IF (inoscr(6) < 0) THEN
  GO TO    50
ELSE
  GO TO    30
END IF
!*****
! TRY AGAIN IF EXECUTE FLAG IS OFF
!*****
50 CALL CLOSE(pool,2)
typecd= andf(inoscr(3),mask)
!*****
! NOW DETERMINE TYPE OF OSCAR FORMAT
!*****
IF(typecd > 2)  GO TO 200
!*****
!*****
! NOW PROCESSING TYPE O AND F
!*****
60 modno= inoscr(2)
fist(2)= fstrst
opntr = 7
ASSIGN 110 TO mm
fistnm=101
!*****
! PROCESS FILES IN OSCAR ENTRY.
!*****
70 j=inoscr(opntr)
opntr=opntr+1
IF(j == 0) GO TO 100
DO  i=1,j
  CALL gnfist(inoscr(opntr),fistnm,modno)
  IF(modno < 0) THEN
    GO TO    60
  ELSE IF (modno == 0) THEN
    GO TO   260
  END IF
  80 opntr= opntr+ 3
  fistnm=fistnm+1
END DO
100 GO TO mm,(110,120)
!*****
! SETUP TO PROCESS OUTPUT FILES
!*****
110 IF(typecd == 2) GO TO 120
ASSIGN 120 TO mm
fistnm=201
GO TO 70
!*****
! PROCESS SCRATCH FILES
!*****
120 j1= inoscr(opntr)
IF(j1 == 0) GO TO 140
fistnm= 301
scrtch(2) = scrtch(3)
ll = 1
l  = 0
DO  j=1,j1
  l = l + 1
  IF ( l == 10 ) scrtch(2) = khrfn1(scrtch(2),3,numbr(ll),1)
  scrtch(2) = khrfn1(scrtch(2),4,numbr(l),1)
  CALL gnfist(scrtch,fistnm,modno)
  IF ( l /= 10 ) GO TO 125
  l  = 0
  ll = ll + 1
  125 IF(modno < 0) THEN
    GO TO    60
  ELSE IF (modno == 0) THEN
    GO TO   260
  ELSE
    GO TO   130
  END IF
  fistnm=fistnm+1
END DO
140 opntr=opntr+1
!*****
! NOW PROCESS PARAMETER LIST IN OSCAR
!  PARMN = NO. OF PARAMETERS TO PROCESS
!*****
parmn=inoscr(opntr)
IF(parmn == 0)  GO TO 200
ii=1
opntr= opntr+ 1
DO  j2=1,parmn
  IF(inoscr(opntr) < 0) THEN
    GO TO   170
  END IF
!*****
! NOW PROCESS CONSTANT PARAMETER
!*****
  150 parml=inoscr(opntr)
  opntr=opntr+1
  DO  j3=1,parml
    param(ii)=inoscr(opntr)
    ii=ii+1
    opntr=opntr+1
  END DO
  CYCLE
!*****
! MOVE VARIABLE INTO COMMON VIA VPS TABLE
!*****
  170 vpsx= andf(inoscr(opntr),mask3)
  opntr=opntr+1
  vparml=vps(vpsx-1)
  DO  j5=1,vparml
    param(ii)=vps(vpsx)
    ii=ii+1
    vpsx=vpsx+1
  END DO
END DO
200 modx = rshift(inoscr(3),16)
!*****
! MODULE IS IN THIS LINK
! PRINT TIME MODULE BEGAN EXECUTION IF FUNCTIONAL MODULE
!*****
245 wordb(2) = inoscr(4)
wordb(3) = inoscr(5)
IF (inoscr(4) /= xequ.AND.inoscr(4) /= xpur) GO TO 250
IF (inoscr(4) /= xequ) GO TO 248
wordb(2) = equiv(1)
wordb(3) = equiv(2)
GO TO 250
248 wordb(2) = purge(1)
wordb(3) = purge(2)
250 CALL tmtogo (ktime)
IF (ktime <= 0.AND.wordb(2) /= EXIT) CALL mesage (-50, 0, wordb(2))
IF (inoscr(4) == xsav.OR.inoscr(4) == ychk) GO TO 1000
wordb(1) = iblnk
wordb(4) = worde(1)

!     EXTRACT DMAP SEQUENCE NUMBER

idin  = andf(inoscr(6),mask)
!WKBIB 5/95
WRITE( wordc, 251 ) idin
251   FORMAT( i4 )
!WKBIE 5/95
!WKBDB 5/95
!      DO 251  I =1,4
!      ICHR  = IDIN -(IDIN/10)*10 +1
!      L = NBPW-NBPC
!      IF (.NOT.LVAX)  WORDB(1) =
!     *    ORF(RSHIFT(WORDB(1),NBPC),LSHIFT(RSHIFT(NUMBR(ICHR),L),L))
!      IF (LVAX)  WORDB(1)=KHRFN1(WORDB(1),5-I,NUMBR(ICHR),1)
!      IDIN = IDIN/10
!      IF(IDIN .EQ. 0)  GO TO 252
!  251 CONTINUE
!WKBDE 5/95
252 CONTINUE
!      CALL CONMSG(WORDB,4,0)
CALL conmsg(wordb,4,111111)
GO TO 1000
!*****
!                   E R R O R   M E S S A G E S
!*****
! MODULE REQUIREMENTS EXCEED AVAILABLE FILES
260 inoscr(6) = andf(inoscr(6),mask)
CALL mesage(-18,inoscr(6),inoscr(4))

! UNEXPECTED ALTERNATE RETURN TAKEN WHILE ATTEMPTING TO OPEN POOL TAPE.
270 CONTINUE
kode = 270
GO TO 990

! OSCAR FILE POSITIONED INCORRECTLY - HIT EOF.
280 CONTINUE
kode = 280
GO TO 990

! OSCAR RECORD TOO LARGE FOR /OSCENT/
290 CONTINUE
kode = 290
GO TO 990

! LINK SPECIFICATIONS INCORRECT FOR THIS MODULE.
940 CONTINUE
WRITE (nout,945) wordb,modx
945 FORMAT (/1X,4A4,i9)
kode = 940
GO TO 990


990 CONTINUE
WRITE(nout,991) kode
991 FORMAT(64H0*** system fatal message 1006, link driver logic error  &
    - code =,i4)
CALL mesage(-37,0,subnam)
!**********************************************************************
! EXECUTE MODULE
1000 CALL sswtch ( 2, ldiag )
!     IF ( LDIAG .NE. 0 .AND. MODX .GT. 14 ) CALL DBMDIA
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    ( 940,  940, 2003,  940, 2005, 2006, 2007, 2008, 2009, 2010),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2031, 2032, 2033, 2034, 2035, 2036, 2037, 2038, 2039, 2040),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2041, 2042, 2043, 2044, 2045, 2046, 2047, 2048, 2049, 2050),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2051, 2052, 2053, 2054, 2055, 2056, 2057, 2058, 2059, 2060),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2061, 2062, 2063, 2064, 2065, 2066, 2067, 2068, 2069, 2070),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2071, 2072, 2073, 2074, 2075, 2076, 2077, 2078, 2079, 2080),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2081, 2082, 2083, 2084, 2085, 2086, 2087, 2088, 2089, 2090),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2091, 2092, 2093, 2094, 2095, 2096, 2097, 2098, 2099, 2100),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2101, 2102, 2103, 2104, 2105, 2106, 2107, 2108, 2109, 2110),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2111, 2112, 2113, 2114, 2115, 2116, 2117, 2118, 2119, 2120),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2121, 2122, 2123, 2124, 2125, 2126, 2127, 2128, 2129, 2130),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2131, 2132, 2133, 2134, 2135, 2136, 2137, 2138, 2139, 2140),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2141, 2142, 2143, 2144, 2145, 2146, 2147, 2148, 2149, 2150),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2151, 2152, 2153, 2154, 2155, 2156, 2157, 2158, 2159, 2160),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2161, 2162, 2163, 2164, 2165, 2166, 2167, 2168, 2169, 2170),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2171, 2172, 2173, 2174, 2175, 2176, 2177, 2178, 2179, 2180),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2181, 2182, 2183, 2184, 2185, 2186, 2187, 2188, 2189, 2190),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2191, 2192, 2193, 2194, 2195, 2196, 2197, 2198, 2199, 2200),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <= 10 ) GO TO  &
    (2201, 2202, 2203, 2204, 2205, 2206, 2207, 2208, 2209, 2210),modx
modx = modx - 10
IF ( modx >=  1  .AND. modx <=  7 ) GO TO  &
    (2211, 2212, 2213, 2214, 2215, 2216, 2217), modx
GO TO 940
2003 CALL xchk
GO TO 10
2005 CALL xcei
GO TO 10
2006 CALL xcei
GO TO 10
2007 CALL xcei
GO TO 10
2008 CALL xsave
GO TO 10
2009 CALL xpurge
GO TO 10
2010 CALL xequiv
GO TO 10
2011 CALL xcei
GO TO 10
2012 CALL xcei
GO TO 10
2013 CALL xcei
GO TO 10
2014 CALL dadd
GO TO 10
2015 CALL dadd5
GO TO 10
2016 CALL amg
GO TO 10
2017 CALL amp
GO TO 10
2018 CALL apd
GO TO 10
2019 CALL bmg
GO TO 10
2020 CALL case
GO TO 10
2021 CALL cyct1
GO TO 10
2022 CALL cyct2
GO TO 10
2023 CALL cead
GO TO 10
2024 CALL curv
GO TO 10
2025 CONTINUE
GO TO 10
2026 CALL ddr
GO TO 10
2027 CALL ddr1
GO TO 10
2028 CALL ddr2
GO TO 10
2029 CALL ddrmm
GO TO 10
2030 CALL ddcomp
GO TO 10
2031 CALL diagon
GO TO 10
2032 CALL dpd
GO TO 10
2033 CALL dschk
GO TO 10
2034 CALL dsmg1
GO TO 10
2035 CALL dsmg2
GO TO 10
2036 CONTINUE
GO TO 10
2037 CALL dumod1
GO TO 10
2038 CALL dumod2
GO TO 10
2039 CALL dumod3
GO TO 10
2040 CALL dumod4
GO TO 10
2041 CONTINUE
GO TO 10
2042 CALL ema1
GO TO 10
!         SET LINKNO TO FLAG SUBROUTINE SMA1B TO CALL EMG1B
2043 linkno = linknm(8)
CALL emg
linkno = linknm(1)
GO TO 10
2044 CALL fa1
GO TO 10
2045 CALL fa2
GO TO 10
2046 CALL dfbs
GO TO 10
2047 CALL frlg
GO TO 10
2048 CALL frrd
GO TO 10
2049 CONTINUE
GO TO 10
2050 CALL gi
GO TO 10
2051 CALL gkad
GO TO 10
2052 CALL gkam
GO TO 10
2053 CALL gp1
GO TO 10
2054 CALL gp2
GO TO 10
2055 CALL gp3
GO TO 10
2056 CALL gp4
GO TO 10
2057 CALL gpcyc
GO TO 10
2058 CALL gpfdr
GO TO 10
2059 CALL dumod5
GO TO 10
2060 CALL gpwg
GO TO 10
2061 CONTINUE
GO TO 10
2062 CALL INPUT
GO TO 10
2063 CALL inptt1
GO TO 10
2064 CALL inptt2
GO TO 10
2065 CALL inptt3
GO TO 10
2066 CALL inptt4
GO TO 10
2067 CALL matgen
GO TO 10
2068 CALL matgpr
GO TO 10
2069 CALL matprn
GO TO 10
2070 CALL prtint
GO TO 10
2071 CALL mce1
GO TO 10
2072 CALL mce2
GO TO 10
2073 CALL merge1
GO TO 10
2074 CONTINUE
GO TO 10
2075 CALL moda
GO TO 10
2076 CALL modacc
GO TO 10
2077 CALL modb
GO TO 10
2078 CALL modc
GO TO 10
2079 CALL dmpyad
GO TO 10
2080 CALL mtrxin
GO TO 10
2081 CALL ofp
GO TO 10
2082 CALL optpr1
GO TO 10
2083 CALL optpr2
GO TO 10
2084 CONTINUE
GO TO 10
2085 CALL outpt
GO TO 10
2086 CALL outpt1
GO TO 10
2087 CALL outpt2
GO TO 10
2088 CALL outpt3
GO TO 10
2089 CALL outpt4
GO TO 10
2090 CALL qparam
GO TO 10
2091 CALL paraml
GO TO 10
2092 CALL qparmr
GO TO 10
2093 CALL partn1
GO TO 10
2094 CONTINUE
GO TO 10
2095 CALL mred1
GO TO 10
2096 CALL mred2
GO TO 10
2097 CALL cmrd2
GO TO 10
2098 CALL pla1
GO TO 10
2099 CALL pla2
GO TO 10
2100 CALL pla3
GO TO 10
2101 CALL pla4
GO TO 10
2102 CONTINUE
GO TO 10
2103 CALL dplot
GO TO 10
2104 CALL dpltst
GO TO 10
2105 CALL plttra
GO TO 10
2106 CALL prtmsg
GO TO 10
2107 CALL prtprm
GO TO 10
2108 CALL random
GO TO 10
2109 CALL rbmg1
GO TO 10
2110 CALL rbmg2
GO TO 10
2111 CALL rbmg3
GO TO 10
2112 CALL rbmg4
GO TO 10
2113 CONTINUE
GO TO 10
2114 CALL reig
GO TO 10
2115 CALL rmg
GO TO 10
2116 CALL scalar
GO TO 10
2117 CALL sce1
GO TO 10
2118 CALL sdr1
GO TO 10
2119 CALL sdr2
GO TO 10
2120 CALL sdr3
GO TO 10
2121 CALL sdrht
GO TO 10
2122 CALL seemat
GO TO 10
2123 CONTINUE
GO TO 10
2124 CALL setval
GO TO 10
2125 CALL sma1
GO TO 10
2126 CALL sma2
GO TO 10
2127 CALL sma3
GO TO 10
2128 CALL smp1
GO TO 10
2129 CALL smp2
GO TO 10
2130 CALL smpyad
GO TO 10
2131 CALL solve
GO TO 10
2132 CONTINUE
GO TO 10
2133 CALL ssg1
GO TO 10
2134 CALL ssg2
GO TO 10
2135 CALL ssg3
GO TO 10
2136 CALL ssg4
GO TO 10
2137 CALL ssght
GO TO 10
2138 CALL ta1
GO TO 10
2139 CALL tabpch
GO TO 10
2140 CONTINUE
GO TO 10
2141 CALL tabfmt
GO TO 10
2142 CALL tabpt
GO TO 10
2143 CONTINUE
GO TO 10
2144 CALL timtst
GO TO 10
2145 CALL trd
GO TO 10
2146 CALL trht
GO TO 10
2147 CALL trlg
GO TO 10
2148 CALL dtranp
GO TO 10
2149 CALL dumerg
GO TO 10
2150 CALL dupart
GO TO 10
2151 CALL vdr
GO TO 10
2152 CALL vec
GO TO 10
2153 CONTINUE
GO TO 10
2154 CALL xyplot
GO TO 10
2155 CALL xyprpt
GO TO 10
2156 CALL xytran
GO TO 10
2157 CONTINUE
GO TO 10
2158 CALL comb1
GO TO 10
2159 CALL comb2
GO TO 10
2160 CALL exio
GO TO 10
2161 CALL rcovr
GO TO 10
2162 CALL emfld
GO TO 10
2163 CONTINUE
GO TO 10
2164 CALL rcovr3
GO TO 10
2165 CALL reduce
GO TO 10
2166 CALL sgen
GO TO 10
2167 CALL sofi
GO TO 10
2168 CALL sofo
GO TO 10
2169 CALL sofut
GO TO 10
2170 CALL subph1
GO TO 10
2171 CALL pltmrg
GO TO 10
2172 CONTINUE
GO TO 10
2173 CALL copy
GO TO 10
2174 CALL switch
GO TO 10
2175 CALL mpy3
GO TO 10
2176 CALL ddcmps
GO TO 10
2177 CALL lodapp
GO TO 10
2178 CALL gpstgn
GO TO 10
2179 CALL eqmck
GO TO 10
2180 CALL adr
GO TO 10
2181 CALL frrd2
GO TO 10
2182 CALL gust
GO TO 10
2183 CALL ift
GO TO 10
2184 CALL lamx
GO TO 10
2185 CALL ema
GO TO 10
2186 CALL anisop
GO TO 10
2187 CONTINUE
GO TO 10
2188 CALL gencos
GO TO 10
2189 CALL ddamat
GO TO 10
2190 CALL ddampg
GO TO 10
2191 CALL nrlsum
GO TO 10
2192 CALL genpar
GO TO 10
2193 CALL casege
GO TO 10
2194 CALL desvel
GO TO 10
2195 CALL prolat
GO TO 10
2196 CALL magbdy
GO TO 10
2197 CALL comugv
GO TO 10
2198 CALL flbmg
GO TO 10
2199 CALL gfsma
GO TO 10
2200 CALL trail
GO TO 10
2201 CALL scan
GO TO 10
2202 CONTINUE
GO TO 10
2203 CALL pthbdy
GO TO 10
2204 CALL varian
GO TO 10
2205 CALL fvrst1
GO TO 10
2206 CALL fvrst2
GO TO 10
2207 CALL alg
GO TO 10
2208 CALL apdb
GO TO 10
2209 CALL prompt
GO TO 10
2210 CALL olplot
GO TO 10
2211 CALL inptt5
GO TO 10
2212 CALL outpt5
GO TO 10
2213 CONTINUE
GO TO 10
2214 CALL qparmd
GO TO 10
2215 CALL ginofl
GO TO 10
2216 CALL dbase
GO TO 10
2217 CALL normal
GO TO 10
END SUBROUTINE xsem00
