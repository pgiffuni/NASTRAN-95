BLOCK DATA sma2bd
!SMA2BD
 
 INTEGER :: clsrw,  clsnrw, eor,   outrw
 COMMON /sma2io/ ifcstm, ifmpt,  ifdit, idum1, ifecpt, igecpt,  &
     ifgpct, iggpct, idum2, idum3, ifmgg,  igmgg,  &
     ifbgg,  igbgg,  idum4, idum5, inrw,   outrw,  &
     clsnrw, clsrw,  neor,  eor, mcbmgg(7),mcbbgg(7)
 
!     SMA2 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON /sma2bk/ icstm,  ncstm,  igpct, ngpct, ipoint, npoint,  &
     i6x6m,  n6x6m,  i6x6b, n6x6b
 
!     SMA2 PROGRAM CONTROL PARAMETERS
 
 COMMON /sma2cl/ iopt4,  bggind, npvt,  left,  frowic, lrowic,  &
     nrowsc, tnrows, jmax,  nlinks,link(10),nogo, dummy(202)
 
!     ECPT COMMON BLOCK
 
 COMMON /sma2et/ ecpt(200)
 
 DATA    nlinks/ 10 /
 DATA    nogo  / 0  /
 DATA    ifcstm, ifmpt,ifecpt,ifgpct,ifdit  / 101,102,103,104,105 /
 DATA    ifmgg , ifbgg       /     201,202  /
 DATA    inrw  , clsrw,clsnrw,eor,neor,outrw/ 0,1,2,1,0,1 /
END BLOCK DATA
