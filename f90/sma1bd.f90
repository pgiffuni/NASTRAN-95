BLOCK DATA sma1bd
 
 INTEGER :: clsrw    ,clsnrw   ,eor      ,outrw
 DOUBLE PRECISION :: dpdum(514)
 COMMON  /sma1io/  ifcstm   ,ifmpt    ,ifdit    ,idum1    ,  &
     ifecpt   ,igecpt   ,ifgpct   ,iggpct   ,  &
     ifgei    ,iggei    ,ifkgg    ,igkgg    ,  &
     if4gg    ,ig4gg    ,ifgpst   ,iggpst   ,  &
     inrw     ,outrw    ,clsnrw   ,clsrw    ,  &
     neor     ,eor      ,mcbkgg(7),mcb4gg(7)
 
!     SMA1 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON  /sma1bk/  icstm    ,ncstm    ,igpct    ,ngpct    ,  &
     ipoint   ,npoint   ,i6x6k    ,n6x6k    , i6x64    ,n6x64
 COMMON  /sma1dp/  dpdum
 
!     SMA1 PROGRAM CONTROL PARAMETERS
 
 COMMON  /sma1cl/  iopt4    ,k4ggsw   ,npvt     ,left     ,  &
     frowic   ,lrowic   ,nrowsc   ,tnrows   ,  &
     jmax     ,nlinks   ,link(10) ,idetck   , dodet    ,nogoo    ,dummy(200)
 
!     ECPT COMMON BLOCK
 
 COMMON  /sma1et/  ecpt(200)
 
 
 DATA     nlinks / 10 /
 DATA     nogoo  /  0 /
 DATA     ifcstm,ifmpt,ifecpt,ifgpct,ifdit / 101,102,103,104,105 /
 DATA     ifkgg,if4gg,ifgpst / 201,202,203 /
 DATA     inrw,clsrw,clsnrw,eor,neor,outrw / 0,1,2,1,0,1 /
 
END
