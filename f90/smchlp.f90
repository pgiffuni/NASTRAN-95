SUBROUTINE smchlp
     
! SMCHLP WRITES THE CONTENTS OF THE KEY PARAMETERS IN THE SYMMETRIC
! DECOMPOSITION COMMON BLOCKS
 
 INCLUDE 'SMCOMX.COM'
 WRITE ( nout, 9001 ) ncol   , ierror  , ivwrds , maxnac  &
     ,   nspill , maxinlop, idbase , idbmax  &
     ,   ibuf1  , ibuf2   , opnscr , ioloop  &
     ,   lascol , krow    , krows  , krown ,   kridx  , kridxn  , jridxn , jrow  &
     ,   jrows  , jrown   , jridx  , jvidx  &
     ,   irow1  , irown   , kfrcol , klscol  &
     ,   klsrow , iol     , iil    , ktype  &
     ,   iskip  , indexv  , kcol   , maxncol  &
     ,   memfre , memcol1 , memlck , memlas ,   memcoln, ispill  , nbandw , nvterm
 9001  FORMAT(//  &
     ,  ' NCOL   =',i9,' IERROR  =',i9,' IVWRDS =',i9,' MAXNAC =',i9  &
     ,/,' NSPILL =',i9,' MAXINLOP=',i9,' IDBASE =',i9,' IDBMAX =',i9  &
     ,/,' IBUF1  =',i9,' IBUF2   =',i9,' OPNSCR =',l9,' IOLOOP =',i9  &
     ,/,' LASCOL =',i9,' KROW    =',i9,' KROWS  =',i9,' KROWN  =',i9  &
     ,/,' KRIDX  =',i9,' KRIDXN  =',i9,' JRIDXN =',i9,' JROW   =',i9  &
     ,/,' JROWS  =',i9,' JROWN   =',i9,' JRIDX  =',i9,' JVIDX  =',i9  &
     ,/,' IROW1  =',i9,' IROWN   =',i9,' KFRCOL =',i9,' KLSCOL =',i9  &
     ,/,' KLSROW =',i9,' IOL     =',i9,' IIL    =',i9,' KTYPE  =',i9  &
     ,/,' ISKIP  =',i9,' INDEXV  =',i9,' KCOL   =',i9,' MAXNCOL=',i9  &
     ,/,' MEMFRE =',i9,' MEMCOL1 =',i9,' MEMLCK =',i9,' MEMLAS =',i9  &
     ,/,' MEMCOLN=',i9,' ISPILL  =',i9,' NBANDW =',i9,' NVTERM =',i9 )
 WRITE ( nout, 9002 ) isysbf, isprec
 9002  FORMAT( /,' ISYSBF =',i9,' ISPREC   =',i9 )
 WRITE ( nout, 9003 ) mblk, moblk
 9003  FORMAT(/,' MBLK (INPUT MATRIX STRING BLOCK)=',/,3(5I10,/)  &
     ,/,' MOBLK (OUTPUT MATRIX STRING BLOCK)=',/,3(5I10,/))
 WRITE ( nout, 9004 ) lcore, power, mindd, chlsky, iscr1
 9004  FORMAT( /,' LCORE  =',i9,' POWER   =',i9,' MINDD   =',e16.8  &
     ,/,' CHLSKY =',i9,' ISCR1   =',i9 )
 WRITE ( nout, 9005 ) mcb, lll
 9005  FORMAT('  INPUT MATRIX MCB=',/,7I8, /,    ' OUTPUT MATRIX MCB=',/,7I8 )
 RETURN
END SUBROUTINE smchlp
