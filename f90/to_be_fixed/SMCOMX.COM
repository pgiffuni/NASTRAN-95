DOUBLE PRECISION :: dajjr   ,dajji   ,ddr     ,ddc     ,mindd
INTEGER :: power   ,sturm   ,chlsky
LOGICAL :: opnscr
COMMON  /smcomx/  ncol    ,ierror  ,ivwrds  ,maxnac  ,nspill,  &
                  maxinlop,idbase  ,idbmax  ,ibuf1   ,ibuf2 ,  &
                  opnscr  ,ioloop  ,iiloop  ,lascol  ,krow  ,  &
                  krows   ,krown   ,kridx   ,kridxn  ,jridxn,  &
                  jrow    ,jrows   ,jrown   ,jridx   ,jvidx ,  &
                  irow1   ,irown   ,kfrcol  ,klscol  ,klsrow,  &
                  iol     ,iil     ,ktype   ,iskip   ,indexv,  &
                  indexvd ,jcol    ,karows  ,mxrecl  ,nvterm,  &
                  kcol    ,maxncol ,memfre  ,memcol1 ,memlck,  &
                  memlas  ,memcoln ,ispill  ,kprec   ,nbandw,  &
                  maxnar  ,mblk(15),moblk(15)
COMMON  /smcomy/  dajjr   ,dajji   ,ajjr    ,ajji
COMMON  /system/  isysbf  ,nout    ,dum1(37),nbpw    ,dum2(14),&
                  isprec
COMMON  /sfact /  mcb(7)  ,lll(7)  ,dbc(7)  ,iscr1   ,iscr2 ,  &
                  lcore   ,ddr     ,ddc     ,power   ,scr3  ,  &
                  mindd   ,chlsky
COMMON  /sturmx/  sturm   ,shftpt  ,KEEP    ,ptshft  ,nr
