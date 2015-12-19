BLOCK DATA dpdcbd
!DPDCBD
! BLOCK DATA PROGRAM FOR THE DYNAMICS POOL DISTRIBUTOR
!*****
 
 INTEGER :: gpl   ,sil   ,uset  ,usetd ,gpld  ,sild  ,dpool  &
     ,dlt   ,frl   ,tfl   ,trl   ,psdl  ,eed   ,scr1  &
     ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  ,buf3  &
     ,buf4  ,eqdyn ,sdt   ,epoint,seqep ,eigc  ,eigb  &
     ,loads ,dload ,freq1 ,freq  ,tic   ,tstep ,tf ,psd   ,eigr
 
 DIMENSION buf(24)   ,epoint(2)    ,seqep(2)     ,mcb(7)  &
     ,nam(2)    ,loads(32)    ,dload(2)     ,freq1(2)  &
     ,freq(2)   ,eigc(2)      ,eigb(2)      ,nolin(21)  &
     ,tic(2)    ,tstep(2)     ,tf(2)        ,psd(2) ,msg(3)    ,eigr(2)
 
 COMMON/dpdcom/dpool ,gpl   ,sil   ,uset  ,gpld  ,sild  ,usetd  &
     ,dlt   ,frl   ,nlft  ,tfl   ,trl   ,psdl  ,eed  &
     ,scr1  ,scr2  ,scr3  ,scr4  ,buf   ,buf1  ,buf2  &
     ,buf3  ,buf4  ,epoint,seqep ,l     ,kn    ,neqdyn  &
     ,loads ,dload ,freq1 ,freq  ,nolin ,nogo  &
     ,msg   ,tic   ,tstep ,tf    ,psd   ,eigr  ,eigb  &
     ,eigc  ,mcb   ,nam   ,eqdyn ,sdt   ,ineq
!*****
! INPUT FILES
!*****
 DATA    dpool/101/  ,gpl/102/     ,sil/103/     ,uset/104/
!*****
! OUTPUT FILES
!*****
 DATA    gpld  /201/ ,sild  /202/  ,usetd /203/  ,tfl   /204/  &
     ,dlt   /205/ ,psdl  /206/  ,frl   /207/  ,nlft  /208/  &
     ,trl   /209/ ,eed   /210/  ,eqdyn /211/  ,sdt   /212/
!*****
! SCRATCH FILES
!*****
 DATA    scr1/301/   ,scr2/302/    ,scr3/303/    ,scr4/304/
!*****
! DATA DEFINING INPUT CARDS
!*****
 DATA epoint  /   707,     7/ ,seqep   /  5707,    57/  &
     ,loads   /    27,    17,     0,     0 ,    37,    18,     0,     0  &
     ,    77,    19,     0,     0 ,  5107,    51,     6,     0  &
     ,  5207,    52,     6,     0 ,  7107,    71,     5,     0  &
     ,  7207,    72,    10,     0 ,     0,     0,     0,     0/  &
     ,dload   /    57,     5/ ,freq1   /  1007,    10/  &
     ,freq    /  1307,    13/
 DATA nolin   /  3107,    31,     8 ,  3207,    32,     8  &
     ,  3307,    33,     8 ,  3407,    34,     8  &
     ,  3507,    35,    16 ,  3607,    36,     5  &
     ,  3707,    37,     8/
 DATA tic     /  6607,    66/ ,tstep   /  8307,    83/  &
     ,tf      /  6207,    62/ ,eigr    /   307,     3/  &
     ,eigb    /   107,     1/ ,eigc    /   207,     2/
!*****
! MISC DATA
!*****
 DATA mcb     /   7*0/ ,nam     /4HDPD ,4H    /

END BLOCK DATA
