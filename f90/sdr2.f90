SUBROUTINE sdr2
     
!     SDR2 IS THE EXECUTIVE CONTROL PROGRAM FOR THE SDR2 MODULE.
 
 INTEGER :: any   ,loads ,displ ,vel   ,acc   ,spcf  ,plots , casecc
COMMON /sdr2x4/ nam(2),END   ,mset  ,icb(7),ocb(7),mcb(7),dtype(8)  &
    ,               icstm ,ncstm ,ivec  ,ivecn ,temp  ,deform,FILE  ,  &
    buf1  ,buf2  ,buf3  ,buf4  ,buf5  ,any   ,all   ,  &
    tloads,eldef ,symflg,branch,ktype ,loads ,spcf  ,  &
    displ ,vel   ,acc   ,stress,force ,kwdest,kwdedt,  &
    kwdgpt,kwdcc ,nrigds,sta(2),rei(2),ds0(2),ds1(2),  &
    frq(2),trn(2),bk0(2),bk1(2),cei(2),pla(22)      ,  &
    nrings,nharms,axic  ,knset ,isopl ,strspt,ddrmm
COMMON /sdr2x2/ casecc
COMMON /system/ sysbuf,opte  ,nogo  ,intap ,mpcn  ,spcn  ,method,  &
    loadnn,symm  ,stftmp,page  ,line  ,tline ,maxlin, date(3),time  ,echo  ,plots

!     EXECUTE THE PHASES OF SDR2.

casecc = 101
CALL sdr2aa
CALL sdr2a
IF (any /= 0) CALL sdr2b
k = loads + spcf + displ + vel + acc + plots
IF (k   /= 0) CALL sdr2c
IF (any /= 0) CALL sdr2d
RETURN
END SUBROUTINE sdr2
