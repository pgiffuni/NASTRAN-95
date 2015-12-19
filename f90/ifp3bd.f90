BLOCK DATA ifp3bd
!IFP3BD
!     B L O C K   D A T A   F O R   I F P 3
 
 
 INTEGER :: FILE          ,iname         ,cdtype
 INTEGER :: axic1         ,cconex        ,forcex
 INTEGER :: force         ,grav          ,load
 INTEGER :: momax         ,moment        ,mpcadd
 INTEGER :: mpcax         ,omitax        ,pointx
 INTEGER :: presax        ,ringax        ,sectax
 INTEGER :: seqgp         ,spcax         ,supax
 INTEGER :: tempax        ,tempd         ,pload
 INTEGER :: mpc           ,spc           ,grid
 INTEGER :: suport        ,neg111        ,t65535
 INTEGER :: temp          ,omit          ,spcadd
 INTEGER :: one           ,zero
 INTEGER :: ctriaa        ,ctrapa
 INTEGER :: rforce
 
 COMMON /ifp3cm /  FILE(6)        ,iname(12)     ,cdtype(50)  &
     ,axic1(3)      ,cconex(3)     ,forcex(3)  &
     ,force(3)      ,grav(3)       ,load(3)  &
     ,momax(3)      ,moment(3)     ,mpcadd(3)  &
     ,mpcax(3)      ,omitax(3)     ,pointx(3)  &
     ,presax(3)     ,ringax(3)     ,sectax(3)  &
     ,seqgp(3)      ,spcax(3)      ,supax(3)  &
     ,tempax(3)     ,tempd(3)      ,pload(3)  &
     ,mpc(3)        ,spc(3)        ,grid(3)  &
     ,suport(3)     ,neg111(3)     ,t65535(3)  &
     ,temp(3)       ,omit(3)       ,spcadd(3)  &
     ,one           ,zero          ,iheadb(96)  &
     ,ctriaa(3)     ,ctrapa(3)     ,iconso ,rforce(3)
 
 DATA one/1/ , zero/0/
 
 DATA FILE  ( 1), FILE  ( 2) / 201   , 208    /
 DATA FILE  ( 3), FILE  ( 4) / 209   , 210    /
 DATA FILE  ( 5), FILE  ( 6) / 301   , 215    /
 
 DATA iname ( 1), iname ( 2) / 4HGEOM, 4H1    /
 DATA iname ( 3), iname ( 4) / 4HGEOM, 4H2    /
 DATA iname ( 5), iname ( 6) / 4HGEOM, 4H3    /
 DATA iname ( 7), iname ( 8) / 4HGEOM, 4H4    /
 DATA iname ( 9), iname (10) / 4HSCRT, 4HCH   /
 DATA iname (11), iname (12) / 4HAXIC, 4H     /
 
 DATA cdtype( 1), cdtype( 2) / 4HAXIC, 4H     /
 DATA cdtype( 3), cdtype( 4) / 4HCCON, 4HEAX  /
 DATA cdtype( 5), cdtype( 6) / 4HFORC, 4HEAX  /
 DATA cdtype( 7), cdtype( 8) / 4HFORC, 4HE    /
 DATA cdtype( 9), cdtype(10) / 4HGRAV, 4H     /
 DATA cdtype(11), cdtype(12) / 4HLOAD, 4H     /
 DATA cdtype(13), cdtype(14) / 4HMOMA, 4HX    /
 DATA cdtype(15), cdtype(16) / 4HMOME, 4HNT   /
 DATA cdtype(17), cdtype(18) / 4HMPCA, 4HDD   /
 DATA cdtype(19), cdtype(20) / 4HMPCA, 4HX    /
 DATA cdtype(21), cdtype(22) / 4HOMIT, 4HAX   /
 DATA cdtype(23), cdtype(24) / 4HPOIN, 4HTAX  /
 DATA cdtype(25), cdtype(26) / 4HPRES, 4HAX   /
 DATA cdtype(27), cdtype(28) / 4HRING, 4HAX   /
 DATA cdtype(29), cdtype(30) / 4HSECT, 4HAX   /
 DATA cdtype(31), cdtype(32) / 4HSEQG, 4HP    /
 DATA cdtype(33), cdtype(34) / 4HSPCA, 4HDD   /
 DATA cdtype(35), cdtype(36) / 4HSPCA, 4HX    /
 DATA cdtype(37), cdtype(38) / 4HSUPA, 4HX    /
 DATA cdtype(39), cdtype(40) / 4HTEMP, 4HAX   /
 DATA cdtype(41), cdtype(42) / 4HTEMP, 4HD    /
 DATA cdtype(43), cdtype(44) / 4HCTRI, 4HAAX  /
 DATA cdtype(45), cdtype(46) / 4HCTRA, 4HPAX  /
 DATA cdtype(47), cdtype(48) / 4HRFOR, 4HCE   /
 
 DATA axic1  /515     ,5      ,0       /
 DATA cconex /8515    ,85     ,0       /
 DATA forcex /2115    ,21     ,0       /
 DATA force  /4201    ,42     ,0       /
 DATA grav   /4401    ,44     ,0       /
 DATA load   /4551    ,61     ,0       /
 DATA momax  /3815    ,38     ,0       /
 DATA moment /4801    ,48     ,0       /
 DATA mpcadd /4891    ,60     ,0       /
 DATA mpcax  /4015    ,40     ,0       /
 DATA omitax /4315    ,43     ,0       /
 DATA pointx /4915    ,49     ,0       /
 DATA presax /5215    ,52     ,0       /
 DATA ringax /5615    ,56     ,0       /
 DATA sectax /6315    ,63     ,0       /
 DATA seqgp  /5301    ,53     ,0       /
 DATA spcax  /6215    ,62     ,0       /
 DATA supax  /6415    ,64     ,0       /
 DATA tempax /6815    ,68     ,0       /
 DATA tempd  /5641    ,65     ,0       /
 DATA pload  /5101    ,51     ,0       /
 DATA mpc    /4901    ,49     ,0       /
 DATA spc    /5501    ,55     ,0       /
 DATA grid   /4501    ,45     ,0       /
 DATA suport /5601    ,56     ,0       /
 DATA temp   /5701    ,57     ,0       /
 DATA omit   /5001    ,50     ,0       /
 DATA spcadd /5491    ,59     ,0       /
 DATA ctriaa /7012    ,70     ,0       /
 DATA ctrapa /7042    ,74     ,0       /
 DATA rforce /5509    ,55     ,0       /
 DATA iconso / 0 /
 DATA neg111 /-1      ,-1     ,-1      /
 DATA t65535/ 65535, 65535, 65535 /
 DATA iheadb / 4HI n ,4HP u ,4HT   ,4HD a ,4HT a ,4H  e  &
     ,4HR r ,4HO r ,4HS   ,4HD e ,4HT e ,4HC t  &
     ,4HE d ,4H  b ,4HY   ,4HI f ,4HP 3 ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H (ax,4HIS-s,4HYMME  &
     ,4HTRIC,4H con,4HICAL,4H she,4HLL d,4HATA  &
     ,4HPROC,4HESSO,4HR-GE,4HNERA,4HTOR),4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H ===  &
     ,4H====,4H====,4H====,4H====,4H====,4H====  &
     ,4H====,4H====,4H====,4H====,4H====,4H====  &
     ,4H====,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H     &
     ,4H    ,4H    ,4H    ,4H    ,4H    ,4H      /
 
END BLOCK DATA 
