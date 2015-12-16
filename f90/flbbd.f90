BLOCK DATA flbbd
!FLBBD
!     FLBBD - BLOCK DATA FOR MODULE FLBMG
 
 INTEGER :: geom2    ,ect      ,bgpdt    ,sil      ,geom3    ,  &
     cstm     ,uset     ,eqexin   ,usetf    ,usets    ,  &
     af       ,dkgg     ,fbelm    ,frelm    ,conect   ,  &
     afmat    ,afdict   ,kgmat    ,kgdict
 
!     GINO FILES
 
 COMMON /flbfil/ geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict
 
!     INPUT DATA BLOCKS
 
 DATA            geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin             /  &
     101      ,102      ,103      ,104      ,105      ,  &
     106      ,107      ,108      ,109                /
 
!     OUTPUT DATA BLOCKS
 
 DATA            usetf    ,usets    ,af       ,dkgg               /  &
     201      ,202      ,203      ,204                /
 
!     INTERNAL SCRATCH FILES
 
 DATA            fbelm    ,frelm    ,conect   ,afmat    ,afdict   ,  &
     kgmat    ,kgdict                                 /  &
     301      ,302      ,303      ,304      ,305      ,  &
     306      ,307                                    /
 
END
