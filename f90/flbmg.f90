SUBROUTINE flbmg
     
!     DRIVER FOR MODULE FLBMG
 
!     COMPUTES THE HYDROELASTIC AREA FACTOR MATRIX AND THE
!     GRAVITATIONAL STIFFNESS MATRIX.
 
!     THE HYDROELASTIC USET VECTOR IA ALSO BUILT.
 
!     DMAP CALL
 
!        FLBMG    GEOM2,ECT,BGPDT,SIL,MPT,GEOM3,CSTM,USET,EQEXIN/
!                 USETF,USETS,AF,DKGG/S,N,NOGRAV/S,N,NOFREE/S,N,TILT $
 
!     INPUT DATA BLOCKS
 
!        GEOM2  - FLUID ELEMENT BOUNDARY DATA
!        ECT    - ELEMENT CONNECTION TABLE
!        BGPDT  - BASIC GRID POINT DEFINITION TABLE
!        SIL    - SCALAR INDEX LIST
!        MPT    - MATERIAL PROPERTIES TABLE
!        GEOM3  - GRAVITY LOAD DATA
!        CSTM   - COORDINATE SYSTEM TRANSFORMATION MATRICES
!        USET   - DISPLACEMENT SET DEFINITION TABLE
!        EQEXIN - EQUIVALENCE BETWEEN EXTERNAL AND INTERNAL GRID POINTS
 
!     OUTPUT DATA BLOCK
 
!        USETF  - FLUID AND STRUCTURAL POINT SET DEFINITION TABLE
!        USETS  - STRUCTURAL POINT SET DEFINITION TABLE
!        AF     - FLUID AREA FACTOR MATRIX
!        DKGG   - STRUCTURAL GRAVITY STIFFNESS AMTRIX
 
!     PARAMETERS
 
!        NOGRAV - INPUT  - FLAG WHICH SPECIFIES WHETHER GRAVITY
!                          EFFECTS ARE TO BE COMPUTED.
!        NOFREE - OUTPUT - FLAG WHICH SPECIFIES WHETHER A FLUID FREE
!                          SURFACE EXISTS.
!        TILT   - OUTPUT - FREE SURFACE TILT VECTOR USED IN PLOTTING
 
!     USER PRINT OPTIONS
 
!        DIAG 32 - PRINTS HYDROELASTIC SET DEFINITION.
!        DIAG 33 - PRINTS HYDROELASTIC DEGREE OF FREEDOM DEFINITION.
 
 
 LOGICAL :: error
 INTEGER :: geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict   , z1       ,z2(1)    ,sysbuf
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim
 COMMON /flbfil/ geom2    ,ect      ,bgpdt    ,sil      ,mpt      ,  &
     geom3    ,cstm     ,uset     ,eqexin   ,usetf    ,  &
     usets    ,af       ,dkgg     ,fbelm    ,frelm    ,  &
     conect   ,afmat    ,afdict   ,kgmat    ,kgdict
 COMMON /flbptr/ error    ,icore    ,lcore    ,ibgpdt   ,nbgpdt   ,  &
     isil     ,nsil     ,igrav    ,ngrav    ,igrid    ,  &
     ngrid    ,ibuf1    ,ibuf2    ,ibuf3    ,ibuf4    , ibuf5
 COMMON /system/ sysbuf   ,nout
 COMMON /BLANK / nograv   ,nofree   ,tilt(2)
 COMMON /zzzzzz/ z1(1)
 EQUIVALENCE     (z2(1),z1(1))
 
 
!     INITILIZE OPEN CORE FOR ELEMENT MATRIX GENERATION PHASE
 
 error =.false.
 lcore = korsz(z1(1))
 icore = 1
 ibuf1 = lcore - sysbuf - 1
 ibuf2 = ibuf1 - sysbuf
 ibuf3 = ibuf2 - sysbuf
 ibuf4 = ibuf3 - sysbuf
 ibuf5 = ibuf4 - sysbuf
 
!     PROCESS FLUID ELEMENTS ON THE FLUID / STRUCTURE BOUNDARY
!     AND THE FREE SURFACE .
 
 CALL flbelm
 IF (error) GO TO 20
 
!     BUILD THE HYDROELASTIC USET VECTOR
 
 CALL flbset
 IF (error) GO TO 20
 
!     GENERATE THE ELEMENT MATRICES
 
 CALL flbemg
 IF (error) GO TO 20
 
!     INITIALIZE CORE FOR THE MATRIX ASSEMBLY PHASE
 
 lcore = korsz(z2(1))
 icore = 1
 ibuf1 = lcore - sysbuf - 1
 ibuf2 = ibuf1 - sysbuf
 
!     ASSEMBLE THE AREA FACTOR MATRIX
 
 CALL flbema (1)
 
!     IF GRAVITY LOADS - ASSEMBLE THE GRAVITY STIFFNESS MATRIX
 
 IF (nograv < 0) GO TO 10
 CALL flbema (2)
 
!     MODULE COMPLETION
 
 10 CONTINUE
 RETURN
 
!     FATAL ERROR OCCURED DURING PROCESSING - TERMINATE RUN
 
 20 WRITE  (nout,30) uim
 30 FORMAT (a29,' 8000, MODULE FLBMG TERMINATED DUE TO ABOVE ERRORS.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE flbmg
