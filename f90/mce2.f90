SUBROUTINE mce2
     
!     MCE2 PARTITIONS KGG INTO KNNB, KMNB AND KMMB THEN COMPUTES
 
!     KNN = KNNB + GM(T)*KMNB + KMNB(T)*GM + GM(T)*KMMB*GM
 
!     SIMILAR OPERATIONS ARE PERFORMED ON MGG, BGG AND K4GG IF THE
!     MATRIX HAS NOT BEEN PURGED.
 
 INTEGER :: scr1  ,scr2  ,scr6  ,uset  ,mcb(7),gm    ,ug    ,  &
     un    ,um    ,bgg   ,bnnb  ,bmnb  ,bnn   ,bmmb
 COMMON /bitpos/ um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua    ,  &
     uf    ,us    ,un    ,ug    ,ue    ,up
 
!     INPUT AND OUTPUT FILES
 
 DATA    uset  , gm , kgg, mgg, bgg, k4gg,   knn, mnn, bnn, k4nn /  &
     101   , 102, 103, 104, 105, 106 ,   201, 202, 203, 204  /
 
!     SCRATCH FILES
 
 DATA    scr1  , scr2, scr6 / 301, 302, 306 /
 DATA    knnb  , kmnb, kmmb / 303, 304, 305 /
 DATA    mnnb  , mmnb, mmmb / 303, 304, 305 /
 DATA    bnnb  , bmnb, bmmb / 303, 304, 305 /
 DATA    k4nnb , k4mnb,k4mmb/ 303, 304, 305 /
 
!     ARITHMETIC TYPES
 
!     PARTITION KGG INTO KNNB, KMNB, AND KMMB
 
 CALL upart (uset,scr1,ug,un,um)
 CALL mpart (kgg,knnb,kmnb,0,kmmb)
 
!     COMPUTE KNN
 
 CALL elim (knnb,kmnb,kmmb,gm,knn,scr1,scr2,scr6)
 
!     TEST TO SEE IF MGG IS PRESENT
 
 mcb(1) = mgg
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 110
 
!     IF MGG PRESENT, PARTITION INTO MNNB, MMNB, AND MMMB
!     THEN COMPUTE MNN
 
 CALL upart (uset,scr1,ug,un,um)
 CALL mpart (mgg,mnnb,mmnb,0,mmmb)
 CALL elim  (mnnb,mmnb,mmmb,gm,mnn,scr1,scr2,scr6)
 
!     TEST TO SEE IF BGG IS PRESENT
 
 110 mcb(1) = bgg
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) GO TO 130
 
!     IF BGG PRESENT, PARTITION INTO BNNB, BMNB, AND BMMB
!     THEN COMPUTE BNN
 
 CALL upart (uset,scr1,ug,un,um)
 CALL mpart (bgg,bnnb,bmnb,0,bmmb)
 CALL elim  (bnnb,bmnb,bmmb,gm,bnn,scr1,scr2,scr6)
 
!     TEST TO SEE IF K4GG IS PRESENT
 
 130 mcb(1) = k4gg
 CALL rdtrl (mcb)
 IF (mcb(1) < 0) RETURN
 
!     IF K4GG IS PRESENT, PARTITION INTO K4NNB, K4MNB, AND K4MMB
!     THEN COMPUTE K4NN
 
 CALL upart (uset,scr1,ug,un,um)
 CALL mpart (k4gg,k4nnb,k4mnb,0,k4mmb)
 CALL elim  (k4nnb,k4mnb,k4mmb,gm,k4nn,scr1,scr2,scr6)
 RETURN
END SUBROUTINE mce2
