SUBROUTINE rbmg1
!*****
! RBMG1 PARTITIONS KAA INTO KLL, KLR AND KRR AND MAA SIMILARLY.
!*****
 
 INTEGER :: uset  ,ua    ,ul    ,ur    ,scr1
 COMMON/bitpos/um    ,uo    ,ur    ,usg   ,usb   ,ul    ,ua  &
     ,uf    ,us    ,un    ,ug    ,ue    ,up
!*****
!     INPUT DATA FILES
!*****
 DATA uset,kaa,maa/101,102,103/
!*****
!     OUTPUT DATA FILES
!*****
 DATA  kll,klr,krr,mll,mlr,mrr/201,202,203,204,205,206/
!*****
!     SCRATCH DATA FILES
!*****
 DATA scr1/301/
!*****
!     PARTITION  KAA INTO KLL,KLR, AND KRR
!     PARTITION  MAA INTO MLL,MLR, AND MRR
!*****
 CALL upart(uset,scr1,ua,ul,ur)
 CALL mpart(kaa,kll,0,klr,krr)
 CALL mpart(maa,mll,0,mlr,mrr)
 RETURN
END SUBROUTINE rbmg1
