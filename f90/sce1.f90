SUBROUTINE sce1
     
!     MODULE 2.6 SCE PARTITIONS KNN,MNN,BNN,AND K4NN
 
!     TO ELIMINATE THE EFFECTS OF SINGLE POINT CONSTRAINTS IF US IS NOT
!     NULL
 
 INTEGER :: us,uset,bnn,bff,un,uf,pvect
 COMMON /patx / n(2),n3,nn(3)
 DATA    un,uf, us / 27,26, 31 /
 DATA    uset , knn,mnn,bnn,k4nn,kff,kfs,kss,mff,bff,k4ff,pvect /  &
     101  , 102,103,104,105 ,201,202,203,204,205,206 ,301   /
 
 CALL upart (uset,pvect,un,uf,us)
 CALL mpart (knn,kff,0,kfs,kss)
 CALL mpart (mnn,mff,0,0,0)
 CALL mpart (bnn,bff,0,0,0)
 CALL mpart (k4nn,k4ff,0,0,0)
 RETURN
END SUBROUTINE sce1
