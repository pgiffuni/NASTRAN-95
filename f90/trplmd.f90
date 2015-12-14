SUBROUTINE trplmd (gmat,dmat,bmat,bmat1,bmat2,mattyp,jcor,wtk)
     
!     ROUTINE TO PERFORM THE TRIPLE MULTIPLY AT EACH INTEGRATION
!     POINT FOR THE QUAD4 ELEMENT.
!     DIFFERENT PATHS ARE TAKEN BASED ON THE FOLLOWING CRITERIA -
!      1- ELEMENT BEING A MEMBRANE ONLY, OR BENDING ONLY, OR BOTH
!         MEMBRANE AND BENDING ELEMENT.
!      2- THE MATERIAL PROPERTIES BEING ISOTROPIC OR NOT.
!      3- THE MACHINE THIS CODE IS RUNNING ON. (TENTATIVE)
 
 
 DOUBLE PRECISION, INTENT(IN)             :: gmat(10,10)
 DOUBLE PRECISION, INTENT(IN)             :: dmat(7,7)
 DOUBLE PRECISION, INTENT(OUT)            :: bmat(240)
 DOUBLE PRECISION, INTENT(IN)             :: bmat1(1)
 DOUBLE PRECISION, INTENT(IN)             :: bmat2(1)
 INTEGER, INTENT(IN OUT)                  :: mattyp
 INTEGER, INTENT(IN OUT)                  :: jcor
 DOUBLE PRECISION, INTENT(IN)             :: wtk
 DOUBLE PRECISION :: akgg
 
 DOUBLE PRECISION :: dbm(240),dmat1(3,3),dmat2(4,4)
 
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth
 
 COMMON /terms / membrn,bendng,shrflx,mbcoup,norpth
 COMMON /zzzzzz/ akgg(1)
 COMMON /trplm / ndof,ibot,iptx1,iptx2,ipty1,ipty2
 
!*****
!     INITIALIZE
!*****
 nd1 = ndof
 nd2 = nd1 * 2
 nd3 = nd1 * 3
 nd4 = nd1 * 4
 nd5 = nd1 * 5
 nd6 = nd1 * 6
 nd7 = nd1 * 7
 nd8 = nd1 * 8
 nd9 = nd1 * 9
 nda = nd1 * 10
 IF (.NOT.norpth) GO TO 500
!*****
!    ALL MIDS ARE THE SAME AND THERE IS NO COUPLING.
!    IF THE MATERIAL IS ISOTROPIC, PERFORM THE 1ST MUTIPLY EXPLICITLY.
!    IF NOT, USE GMMATD. IN EITHER CASE, THE 2ND MULTIPLY USES GMMATD.
!*****
 DO  i=1,nd1
   bmat(i+nd2) = bmat2(i+ibot     )
   bmat(i+nd3) = bmat1(i+ipty1    )
   bmat(i+nd4) = bmat1(i+ipty2    )
   bmat(i+nd5) = bmat1(i+iptx1+nd1)
   bmat(i+nd6) = bmat1(i+iptx2+nd1)
 END DO
 
 IF (mattyp /= 1) GO TO 300
 DO  i=1,nd1
   dbm (i    ) = dmat(1,1)*bmat(i    ) + dmat(1,2)*bmat(i+nd1)
   dbm (i+nd1) = dmat(2,1)*bmat(i    ) + dmat(2,2)*bmat(i+nd1)
   dbm (i+nd2) = dmat(3,3)*bmat(i+nd2)
   dbm (i+nd3) = dmat(4,4)*bmat(i+nd3) + dmat(4,5)*bmat(i+nd4)
   dbm (i+nd4) = dmat(5,4)*bmat(i+nd3) + dmat(5,5)*bmat(i+nd4)
   dbm (i+nd5) = dmat(6,6)*bmat(i+nd5) + dmat(6,7)*bmat(i+nd6)
   dbm (i+nd6) = dmat(7,6)*bmat(i+nd5) + dmat(7,7)*bmat(i+nd6)
 END DO
 GO TO 400
 
 300 CALL gmmatd (dmat,7,7,0,bmat,7,nd1,0,dbm)
 
 400 DO  i=1,nd7
   bmat(i) = bmat(i)*wtk
 END DO
 CALL gmmatd (bmat,7,nd1,-1,dbm,7,nd1,0,akgg(jcor))
 RETURN
!*****
!     MIDS ARE NOT THE SAME. CHECK FOR MEMBRANE ONLY AND BENDING ONLY
!     CASES AND BRANCH APPROPRIATELY. IF BOTH ARE THERE, CONTINUE.
!*****
 500 IF (.NOT.bendng) GO TO 800
 IF (.NOT.membrn) GO TO 1200
 DO  i=1,nd1
   bmat(i+nd2) = bmat2(i+ibot     )
   bmat(i+nd5) = bmat2(i+ibot+nd1 )
   bmat(i+nd6) = bmat1(i+ipty1    )
   bmat(i+nd7) = bmat1(i+ipty2    )
   bmat(i+nd8) = bmat1(i+iptx1+nd1)
   bmat(i+nd9) = bmat1(i+iptx2+nd1)
 END DO
 
 CALL gmmatd (gmat,10,10,0,bmat,10,nd1,0,dbm)
 
 DO  i=1,nda
   bmat(i) = bmat(i)*wtk
 END DO
 CALL gmmatd (bmat,10,nd1,-1,dbm,10,nd1,0,akgg(jcor))
 RETURN
!*****
!     MEMBRANE ONLY ELEMENT. ONLY THE FIRST 3X3 OF GMAT AND THE FIRST
!     3 ROWS OF BMAT ARE MULTIPLIED.
!*****
 800 DO  i=1,nd1
   bmat(i+nd2) = bmat2(i+ibot)
 END DO
 
 IF (mattyp /= 1) GO TO 950
 DO  i=1,nd1
   dbm (i    ) = gmat(1,1)*bmat(i    ) + gmat(1,2)*bmat(i+nd1)
   dbm (i+nd1) = gmat(2,1)*bmat(i    ) + gmat(2,2)*bmat(i+nd1)
   dbm (i+nd2) = gmat(3,3)*bmat(i+nd2)
 END DO
 GO TO 1050
 
 950 DO  i=1,3
   DO  j=1,3
     dmat1(i,j) = gmat(i,j)
   END DO
 END DO
 CALL gmmatd (dmat1,3,3,0,bmat(1),3,nd1,0,dbm(1))
 
 1050 DO  i=1,nd3
   bmat(i) = bmat(i)*wtk
 END DO
 CALL gmmatd (bmat,3,nd1,-1,dbm,3,nd1,0,akgg(jcor))
 RETURN
!*****
!     BENDING ONLY ELEMENT. THE FIRST 3 ROWS AND COLUMNS OF GMAT AND
!     THE FIRST 3 ROWS OF BMAT WILL BE EXCLUDED FROM MULTIPLICATIONS.
!*****
 1200 DO  i=1,nd1
   bmat(i+nd6) = bmat1(i+ipty1    )
   bmat(i+nd7) = bmat1(i+ipty2    )
   bmat(i+nd8) = bmat1(i+iptx1+nd1)
   bmat(i+nd9) = bmat1(i+iptx2+nd1)
 END DO
 
 DO  i=1,3
   DO  j=1,3
     dmat1(i,j) = gmat(i+3,j+3)
   END DO
 END DO
 DO  i=1,4
   DO  j=1,4
     dmat2(i,j) = gmat(i+6,j+6)
   END DO
 END DO
 
 CALL gmmatd (dmat1,3,3,0,bmat(nd3+1),3,nd1,0,dbm(1    ))
 CALL gmmatd (dmat2,4,4,0,bmat(nd6+1),4,nd1,0,dbm(nd3+1))
 
 DO  i=nd3+1,nda
   bmat(i) = bmat(i)*wtk
 END DO
 CALL gmmatd (bmat(nd3+1),7,nd1,-1,dbm,7,nd1,0,akgg(jcor))
 RETURN
 
END SUBROUTINE trplmd
