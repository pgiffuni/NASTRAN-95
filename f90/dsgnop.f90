SUBROUTINE dsgnop
     INCLUDE 'DSIOF.COM'
 INCLUDE 'GINOX.COM'
 INCLUDE 'XNSTRN.COM'
 CHARACTER (LEN=4) :: cbuff( 3 )
 EQUIVALENCE     (cbuff,ibase)
 idsn    = ifilex
 IF ( iocode >= 2 ) GO TO 10
 IF ( iocode /= 1 ) GO TO 30
 CALL dsrlse
 GO TO 30
 10   inext   = IAND( mdsfcb( 3, idsn ), maskh2 )
 IF ( inext == 0 ) GO TO 30
 itest =  fcb( 6, idsn )
 IF ( nblock <= itest ) GO TO 30
 idsn    = inext
 GO TO 10
 30   iop = MOD ( iocode,2 )
 mdsfcb( 2,ifilex ) = idsn
 mdsfcb( 1,idsn )   = ior( mdsfcb( 1,idsn ), maskh2 )
 40   CONTINUE
 CALL dsopen( mdsnam( idsn ), idsn, iocode)
 RETURN
END SUBROUTINE dsgnop
