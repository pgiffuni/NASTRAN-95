SUBROUTINE CLOSE ( FILE, iop )
!***************************************************************
!                          NOTICE
 
!     THIS PROGRAM BELONGS TO RPK CORPORATION.  IT IS CONSIDERED
! A TRADE SECRET AND IS NOT TO BE DIVULGED OR USED BY PARTIES
! WHO HAVE NOT RECEIVED WRITTEN AUTHORIZATION FROM RPK.
!***************************************************************
 INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: iop
 INTEGER*2         iunit
 COMMON / dsunit / iunit( 220 )
 COMMON / system / isysbf, dum1( 77 ), idiag, dum2( 21 )
 
 
 NAME   = FILE
 iocode = iop
 iretrn = 77
 CALL dsgefl
 IF ( ifilex == 0 ) GO TO 50
 iretrn = 0
 IF ( IAND( idiag,2**14 ) /= 0 ) CALL dsmsg( 2 )
 IF ( iocode /= 1 ) GO TO 20
 IF ( iprvop == 0 ) GO TO 10
 CALL dsefwr
 IF ( ( indclr-indbas ) == 5 ) GO TO 5
 ibase( indbas+4 ) = indclr - indbas + 1
 CALL dbmmgr( 4 )
 5     CALL dsxfsz
 10    CONTINUE
 CALL dbmmgr( 2 )
 nblock  = 1
 indclr  = indbas + 5
 indcbp  = indclr
 GO TO 40
 20    IF ( iprvop == 0 ) GO TO 30
 CALL dsefwr
 ibase( indbas+4 ) = indclr - indbas + 1
! SAVE INDBAS TO ALLOW DSBRC1 TO CORRECTLY BACKSPACE FILE OPENNED FOR WRITE
 isave = indbas
 CALL dbmmgr( 4 )
 CALL dsxfsz
 indbas = isave
 IF ( iocode /= -2 ) CALL dsbrc1
!      CALL DSGNCL
 CALL dbmmgr( 2 )
 GO TO 40
 30    IF ( indcbp == indclr ) GO TO 35
 CALL dsskrc
 35    CONTINUE
 CALL dbmmgr( 2 )
 40    CALL dssdcb
 fcb( 2,ifilex ) = 0
 fcb(12,ifilex ) = 0
 IF ( NAME < 101 .OR. NAME > 320 ) GO TO 50
 iunit( NAME-100 ) = 0
 50    RETURN
!***************************************************************
!                          NOTICE
 
!     THIS PROGRAM BELONGS TO RPK CORPORATION.  IT IS CONSIDERED
! A TRADE SECRET AND IS NOT TO BE DIVULGED OR USED BY PARTIES
! WHO HAVE NOT RECEIVED WRITTEN AUTHORIZATION FROM RPK.
!***************************************************************
END SUBROUTINE CLOSE
