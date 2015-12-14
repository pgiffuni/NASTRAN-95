SUBROUTINE fndgrd( isub , icomp , igrid , ip , ic , n )
     
 
 INTEGER, INTENT(IN)                      :: isub
 INTEGER, INTENT(IN)                      :: icomp
 INTEGER, INTENT(IN OUT)                  :: igrid
 INTEGER, INTENT(IN OUT)                  :: ip(6)
 INTEGER, INTENT(IN OUT)                  :: ic(6)
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER :: aaa(2),scsfil,buf3,z,score
 COMMON/cmb001/junk(3),scsfil
 COMMON/cmb002/junk1(2),buf3,junk2(2),score,lcore
 COMMON/cmbfnd/  inam(2),ierr
 COMMON/zzzzzz/z(1)
 DATA aaa/ 4HFNDG,4HRD   /
 
 CALL OPEN(*2001,scsfil,z(buf3),0)
 nfil = isub-1
 CALL skpfil( scsfil , nfil )
 nrec = icomp - 1
 IF( nrec == 0 ) GO TO 3
 DO  i=1,nrec
   CALL fwdrec(*2002,scsfil)
 END DO
 3     CALL READ(*2002,*2,scsfil,z(score),lcore,1,nwd)
 GO TO 2004
 2     CONTINUE
 CALL gridip( igrid , score , nwd , ip , ic , n , z , lloc )
 CALL CLOSE( scsfil , 1 )
 RETURN
 2001  CALL mesage( -1 , scsfil , aaa )
 2002  CALL mesage( -2 , scsfil , aaa )
 2004  CALL mesage( -8 , scsfil , aaa )
 RETURN
END SUBROUTINE fndgrd
