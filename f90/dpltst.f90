SUBROUTINE dpltst
     
 IMPLICIT INTEGER (a-z)
 INTEGER :: errttl(32)
 COMMON /BLANK / ngp,nsets,skp1(8),pcdb,eqexin,ect,skp2(7),  &
     merr,parm,gpset,elset,skp3(6),mset,pect
 COMMON /system/ bufsiz
 COMMON /zzzzzz/ x(1)
 DATA    outrew, rew / 1,1 /
 DATA    errttl/ 8*2H  ,4HERRO,4HR me,4HSSAG,4HES f,4HROM ,  &
     4HTHE ,4HPLOT,4H set,4H def,4HINIT,4HION ,  &
     4HMODU,4HLE (,4HPLTS,4HET) ,9*1H  /
 
 nsets  = 0
 pcdb   = 101
 eqexin = 102
 ect    = 103
 ept    = 104
 merr   = 201
 parm   = 202
 gpset  = 203
 elset  = 204
 mset   = 301
 pect   = 302
 CALL totape (1,x(1))
 
 x(1) = eqexin
 CALL rdtrl (x)
 i2   = 2
 i3   = 3
 ngp  = x(i2) - x(i3)
 i1   = korsz(x) - bufsiz + 1
 CALL gopen (merr,x(i1),outrew)
 CALL WRITE (merr,-4,1,0)
 CALL WRITE (merr,errttl,32,0)
 CALL setinp
 IF (nsets /= 0) GO TO 150
 nsets = -1
 GO TO 200
 150 i1 = nsets + 1
 i2 = i1 + ngp
 i3 = korsz(x) - 4*bufsiz + 1
 CALL comect (x(i2),i3-i2)
 CALL cnstrc (x(i1),x(i2),x(i3),i3-i2)
 
 200 CALL clstab (merr,rew)
 RETURN
END SUBROUTINE dpltst
