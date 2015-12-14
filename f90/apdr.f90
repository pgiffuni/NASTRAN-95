SUBROUTINE apdr (FILE,z,core,in,out,wr,buf,TYPE)
     
 
    INTEGER, INTENT(IN OUT)                  :: FILE
    INTEGER, INTENT(IN OUT)                  :: z(1)
    INTEGER, INTENT(OUT)                     :: core
    INTEGER, INTENT(OUT)                     :: in
    INTEGER, INTENT(IN OUT)                  :: out
    INTEGER, INTENT(OUT)                     :: wr
    INTEGER, INTENT(IN OUT)                  :: buf
    INTEGER, INTENT(IN OUT)                  :: TYPE(3)
    INTEGER :: flag, NAME(2)
    DATA    NAME /4HAPD ,4HR   /
 
    wr = 0
    in = 0
    CALL locate (*20,z(buf),TYPE,flag)
    in = out + 1
    CALL READ (*90,*10,FILE,z(in),core,0,wr)
    GO TO 80
10  out  = in + wr - 1
20  core = core - wr
    RETURN
 
80  CALL mesage (-3,FILE,NAME)
90  CALL mesage (-2,FILE,NAME)
    GO TO 20

END SUBROUTINE apdr
