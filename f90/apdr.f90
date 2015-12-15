SUBROUTINE apdr (FILE,z,core,in,out,wr,buf,TYPE)
     
 
    INTEGER, INTENT(IN OUT)                  :: file
    INTEGER, INTENT(IN OUT)                  :: z(1)
    INTEGER, INTENT(OUT)                     :: core
    INTEGER, INTENT(OUT)                     :: in
    INTEGER, INTENT(IN OUT)                  :: out
    INTEGER, INTENT(OUT)                     :: wr
    INTEGER, INTENT(IN OUT)                  :: buf
    INTEGER, INTENT(IN OUT)                  :: type(3)
    INTEGER :: flag, NAME(2)
    DATA    NAME /4HAPD ,4HR   /
 
    wr = 0
    in = 0
    CALL locate (*20,z(buf),type,flag)
    in = out + 1
    CALL READ (*90,*10,file,z(in),core,0,wr)
    GO TO 80
10  out  = in + wr - 1
20  core = core - wr
    RETURN
 
80  CALL mesage (-3,file,name)
90  CALL mesage (-2,file,name)
    GO TO 20

END SUBROUTINE apdr
