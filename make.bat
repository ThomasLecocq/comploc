@ECHO OFF
if "%1" == "datetime" (
    gcc datetime.f -c
)

if "%1" == "vzfillin" (
    gcc vzfillin.f -lm -lgfortran  -o vzfillin
)
if "%1" == "deptable" (
    gcc deptable.f -lm -lgfortran  -o deptable
)
if "%1" == "getstlist" (
    gcc getstlist.f -lm -lgfortran  -o getstlist
)
if "%1" == "phase2bed3" (
    gcc phase2bed3.f datetime.o -lm -lgfortran  -o phase2bed3
)
if "%1" == "listbed3" (
    gcc listbed3.f datetime.o -lm -lgfortran  -o listbed3
)
if "%1" == "comploc" (
    gcc comploc.f datetime.o -lm -lgfortran  -o comploc
)


if "%1" == "all" (
    gcc datetime.f -c
    gcc vzfillin.f -lm -lgfortran  -o vzfillin
    gcc deptable.f -lm -lgfortran  -o deptable
    gcc getstlist.f -lm -lgfortran  -o getstlist
    gcc phase2bed3.f datetime.o -lm -lgfortran  -o phase2bed3
    gcc listbed3.f datetime.o -lm -lgfortran  -o listbed3
    gcc comploc.f datetime.o -lm -lgfortran  -o comploc
)

if "%1" == "clean" (
    erase datetime.o
    erase vzfillin.exe
    erase deptable.exe
    erase getstlist.exe
    erase phase2bed3.exe
    erase listbed3.exe
    erase comploc.exe
)








