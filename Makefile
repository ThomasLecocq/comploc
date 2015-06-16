OBJS1= datetime.o


vzfillin: Makefile vzfillin.f
	f77 -C vzfillin.f -o vzfillin

deptable: Makefile deptable.f
	f77 -C deptable.f -o deptable

getstlist: Makefile getstlist.f
	f77 -C getstlist.f -o getstlist

phase2bed3: Makefile phase2bed3.f $(OBJS1)
	f77 -O phase2bed3.f $(OBJS1) -o phase2bed3

listbed3: Makefile listbed3.f $(OBJS1)
	f77 -O listbed3.f $(OBJS1) -o listbed3

comploc: Makefile comploc.f $(OBJS1)
	f77 -O comploc.f $(OBJS1) -o comploc

