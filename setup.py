import os
import sys
#
memc=sys.path[0]
# memc=memc.replace("/utils","/")
memc=memc+"/"
#
if len(sys.argv)==2:
	if sys.argv[1]=='-d' or sys.argv[1]=='-debug'\
		or sys.argv[1]=='debug' or sys.argv[1]=='d':
		# os.system("ln -s " + memc + "main.c .")
		# os.system("ln -s " + memc + "main_mpi.c .")
		# os.system("cp "+ memc +"Makefile .")
		os.system("mkdir utils")
		os.system("mkdir hosts")
		os.system("mkdir src")
		os.system("mkdir includes")
		os.system("mkdir mains")
		#
		os.system("ln -s "+memc+"src/*.cpp src/")
		os.system("ln -s "+memc+"src/*.f90 src/")
		os.system("ln -s "+memc+"includes/*.h includes/")
		os.system("ln -s "+memc+"hosts/* hosts/")
		os.system("ln -s "+memc+"utils/*.py utils/")
		os.system("ln -s "+memc+"utils/*.cpp utils/")
		os.system("ln -s "+memc+"mains/*.cpp mains/")
		# os.system("cp "+memc+"input/*.h5 input/")
		os.system("cp "+memc+"para_file.in .")
		os.system("cp "+memc+"Makefile .")
		os.system("cp -r"+memc+"/conf/ .")
else:
	curr_dir = os.getcwd()
	os.chdir(memc)
	os.system("make")
	os.chdir(curr_dir)
	os.system("cp "+memc+"bin/* .")
	# os.system("cp "+memc+"exe_memc .")
	os.system("cp -r"+memc+"conf/ .")
	os.system("cp "+memc+"para_file.in .")
	os.system("mkdir utils")
	os.system("ln -s "+memc+"utils/*.py utils/")
	os.system("ln -s "+memc+"utils/*.cpp utils/")