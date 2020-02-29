# SRwithMklProj
 Project to create executable using the SRlibSimple library and the Mkl pardiso solver
 
 linux version:
 go into either linuxDebug or linuxRelease directory and make all there.
 note that both the makefile and all ".mk" files are needed because the makefile references them
 Please doublecheck the -L directive in the makefile, and the -I directive in subdir.mk, there may be relative paths.
 -L should have the path to your SRlibSimple library as well as the path to the intel mkl library.
 -I should point to the corresponding directories containing the .h files.
 I was unsuccesful at linking to the Intel Mkl Solver. I used the same mkl libraries that work on windows:
 (mkl_core, mkl_sequential, and mkl_intel_lp64). These were also recommended by Intel's link line advisor
 (https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor). On windows this works fine, on linux there are lots of undefineds
 referenced from the 3 libraries liste.
 If anyone figures this out please post a comment.
 The executable will link with the solver suppressed.
 To try out mkl, 
 a. comment out the line "#define NOSOLVER" in SRlinux.h
 b. download the free version of intel mkl to your machine (https://software.intel.com/en-us/mkl/choose-download)
 


 Windows Visual Studio version that was last modified using VS2013. That and any newer version should work to compile this project.
 It is in C++. Doubleclick on the ".sln" file to open the project.
 If using a newer version than 2013, you may get a message about updating to the latest VS version, to which you should say yes.
 Choose a configuration at the top (debug or release). Next to that it should say x64. You can change that to win32 but executation
 is faster in x64 mode. Now hit build to compile. A subfolder is created x64/debug or x64/release depending on the configuration.
 Inside that folder will be the executable SRwithMkl.exe
