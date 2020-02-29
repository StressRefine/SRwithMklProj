################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../SRPostfF06.cpp \
../SRanalysis.cpp \
../SRerrorCheck.cpp \
../SRinput.cpp \
../SRinputUtilities.cpp \
../SRmklUtil.cpp \
../SRpardiso.cpp \
../SRpardisoSimple.cpp \
../SRpostProcess.cpp \
../SRsolver.cpp \
../SRwithMkl.cpp 

OBJS += \
./SRPostfF06.o \
./SRanalysis.o \
./SRerrorCheck.o \
./SRinput.o \
./SRinputUtilities.o \
./SRmklUtil.o \
./SRpardiso.o \
./SRpardisoSimple.o \
./SRpostProcess.o \
./SRsolver.o \
./SRwithMkl.o 

CPP_DEPS += \
./SRPostfF06.d \
./SRanalysis.d \
./SRerrorCheck.d \
./SRinput.d \
./SRinputUtilities.d \
./SRmklUtil.d \
./SRpardiso.d \
./SRpardisoSimple.d \
./SRpostProcess.d \
./SRsolver.d \
./SRwithMkl.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/rich/opt/newEclipseWS/SRlibSimple -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


