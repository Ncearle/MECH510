################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Prog_asgn_3/Prog_asgn_3/Source.cpp 

OBJS += \
./Prog_asgn_3/Prog_asgn_3/Source.o 

CPP_DEPS += \
./Prog_asgn_3/Prog_asgn_3/Source.d 


# Each subdirectory must supply rules for building sources it contributes
Prog_asgn_3/Prog_asgn_3/%.o: ../Prog_asgn_3/Prog_asgn_3/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


