################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/LSST/M1M3/SS/Gyro/Gyro.cpp \
../src/LSST/M1M3/SS/Gyro/GyroSettings.cpp 

OBJS += \
./src/LSST/M1M3/SS/Gyro/Gyro.o \
./src/LSST/M1M3/SS/Gyro/GyroSettings.o 

CPP_DEPS += \
./src/LSST/M1M3/SS/Gyro/Gyro.d \
./src/LSST/M1M3/SS/Gyro/GyroSettings.d 


# Each subdirectory must supply rules for building sources it contributes
src/LSST/M1M3/SS/Gyro/%.o: ../src/LSST/M1M3/SS/Gyro/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"../../../boost_1_65_1" -I"../../../ts_opensplice/OpenSpliceDDS/V6.4.1/HDE/x86_64.linux/examples/include" -I"../../../ts_opensplice/OpenSpliceDDS/V6.4.1/HDE/x86_64.linux/examples" -I"../../../ts_opensplice/OpenSpliceDDS/V6.4.1/HDE/x86_64.linux/include" -I"../../../ts_opensplice/OpenSpliceDDS/V6.4.1/HDE/x86_64.linux/include/sys" -I"../../../ts_opensplice/OpenSpliceDDS/V6.4.1/HDE/x86_64.linux/include/dcps/C++/SACPP" -I"../../../ts_sal/test/MTM1M3/cpp/src" -I"../../../ts_sal/test/MTM1M3/cpp" -I"../../../ts_sal/test/MTMount/cpp/src" -I"../../../ts_sal/lsstsal/include" -I"../../../ts_sal/test/MTMount/cpp" -I"../../../ts_sal/test" -I../../../ts_sal/test/include -I"../src/LSST/M1M3/SS/DigitalInputOutput" -I"../src/LSST/M1M3/SS/FirmwareUpdate" -I"../src/LSST/M1M3/SS/Accelerometer" -I"../src/LSST/M1M3/SS/AutomaticOperationsController" -I"../src/LSST/M1M3/SS/BusLists" -I"../src/LSST/M1M3/SS/CommandFactory" -I"../src/LSST/M1M3/SS/Displacement" -I"../src/LSST/M1M3/SS/Inclinometer" -I"../src/LSST/M1M3/SS/ForceController" -I"../src/LSST/M1M3/SS/Commands" -I"../src/LSST/M1M3/SS/Context" -I"../src/LSST/M1M3/SS/Controller" -I"../src/LSST/M1M3/SS/Domain" -I"../src/LSST/M1M3/SS/FPGA" -I"../src/LSST/M1M3/SS/Gyro" -I"../src/LSST/M1M3/SS/ILC" -I"../src/LSST/M1M3/SS/Include" -I"../src/LSST/M1M3/SS/Logging" -I"../src/LSST/M1M3/SS/Modbus" -I"../src/LSST/M1M3/SS/Model" -I"../src/LSST/M1M3/SS/PID" -I"../src/LSST/M1M3/SS/PositionController" -I"../src/LSST/M1M3/SS/PowerController" -I"../src/LSST/M1M3/SS/Publisher" -I"../src/LSST/M1M3/SS/SafetyController" -I"../src/LSST/M1M3/SS/Settings" -I"../src/LSST/M1M3/SS/StateFactory" -I"../src/LSST/M1M3/SS/States" -I"../src/LSST/M1M3/SS/Subscriber" -I"../src/LSST/M1M3/SS/Threads" -I"../src/LSST/M1M3/SS/Utility" -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


