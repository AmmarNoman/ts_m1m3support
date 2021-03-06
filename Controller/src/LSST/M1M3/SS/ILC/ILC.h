/*
 * ILC.h
 *
 *  Created on: Oct 1, 2017
 *      Author: ccontaxis
 */

#ifndef ILC_H_
#define ILC_H_

#include <ILCDataTypes.h>
#include <ILCSubnetData.h>
#include <ModbusBuffer.h>
#include <ILCMessageFactory.h>
#include <SetADCChanneOffsetAndSensitivityBusList.h>
#include <ChangeILCModeBusList.h>
#include <ReadBoostValveDCAGainBusList.h>
#include <ReadCalibrationBusList.h>
#include <ReportADCScanRateBusList.h>
#include <ReportDCAIDBusList.h>
#include <ReportDCAStatusBusList.h>
#include <ReportServerIDBusList.h>
#include <ReportServerStatusBusList.h>
#include <ResetBustList.h>
#include <SetADCScanRateBusList.h>
#include <SetBoostValveDCAGainBusList.h>
#include <FreezeSensorBusList.h>
#include <RaisedBusList.h>
#include <ActiveBusList.h>
#include <ILCResponseParser.h>
#include <SAL_m1m3C.h>
#include <FirmwareUpdate.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class M1M3SSPublisher;
class FPGA;
class IBusList;
class ILCApplicationSettings;
class ForceActuatorApplicationSettings;
class ForceActuatorSettings;
class HardpointActuatorApplicationSettings;
class HardpointActuatorSettings;
class HardpointMonitorApplicationSettings;
class PositionController;
class SafetyController;

/*!
 * The ILC class used to communicate with the M1M3's 5 subnets.
 */
class ILC {
private:
	M1M3SSPublisher* publisher;
	FPGA* fpga;
	SafetyController* safetyController;
	ILCSubnetData subnetData;
	ILCMessageFactory ilcMessageFactory;
	ILCResponseParser responseParser;

	SetADCChanneOffsetAndSensitivityBusList busListSetADCChannelOffsetAndSensitivity;
	SetADCScanRateBusList busListSetADCScanRate;
	SetBoostValveDCAGainBusList busListSetBoostValveDCAGains;
	ResetBustList busListReset;
	ReportServerIDBusList busListReportServerID;
	ReportServerStatusBusList busListReportServerStatus;
	ReportADCScanRateBusList busListReportADCScanRate;
	ReadCalibrationBusList busListReadCalibration;
	ReadBoostValveDCAGainBusList busListReadBoostValveDCAGains;
	ReportDCAIDBusList busListReportDCAID;
	ReportDCAStatusBusList busListReportDCAStatus;
	ChangeILCModeBusList busListChangeILCModeDisabled;
	ChangeILCModeBusList busListChangeILCModeEnabled;
	ChangeILCModeBusList busListChangeILCModeStandby;
	ChangeILCModeBusList busListChangeILCModeClearFaults;
	FreezeSensorBusList busListFreezeSensor;
	RaisedBusList busListRaised;
	ActiveBusList busListActive;
	FirmwareUpdate firmwareUpdate;

	HardpointActuatorSettings* hardpointActuatorSettings;
	m1m3_HardpointActuatorDataC* hardpointActuatorData;
	ForceActuatorApplicationSettings* forceActuatorApplicationSettings;
	ForceActuatorSettings* forceActuatorSettings;
	m1m3_ForceActuatorDataC* forceActuatorData;
	m1m3_logevent_HardpointActuatorInfoC* hardpointActuatorInfo;
	PositionController* positionController;

	int8_t hpStepCommand[6];
	int32_t hpSteps[6];

	uint16_t u16Buffer[1];
	ModbusBuffer rxBuffer;

	int32_t controlListToggle;

public:
	ILC(M1M3SSPublisher* publisher, FPGA* fpga, PositionController* positionController, ILCApplicationSettings* ilcApplicationSettings, ForceActuatorApplicationSettings* forceActuatorApplicationSettings, ForceActuatorSettings* forceActuatorSettings, HardpointActuatorApplicationSettings* hardpointActuatorApplicationSettings, HardpointActuatorSettings* hardpointActuatorSettings, HardpointMonitorApplicationSettings* hardpointMonitorApplicationSettings, SafetyController* safetyController);
	virtual ~ILC();

	void programILC(int32_t actuatorId, std::string filePath);
	void modbusTransmit(int32_t actuatorId, int32_t functionCode, int32_t dataLength, int16_t* data);

	void writeCalibrationDataBuffer();
	void writeSetADCScanRateBuffer();
	void writeSetBoostValveDCAGainBuffer();
	void writeResetBuffer();
	void writeReportServerIDBuffer();
	void writeReportServerStatusBuffer();
	void writeReportADCScanRateBuffer();
	void writeReadCalibrationDataBuffer();
	void writeReadBoostValveDCAGainBuffer();
	void writeReportDCAIDBuffer();
	void writeReportDCAStatusBuffer();
	void writeSetModeDisableBuffer();
	void writeSetModeEnableBuffer();
	void writeSetModeStandbyBuffer();
	void writeSetModeClearFaultsBuffer();
	void writeFreezeSensorListBuffer();
	void writeRaisedListBuffer();
	void writeActiveListBuffer();
	void writeControlListBuffer();

	void triggerModbus();

	void waitForSubnet(int32_t subnet, int32_t timeout);
	void waitForAllSubnets(int32_t timeout);

	void read(uint8_t subnet);
	void readAll();
	void flush(uint8_t subnet);
	void flushAll();

	void calculateHPPostion();
	void calculateHPMirrorForces();
	void calculateFAMirrorForces();
	void clearResponses();
	void verifyResponses();

	void publishForceActuatorInfo();
	void publishForceActuatorStatus();
	void publishForceActuatorData();
	void publishHardpointActuatorInfo();
	void publishHardpointStatus();
	void publishHardpointData();
	void publishHardpointMonitorInfo();
	void publishHardpointMonitorStatus();
	void publishHardpointMonitorData();

private:
	uint8_t subnetToRxAddress(uint8_t subnet);
	uint8_t subnetToTxAddress(uint8_t subnet);

	void writeBusList(IBusList* busList);

	void updateHPSteps();
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* ILC_H_ */
