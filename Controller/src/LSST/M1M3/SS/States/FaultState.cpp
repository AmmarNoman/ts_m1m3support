/*
 * FaultState.cpp
 *
 *  Created on: Oct 26, 2017
 *      Author: ccontaxis
 */

#include <FaultState.h>
#include <InterlockController.h>
#include <Model.h>
#include <ILC.h>
#include <Displacement.h>
#include <Inclinometer.h>
#include <PowerController.h>
#include <SafetyController.h>
#include <ForceController.h>
#include <unistd.h>
#include <M1M3SSPublisher.h>
#include <Accelerometer.h>
#include <Log.h>

namespace LSST {
namespace M1M3 {
namespace SS {

FaultState::FaultState(M1M3SSPublisher* publisher) : State(publisher, "FaultState") { }
FaultState::FaultState(M1M3SSPublisher* publisher, std::string name) : State(publisher, name) { }

States::Type FaultState::update(UpdateCommand* command, Model* model) {
	Log.Trace("FaultState: update()");
	Accelerometer* accelerometer = model->getAccelerometer();
	M1M3SSPublisher* publisher = model->getPublisher();
	model->getILC()->writeFreezeSensorListBuffer();
	model->getILC()->triggerModbus();
	model->getDisplacement()->writeDataRequest();
	model->getInclinometer()->writeDataRequest();
	model->getAccelerometer()->sampleData();
	model->getILC()->waitForAllSubnets(5000);
	model->getILC()->readAll();
	model->getDisplacement()->readDataResponse();
	model->getInclinometer()->readDataResponse();
	model->getILC()->calculateHPPostion();
	model->getILC()->calculateHPMirrorForces();
	model->getILC()->calculateFAMirrorForces();
	model->getILC()->verifyResponses();
	usleep(50000);
	model->queryFPGAData();
	usleep(10000);
	model->publishFPGAData();
	model->getILC()->publishForceActuatorStatus();
	model->getILC()->publishForceActuatorData();
	model->getILC()->publishHardpointStatus();
	model->getILC()->publishHardpointData();
	//model->getAirController()->checkStatus();
	model->getInterlockController()->tryToggleHeartbeat();
	return States::NoStateTransition;
}

States::Type FaultState::standby(StandbyCommand* command, Model* model) {
	Log.Trace("FaultState: standby()");
	States::Type newState = States::StandbyState;
	model->getILC()->writeSetModeStandbyBuffer();
	model->getILC()->triggerModbus();
	model->getILC()->waitForAllSubnets(5000);
	model->getILC()->readAll();
	model->getILC()->verifyResponses();
	model->getPowerController()->setAllPowerNetworks(false);
	model->getSafetyController()->clearErrorCode();
	return newState;
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
