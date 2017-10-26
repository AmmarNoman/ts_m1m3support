/*
 * DisabledState.cpp
 *
 *  Created on: Sep 27, 2017
 *      Author: ccontaxis
 */

#include <DisabledState.h>
#include <SAL_m1m3C.h>
#include <IModel.h>
#include <IDisplacement.h>
#include <IInclinometer.h>
#include <unistd.h>
#include <IPublisher.h>
#include <IILC.h>
#include <IFPGA.h>
#include <AirController.h>
#include <iostream>

using namespace std;

namespace LSST {
namespace M1M3 {
namespace SS {

States::Type DisabledState::enable(EnableCommand* command, IModel* model) {
	States::Type newState = States::ParkedEngineeringState;
	model->publishStateChange(newState);
	return newState;
}

States::Type DisabledState::standby(StandbyCommand* command, IModel* model) {
	States::Type newState = States::StandbyState;
	model->getILC()->writeSetModeStandbyBuffer();
	model->getILC()->triggerModbus();
	model->getILC()->waitForAllSubnets(5000);
	model->getILC()->readAll();
	model->getILC()->verifyResponses();
	model->publishStateChange(newState);
	return newState;
}

States::Type DisabledState::update(UpdateCommand* command, IModel* model) {
	model->getILC()->writeFreezeSensorListBuffer();
	model->getILC()->triggerModbus();
	model->getDisplacement()->writeDataRequest();
	model->getInclinometer()->writeDataRequest();
	model->getILC()->waitForAllSubnets(5000);
	model->getILC()->readAll();
	model->getDisplacement()->readDataResponse();
	model->getInclinometer()->readDataResponse();
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
	return States::Ignore;
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
