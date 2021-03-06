/*
 * EngineeringState.cpp
 *
 *  Created on: Nov 3, 2017
 *      Author: ccontaxis
 */

#include <EngineeringState.h>
#include <Model.h>
#include <DigitalInputOutput.h>
#include <PositionController.h>
#include <ForceController.h>
#include <PowerController.h>
#include <MoveHardpointActuatorsCommand.h>
#include <EnableHardpointChaseCommand.h>
#include <DisableHardpointChaseCommand.h>
#include <ApplyOffsetForcesCommand.h>
#include <ApplyOffsetForcesByMirrorForceCommand.h>
#include <TurnPowerOnCommand.h>
#include <TurnPowerOffCommand.h>
#include <SafetyController.h>
#include <M1M3SSPublisher.h>
#include <ILC.h>
#include <ModbusTransmitCommand.h>
#include <Log.h>

namespace LSST {
namespace M1M3 {
namespace SS {

EngineeringState::EngineeringState(M1M3SSPublisher* publisher) : EnabledState(publisher, "EngineeringState") { }
EngineeringState::EngineeringState(M1M3SSPublisher* publisher, std::string name) : EnabledState(publisher, name) { }

States::Type EngineeringState::turnAirOn(TurnAirOnCommand* command, Model* model) {
	Log.Info("%s: turnAirOn()", this->name.c_str());
	model->getDigitalInputOutput()->turnAirOn();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::turnAirOff(TurnAirOffCommand* command, Model* model) {
	Log.Info("%s: turnAirOff()", this->name.c_str());
	model->getDigitalInputOutput()->turnAirOff();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::stopHardpointMotion(StopHardpointMotionCommand* command, Model* model) {
	Log.Info("%s: stopHardpointMotion()", this->name.c_str());
	model->getPositionController()->stopMotion();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::moveHardpointActuators(MoveHardpointActuatorsCommand* command, Model* model) {
	Log.Info("%s: moveHardpointActuators()", this->name.c_str());
	if (!model->getPositionController()->move(command->getData()->Steps)) {
		model->getPublisher()->logCommandRejectionWarning("MoveHardpointActuators", "At least one hardpoint actuator commanded to move is already MOVING or CHASING.");
	}
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::enableHardpointChase(EnableHardpointChaseCommand* command, Model* model) {
	Log.Info("%s: enableHardpointChase()", this->name.c_str());
	if (!model->getPositionController()->enableChase(command->getData()->HardpointActuator)) {
		model->getPublisher()->logCommandRejectionWarning("EnableHardpointChase", "At least one hardpoint actuator commanded to chase is already MOVING or CHASING.");
	}
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::disableHardpointChase(DisableHardpointChaseCommand* command, Model* model) {
	Log.Info("%s: disableHardpointChase()", this->name.c_str());
	model->getPositionController()->disableChase(command->getData()->HardpointActuator);
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::applyOffsetForces(ApplyOffsetForcesCommand* command, Model* model) {
	Log.Info("%s: applyOffsetForces()", this->name.c_str());
	model->getForceController()->applyOffsetForces(command->getData()->XForces, command->getData()->YForces, command->getData()->ZForces);
	model->getForceController()->processAppliedForces();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::applyOffsetForcesByMirrorForce(ApplyOffsetForcesByMirrorForceCommand* command, Model* model) {
	Log.Info("%s: applyOffsetForcesByMirrorForce()", this->name.c_str());
	model->getForceController()->applyOffsetForcesByMirrorForces(command->getData()->XForce, command->getData()->YForce, command->getData()->ZForce, command->getData()->XMoment, command->getData()->YMoment, command->getData()->ZMoment);
	model->getForceController()->processAppliedForces();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::clearOffsetForces(ClearOffsetForcesCommand* command, Model* model) {
	Log.Info("%s: clearOffsetForces()", this->name.c_str());
	model->getForceController()->zeroOffsetForces();
	model->getForceController()->processAppliedForces();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::turnLightsOn(TurnLightsOnCommand* command, Model* model) {
	Log.Info("%s: turnLightsOn()", this->name.c_str());
	model->getDigitalInputOutput()->turnCellLightsOn();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::turnLightsOff(TurnLightsOffCommand* command, Model* model) {
	Log.Info("%s: turnLightsOff()", this->name.c_str());
	model->getDigitalInputOutput()->turnCellLightsOff();
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::turnPowerOn(TurnPowerOnCommand* command, Model* model) {
	Log.Info("%s: turnPowerOn()", this->name.c_str());
	if (command->getData()->TurnPowerNetworkAOn) {
		model->getPowerController()->setPowerNetworkA(true);
	}
	if (command->getData()->TurnPowerNetworkBOn) {
		model->getPowerController()->setPowerNetworkB(true);
	}
	if (command->getData()->TurnPowerNetworkCOn) {
		model->getPowerController()->setPowerNetworkC(true);
	}
	if (command->getData()->TurnPowerNetworkDOn) {
		model->getPowerController()->setPowerNetworkD(true);
	}
	if (command->getData()->TurnAuxPowerNetworkAOn) {
		model->getPowerController()->setAuxPowerNetworkA(true);
	}
	if (command->getData()->TurnAuxPowerNetworkBOn) {
		model->getPowerController()->setAuxPowerNetworkB(true);
	}
	if (command->getData()->TurnAuxPowerNetworkCOn) {
		model->getPowerController()->setAuxPowerNetworkC(true);
	}
	if (command->getData()->TurnAuxPowerNetworkDOn) {
		model->getPowerController()->setAuxPowerNetworkD(true);
	}
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::turnPowerOff(TurnPowerOffCommand* command, Model* model) {
	Log.Info("%s: turnPowerOff()", this->name.c_str());
	if (command->getData()->TurnPowerNetworkAOff) {
		model->getPowerController()->setPowerNetworkA(false);
	}
	if (command->getData()->TurnPowerNetworkBOff) {
		model->getPowerController()->setPowerNetworkB(false);
	}
	if (command->getData()->TurnPowerNetworkCOff) {
		model->getPowerController()->setPowerNetworkC(false);
	}
	if (command->getData()->TurnPowerNetworkDOff) {
		model->getPowerController()->setPowerNetworkD(false);
	}
	if (command->getData()->TurnAuxPowerNetworkAOff) {
		model->getPowerController()->setAuxPowerNetworkA(false);
	}
	if (command->getData()->TurnAuxPowerNetworkBOff) {
		model->getPowerController()->setAuxPowerNetworkB(false);
	}
	if (command->getData()->TurnAuxPowerNetworkCOff) {
		model->getPowerController()->setAuxPowerNetworkC(false);
	}
	if (command->getData()->TurnAuxPowerNetworkDOff) {
		model->getPowerController()->setAuxPowerNetworkD(false);
	}
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

States::Type EngineeringState::modbusTransmit(ModbusTransmitCommand* command, Model* model) {
	Log.Info("%s: modbusTransmit()", this->name.c_str());
	model->getILC()->modbusTransmit(command->getData()->ActuatorId, command->getData()->FunctionCode, command->getData()->DataLength, command->getData()->Data);
	return model->getSafetyController()->checkSafety(States::NoStateTransition);
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
