/*
 * BalanceForceComponent.h
 *
 *  Created on: Jul 9, 2018
 *      Author: ccontaxis
 */

#ifndef LSST_M1M3_SS_FORCECONTROLLER_BALANCEFORCECOMPONENT_H_
#define LSST_M1M3_SS_FORCECONTROLLER_BALANCEFORCECOMPONENT_H_

#include <ForceComponent.h>
#include <PID.h>
#include <SAL_m1m3C.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class M1M3SSPublisher;
class SafetyController;
class ForceActuatorApplicationSettings;
class ForceActuatorSettings;
class PIDSettings;

class BalanceForceComponent: public ForceComponent {
private:
	M1M3SSPublisher* publisher;
	SafetyController* safetyController;
	ForceActuatorApplicationSettings* forceActuatorApplicationSettings;
	ForceActuatorSettings* forceActuatorSettings;
	PIDSettings* pidSettings;

	PID fx;
	PID fy;
	PID fz;
	PID mx;
	PID my;
	PID mz;

	m1m3_logevent_ForceActuatorStateC* forceActuatorState;
	m1m3_logevent_ForceSetpointWarningC* forceSetpointWarning;
	m1m3_logevent_AppliedBalanceForcesC* appliedBalanceForces;
	m1m3_logevent_RejectedBalanceForcesC* rejectedBalanceForces;

public:
	BalanceForceComponent(M1M3SSPublisher* publisher, SafetyController* safetyController, ForceActuatorApplicationSettings* forceActuatorApplicationSettings, ForceActuatorSettings* forceActuatorSettings, PIDSettings* pidSettings);

	void applyBalanceForces(float* x, float* y, float* z);
	void applyBalanceForcesByMirrorForces(float xForce, float yForce, float zForce, float xMoment, float yMoment, float zMoment);

	void updatePID(int id, PIDParameters parameters);
	void resetPID(int id);
	void resetPIDs();

protected:
	void postEnableDisableActions();
	void postUpdateActions();

private:
	PID* idToPID(int id);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* LSST_M1M3_SS_FORCECONTROLLER_BALANCEFORCECOMPONENT_H_ */
