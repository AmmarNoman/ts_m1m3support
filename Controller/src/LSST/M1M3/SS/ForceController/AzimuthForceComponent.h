/*
 * AzimuthForceComponent.h
 *
 *  Created on: Jul 9, 2018
 *      Author: ccontaxis
 */

#ifndef LSST_M1M3_SS_FORCECONTROLLER_AZIMUTHFORCECOMPONENT_H_
#define LSST_M1M3_SS_FORCECONTROLLER_AZIMUTHFORCECOMPONENT_H_

#include <ForceComponent.h>
#include <SAL_m1m3C.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class M1M3SSPublisher;
class SafetyController;
class ForceActuatorApplicationSettings;
class ForceActuatorSettings;

class AzimuthForceComponent: public ForceComponent {
private:
	M1M3SSPublisher* publisher;
	SafetyController* safetyController;
	ForceActuatorApplicationSettings* forceActuatorApplicationSettings;
	ForceActuatorSettings* forceActuatorSettings;

	m1m3_logevent_ForceActuatorStateC* forceActuatorState;
	m1m3_logevent_ForceSetpointWarningC* forceSetpointWarning;
	m1m3_logevent_AppliedAzimuthForcesC* appliedAzimuthForces;
	m1m3_logevent_RejectedAzimuthForcesC* rejectedAzimuthForces;

public:
	AzimuthForceComponent(M1M3SSPublisher* publisher, SafetyController* safetyController, ForceActuatorApplicationSettings* forceActuatorApplicationSettings, ForceActuatorSettings* forceActuatorSettings);

	void applyAzimuthForces(float* x, float* y, float* z);
	void applyAzimuthForcesByAzimuthAngle(float azimuthAngle);

protected:
	void postEnableDisableActions();
	void postUpdateActions();
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* LSST_M1M3_SS_FORCECONTROLLER_AZIMUTHFORCECOMPONENT_H_ */
