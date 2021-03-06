/*
 * ElevationForceComponent.cpp
 *
 *  Created on: Jul 9, 2018
 *      Author: ccontaxis
 */

#include <ElevationForceComponent.h>
#include <M1M3SSPublisher.h>
#include <SafetyController.h>
#include <ForceActuatorApplicationSettings.h>
#include <ForceActuatorSettings.h>
#include <Range.h>
#include <ForcesAndMoments.h>
#include <ForceConverter.h>
#include <DistributedForces.h>
#include <Log.h>

namespace LSST {
namespace M1M3 {
namespace SS {

ElevationForceComponent::ElevationForceComponent(M1M3SSPublisher* publisher, SafetyController* safetyController, ForceActuatorApplicationSettings* forceActuatorApplicationSettings, ForceActuatorSettings* forceActuatorSettings) {
	this->name = "Elevation";

	this->publisher = publisher;
	this->safetyController = safetyController;
	this->forceActuatorApplicationSettings = forceActuatorApplicationSettings;
	this->forceActuatorSettings = forceActuatorSettings;
	this->forceActuatorState = this->publisher->getEventForceActuatorState();
	this->forceSetpointWarning = this->publisher->getEventForceSetpointWarning();
	this->appliedElevationForces = this->publisher->getEventAppliedElevationForces();
	this->rejectedElevationForces = this->publisher->getEventRejectedElevationForces();
	this->maxRateOfChange = this->forceActuatorSettings->ElevationComponentSettings.MaxRateOfChange;
	this->nearZeroValue = this->forceActuatorSettings->ElevationComponentSettings.NearZeroValue;
}

void ElevationForceComponent::applyElevationForces(float* x, float* y, float* z) {
	Log.Trace("ElevationForceComponent: applyElevationForces()");
	if (!this->enabled) {
		Log.Error("ElevationForceComponent: applyElevationForces() called when the component is not applied");
		return;
	}
	if (this->disabling) {
		Log.Warn("ElevationForceComponent: applyElevationForces() called when the component is disabling");
		return;
	}
	for(int i = 0; i < 156; ++i) {
		if (i < 12) {
			this->xTarget[i] = x[i] * this->forceActuatorState->SupportPercentage;
		}

		if (i < 100) {
			this->yTarget[i] = y[i] * this->forceActuatorState->SupportPercentage;
		}

		this->zTarget[i] = z[i] * this->forceActuatorState->SupportPercentage;
	}
}

void ElevationForceComponent::applyElevationForcesByElevationAngle(float elevationAngle) {
	Log.Trace("ElevationForceComponent: applyElevationForcesByMirrorForces(%0.1f)", elevationAngle);
	DistributedForces forces = ForceConverter::calculateForceFromElevationAngle(this->forceActuatorSettings, elevationAngle);
	float xForces[12];
	float yForces[100];
	float zForces[156];
	for(int zIndex = 0; zIndex < FA_COUNT; ++zIndex) {
		int xIndex = this->forceActuatorApplicationSettings->ZIndexToXIndex[zIndex];
		int yIndex = this->forceActuatorApplicationSettings->ZIndexToYIndex[zIndex];

		if (xIndex != -1) {
			xForces[xIndex] = forces.XForces[zIndex];
		}
		if (yIndex != -1) {
			yForces[yIndex] = forces.YForces[zIndex];
		}
		zForces[zIndex] = forces.ZForces[zIndex];
	}
	this->applyElevationForces(xForces, yForces, zForces);
}

void ElevationForceComponent::postEnableDisableActions() {
	Log.Debug("ElevationForceComponent: postEnableDisableActions()");

	this->forceActuatorState->Timestamp = this->publisher->getTimestamp();
	this->forceActuatorState->ElevationForcesApplied = this->enabled;
	this->publisher->tryLogForceActuatorState();
}

void ElevationForceComponent::postUpdateActions() {
	Log.Trace("ElevationForceController: postUpdateActions()");

	bool notInRange = false;
	bool rejectionRequired = false;
	this->appliedElevationForces->Timestamp = this->publisher->getTimestamp();
	this->rejectedElevationForces->Timestamp = this->appliedElevationForces->Timestamp;
	for(int zIndex = 0; zIndex < 156; ++zIndex) {
		int xIndex = this->forceActuatorApplicationSettings->ZIndexToXIndex[zIndex];
		int yIndex = this->forceActuatorApplicationSettings->ZIndexToYIndex[zIndex];

		this->forceSetpointWarning->ElevationForceWarning[zIndex] = false;

		if (xIndex != -1) {
			float xLowFault = this->forceActuatorSettings->ElevationLimitXTable[xIndex].LowFault;
			float xHighFault = this->forceActuatorSettings->ElevationLimitXTable[xIndex].HighFault;
			this->rejectedElevationForces->XForces[xIndex] = this->xCurrent[xIndex];
			notInRange = !Range::InRangeAndCoerce(xLowFault, xHighFault, this->rejectedElevationForces->XForces[xIndex], this->appliedElevationForces->XForces + xIndex);
			this->forceSetpointWarning->ElevationForceWarning[zIndex] = this->forceSetpointWarning->ElevationForceWarning[zIndex] || notInRange;
		}

		if (yIndex != -1) {
			float yLowFault = this->forceActuatorSettings->ElevationLimitYTable[yIndex].LowFault;
			float yHighFault = this->forceActuatorSettings->ElevationLimitYTable[yIndex].HighFault;
			this->rejectedElevationForces->YForces[yIndex] = this->yCurrent[yIndex];
			notInRange = !Range::InRangeAndCoerce(yLowFault, yHighFault, this->rejectedElevationForces->YForces[yIndex], this->appliedElevationForces->YForces + yIndex);
			this->forceSetpointWarning->ElevationForceWarning[zIndex] = this->forceSetpointWarning->ElevationForceWarning[zIndex] || notInRange;
		}

		float zLowFault = this->forceActuatorSettings->ElevationLimitZTable[zIndex].LowFault;
		float zHighFault = this->forceActuatorSettings->ElevationLimitZTable[zIndex].HighFault;
		this->rejectedElevationForces->ZForces[zIndex] = this->zCurrent[zIndex];
		notInRange = !Range::InRangeAndCoerce(zLowFault, zHighFault, this->rejectedElevationForces->ZForces[zIndex], this->appliedElevationForces->ZForces + zIndex);
		this->forceSetpointWarning->ElevationForceWarning[zIndex] = this->forceSetpointWarning->ElevationForceWarning[zIndex] || notInRange;
		rejectionRequired = rejectionRequired || this->forceSetpointWarning->ElevationForceWarning[zIndex];
	}

	ForcesAndMoments fm = ForceConverter::calculateForcesAndMoments(this->forceActuatorApplicationSettings, this->forceActuatorSettings, this->appliedElevationForces->XForces, this->appliedElevationForces->YForces, this->appliedElevationForces->ZForces);
	this->appliedElevationForces->Fx = fm.Fx;
	this->appliedElevationForces->Fy = fm.Fy;
	this->appliedElevationForces->Fz = fm.Fz;
	this->appliedElevationForces->Mx = fm.Mx;
	this->appliedElevationForces->My = fm.My;
	this->appliedElevationForces->Mz = fm.Mz;
	this->appliedElevationForces->ForceMagnitude = fm.ForceMagnitude;

	fm = ForceConverter::calculateForcesAndMoments(this->forceActuatorApplicationSettings, this->forceActuatorSettings, this->rejectedElevationForces->XForces, this->rejectedElevationForces->YForces, this->rejectedElevationForces->ZForces);
	this->rejectedElevationForces->Fx = fm.Fx;
	this->rejectedElevationForces->Fy = fm.Fy;
	this->rejectedElevationForces->Fz = fm.Fz;
	this->rejectedElevationForces->Mx = fm.Mx;
	this->rejectedElevationForces->My = fm.My;
	this->rejectedElevationForces->Mz = fm.Mz;
	this->rejectedElevationForces->ForceMagnitude = fm.ForceMagnitude;

	this->safetyController->forceControllerNotifyElevationForceClipping(rejectionRequired);

	this->publisher->tryLogForceSetpointWarning();
	if (rejectionRequired) {
		this->publisher->logRejectedElevationForces();
	}
	this->publisher->logAppliedElevationForces();
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
