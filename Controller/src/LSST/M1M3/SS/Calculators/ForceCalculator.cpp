/*
 * ForceCalculator.cpp
 *
 *  Created on: Oct 23, 2017
 *      Author: ccontaxis
 */

#include <ForceCalculator.h>
#include <ForceActuatorApplicationSettings.h>
#include <SAL_m1m3C.h>
#include <cmath>
#include <vector>

namespace LSST {
namespace M1M3 {
namespace SS {

ForceCalculator::ForceCalculator(ForceActuatorApplicationSettings* forceActuatorSettings, m1m3_logevent_ForceActuatorInfoC* forceInfo, m1m3_InclinometerDataC* inclinometerData, m1m3_ForceActuatorDataC* forceData) {
	this->forceActuatorSettings = forceActuatorSettings;
	this->forceInfo = forceInfo;
	this->inclinometerData = inclinometerData;
	this->forceData = forceData;
}

void ForceCalculator::calculate() {
	this->calculateElevationForces();
	this->calculateForces();
	this->calculateSetpoints();
}

void ForceCalculator::calculateElevationForces() {
	double elevationAngle = this->inclinometerData->InclinometerAngle;
	double elevationMatrix[] = { std::pow(elevationAngle, 5.0), std::pow(elevationAngle, 4.0), std::pow(elevationAngle, 3.0), std::pow(elevationAngle, 2.0), elevationAngle, 1 };
	for(int i = 0; i < FA_COUNT; ++i) {
		int mIndex = i * 6;
		this->forceData->ElevationXSetpoint[i] =
				this->forceActuatorSettings->ElevationXAxisCoefficients[mIndex + 0] * elevationMatrix[0] +
				this->forceActuatorSettings->ElevationXAxisCoefficients[mIndex + 1] * elevationMatrix[1] +
				this->forceActuatorSettings->ElevationXAxisCoefficients[mIndex + 2] * elevationMatrix[2] +
				this->forceActuatorSettings->ElevationXAxisCoefficients[mIndex + 3] * elevationMatrix[3] +
				this->forceActuatorSettings->ElevationXAxisCoefficients[mIndex + 4] * elevationMatrix[4] +
				this->forceActuatorSettings->ElevationXAxisCoefficients[mIndex + 5];
		this->forceData->ElevationYSetpoint[i] =
				this->forceActuatorSettings->ElevationYAxisCoefficients[mIndex + 0] * elevationMatrix[0] +
				this->forceActuatorSettings->ElevationYAxisCoefficients[mIndex + 1] * elevationMatrix[1] +
				this->forceActuatorSettings->ElevationYAxisCoefficients[mIndex + 2] * elevationMatrix[2] +
				this->forceActuatorSettings->ElevationYAxisCoefficients[mIndex + 3] * elevationMatrix[3] +
				this->forceActuatorSettings->ElevationYAxisCoefficients[mIndex + 4] * elevationMatrix[4] +
				this->forceActuatorSettings->ElevationYAxisCoefficients[mIndex + 5];
		this->forceData->ElevationZSetpoint[i] =
				this->forceActuatorSettings->ElevationZAxisCoefficients[mIndex + 0] * elevationMatrix[0] +
				this->forceActuatorSettings->ElevationZAxisCoefficients[mIndex + 1] * elevationMatrix[1] +
				this->forceActuatorSettings->ElevationZAxisCoefficients[mIndex + 2] * elevationMatrix[2] +
				this->forceActuatorSettings->ElevationZAxisCoefficients[mIndex + 3] * elevationMatrix[3] +
				this->forceActuatorSettings->ElevationZAxisCoefficients[mIndex + 4] * elevationMatrix[4] +
				this->forceActuatorSettings->ElevationZAxisCoefficients[mIndex + 5];
	}
}

void ForceCalculator::calculateForces() {
	for(int i = 0; i < FA_COUNT; ++i) {
		this->forceData->XSetpoint[i] = this->forceData->ElevationXSetpoint[i];
		this->forceData->YSetpoint[i] = this->forceData->ElevationYSetpoint[i];
		this->forceData->ZSetpoint[i] = this->forceData->ElevationZSetpoint[i];
	}
}

void ForceCalculator::calculateSetpoints() {
	for(int i = 0; i < FA_COUNT; i++) {
		switch(this->forceInfo->ActuatorOrientation[i]) {
		case ForceActuatorOrientations::PositiveY:
			this->forceData->PrimaryCylinderSetpointCommanded[i] = toInt24(this->forceData->ZSetpoint[i] - this->forceData->YSetpoint[i]);
			this->forceData->SecondaryCylinderSetpointCommanded[i] = toInt24(this->forceData->YSetpoint[i] * sqrt2);
			break;
		case ForceActuatorOrientations::NA:
			this->forceData->PrimaryCylinderSetpointCommanded[i] = toInt24(this->forceData->ZSetpoint[i]);
			break;
		case ForceActuatorOrientations::PositiveX:
			this->forceData->PrimaryCylinderSetpointCommanded[i] = toInt24(this->forceData->ZSetpoint[i] - this->forceData->XSetpoint[i]);
			this->forceData->SecondaryCylinderSetpointCommanded[i] = toInt24(this->forceData->XSetpoint[i] * sqrt2);
			break;
		case ForceActuatorOrientations::NegativeX:
			this->forceData->PrimaryCylinderSetpointCommanded[i] = toInt24(this->forceData->ZSetpoint[i] - -this->forceData->XSetpoint[i]);
			this->forceData->SecondaryCylinderSetpointCommanded[i] = toInt24(-this->forceData->XSetpoint[i] * sqrt2);
			break;
		case ForceActuatorOrientations::NegativeY:
			this->forceData->PrimaryCylinderSetpointCommanded[i] = toInt24(this->forceData->ZSetpoint[i] - -this->forceData->YSetpoint[i]);
			this->forceData->SecondaryCylinderSetpointCommanded[i] = toInt24(-this->forceData->YSetpoint[i] * sqrt2);
			break;
		}
		this->forceData->SecondaryCylinderSetpointCommanded[i] = toInt24(this->forceData->XSetpoint[i] + this->forceData->YSetpoint[i]);
	}
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
