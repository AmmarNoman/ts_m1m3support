/*
 * ProfileController.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: ccontaxis
 */

#include <ProfileController.h>
#include <Log.h>

namespace LSST {
namespace M1M3 {
namespace SS {

ProfileController::ProfileController() {
	Log.Debug("ProfileController: ProfileController()");
}

void ProfileController::setupMirrorForceProfile(float* xForce, float* yForce, float* zForce, float* xMoment, float* yMoment, float* zMoment) {
	Log.Info("ProfileController: setupMirrorForceProfile()");
	this->mirrorForceProfile.Index = 0;
	for (int i = 0; i < 1000; ++i) {
		this->mirrorForceProfile.Records[i].XForce = xForce[i];
		this->mirrorForceProfile.Records[i].YForce = yForce[i];
		this->mirrorForceProfile.Records[i].ZForce = zForce[i];
		this->mirrorForceProfile.Records[i].XMoment = xMoment[i];
		this->mirrorForceProfile.Records[i].YMoment = yMoment[i];
		this->mirrorForceProfile.Records[i].ZMoment = zMoment[i];
	}
}

bool ProfileController::incMirrorForceProfile() {
	Log.Trace("ProfileController: incMirrorForceProfile()");
	++this->mirrorForceProfile.Index;
	return this->mirrorForceProfile.Index >= 1000;
}

MirrorForceProfileRecord ProfileController::getMirrorForceProfileData() {
	return this->mirrorForceProfile.Records[this->mirrorForceProfile.Index];
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
