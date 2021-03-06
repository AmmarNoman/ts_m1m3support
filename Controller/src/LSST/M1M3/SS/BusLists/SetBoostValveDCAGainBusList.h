/*
 * SetBoostValveDCAGainBusList.h
 *
 *  Created on: Oct 19, 2017
 *      Author: ccontaxis
 */

#ifndef SETBOOSTVALVEDCAGAINBUSLIST_H_
#define SETBOOSTVALVEDCAGAINBUSLIST_H_

#include <BusList.h>
#include <SAL_m1m3C.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class SetBoostValveDCAGainBusList: public BusList {
private:
	m1m3_logevent_ForceActuatorInfoC* forceInfo;

public:
	SetBoostValveDCAGainBusList(ILCSubnetData* subnetData, ILCMessageFactory* ilcMessageFactory, m1m3_logevent_ForceActuatorInfoC* forceInfo);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* SETBOOSTVALVEDCAGAINBUSLIST_H_ */
