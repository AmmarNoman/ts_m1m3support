/*
 * EnabledState.h
 *
 *  Created on: Sep 27, 2017
 *      Author: ccontaxis
 */

#ifndef ENABLEDSTATE_H_
#define ENABLEDSTATE_H_

#include <State.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class EnabledState: public State {
public:
	States::Type storeTMAAzimuthSample(TMAAzimuthSampleCommand* command, IModel* model);
	States::Type storeTMAElevationSample(TMAElevationSampleCommand* command, IModel* model);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* ENABLEDSTATE_H_ */
