/*
 * TestForceActuatorCommand.h
 *
 *  Created on: Oct 30, 2017
 *      Author: ccontaxis
 */

#ifndef TESTFORCEACTUATORCOMMAND_H_
#define TESTFORCEACTUATORCOMMAND_H_

#include <Command.h>
#include <SAL_m1m3C.h>
#include <DataTypes.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class TestForceActuatorCommand: public Command {
private:
	Context* context;
	M1M3SSPublisher* publisher;
	m1m3_command_TestForceActuatorC data;

public:
	TestForceActuatorCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_TestForceActuatorC* data);

	m1m3_command_TestForceActuatorC* getData() { return &this->data; }

	bool validate();
	void execute();
	void ackInProgress();
	void ackComplete();
	void ackFailed(std::string reason);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* TESTFORCEACTUATORCOMMAND_H_ */
