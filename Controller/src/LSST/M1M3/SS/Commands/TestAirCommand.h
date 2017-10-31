/*
 * TestAirCommand.h
 *
 *  Created on: Oct 30, 2017
 *      Author: ccontaxis
 */

#ifndef TESTAIRCOMMAND_H_
#define TESTAIRCOMMAND_H_

#include <Command.h>
#include <SAL_m1m3C.h>
#include <DataTypes.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class TestAirCommand: public Command {
private:
	IContext* context;
	IPublisher* publisher;
	int32_t commandID;
	m1m3_command_TestAirC data;

public:
	TestAirCommand(IContext* context, IPublisher* publisher, int32_t commandID, m1m3_command_TestAirC* data);

	int32_t getCommandID() { return this->commandID; }
	m1m3_command_TestAirC* getData() { return &this->data; }

	bool validate();
	void execute();
	void ackInProgress();
	void ackComplete();
	void ackFailed(std::string reason);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* TESTAIRCOMMAND_H_ */
