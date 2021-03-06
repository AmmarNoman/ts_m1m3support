/*
 * TurnPowerOnCommand.h
 *
 *  Created on: Dec 7, 2017
 *      Author: ccontaxis
 */

#ifndef TURNPOWERONCOMMAND_H_
#define TURNPOWERONCOMMAND_H_

#include <Command.h>
#include <SAL_m1m3C.h>
#include <DataTypes.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class TurnPowerOnCommand: public Command {
private:
	Context* context;
	M1M3SSPublisher* publisher;
	m1m3_command_TurnPowerOnC data;

public:
	TurnPowerOnCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_TurnPowerOnC* data);

	m1m3_command_TurnPowerOnC* getData() { return &this->data; }

	bool validate();
	void execute();
	void ackInProgress();
	void ackComplete();
	void ackFailed(std::string reason);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* TURNPOWERONCOMMAND_H_ */
