/*
 * MoveHardpointActuatorsCommand.h
 *
 *  Created on: Oct 30, 2017
 *      Author: ccontaxis
 */

#ifndef MOVEHARDPOINTACTUATORSCOMMAND_H_
#define MOVEHARDPOINTACTUATORSCOMMAND_H_

#include <Command.h>
#include <SAL_m1m3C.h>
#include <DataTypes.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class MoveHardpointActuatorsCommand: public Command {
private:
	Context* context;
	M1M3SSPublisher* publisher;
	m1m3_command_MoveHardpointActuatorsC data;

public:
	MoveHardpointActuatorsCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_MoveHardpointActuatorsC* data);

	m1m3_command_MoveHardpointActuatorsC* getData() { return &this->data; }

	bool validate();
	void execute();
	void ackInProgress();
	void ackComplete();
	void ackFailed(std::string reason);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* MOVEHARDPOINTACTUATORSCOMMAND_H_ */
