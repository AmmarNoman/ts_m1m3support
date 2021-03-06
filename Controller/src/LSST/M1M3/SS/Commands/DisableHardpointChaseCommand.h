/*
 * DisableHardpointChaseCommand.h
 *
 *  Created on: Oct 30, 2017
 *      Author: ccontaxis
 */

#ifndef DISABLEHARDPOINTCHASECOMMAND_H_
#define DISABLEHARDPOINTCHASECOMMAND_H_

#include <Command.h>
#include <SAL_m1m3C.h>
#include <DataTypes.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class DisableHardpointChaseCommand: public Command {
private:
	Context* context;
	M1M3SSPublisher* publisher;
	m1m3_command_DisableHardpointChaseC data;

public:
	DisableHardpointChaseCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_DisableHardpointChaseC* data);

	m1m3_command_DisableHardpointChaseC* getData() { return &this->data; }

	bool validate();
	void execute();
	void ackInProgress();
	void ackComplete();
	void ackFailed(std::string reason);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* DISABLEHARDPOINTCHASECOMMAND_H_ */
