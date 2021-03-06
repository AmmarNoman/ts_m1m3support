/*
 * EnterEngineeringCommand.h
 *
 *  Created on: Oct 30, 2017
 *      Author: ccontaxis
 */

#ifndef ENTERENGINEERINGCOMMAND_H_
#define ENTERENGINEERINGCOMMAND_H_

#include <Command.h>
#include <SAL_m1m3C.h>
#include <DataTypes.h>

namespace LSST {
namespace M1M3 {
namespace SS {

class EnterEngineeringCommand: public Command {
private:
	Context* context;
	M1M3SSPublisher* publisher;
	m1m3_command_EnterEngineeringC data;

public:
	EnterEngineeringCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_EnterEngineeringC* data);

	m1m3_command_EnterEngineeringC* getData() { return &this->data; }

	bool validate();
	void execute();
	void ackInProgress();
	void ackComplete();
	void ackFailed(std::string reason);
};

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */

#endif /* ENTERENGINEERINGCOMMAND_H_ */
