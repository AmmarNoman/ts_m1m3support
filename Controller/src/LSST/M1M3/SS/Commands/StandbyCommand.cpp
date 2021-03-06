/*
 * StandbyCommand.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: ccontaxis
 */

#include <StandbyCommand.h>
#include <Context.h>
#include <M1M3SSPublisher.h>

namespace LSST {
namespace M1M3 {
namespace SS {

StandbyCommand::StandbyCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_StandbyC* data) {
	this->context = context;
	this->publisher = publisher;
	this->commandID = commandID;
	this->data.Standby = data->Standby;
}

bool StandbyCommand::validate() {
	if (!this->data.Standby) {
		this->publisher->logCommandRejectionWarning("Standby", "The field Standby is not TRUE.");
	}
	return this->data.Standby;
}

void StandbyCommand::execute() {
	this->context->standby(this);
}

void StandbyCommand::ackInProgress() {
	this->publisher->ackCommandStandby(this->commandID, ACK_INPROGRESS, "In-Progress");
}

void StandbyCommand::ackComplete() {
	this->publisher->ackCommandStandby(this->commandID, ACK_COMPLETE, "Complete");
}

void StandbyCommand::ackFailed(std::string reason) {
	this->publisher->ackCommandStandby(this->commandID, ACK_COMPLETE, "Failed: " + reason);
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
