/*
 * TurnLightsOffCommand.cpp
 *
 *  Created on: Nov 30, 2017
 *      Author: ccontaxis
 */

#include <TurnLightsOffCommand.h>
#include <Context.h>
#include <M1M3SSPublisher.h>

namespace LSST {
namespace M1M3 {
namespace SS {

TurnLightsOffCommand::TurnLightsOffCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_TurnLightsOffC* data) {
	this->context = context;
	this->publisher = publisher;
	this->commandID = commandID;
	this->data.TurnLightsOff = data->TurnLightsOff;
}

bool TurnLightsOffCommand::validate() {
	if (!this->data.TurnLightsOff) {
		this->publisher->logCommandRejectionWarning("TurnLightsOff", "The field TurnLightsOff is not TRUE.");
	}
	return this->data.TurnLightsOff;
}

void TurnLightsOffCommand::execute() {
	this->context->turnLightsOff(this);
}

void TurnLightsOffCommand::ackInProgress() {
	this->publisher->ackCommandTurnLightsOff(this->commandID, ACK_INPROGRESS, "In-Progress");
}

void TurnLightsOffCommand::ackComplete() {
	this->publisher->ackCommandTurnLightsOff(this->commandID, ACK_COMPLETE, "Complete");
}

void TurnLightsOffCommand::ackFailed(std::string reason) {
	this->publisher->ackCommandTurnLightsOff(this->commandID, ACK_COMPLETE, "Failed: " + reason);
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
