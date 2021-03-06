/*
 * StartCommand.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: ccontaxis
 */

#include <StartCommand.h>
#include <Context.h>
#include <M1M3SSPublisher.h>

namespace LSST {
namespace M1M3 {
namespace SS {

StartCommand::StartCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_StartC* data) {
	this->context = context;
	this->publisher = publisher;
	this->commandID = commandID;
	this->data.Start = data->Start;
	this->data.SettingsToApply = data->SettingsToApply;
}

bool StartCommand::validate() {
	if (!this->data.Start) {
		this->publisher->logCommandRejectionWarning("Start", "The field Start is not TRUE.");
	}
	return this->data.Start;
}

void StartCommand::execute() {
	this->context->start(this);
}

void StartCommand::ackInProgress() {
	this->publisher->ackCommandStart(this->commandID, ACK_INPROGRESS, "In-Progress");
}

void StartCommand::ackComplete() {
	this->publisher->ackCommandStart(this->commandID, ACK_COMPLETE, "Complete");
}

void StartCommand::ackFailed(std::string reason) {
	this->publisher->ackCommandStart(this->commandID, ACK_COMPLETE, "Failed: " + reason);
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
