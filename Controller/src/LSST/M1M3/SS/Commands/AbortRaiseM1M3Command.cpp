/*
 * AbortRaiseM1M3Command.cpp
 *
 *  Created on: Oct 30, 2017
 *      Author: ccontaxis
 */

#include <AbortRaiseM1M3Command.h>
#include <Context.h>
#include <M1M3SSPublisher.h>

namespace LSST {
namespace M1M3 {
namespace SS {

AbortRaiseM1M3Command::AbortRaiseM1M3Command(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_AbortRaiseM1M3C* data) {
	this->context = context;
	this->publisher = publisher;
	this->commandID = commandID;
	this->data.AbortRaiseM1M3 = data->AbortRaiseM1M3;
}

bool AbortRaiseM1M3Command::validate() {
	if (!this->data.AbortRaiseM1M3) {
		this->publisher->logCommandRejectionWarning("AbortRaiseM1M3", "The field AbortRaiseM1M3 is not TRUE.");
	}
	return this->data.AbortRaiseM1M3;
}

void AbortRaiseM1M3Command::execute() {
	this->context->abortRaiseM1M3(this);
}

void AbortRaiseM1M3Command::ackInProgress() {
	this->publisher->ackCommandAbortRaiseM1M3(this->commandID, ACK_INPROGRESS, "In-Progress");
}

void AbortRaiseM1M3Command::ackComplete() {
	this->publisher->ackCommandAbortRaiseM1M3(this->commandID, ACK_COMPLETE, "Completed");
}

void AbortRaiseM1M3Command::ackFailed(std::string reason) {
	this->publisher->ackCommandAbortRaiseM1M3(this->commandID, ACK_COMPLETE, "Failed: " + reason);
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
