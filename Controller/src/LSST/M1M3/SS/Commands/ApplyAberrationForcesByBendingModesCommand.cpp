/*
 * ApplyAberrationForcesByBendingModesCommand.cpp
 *
 *  Created on: Oct 26, 2017
 *      Author: ccontaxis
 */

#include <ApplyAberrationForcesByBendingModesCommand.h>
#include <Context.h>
#include <M1M3SSPublisher.h>

namespace LSST {
namespace M1M3 {
namespace SS {

ApplyAberrationForcesByBendingModesCommand::ApplyAberrationForcesByBendingModesCommand(Context* context, M1M3SSPublisher* publisher, int32_t commandID, m1m3_command_ApplyAberrationForcesByBendingModesC* data) {
	this->context = context;
	this->publisher = publisher;
	this->commandID = commandID;
	for(int i = 0; i < BENDING_MODES; i++) {
		this->data.Coefficients[i] = data->Coefficients[i];
	}
}

bool ApplyAberrationForcesByBendingModesCommand::validate() {
	return true;
}

void ApplyAberrationForcesByBendingModesCommand::execute() {
	this->context->applyAberrationForcesByBendingModes(this);
}

void ApplyAberrationForcesByBendingModesCommand::ackInProgress() {
	this->publisher->ackCommandApplyAberrationForcesByBendingModes(this->commandID, ACK_INPROGRESS, "In-Progress");
}

void ApplyAberrationForcesByBendingModesCommand::ackComplete() {
	this->publisher->ackCommandApplyAberrationForcesByBendingModes(this->commandID, ACK_COMPLETE, "Complete");
}

void ApplyAberrationForcesByBendingModesCommand::ackFailed(std::string reason) {
	this->publisher->ackCommandApplyAberrationForcesByBendingModes(this->commandID, ACK_COMPLETE, "Failed: " + reason);
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
