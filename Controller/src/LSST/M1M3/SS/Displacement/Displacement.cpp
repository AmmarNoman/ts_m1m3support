/*
 * Displacement.cpp
 *
 *  Created on: Oct 6, 2017
 *      Author: ccontaxis
 */

#include <Displacement.h>
#include <IPublisher.h>
#include <ISafetyController.h>
#include <IFPGA.h>
#include <FPGAAddresses.h>
#include <U8ArrayUtilities.h>
#include <Timestamp.h>
#include <SAL_m1m3C.h>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <string>
#include <cstring>

namespace LSST {
namespace M1M3 {
namespace SS {

Displacement::Displacement(IPublisher* publisher, ISafetyController* safetyController, IFPGA* fpga) {
	this->publisher = publisher;
	this->safetyController = safetyController;
	this->fpga = fpga;
	this->imsData = this->publisher->getIMSData();
	this->displacementWarning = this->publisher->getEventDisplacementSensorWarning();

	memset(this->imsData, 0, sizeof(m1m3_IMSDataC));
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));

	this->createTxBuffer();
}

void Displacement::writeDataRequest() {
	this->fpga->writeCommandFIFO(this->txBuffer, 6, 0);
}

void Displacement::readDataResponse() {
	this->fpga->writeRequestFIFO(FPGAAddresses::Displacement, 0);
	this->fpga->readU8ResponseFIFO(this->rxBuffer, 2, 10);
	uint16_t length = U8ArrayUtilities::u16(this->rxBuffer, 0);
	if (length != 0) {
		this->fpga->readU8ResponseFIFO(this->rxBuffer, length, 10);
		double timestamp = Timestamp::fromRaw(U8ArrayUtilities::u64(this->rxBuffer, 0));
		std::string response = U8ArrayUtilities::toString(this->rxBuffer, 8, length - 10);
		typedef boost::tokenizer< boost::escaped_list_separator<char> > tokenizer;
		tokenizer tokenize(response);
		tokenizer::iterator token = tokenize.begin();
		if ((*token) == "M0") {
			this->imsData->Timestamp = timestamp;
			++token;
			this->imsData->RawSensorData[0] = boost::lexical_cast<double>(*token);
			++token;
			this->imsData->RawSensorData[1] = boost::lexical_cast<double>(*token);
			++token;
			this->imsData->RawSensorData[2] = boost::lexical_cast<double>(*token);
			++token;
			this->imsData->RawSensorData[3] = boost::lexical_cast<double>(*token);
			++token;
			this->imsData->RawSensorData[4] = boost::lexical_cast<double>(*token);
			++token;
			this->imsData->RawSensorData[5] = boost::lexical_cast<double>(*token);
			this->clearWarning(timestamp);
			this->publisher->putIMSData();
		}
		else if ((*token) == "ER") {
			++token;
			++token;
			switch(boost::lexical_cast<int32_t>(*token)) {
			case 0: this->warnSensorReportsInvalidCommand(timestamp); break;
			case 2: this->warnSensorReportsCommunicationTimeoutError(timestamp); break;
			case 20: this->warnSensorReportsDataLengthError(timestamp); break;
			case 21: this->warnSensorReportsNumberOfParametersError(timestamp); break;
			case 22: this->warnSensorReportsParameterError(timestamp); break;
			case 29: this->warnSensorReportsCommunicationError(timestamp); break;
			case 65: this->warnSensorReportsIDNumberError(timestamp); break;
			case 66: this->warnSensorReportsExpansionLineError(timestamp); break;
			case 67: this->warnSensorReportsWriteControlError(timestamp); break;
			default: this->warnUnknownProblem(timestamp); break;
			}
		}
		else {
			this->warnUnknownCommand(timestamp);
		}
	}
	else {
		this->warnResponseTimeoutError(this->publisher->getTimestamp());
	}
}

void Displacement::createTxBuffer() {
	this->txBuffer[0] = FPGAAddresses::Displacement;
	this->txBuffer[1] = 4;
	this->txBuffer[2] = (uint8_t)'M';
	this->txBuffer[3] = (uint8_t)'0';
	this->txBuffer[4] = (uint8_t)'\r';
	this->txBuffer[5] = (uint8_t)'\n';
}

void Displacement::clearWarning(double timestamp) {
	if (this->displacementWarning->AnyWarning) {
		memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
		this->displacementWarning->Timestamp = timestamp;
		this->safetyController->displacementNotifySensorReportsInvalidCommand(this->displacementWarning->SensorReportsInvalidCommand);
		this->safetyController->displacementNotifySensorReportsCommunicationTimeoutError(this->displacementWarning->SensorReportsCommunicationTimeoutError);
		this->safetyController->displacementNotifySensorReportsNumberOfParametersError(this->displacementWarning->SensorReportsNumberOfParametersError);
		this->safetyController->displacementNotifySensorReportsParameterError(this->displacementWarning->SensorReportsParameterError);
		this->safetyController->displacementNotifySensorReportsCommunicationError(this->displacementWarning->SensorReportsCommunicationError);
		this->safetyController->displacementNotifySensorReportsIDNumberError(this->displacementWarning->SensorReportsIDNumberError);
		this->safetyController->displacementNotifySensorReportsExpansionLineError(this->displacementWarning->SensorReportsExpansionLineError);
		this->safetyController->displacementNotifySensorReportsWriteControlError(this->displacementWarning->SensorReportsWriteControlError);
		this->safetyController->displacementNotifyResponseTimeoutError(this->displacementWarning->ResponseTimeout);
		this->safetyController->displacementNotifyInvalidLength(this->displacementWarning->InvalidLength);
		this->safetyController->displacementNotifyUnknownCommand(this->displacementWarning->UnknownCommand);
		this->safetyController->displacementNotifyUnknownProblem(this->displacementWarning->UnknownProblem);
		this->publisher->logDisplacementSensorWarning();
	}
}

void Displacement::warnSensorReportsInvalidCommand(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsInvalidCommand = true;
	this->safetyController->displacementNotifySensorReportsInvalidCommand(this->displacementWarning->SensorReportsInvalidCommand);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsCommunicationTimeoutError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsCommunicationTimeoutError = true;
	this->safetyController->displacementNotifySensorReportsCommunicationTimeoutError(this->displacementWarning->SensorReportsCommunicationTimeoutError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsDataLengthError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsDataLengthError = true;
	this->safetyController->displacementNotifySensorReportsDataLengthError(this->displacementWarning->SensorReportsDataLengthError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsNumberOfParametersError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsNumberOfParametersError = true;
	this->safetyController->displacementNotifySensorReportsNumberOfParametersError(this->displacementWarning->SensorReportsNumberOfParametersError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsParameterError(double timestamp) {
	memset(&this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsParameterError = true;
	this->safetyController->displacementNotifySensorReportsParameterError(this->displacementWarning->SensorReportsParameterError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsCommunicationError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsCommunicationError = true;
	this->safetyController->displacementNotifySensorReportsCommunicationError(this->displacementWarning->SensorReportsCommunicationError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsIDNumberError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsIDNumberError = true;
	this->safetyController->displacementNotifySensorReportsIDNumberError(this->displacementWarning->SensorReportsIDNumberError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsExpansionLineError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsExpansionLineError = true;
	this->safetyController->displacementNotifySensorReportsExpansionLineError(this->displacementWarning->SensorReportsExpansionLineError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnSensorReportsWriteControlError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->SensorReportsWriteControlError = true;
	this->safetyController->displacementNotifySensorReportsWriteControlError(this->displacementWarning->SensorReportsWriteControlError);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnResponseTimeoutError(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->ResponseTimeout = true;
	this->safetyController->displacementNotifyResponseTimeoutError(this->displacementWarning->ResponseTimeout);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnInvalidLength(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->InvalidLength = true;
	this->safetyController->displacementNotifyInvalidLength(this->displacementWarning->InvalidLength);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnUnknownCommand(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->UnknownCommand = true;
	this->safetyController->displacementNotifyUnknownCommand(this->displacementWarning->UnknownCommand);
	this->publisher->logDisplacementSensorWarning();
}

void Displacement::warnUnknownProblem(double timestamp) {
	memset(this->displacementWarning, 0, sizeof(m1m3_logevent_DisplacementSensorWarningC));
	this->displacementWarning->Timestamp = timestamp;
	this->displacementWarning->AnyWarning = true;
	this->displacementWarning->UnknownProblem = true;
	this->safetyController->displacementNotifyUnknownProblem(this->displacementWarning->UnknownProblem);
	this->publisher->logDisplacementSensorWarning();
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
