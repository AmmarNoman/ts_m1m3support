/*
 * SetADCChanneOffsetAndSensitivityBusList.cpp
 *
 *  Created on: Oct 19, 2017
 *      Author: ccontaxis
 */

#include <SetADCChanneOffsetAndSensitivityBusList.h>
#include <ILCSubnetData.h>
#include <ILCMessageFactory.h>
#include <SAL_m1m3C.h>
#include <iostream>

using namespace std;

namespace LSST {
namespace M1M3 {
namespace SS {

SetADCChanneOffsetAndSensitivityBusList::SetADCChanneOffsetAndSensitivityBusList() : BusList() {
	this->forceInfo = 0;
	this->hardpointInfo = 0;
}

SetADCChanneOffsetAndSensitivityBusList::SetADCChanneOffsetAndSensitivityBusList(ILCSubnetData* subnetData, ILCMessageFactory* ilcMessageFactory, m1m3_logevent_ForceActuatorInfoC* forceInfo, m1m3_logevent_HardpointActuatorInfoC* hardpointInfo)
: BusList(subnetData, ilcMessageFactory) {
	cout << "Start" << endl;
	this->forceInfo = forceInfo;
	this->hardpointInfo = hardpointInfo;
	for(int subnetIndex = 0; subnetIndex < SUBNET_COUNT; ++subnetIndex) {
		cout << "SN" << subnetIndex << " " << this->subnetData->getFACount(subnetIndex) << endl;
		this->startSubnet(subnetIndex);
		for(int faIndex = 0; faIndex < this->subnetData->getFACount(subnetIndex); ++faIndex) {
			cout << "FA" << faIndex << endl;
			uint8_t address = this->subnetData->getFAIndex(subnetIndex, faIndex).Address;
			int32_t dataIndex = this->subnetData->getFAIndex(subnetIndex, faIndex).DataIndex;
			this->ilcMessageFactory->setADCChannelOffsetAndSensitivity(&this->buffer, address, 1, this->forceInfo->PrimaryCylinderSensorOffset[dataIndex], this->forceInfo->PrimaryCylinderSensorSensitivity[dataIndex]);
			this->ilcMessageFactory->setADCChannelOffsetAndSensitivity(&this->buffer, address, 2, this->forceInfo->SecondaryCylinderSensorOffset[dataIndex], this->forceInfo->SecondaryCylinderSensorSensitivity[dataIndex]);
			this->expectedFAResponses[dataIndex] = 2;
		}
		for(int hpIndex = 0; hpIndex < this->subnetData->getHPCount(subnetIndex); ++hpIndex) {
			cout << "HP" << hpIndex << endl;
			uint8_t address = this->subnetData->getHPIndex(subnetIndex, hpIndex).Address;
			int32_t dataIndex = this->subnetData->getHPIndex(subnetIndex, hpIndex).DataIndex;
			this->ilcMessageFactory->setADCChannelOffsetAndSensitivity(&this->buffer, address, 1, this->hardpointInfo->SensorOffset[dataIndex], this->hardpointInfo->SensorSensitivity[dataIndex]);
			this->expectedHPResponses[dataIndex] = 1;
		}
		this->endSubnet();
	}
	this->buffer.setLength(this->buffer.getIndex());
	cout << "Done" << endl;
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
