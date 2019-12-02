/*
 * InclinometerSettings.cpp
 *
 *  Created on: May 11, 2018
 *      Author: ccontaxis
 */

#include <InclinometerSettings.h>
#include <boost/lexical_cast.hpp>
#include <pugixml.hpp>
#include <Log.h>

namespace LSST {
namespace M1M3 {
namespace SS {

void InclinometerSettings::load(const std::string &filename) {
	pugi::xml_document doc;
	pugi::xml_parse_result load_file_xml_parse_result = doc.load_file(filename.c_str());
	if (!load_file_xml_parse_result) {
		Log.Fatal("Settings file %s could not be loaded", filename.c_str());
		Log.Fatal("Error description: %s", load_file_xml_parse_result.description());
	}
	this->Offset = boost::lexical_cast<float>(doc.select_node("//InclinometerSettings/Offset").node().child_value());
	this->FaultLow = boost::lexical_cast<float>(doc.select_node("//InclinometerSettings/Limits/FaultLow").node().child_value());
	this->WarningLow = boost::lexical_cast<float>(doc.select_node("//InclinometerSettings/Limits/WarningLow").node().child_value());
	this->WarningHigh = boost::lexical_cast<float>(doc.select_node("//InclinometerSettings/Limits/WarningHigh").node().child_value());
	this->FaultHigh = boost::lexical_cast<float>(doc.select_node("//InclinometerSettings/Limits/FaultHigh").node().child_value());
}

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
