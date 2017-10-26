/*
 * IState.cpp
 *
 *  Created on: Sep 26, 2017
 *      Author: ccontaxis
 */

#include <IState.h>

namespace LSST {
namespace M1M3 {
namespace SS {

IState::~IState() { }

States::Type IState::boot(BootCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::start(StartCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::enable(EnableCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::disable(DisableCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::standby(StandbyCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::shutdown(ShutdownCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::update(UpdateCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::turnAirOn(TurnAirOnCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::turnAirOff(TurnAirOffCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::applyOffsetForces(ApplyOffsetForcesCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::clearOffsetForces(ClearOffsetForcesCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::raiseM1M3(RaiseM1M3Command* command, IModel* model) { return States::Ignore; }
States::Type IState::lowerM1M3(LowerM1M3Command* command, IModel* model) { return States::Ignore; }
States::Type IState::applyAberrationByBendingModes(ApplyAberrationByBendingModesCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::applyAberrationByForces(ApplyAberrationByForcesCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::clearAberration(ClearAberrationCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::applyAOSCorrectionByBendingModes(ApplyAOSCorrectionByBendingModesCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::applyAOSCorrectionByForces(ApplyAOSCorrectionByForcesCommand* command, IModel* model) { return States::Ignore; }
States::Type IState::clearAOSCorrection(ClearAOSCorrectionCommand* command, IModel* model) { return States::Ignore; }

} /* namespace SS */
} /* namespace M1M3 */
} /* namespace LSST */
