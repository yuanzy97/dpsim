/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <cps/Definitions.h>

namespace CPS {
namespace Base {
namespace Ph1 {
	class Transformer {
	protected:
		/// Nominal voltage of primary side
		/// FIXME: SP wants this to be an attribute, DP not
		Real mNominalVoltageEnd1;
		/// Nominal voltage of secondary side
		/// FIXME: SP wants this to be an attribute, DP not
		Real mNominalVoltageEnd2;
		/// Rated Apparent Power [VA]
		/// FIXME: SP wants this to be an attribute, DP not
		Real mRatedPower;
	public:
		/// Complex transformer ratio
		Attribute<Complex>::Ptr mRatio;
		/// Resistance [Ohm]
		Attribute<Real>::Ptr mResistance;
		/// Inductance [H]
		Attribute<Real>::Ptr mInductance;
		///
		void setParameters(Real nomVoltageEnd1, Real nomVoltageEnd2, Real ratioAbs, Real ratioPhase, Real resistance, Real inductance) {
			mNominalVoltageEnd1 = nomVoltageEnd1;
			mNominalVoltageEnd2 = nomVoltageEnd2;
			**mRatio = std::polar<Real>(ratioAbs, ratioPhase);
			**mResistance = resistance;
			**mInductance = inductance;
		}
	};
}
}
}
