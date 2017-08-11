#ifndef INDUCTOREMT_H
#define INDUCTOREMT_H

#include "BaseComponent.h"

namespace DPsim {

	class InductorEMT : public BaseComponent {
	protected:
		double mInductance;
		double mDeltav;		
		double mCurr;		
		double mCureq;
		double mGl;
		double mP;

	public:
		InductorEMT() { };
		InductorEMT(std::string name, int src, int dest, double inductance);

		void init(Real om, Real dt);
		void applySystemMatrixStamp(SystemModel& system);
		void applyRightSideVectorStamp(SystemModel& system) { }
		void step(SystemModel& system, Real time);
		void postStep(SystemModel& system);
	};
}
#endif