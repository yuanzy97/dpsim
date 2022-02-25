/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
*                     EONERC, RWTH Aachen University
*
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, You can obtain one at https://mozilla.org/MPL/2.0/.
*********************************************************************************/

#pragma once

#include <dpsim/MNASolver.h>
#include <dpsim/DataLogger.h>
#include <dpsim/MNASolverEigenDense.h>
#ifdef WITH_SPARSE
#include <dpsim/MNASolverEigenSparse.h>
#ifdef WITH_KLU
#include <dpsim/MNASolverEigenKLU.h>
#include <dpsim/MNASolverEigenPartialKLU.h>
#endif
#ifdef WITH_NICSLU
#include <dpsim/MNASolverEigenNICSLU.h>
#include <dpsim/MNASolverEigenPartialNICSLU.h>
#endif
#endif
#ifdef WITH_CUDA
	#include <dpsim/MNASolverGpuDense.h>
#ifdef WITH_SPARSE
	#include <dpsim/MNASolverGpuSparse.h>
#ifdef WITH_MAGMA
	#include <dpsim/MNASolverGpuMagma.h>
#endif
#endif
#endif

namespace DPsim {

class MnaSolverFactory {
	public:
	/// \brief The implementations of the MNA solvers MnaSolver can support.
	///
	enum MnaSolverImpl {
		Undef = 0,
		EigenDense,
		EigenSparse,
		EigenKLU,
		EigenPartialKLU,
		EigenNICSLU,
		EigenPartialNICSLU,
		CUDADense,
		CUDASparse,
		CUDAMagma,
	};

	/// MNA implementations supported by this compilation
	static const std::vector<MnaSolverImpl> mSupportedSolverImpls(void) {
		static std::vector<MnaSolverImpl> ret = {
			EigenDense,
#ifdef WITH_SPARSE
			EigenSparse,
#ifdef WITH_KLU
			EigenKLU,
			EigenPartialKLU,
#endif //WITH_KLU
#ifdef WITH_NICSLU
			EigenNICSLU,
			EigenPartialNICSLU,
#endif //WITH_NICSLU
#endif //WITH_SPARSE
#ifdef WITH_CUDA
			CUDADense,
#ifdef WITH_SPARSE
			/* CUDASparse is not currently wokring correctly, DO NOT USE IT! */
			//CUDASparse,
#ifdef WITH_MAGMA
			CUDAMagma,
#endif //WITH_MAGMA
#endif //WITH_SPARSE
#endif //WITH_CUDA
		};
		return ret;
	}

	/// sovlerImpl: choose the most advanced solver implementation available by default
	template <typename VarType>
	static std::shared_ptr<MnaSolver<VarType>> factory(String name,
		CPS::Domain domain = CPS::Domain::DP,
		CPS::Logger::Level logLevel = CPS::Logger::Level::info,
		MnaSolverImpl implementation = mSupportedSolverImpls().back())
	{
		//To avoid regression we use EigenDense in case of undefined implementation
		if (implementation == MnaSolverImpl::Undef) {
			implementation = MnaSolverImpl::EigenDense;
		}
		CPS::Logger::Log log = CPS::Logger::get("MnaSolverFactory", CPS::Logger::Level::info, CPS::Logger::Level::info);

		switch(implementation) {
		case MnaSolverImpl::EigenDense:
			log->info("creating EigenDense solver implementation");
			return std::make_shared<MnaSolverEigenDense<VarType>>(name, domain, logLevel);
#ifdef WITH_SPARSE
		case MnaSolverImpl::EigenSparse:
			log->info("creating EigenSparse solver implementation");
			return std::make_shared<MnaSolverEigenSparse<VarType>>(name, domain, logLevel);
#ifdef WITH_KLU
		case MnaSolverImpl::EigenKLU:
			log->info("creating EigenKLU solver implementation");
			return std::make_shared<MnaSolverEigenKLU<VarType>>(name, domain, logLevel);
		case MnaSolverImpl::EigenPartialKLU:
			log->info("creating EigenPartialKLU solver implementation");
			return std::make_shared<MnaSolverEigenPartialKLU<VarType>>(name, domain, logLevel);
#endif
#ifdef WITH_NICSLU
		case MnaSolverImpl::EigenNICSLU:
			log->info("creating EigenNICSLU solver implementation");
			return std::make_shared<MnaSolverEigenNICSLU<VarType>>(name, domain, logLevel);
		case MnaSolverImpl::EigenPartialNICSLU:
			log->info("creating EigenPartialNICSLU solver implementation");
			return std::make_shared<MnaSolverEigenPartialNICSLU<VarType>>(name, domain, logLevel);
#endif
#endif
#ifdef WITH_CUDA
		case MnaSolverImpl::CUDADense:
			log->info("creating CUDADense solver implementation");
			return std::make_shared<MnaSolverGpuDense<VarType>>(name, domain, logLevel);
#ifdef WITH_SPARSE
		case MnaSolverImpl::CUDASparse:
			log->info("creating CUDASparse solver implementation");
			return std::make_shared<MnaSolverGpuSparse<VarType>>(name, domain, logLevel);
#ifdef WITH_MAGMA
		case MnaSolverImpl::CUDAMagma:
			log->info("creating CUDAMagma solver implementation");
			return std::make_shared<MnaSolverGpuMagma<VarType>>(name, domain, logLevel);
#endif
#endif
#endif
		default:
			throw CPS::SystemError("unsupported MNA implementation.");

		}
	}
};
}
