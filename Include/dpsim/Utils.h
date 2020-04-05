/* Copyright 2017-2020 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#if defined(__GNUC__) && !defined(__clang__)
  #include <cxxabi.h>
#endif

#include <cstdlib>
#include <list>
#include <vector>
#include <experimental/filesystem>

#include <dpsim/Timer.h>
#include <dpsim/Solver.h>
#include <cps/Logger.h>

namespace fs = std::experimental::filesystem;

namespace DPsim {

class CommandLineArgs {

protected:
	struct Argument {
		const char *name;
		int has_arg;
		int *flag;
		int val;
		const char *valdesc;
		const char *desc;
	};

	String mProgramName;
	std::vector<Argument> mArguments;

public:
	CommandLineArgs(int argc, char *argv[],
		/* Default settings */
		String name = "dpsim",
		Real dt = 0.001,
		Real d = 1,
		Real sf = 50,
		Int s = -1,
		CPS::Logger::Level ll = CPS::Logger::Level::info,
		Bool ss = false,
		Bool b = false,
		CPS::Domain sd = CPS::Domain::DP,
		Solver::Type st = Solver::Type::MNA
	);

	void showUsage();
	void showCopyright();

	double timeStep;
	double duration;
	double sysFreq;
	int scenario;

	CPS::Logger::Level logLevel;
	String name;

	bool startSynch;
	bool blocking;

	struct {
		CPS::Domain domain;
		Solver::Type type;
	} solver;

	DPsim::Timer::StartClock::time_point startTime;

	std::list<String> positional;
	std::list<fs::path> positionalPaths() const;

	std::map<String, Real> options;
};

namespace Utils {

String encodeXml(String& data);

template<typename T>
static CPS::String type(const CPS::String &stripPrefix = "CPS::") {
	Int status = 1;
	const char *mangled, *unmangled;

	mangled = typeid(T).name();

#ifdef _MSC_VER
	return CPS::String(mangled);
#else
	unmangled = abi::__cxa_demangle(mangled, NULL, NULL, &status);

	if (status)
		return mangled;
	else {
		CPS::String type = unmangled;

		delete unmangled;

		if (type.find(stripPrefix) == 0)
			type = type.substr(stripPrefix.size());

		return type;
	}
#endif
}

std::vector<std::string> tokenize(std::string s, char delimiter);

fs::path findFile(const fs::path &name,
	const fs::path &hint = fs::path(), const std::string &useEnv = std::string());

std::list<fs::path> findFiles(std::list<fs::path> filennames,
	const fs::path &hint, const std::string &useEnv = std::string());

}
}
