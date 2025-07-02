//© 2024 Colorado State University
#include "Diagnostic_resource.h"

Diagnostic_resource::Diagnostic_resource(int startLevel_in)
	:
	diagnosticCapacity(startLevel_in)
{
}

Diagnostic_resource::Diagnostic_resource()
	:
	diagnosticCapacity(0)
{
}
Diagnostic_resource::~Diagnostic_resource()
{
}

int Diagnostic_resource::get_diagnosticDailyLimit() const
{
// Can make this more complex to check if fixed or value from file should be used
// or generate value from normal distribution. Hence no inlining. But for now...
	return std::get<0>(diagnosticFixedDailyLimit);
}

