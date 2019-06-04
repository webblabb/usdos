#include "Control_resource.h"

Control_resource::Control_resource(resourceLocation locInfo, 
	int startLevel_in)
	:
	x_coordinate(locInfo.x),
	y_coordinate(locInfo.y),
	fips(locInfo.fips),
	state(locInfo.stateAbbrev),
	capacity(startLevel_in)
{
}

Control_resource::Control_resource()
	:
	x_coordinate(0),
	y_coordinate(0),
	fips(""),
	state(""),
	capacity(0)
{
}

Control_resource::~Control_resource()
{
}

int Control_resource::get_dailyLimit() const
{
// Can make this more complex to check if fixed or value from file should be used
// or generate value from normal distribution. Hence no inlining. But for now...
	return std::get<0>(fixedDailyLimit);
}
