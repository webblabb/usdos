//© 2024 Colorado State University
/* Base class for county, state and grid cells. */

#ifndef REGION_H
#define REGION_H

#include <algorithm>
#include <vector>
#include <string>
#include <unordered_map>
#include <Point.h>

class Farm;

class Region
{
public:
    Region(std::string id);
    Region(std::string id, double x, double y);
    virtual ~Region();
    double distance_to(Region* target) const;
    double measure_distance(const Point* p1, const Point* p2) const;

    void set_position(double x, double y);

    std::string get_id(); //Inlined
    const Point* get_centroid(); //Inlined
    const std::string get_type() const; //inlined

protected:
    std::string type = "region/unchanged";
    std::string id = "unknown";
    Point centroid;
    bool is_set_id = false;
    bool is_set_position = false;
    bool region_initialized = false;

    virtual void set_initialized(bool& parameter);
    virtual void all_initialized(); //If all parameters that need to be initialized
                                    //have been init., set x_initialized to true.
    void not_initialized();
};

inline std::string Region::get_id()
{
// Region does not need to be completely initialized in order to return an id
    return id;
}

inline const std::string Region::get_type() const
{
   return type;
}

inline const Point* Region::get_centroid()
{
    if(!region_initialized)
        not_initialized();

    return &centroid;
}

class Region_status: public Region
{
	private:
		bool reported; // at least one reported premises
		std::vector<std::string> controlStatus; /// Vector of control statuses ever applied
		std::vector<double> probPreventExposure;
		std::vector<double> probPreventTransmission;

		std::unordered_map<std::string, int> start; /// Start times for each status. Assumes unique names among file, disease, and control statuses due to mapping by status name.
		std::unordered_map<std::string, int> end; /// End times for each status. Assumes unique names among file, disease, and control statuses due to mapping by status name.

	public:
		Region_status(Region*);
		~Region_status();
		void report(); //inlined
		void add_controlStatus(std::string); //inlined

		double get_probPreventExposure() const;
		void add_probPreventExposure(double); //inlined
		void rem_probPreventExposure(double);
		double get_probPreventTransmission() const;
		void add_probPreventTransmission(double); //inlined
		void rem_probPreventTransmission(double);

 		void set_start(const std::string, const int); //inlined - set start time for status
 		void set_end(const std::string, const int); //inlined - set end time for status
		int get_start(std::string) const;
 		int get_end(std::string) const;
		bool is_reported() const; // inlined
};

inline void Region_status::report()
{
	reported = true;
}
inline void Region_status::add_controlStatus(std::string status)
{
	controlStatus.emplace_back(status);
}
inline bool Region_status::is_reported() const
{
	return reported;
}
inline void Region_status::add_probPreventExposure(double eff)
{
	probPreventExposure.emplace_back(eff);
}
inline void Region_status::add_probPreventTransmission(double eff)
{
	probPreventTransmission.emplace_back(eff);
}
inline void Region_status::set_start(const std::string status, const int t)
{
	start[status] = t;
}
inline void Region_status::set_end(const std::string status, const int t)
{
	end[status] = t;
}
#endif // REGION_H
