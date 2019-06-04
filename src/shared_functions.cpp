#include "shared_functions.h"
#include "Farm.h"
#include <iterator>

double uniform_rand()
{
	static std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
	static unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return unif_dist(generator);
}

double normal_rand()
{
	static std::normal_distribution<double> norm_dist(0,1);
	static unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return norm_dist(generator);
}

int rand_int(int lo, int hi)
{
	std::uniform_int_distribution<int> unif_dist(lo, hi);
	unsigned long long int seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return unif_dist(generator);
}

/// Used in gridding (binomial method) to determing # of infected farms
/// \param[in]	N	Number of trials (premises in cell)
///	\param[in]	prob	Probability of success for all trials (pmax for all premises in cell)
int draw_binom(int N, double prob)
// draw from a binomial distribution based on N farms and prob (calc with focalInf & gridKern)
{
	std::binomial_distribution<int> binom_dist(N,prob);
	static unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return binom_dist(generator);
}

/// Used to generate the number of shipments that originate from a state in a given timestep.
/// \param[in] lambda Rate of distribution.
int draw_poisson(double lambda)
{
    std::poisson_distribution<int> p_dist(lambda);
	static unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
	static std::mt19937 generator(seed); //Mersenne Twister pseudo-random number generator. Generally considered research-grade.
	return p_dist(generator);
}


unsigned int generate_distribution_seed()
{
    return std::chrono::system_clock::now().time_since_epoch().count();
}

size_t get_day_of_year(size_t current_timestep, size_t start_day)
{
    return ((current_timestep-2 + start_day) % 365) + 1;
}

/// Based on algorithm described at http://www.johndcook.com/blog/cpp_expm1/
/// \param[in]	x	Exponent value
double oneMinusExp(double x)
{
	if (x == 0){
		return 0;
	} else if (std::abs(x) < 1e-5){
		return -(x + 0.5*x*x);
	} else {
		return -(exp(x) - 1.0);
	}
}

/// Used by Status_manager to determine status level progression times
///	\param[in]	params	Tuple of two doubles, the mean and variance
int normDelay(std::tuple<double, double>& params)
{
	double mean = std::get<0>(params);
	double var = std::get<1>(params);
	double normDraw = normal_rand()*var + mean; // scale # drawn from N(0,1) to mean and variance
	int draw = (int)(normDraw+0.5); // round up to nearest day
	if(draw<1){draw = 1;}
	return draw;
}

// used in reading in files
std::vector<std::string>
	split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim))
	{
		if (!item.empty())
		{
			elems.emplace_back(item);
		}
    }
    return elems;
}

std::vector<std::string>
	split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// Skips the Byte Order Mark (BOM) that defines UTF-8 in some text files.
// Credit to user 'Contango' at stackoverflow.com
void skipBOM(std::ifstream &in)
{
    char test[3] = {0};
    in.read(test, 3);
    if ((unsigned char)test[0] == 0xEF &&
        (unsigned char)test[1] == 0xBB &&
        (unsigned char)test[2] == 0xBF)
    {
        return;
    }
    in.seekg(0);
}

std::vector<double> stringToNumVec(std::string& toConvert)//, std::string delimiter/* = "," */)
{
	std::vector<double> output;
    std::string delim = ",";//delimiter; // defaults to ","
    std::string substring;
    double temp;

    auto start = 0;
    auto end = toConvert.find(delim); // first delimiter occurrence
    while (end != std::string::npos)
    {
        substring = toConvert.substr(start, end - start);
        substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
        temp = stringToNum<double>(substring); // convert substring to double
        output.emplace_back(temp); // add double to vector
        start = end + delim.length(); // set end(+delimiter) as new start point
        end = toConvert.find(delim, start); // find new endpoint
    }
    substring = toConvert.substr(start, end); // last substring
    substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
	temp = stringToNum<double>(substring); // convert substring to double
    output.emplace_back(temp); // add double to vector

    return output;
}

std::vector<std::string> semicolonStringToStringVec(std::string& toConvert)
{
	std::vector<std::string> output;
    std::string delim = ";";
    std::string substring;

    auto start = 0;
    auto end = toConvert.find(delim); // first delimiter occurrence
    while (end != std::string::npos)
    {
        substring = toConvert.substr(start, end - start);
        substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
        output.emplace_back(substring); // add string to vector
        start = end + delim.length(); // set end(+delimiter) as new start point
        end = toConvert.find(delim, start); // find new endpoint
    }
    substring = toConvert.substr(start, end); // last substring
    substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
    output.emplace_back(substring); // add string to vector

    return output;
}


std::vector<int> stringToIntVec(std::string& toConvert)//, std::string delimiter/* = "," */)
{
	std::vector<int> output;
    std::string delim = ",";//delimiter;
    std::string substring;
    int temp;

    auto start = 0;
    auto end = toConvert.find(delim); // first delimiter occurrence
    while (end != std::string::npos)
    {
        substring = toConvert.substr(start, end - start);
        substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
        temp = stringToNum<int>(substring); // convert substring to double
        output.emplace_back(temp); // add double to vector
        start = end + delim.length(); // set end(+delimiter) as new start point
        end = toConvert.find(delim, start); // find new endpoint
    }
    substring = toConvert.substr(start, end); // last substring
    substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
	temp = stringToNum<int>(substring); // convert substring to double
    output.emplace_back(temp); // add double to vector

    return output;
}

std::vector<std::string> stringToStringVec(std::string& toConvert)
{
	std::vector<std::string> output;
    std::string delim = ",";
    std::string substring;

    auto start = 0;
    auto end = toConvert.find(delim); // first delimiter occurrence
    while (end != std::string::npos)
    {
        substring = toConvert.substr(start, end - start);
        substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
        output.emplace_back(substring); // add string to vector
        start = end + delim.length(); // set end(+delimiter) as new start point
        end = toConvert.find(delim, start); // find new endpoint
    }
    substring = toConvert.substr(start, end); // last substring
    substring.erase(std::remove_if(substring.begin(), substring.end(), isspace), substring.end()); // remove whitespace
    output.emplace_back(substring); // add string to vector

    return output;
}

std::string vecToCommaSepString(const std::vector<int> vecToPaste)
{
	std::string output;
	char temp[10];
	for (auto& v:vecToPaste){
		sprintf(temp, "%d,", v);
		output += temp;
	}
	output.pop_back(); // remove last comma
	return output;
}

std::string vecToCommaSepString(const std::vector<std::string> vecToPaste)
{
	std::string output;
	for (auto& v:vecToPaste){
		output += v;
		output += ",";
	}
	output.pop_back();  // remove last comma
	return output;
}

void addItemTab(std::string& outString, int toAdd){
	char temp[20];
	sprintf(temp, "%d\t", toAdd);
	outString += temp;
}

void addItemTab(std::string& outString, double toAdd){
	char temp[50];
	sprintf(temp, "%.2f\t", toAdd);
	outString += temp;
}

void addItemTab(std::string& outString, std::string toAdd){
	outString += toAdd;
	outString +="\t";
}

/// Adds a formatted string (including tabs, newline) to an output file
void printLine(std::string& outputFile, std::string& printString)
{
	std::ofstream outfile;
	outfile.open(outputFile, std::ios::app); // append to existing file
	if(!outfile){std::cout<<"File "<<outputFile<<" not open."<<std::endl;}
	outfile << printString;
	std::flush(outfile);
	outfile.close();
}

unsigned int get_n_lines(std::ifstream& f)
{
//    unsigned int n_lines = 0;
//    std::string tempstring;
//    f.seekg(0, f.beg);
//    while(!f.eof())
//    {
//        f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
//        n_lines++;
//    }
//    f.clear();
//    f.seekg(0, f.beg);

    f.unsetf(std::ios_base::skipws);
    unsigned n_lines = std::count(
        std::istream_iterator<char>(f),
        std::istream_iterator<char>(),
        '\n');
    f.clear();
    f.seekg(0, f.beg);
    return n_lines;
}
