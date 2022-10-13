#include <stdio.h>
#include "shared_functions.h"
#include "Farm.h"
#include <iterator>
#include <gsl/gsl_randist.h>

class Gsl_rng_wrapper
{
    gsl_rng* r;
    public:
        Gsl_rng_wrapper()
        {
            unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
            const gsl_rng_type* rng_type = gsl_rng_default;
            r = gsl_rng_alloc(rng_type);
            gsl_rng_set(r, seed);
        }
        ~Gsl_rng_wrapper() { gsl_rng_free(r); }
        gsl_rng* get_r() { return r; }
};

double uniform_rand()
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_flat(r, 0.0, 1.0);
}

int rand_int(int lo, int hi)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return int( gsl_ran_flat(r, double(lo), double(hi)) );
}

double draw_norm(double mean, double stdev)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return mean + gsl_ran_gaussian(r, stdev);
}

/// Used in gridding (binomial method) to determing # of infected farms
/// \param[in]	N	Number of trials (premises in cell)
///	\param[in]	prob	Probability of success for all trials (pmax for all premises in cell)
int draw_binom(int N, double prob)
// draw from a binomial distribution based on N farms and prob (calc with focalInf & gridKern)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_binomial(r, prob, N);
}

/// Used to generate the number of shipments that originate from a state in a given timestep.
/// \param[in] lambda Rate of distribution.
int draw_poisson(double lambda)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_poisson(r, lambda);
}

double draw_gamma(double shape, double scale)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    return gsl_ran_gamma(r, shape, scale);
}

void draw_multinomial(int n, const std::vector<double>& weights, std::vector<unsigned int>& outcome)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
    outcome.clear();
    outcome.resize(weights.size(), 0);
    gsl_ran_multinomial(r, weights.size(), n, weights.data(), outcome.data());
}

void draw_multivariate_hypergeometric(int k, const std::vector<unsigned int>& bins, std::vector<unsigned int>& outcome)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }

    size_t nbins = bins.size();
    outcome.clear();
    outcome.resize(nbins, 0);
    int M = std::accumulate(bins.begin(), bins.end(), 0); //Total number of objects.
    int M_left = M;
    for(size_t i=0; i<nbins-1; ++i)
    {
        M_left = M_left - bins[i]; //Total number of objects in the bins that have not been sampled and which are not in the bin we are about to sample from.
        outcome[i] = gsl_ran_hypergeometric(r, bins[i], M_left, k);
        k -= outcome[i];
        if(k == 0)
        {
            break;
        }
    }
    outcome[nbins-1] = k; //No need to do last bin, it always gets the remaining.
}

/// Used to generate the success of the diagnostic tests
/// \param[in] alpha number of successful (positive) tests +1.
/// \param[in] beta number of unsuccessful (negative) tests +1.
double draw_beta(double alpha, double beta)
{
    static Gsl_rng_wrapper rng_w;
    static gsl_rng* r;
    if(r == nullptr){ r = rng_w.get_r(); }
	return gsl_ran_beta(r, alpha, beta);
}

unsigned int generate_distribution_seed()
{
    return std::chrono::system_clock::now().time_since_epoch().count();
}

size_t get_day_of_year(size_t days_after_start, size_t start_day)
{
    return ((days_after_start-2 + start_day) % 365) + 1;
}

size_t get_month_of_year(size_t day_of_year)
{
    if(day_of_year > 0 and day_of_year <= 31) { return 0; }
    else if(day_of_year > 31 and day_of_year <= 59) { return 1; }
    else if(day_of_year > 59 and day_of_year <= 90) { return 2; }
    else if(day_of_year > 90 and day_of_year <= 120) { return 3; }
    else if(day_of_year > 120 and day_of_year <= 151) { return 4; }
    else if(day_of_year > 151 and day_of_year <= 181) { return 5; }
    else if(day_of_year > 181 and day_of_year <= 212) { return 6; }
    else if(day_of_year > 212 and day_of_year <= 243) { return 7; }
    else if(day_of_year > 243 and day_of_year <= 273) { return 8; }
    else if(day_of_year > 273 and day_of_year <= 304) { return 9; }
    else if(day_of_year > 304 and day_of_year <= 334) { return 10; }
    else if(day_of_year > 334 and day_of_year <= 365) { return 11; }
    else
    {
        std::cout << "Day " << day_of_year << " is outside of the year." << std::endl;
        exit(EXIT_FAILURE);
    }

}

int get_quarter(size_t current_day_of_year)
{
    if(current_day_of_year <= 90) { return 0; }
    else if(current_day_of_year > 90 and current_day_of_year <= 181) { return 1; }
    else if(current_day_of_year > 181 and current_day_of_year <= 273) { return 2; }
    return 3;
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
	double stdev = std::get<1>(params);
	double normDraw = draw_norm(mean, stdev);
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
	if(!outfile)
    {
        std::cout<<"File "<<outputFile<<" not open."<<std::endl;
    }
	else
    {
        outfile << printString;
        std::flush(outfile);
        outfile.close();
    }
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

int read_table(std::string fname, char delim, bool has_header,
               std::vector<std::vector<std::string>>& out,
               char comment_char)
{
    std::ifstream f(fname, std::ifstream::in);
    std::vector<std::vector<std::string>> lines;
    if(f.is_open())
    {
        skipBOM(f);
        unsigned int n_lines = get_n_lines(f);
        lines.reserve(n_lines);

        std::string line;
        std::string header = "";
        if(has_header)
        {
            std::getline(f, header);
        }

        while(std::getline(f, line))
        {
            if(!line.empty())
            {
                //Remove comments
                std::stringstream ss;
                ss << line;
                std::getline(ss, line, comment_char); //Read up until comment char is found.

                //Split the rest into a vector.
                std::vector<std::string> line_vector = split(line, delim);

                //Remove whitespace from string.
                for(size_t i=0; i<line_vector.size(); i++)
                {
                    line_vector.at(i) = trim(line_vector.at(i));
                }

                if(!line_vector.empty())
                {
                    lines.push_back(line_vector);
                }
            }
        }
        out.swap(lines);
    }
    else
    {
        return 1;
    }
    return 0;
}


std::string trim(const std::string &s) //Removes whitespace around string.
{
    auto wsfront=std::find_if_not(s.begin(),s.end(),[](int c){return std::isspace(c);});
    auto wsback=std::find_if_not(s.rbegin(),s.rend(),[](int c){return std::isspace(c);}).base();
    return (wsback<=wsfront ? std::string() : std::string(wsfront,wsback));
}
