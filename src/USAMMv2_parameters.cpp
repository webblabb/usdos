#include "USAMMv2_parameters.h"

#include <string>
#include <typeinfo>
#include <math.h>

#include "File_manager.h"
#include "shared_functions.h"
#include "Farm.h"
#include "County.h"

USAMMv2_parameters::USAMMv2_parameters()
{

}

USAMMv2_parameters::USAMMv2_parameters(const Parameters* p, Farm_type* ft) :
    p(p), farm_type(ft)
{
    species_index = farm_type->get_index();
    species_name = farm_type->get_species();
    for(std::string s : p->USAMM_temporal_order)
    {
        config_time_periods.insert(s);
    }
    o_covariates_loaded = false;
    d_covariates_loaded = false;

    //Read and store county-level origin covariates
    initialize_covariates(p->USAMM_ocov_files.at(species_index),
                          county_ocov_map, ocov_county_par_names,
                          o_covariates_loaded);
    //...and destination covariates
    initialize_covariates(p->USAMM_dcov_files.at(species_index),
                          county_dcov_map, dcov_county_par_names,
                          d_covariates_loaded);
    //Read and store usamm parameter values.
    initialize_parameters(p->USAMM_parameter_files.at(species_index));

    //Read and store county-level superspreader /-shipper covariates.
    initialize_supernode_cov(p->USAMM_supernode_files.at(species_index),
                             county_ocov_map, county_dcov_map);

    use_raw_parameters = false;
}

USAMMv2_parameters::~USAMMv2_parameters()
{
    //dtor
}

void USAMMv2_parameters::initialize_parameters(std::string fname)
{
    o_frm_active = false;
    o_feedl_active = false;
    o_mkt_active = false;
    d_frm_active = false;
    d_feedl_active = false;
    d_mkt_active = false;
    std::vector<std::string> header_vector;
    std::vector<std::string> data_vector;
    std::ifstream f(fname);
    std::stringstream generation_data_ss;
    if(f.is_open())
    {
        n_lines_in_posterior_file = p->USAMM_post_lengths[species_index]; //Member var.
        skipBOM(f);
        std::string header;
        std::getline(f, header);
        generation_data_ss << header << std::endl;
        header_vector = split(header, '\t');
        std::string data;
        do
        {
            get_one_random_line(f, data, true);
        } while(data.compare("") == 0);

        generation_data_ss << data;
        generation_string = generation_data_ss.str();
        data_vector = split(data, '\t');
        f.close();
        for(size_t col_i = 0; col_i < header_vector.size(); col_i++)
        {
            double this_val = std::stod(data_vector.at(col_i));
            std::vector<std::string> colname_vector = split(header_vector.at(col_i), '_');

            if(colname_vector.size() == 3 and colname_vector.at(2).compare("like") != 0)
            {
                //std, kurt, n or s variable (state, time, parameter name)
                // or covariate variable (type, time, parameter name)
                std::string temporal_name = colname_vector.at(1);
                if(!temporal_name_exists(temporal_name))
                {
                    std::cout << "The temporal identifier \"" << temporal_name << "\" in"
                              << p->USAMM_parameter_files.at(species_index) << " was not found "
                              << "among the temporal ordering in the config file (option 45). "
                              << "Exiting..." << std::endl;
                    exit(EXIT_FAILURE);
                }

                USAMM_time_periods.insert(temporal_name);
                if(colname_vector.at(0).compare("o") == 0)
                {
                    //Origin covariate
                    if(colname_vector.at(2).compare("FarmExp") == 0 or
                       colname_vector.at(2).compare("FeedlotExp") == 0 or
                       colname_vector.at(2).compare("MarketExp") == 0 or
                       colname_vector.at(2).compare("FarmCoeff") == 0 or
                       colname_vector.at(2).compare("FeedlotCoeff") == 0 or
                       colname_vector.at(2).compare("MarketCoeff") == 0)
                    {
                        ocov_prem_par_map[temporal_name][colname_vector.at(2)] = this_val;

                        if(colname_vector.at(2).compare("FarmExp") == 0)
                        {
                            o_frm_active = true;
                        }
                        else if(colname_vector.at(2).compare("FeedlotExp") == 0 or
                                colname_vector.at(2).compare("FeedlotCoeff") == 0)
                        {
                            o_feedl_active = true;
                        }
                        else if(colname_vector.at(2).compare("MarketExp") == 0 or
                                colname_vector.at(2).compare("MarketCoeff") == 0)
                        {
                            o_mkt_active = true;
                        }
                    }
                    else
                    {
                        ocov_county_par_map[temporal_name][colname_vector.at(2)] = this_val;
                    }
                }
                else if(colname_vector.at(0).compare("i") == 0)
                {
                    //Destination covariate
                    if(colname_vector.at(2).compare("FarmExp") == 0 or
                       colname_vector.at(2).compare("FeedlotExp") == 0 or
                       colname_vector.at(2).compare("MarketExp") == 0 or
                       colname_vector.at(2).compare("FarmCoeff") == 0 or
                       colname_vector.at(2).compare("FeedlotCoeff") == 0 or
                       colname_vector.at(2).compare("MarketCoeff") == 0)
                    {
                        dcov_prem_par_map[temporal_name][colname_vector.at(2)] = this_val;

                        if(colname_vector.at(2).compare("FarmExp") == 0)
                        {
                            d_frm_active = true;
                        }
                        else if(colname_vector.at(2).compare("FeedlotExp") == 0 or
                                colname_vector.at(2).compare("FeedlotCoeff") == 0)
                        {
                            d_feedl_active = true;
                        }
                        else if(colname_vector.at(2).compare("MarketExp") == 0 or
                                colname_vector.at(2).compare("MarketCoeff") == 0)
                        {
                            d_mkt_active = true;
                        }
                    }
                    else
                    {
                        dcov_county_par_map[temporal_name][colname_vector.at(2)] = this_val;
                    }
                }
                else
                {
                    //Some state specific variable.
                    int state_id = std::stoi(colname_vector.at(0));
                    if(colname_vector.at(2).compare("std") == 0)
                    {
                        std_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("kurt") == 0)
                    {
                        kurt_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("a") == 0)
                    {
                        a_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("b") == 0)
                    {
                        b_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("N") == 0)
                    {
                        N_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("lambda") == 0)
                    {
                        lambda_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("beta") == 0)
                    {
                        s_map[state_id][temporal_name] = this_val;
                    }
                    else if(colname_vector.at(2).compare("like") == 0)
                    {
                        //Ignore this case
                        continue;
                    }
                    else
                    {
                        //Unknown parameter name.
                        std::cout << "Unknown parameter name " << colname_vector.at(2) << " in column "
                                  << header_vector.at(col_i) << " in " << fname << ". Exiting..."
                                  << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
            else
            {
                //Uninteresting parameter, skip.
                continue;
            }
        }
    }
    else
    {
        std::cout << "Failed to open " << fname << ". Exting..." << std::endl;
        exit(EXIT_FAILURE);
    }
    if(USAMM_time_periods != config_time_periods)
    {
        std::cout << "The temporal components of the parameters in "
                  << p->USAMM_parameter_files.at(species_index) << ":" << std::endl;
        for(auto& s : USAMM_time_periods)
        {
          std::cout << "\t" << s;
        }
        std::cout << std::endl << "...does not match those given in the config file"
                  << " (option 45):";
        for(auto& s : config_time_periods)
        {
          std::cout << "\t" << s;
        }
        std::cout << std::endl << "Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

void USAMMv2_parameters::initialize_covariates(std::string fname, str_vec_map& cov_map,
                                             std::vector<std::string>& header, bool& covariates_loaded)
{
    if(fname.compare("*") != 0)
    {
        std::ifstream f(fname);
        if(f.is_open())
        {
            skipBOM(f);
            std::string header_str;
            std::getline(f, header_str);
            header = split(header_str, '\t');
            header.erase(header.begin());
            std::string line;
            std::vector<std::string> line_vector;
            while(std::getline(f, line))
            {
                line_vector = split(line, '\t');
                if(!line_vector.empty())
                {
                    std::string fips_id = line_vector.at(0);
                    cov_map[fips_id].reserve(line_vector.size());
                    for(size_t i = 1; i < line_vector.size(); i++)
                    {
                        cov_map[fips_id].push_back(std::stod(line_vector.at(i)));
                    }
                }
            }
            covariates_loaded = true;
        }
        else
        {
            std::cout << "Failed to open " << fname << ". Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
}

void USAMMv2_parameters::initialize_supernode_cov(std::string fname, str_vec_map& ocov_map,
                                                str_vec_map& dcov_map)
{
    if(fname.compare("*") != 0)
    {
        supernodes_on = true;
        std::ifstream f(fname);
        if(f.is_open())
        {
            skipBOM(f);
            std::string header_str;
            std::getline(f, header_str);
            std::string line;
            std::vector<std::string> line_vector;
            ocov_county_par_names.push_back(std::string("SuperShipper"));
            dcov_county_par_names.push_back(std::string("SuperReceiver"));
            while(std::getline(f, line))
            {
                line_vector = split(line, '\t');
                if(!line_vector.empty())
                {
                    std::string fips_id = line_vector.at(0);
                    ocov_map[fips_id].push_back(std::stod(line_vector.at(1)));
                    dcov_map[fips_id].push_back(std::stod(line_vector.at(2)));
                }
            }
        }
        else
        {
            std::cout << "Failed to open " << fname << ". Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        supernodes_on = false;
    }
}

Farm_type* USAMMv2_parameters::get_farm_type()
{
    return farm_type;
}

double USAMMv2_parameters::get_std(int state_code, std::string time_period)
{
    return std_map.at(state_code).at(time_period);
}

double USAMMv2_parameters::get_kurt(int state_code, std::string time_period)
{
    return kurt_map.at(state_code).at(time_period);
}

double USAMMv2_parameters::get_a(int state_code, std::string time_period)
{
    //If the value have been calculated, return that.
    if(a_map.find(state_code) != a_map.end())
    {
        if(a_map.at(state_code).find(time_period) != a_map.at(state_code).end())
        {
            return a_map.at(state_code).at(time_period);
        }
    }
    //Can we skip the transformation from std and kurt to a and b?
    if(use_raw_parameters)
    {
        return get_std(state_code, time_period);
    }
    //Otherwise calculate a & b and save.
    std::vector<double> a_b_vec = calculate_a_and_b(get_std(state_code, time_period),
                                                get_kurt(state_code, time_period));
    a_map[state_code][time_period] = a_b_vec.at(0);
    b_map[state_code][time_period] = a_b_vec.at(1);
    return a_b_vec.at(0);
}

double USAMMv2_parameters::get_b(int state_code, std::string time_period)
{
    //If the value have been calculated, return that.
    if(b_map.find(state_code) != b_map.end())
    {
        if(b_map.at(state_code).find(time_period) != b_map.at(state_code).end())
        {
            return b_map.at(state_code).at(time_period);
        }
    }
    //Can we skip the transformation from std and kurt to a and b?
    if(use_raw_parameters)
    {
        return get_kurt(state_code, time_period);
    }
    //Otherwise calculate a & b and save.
    std::vector<double> a_b_vec = calculate_a_and_b(get_std(state_code, time_period),
                                                    get_kurt(state_code, time_period));
    a_map[state_code][time_period] = a_b_vec.at(0);
    b_map[state_code][time_period] = a_b_vec.at(1);
    return a_b_vec.at(1);
}

double USAMMv2_parameters::get_N(int state_code, std::string time_period)
{
    double temp_N;
    try
    {
        temp_N = N_map.at(state_code).at(time_period);
    }
    catch(const std::out_of_range& e)
    {
        std::cout << "Retrieveing parameter N from the USAMM parameters object failed." << std::endl
                  << "Maybe you are trying to run USDOS with USAMM v2 parameters," << std::endl
                  << "but with USDOS config setting no 41 set to 1?" << std::endl;
        exit(EXIT_FAILURE);
    }

    return temp_N * p->shipment_rate_factor;
}

double USAMMv2_parameters::get_lambda(int state_code, std::string time_period)
{
    double temp_lambda;
    try
    {
        temp_lambda = lambda_map.at(state_code).at(time_period);
    }
    catch(const std::out_of_range& e)
    {
        std::cout << "Retrieveing parameter lambda from the USAMM parameters object failed." << std::endl
                  << "Maybe you are trying to run USDOS with USAMM v1 parameters," << std::endl
                  << "but with USDOS config setting no 41 set to something other than 1?" << std::endl;
        exit(EXIT_FAILURE);
    }
    return temp_lambda * p->shipment_rate_factor;
}

double USAMMv2_parameters::get_s(int state_code, std::string time_period)
{
    return s_map.at(state_code).at(time_period);
}

std::vector<double> USAMMv2_parameters::get_county_o_covs(County* c)
{
    if(o_covariates_loaded)
    {
        std::string county_id = c->get_id();
        bool has_farms = c->get_n_premises(this->farm_type) > 0;
        bool has_covariates = county_ocov_map.find(county_id) != county_ocov_map.end();
        if(has_farms and has_covariates)
        {
            if(county_ocov_map.at(county_id).size() == ocov_county_par_names.size())
            {
                return county_ocov_map.at(county_id);
            }
            else
            {
                std::cout << "The county " << county_id << " has an incorrect number of origin " <<
                          "covariates (" << county_ocov_map.at(county_id).size() << " found, should be " <<
                          ocov_county_par_names.size() << ". This county is probably missing from either the " <<
                          "covariate or supernode files. All counties must be present in all the files " <<
                          "although their covariate values may be zero." << std::endl;
                exit(1);
            }
        }
        else if(has_farms and !has_covariates)
        {
            std::cout << "The county " << county_id << " has farms of type "
                      << this->farm_type->get_species() << " but no origin covariates " <<
                      "associated with that type. Setting those covariates to 0.0"
                      << std::endl;
            return std::vector<double>(ocov_county_par_names.size(), 0.0);
        }
    }
    return std::vector<double>(ocov_county_par_names.size(), 0.0); //If covariates not loaded, or no farms.
}

std::vector<double> USAMMv2_parameters::get_county_d_covs(County* c)
{
    if(d_covariates_loaded)
    {
        std::string county_id = c->get_id();
        bool has_farms = c->get_n_premises(this->farm_type) > 0;
        bool has_covariates = county_dcov_map.find(county_id) != county_dcov_map.end();
        if(has_farms and has_covariates)
        {
            if(county_dcov_map.at(county_id).size() == dcov_county_par_names.size())
            {
                return county_dcov_map.at(county_id);
            }
            else
            {
                std::cout << "The county " << county_id << " has an incorrect number of origin " <<
                          "covariates (" << county_dcov_map.at(county_id).size() << " found, should be " <<
                          dcov_county_par_names.size() << "). This county is probably missing from either the " <<
                          "covariate or supernode files. All counties must be present in all the files " <<
                          "although their covariate values may be zero." << std::endl;
                exit(1);
            }
        }
        else if(has_farms and !has_covariates)
        {
            std::cout << "The county " << county_id << " has farms of type "
                      << this->farm_type->get_species() << " but no destination covariates " <<
                      "associated with that type. Setting those covariates to 0.0"
                      << std::endl;
            return std::vector<double>(dcov_county_par_names.size(), 0.0);
        }
    }
    return std::vector<double>(dcov_county_par_names.size(), 0.0); //If covariates not loaded, or no farms.
}

std::vector<std::string> USAMMv2_parameters::get_county_ocov_names()
{
    return ocov_county_par_names;
}

std::vector<std::string> USAMMv2_parameters::get_county_dcov_names()
{
    return dcov_county_par_names;
}

double USAMMv2_parameters::get_ocov_county_par(std::string par_name, std::string time_period)
{
    double temp_cov;
    try
    {
        temp_cov = ocov_county_par_map.at(time_period).at(par_name);
    }
    catch(const std::out_of_range& e)
    {
        std::cout << "The covariate with the name " << par_name << " (" << time_period << ") "
                  << "could not be found in the USAMM parameter file. If you are running without "
                  << "origin covariates please set option 47 to * to turn them off." << std::endl;
        exit(EXIT_FAILURE);
    }
    return temp_cov;
}

double USAMMv2_parameters::get_dcov_county_par(std::string par_name, std::string time_period)
{
    double temp_cov;
    try
    {
        temp_cov = dcov_county_par_map.at(time_period).at(par_name);
    }
    catch(const std::out_of_range& e)
    {
        std::cout << "The covariate with the name " << par_name << " (" << time_period << ") "
                  << "could not be found in the USAMM parameter file. If you are running without "
                  << "destination covariates please set option 48 to * to turn them off." << std::endl;
        exit(EXIT_FAILURE);
    }
    return temp_cov;
}

double USAMMv2_parameters::get_ocov_prem_par(std::string par_name, std::string time_period)
{
    double temp_cov;
    try
    {
        temp_cov = ocov_prem_par_map.at(time_period).at(par_name);
    }
    catch(const std::out_of_range& e)
    {
        std::cout << "The covariate with the name " << par_name << " (" << time_period << ") "
                  << "could not be found in the USAMM parameter file. If you are running without "
                  << "origin covariates please set option 47 to * to turn them off." << std::endl;
        exit(EXIT_FAILURE);
    }
    return temp_cov;
}

double USAMMv2_parameters::get_dcov_prem_par(std::string par_name, std::string time_period)
{
    double temp_cov;
    try
    {
        temp_cov = dcov_prem_par_map.at(time_period).at(par_name);
    }
    catch(const std::out_of_range& e)
    {
        std::cout << "The covariate with the name " << par_name << " (" << time_period << ") "
                  << "could not be found in the USAMM parameter file. If you are running without "
                  << "destination covariates please set option 48 to * to turn them off." << std::endl;
        exit(EXIT_FAILURE);
    }
    return temp_cov;
}
bool USAMMv2_parameters::has_ocov_frm_par()
{
    return o_frm_active;
}

bool USAMMv2_parameters::has_ocov_feedl_par()
{
    return o_feedl_active;
}

bool USAMMv2_parameters::has_ocov_mkt_par()
{
    return o_mkt_active;
}

bool USAMMv2_parameters::has_dcov_frm_par()
{
    return d_frm_active;
}

bool USAMMv2_parameters::has_dcov_feedl_par()
{
    return d_feedl_active;
}

bool USAMMv2_parameters::has_dcov_mkt_par()
{
    return d_mkt_active;
}


std::string USAMMv2_parameters::get_generation_string()
{
    return generation_string;
}

std::vector<double> USAMMv2_parameters::calculate_a_and_b(double dx1, double R,
                                                        double x1, double x2) // x1, x2 has default values 0.5 and 0.05
{

    std::vector<double> a_b_vec =  {0.0, 0.0};
    if((p->shipment_kernel).compare("power_exp")==0){
        a_b_vec[1] = -std::log(std::log(1/x1) / std::log(1/x2)) / std::log(R);
        a_b_vec[0] = dx1 * std::pow(std::log(1/x1), -1/a_b_vec[1]);
    } else if ((p->shipment_kernel).compare("one_minus_exp")==0){
        a_b_vec[1] = std::log(std::log(1 - x2) / std::log(1 - x1)) / std::log(R);
        a_b_vec[0] = dx1 * std::pow(-std::log(1 - x1), -(1/a_b_vec[1]));
    } else if ((p->shipment_kernel).compare("local")==0){
        a_b_vec[1] = (std::log(1/x2 - 1) - std::log(1/x1 - 1)) * (1/std::log(R));
        a_b_vec[0] = dx1 / std::pow(1/x1 - 1, 1/a_b_vec[1]);
    } else {
        std::cout << "calculate_a_and_b does not work for kernel = " <<
            p->shipment_kernel << std::endl;
        std::cout <<"file: "<<__FILE__ << "\t line: " << __LINE__<< std::endl;
        exit(EXIT_FAILURE);
    }
    return a_b_vec;
}

bool USAMMv2_parameters::get_parameters(double stdev, double kurt, double const& epsilon,
                                      std::vector<double>& ab_vec)
{
    double log_right_limit = 100;
    double log_left_limit = 0.00001;
    bool status = true; // changes if anything goes wrong

    ab_vec[1] = log_int_hal(log_left_limit, log_right_limit, epsilon,
                            std::log(kurt), status);  // this solves the equations for a and b... L65... integration by halves.
    ab_vec[0] = std::exp(return_log_a(ab_vec[1], stdev));

    return status;
}

double USAMMv2_parameters::log_int_hal(double n0,double n1, double epsilon, double log_kurtosis,
                                     bool& status)
{
    double check;
    status = true;
    if (log_kurtosisf(n0, log_kurtosis, status) < 0 &&
        log_kurtosisf(n1, log_kurtosis, status) > 0)
    {
        while(std::fabs(n1-n0) > epsilon && status)
        {
            check = log_kurtosisf(n0 + (n1-n0)/2, log_kurtosis, status);
            if(check > 0)
            {
                n1 = n0 + (n1 - n0)/2;
            }
            else
            {
                n0 = n0 + (n1 - n0)/2;
            }
        }

        if(status)
        {
            return 1/n0;
        }
        else
        {
            std::cout << "limits: n0 = " << n0 << ", n1 = " << n1 << std::endl;
            std::cout << "log(kurtosis) = " << log_kurtosis << std::endl;
            status = false;
            return -1;
        }
    }
    else
    {
        //error
        std::cout << "limits: n0 = " << n0 << ", n1 = " << n1 << std::endl;
        std::cout << "log(kurtosis) = " << log_kurtosis << std::endl;
        status = false;
        return -1;
    }
}

double USAMMv2_parameters::log_kurtosisf(double n, double log_kurtosis, bool& status)
{
    double res = 0.0;
    try
    {
        res =  std::lgamma(6*n)+(std::lgamma(2*n))-2*(std::lgamma(4*n))-log_kurtosis;
    }
    catch(std::exception &e)
    {
        std::cout << "Error in line = " << __LINE__ << " file = " << __FILE__ << std::endl;
    	std::cout << typeid(e).name() << ":" << e.what() << std::endl;
        status = false;
    }
    return res;
}

double USAMMv2_parameters::return_log_a(double n, double sd)
{
    return std::log(sd) + (std::lgamma(2/n) - std::lgamma(4/n))/2;
}

bool USAMMv2_parameters::get_one_random_line(std::ifstream& f, std::string& res_line, bool header)
{
    res_line.clear();
    f.seekg(0, f.beg);
    size_t get_this_line = 0;
    int start_at = 0;
    if(header)
    {
        start_at += 1;
    }

    if(f.is_open())
    {
        while(res_line.empty())
        {
            get_this_line = rand_int(start_at, n_lines_in_posterior_file);
            f.unsetf(std::ios_base::skipws);
            for(size_t i = 0; i < get_this_line; i++)
            {
             f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            std::getline(f, res_line);
            f.clear();
            f.seekg(0, f.beg);
        }
        return true;
    }
    return false;
}

bool USAMMv2_parameters::temporal_name_exists(std::string temporal_name)
{
    bool exists = false;
    for(auto& this_name : p->USAMM_temporal_order)
    {
        if(temporal_name.compare(this_name) == 0)
        {
            exists = true;
            break;
        }
    }
    return exists;
}
