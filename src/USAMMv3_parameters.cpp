#include "USAMMv3_parameters.h"

#include "County.h"
#include "Farm.h"
#include "File_manager.h"

#include <math.h>
#include <string>
#include <vector>

std::map<Farm_type*, std::map<Prem_class*, std::set<int>>> USAMMv3_parameters::size_bins;

USAMMv3_parameters::USAMMv3_parameters(const Parameters* parameters, Farm_type* fty,
                                       std::set<Prem_class*> pclasses,
                                       std::map<Prem_class*, double> avg_prem_sizes) :
    config(parameters), fty(fty), avg_prem_sizes(avg_prem_sizes)
{
    prem_classes.resize(pclasses.size(), nullptr);
    for(Prem_class* pcl : pclasses)
    {
        prem_classes[pcl->idx] = pcl;
    }

    for(std::string s : config->USAMM_temporal_order)
    {
        config_time_periods.insert(s);
    }

    for(Prem_class* pcl : prem_classes)
    {
        pclass_tag_to_idx_map[pcl->tag] = pcl->idx;
    }
    fty_idx = fty->get_index();
    fty_name = fty->get_species();

    o_cov_loaded = false;
    d_cov_loaded = false;
    //Read and store county-level origin covariates
    initialize_covariates(config->USAMM_ocov_files.at(fty_idx),
                          ocov_values, ocov_par_names,
                          o_cov_loaded);
    //...and destination covariates
    initialize_covariates(config->USAMM_dcov_files.at(fty_idx),
                          dcov_values, dcov_par_names,
                          d_cov_loaded);

    //Read parameters.
    post_fname = config->USAMM_parameter_files.at(fty_idx);
    n_lines_in_posterior_file = config->USAMM_post_lengths.at(fty_idx);
    sample_post(post_fname);

    //If this is dairy, the bin sizes needs to be copied from beef.
    if(fty->get_species() == "dairy")
    {
        Farm_type* beef_ft = nullptr;
        for(auto& fty_size_map_pair : USAMMv3_parameters::size_bins)
        {
            if(fty_size_map_pair.first->get_species() == "beef")
            {
                beef_ft = fty_size_map_pair.first;
                break;
            }
        }
        Prem_class* mkt_pcl = prem_classes[pclass_tag_to_idx_map["Mkt"]];
        USAMMv3_parameters::size_bins[fty][mkt_pcl] = USAMMv3_parameters::size_bins[beef_ft][mkt_pcl];
    }
    evaluate_size_bins();
}

USAMMv3_parameters::~USAMMv3_parameters()
{
}

std::vector<double> USAMMv3_parameters::get_county_o_covs(County* c)
{
    if(o_cov_loaded)
    {
        std::string county_id = c->get_id();
        bool has_farms = c->get_n_premises(this->fty) > 0;
        bool has_covariates = ocov_values.find(county_id) != ocov_values.end();
        if(has_farms and has_covariates)
        {
            if(ocov_values.at(county_id).size() == ocov_par_names.size())
            {
                return ocov_values.at(county_id);
            }
            else
            {
                std::cout << "The county " << county_id << " has an incorrect number of origin " <<
                          "covariates (" << ocov_values.at(county_id).size() << " found, should be " <<
                          ocov_par_names.size() << ". This county is probably missing from either the " <<
                          "covariate or supernode files. All counties must be present in all the files " <<
                          "although their covariate values may be zero." << std::endl;
                exit(1);
            }
        }
        else if(has_farms and !has_covariates)
        {
            std::cout << "The county " << county_id << " has farms of type "
                      << this->fty->get_species() << " but no origin covariates " <<
                      "associated with that type. Setting those covariates to 0.0"
                      << std::endl;
            return std::vector<double>(ocov_par_names.size(), 0.0);
        }
    }
    return std::vector<double>(ocov_par_names.size(), 0.0); //If covariates not loaded, or no farms.
}

std::vector<double> USAMMv3_parameters::get_county_d_covs(County* c)
{
    if(d_cov_loaded)
    {
        std::string county_id = c->get_id();
        bool has_farms = c->get_n_premises(this->fty) > 0;
        bool has_covariates = dcov_values.find(county_id) != dcov_values.end();
        if(has_farms and has_covariates)
        {
            if(dcov_values.at(county_id).size() == dcov_par_names.size())
            {
                return dcov_values.at(county_id);
            }
            else
            {
                std::cout << "The county " << county_id << " has an incorrect number of origin " <<
                          "covariates (" << dcov_values.at(county_id).size() << " found, should be " <<
                          dcov_par_names.size() << "). This county is probably missing from either the " <<
                          "covariate or supernode files. All counties must be present in all the files " <<
                          "although their covariate values may be zero." << std::endl;
                exit(1);
            }
        }
        else if(has_farms and !has_covariates)
        {
            std::cout << "The county " << county_id << " has farms of type "
                      << this->fty->get_species() << " but no destination covariates " <<
                      "associated with that type. Setting those covariates to 0.0"
                      << std::endl;
            return std::vector<double>(dcov_par_names.size(), 0.0);
        }
    }
    return std::vector<double>(dcov_par_names.size(), 0.0); //If covariates not loaded, or no farms.
}


double USAMMv3_parameters::get_c(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr.c_map.at(pcl_pair) * config->shipment_rate_factor;
}

double USAMMv3_parameters::get_kMu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr.kMu_map.at(pcl_pair);
}

double USAMMv3_parameters::get_kNu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr.kNu_map.at(pcl_pair);
}

double USAMMv3_parameters::get_phi_O(Prem_class* pcl_O, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return parstr.phiO_map.at(pcl_O->idx);
}

double USAMMv3_parameters::get_phi_D(Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return parstr.phiD_map.at(pcl_D->idx);
}

std::pair<double, double> USAMMv3_parameters::get_ab(std::string state, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return std::make_pair(parstr.a_map.at(state), parstr.b_map.at(state));
}

double USAMMv3_parameters::get_outflow(std::string state, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return parstr.outflow_map.at(state);
}

double USAMMv3_parameters::get_inflow(std::string state, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return parstr.inflow_map.at(state);
}

double USAMMv3_parameters::get_ocov(std::string ocov_name, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return parstr.ocov_map.at(ocov_name);
}

double USAMMv3_parameters::get_dcov(std::string dcov_name, std::string time_period)
{
    USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(time_period);
    return parstr.dcov_map.at(dcov_name);
}

void USAMMv3_parameters::initialize_covariates(std::string fname, str_vec_map& cov_map,
                                               std::vector<std::string>& header,
                                               bool& covariates_loaded)
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

void USAMMv3_parameters::evaluate_size_bins()
{
    for(std::string period : USAMM_time_periods)
    {
        origin_size_weights[period].resize(prem_classes.size(), std::vector<double>());
        destination_size_weights[period].resize(prem_classes.size(), std::vector<double>());
        for(Prem_class* pcl : prem_classes)
        {
            std::vector<int> bins(USAMMv3_parameters::size_bins.at(fty).at(pcl).begin(),
                                  USAMMv3_parameters::size_bins.at(fty).at(pcl).end());
            double avg_prem_size = avg_prem_sizes.at(pcl);
            for(size_t i=0; i<bins.size(); ++i)
            {
                double binned_size = double(bins[i]);
                double origin_size_weight = std::pow(binned_size / avg_prem_size, get_phi_O(pcl, period));
                double destination_size_weight = std::pow(binned_size / avg_prem_size, get_phi_D(pcl, period));

                origin_size_weights[period][pcl->idx].push_back(origin_size_weight);
                destination_size_weights[period][pcl->idx].push_back(destination_size_weight);

                origin_size_weights_lookup[period][pcl][bins[i]] = origin_size_weights[period][pcl->idx][i];
                dest_size_weights_lookup[period][pcl][bins[i]] = destination_size_weights[period][pcl->idx][i];
            }
        }
    }
}

size_t USAMMv3_parameters::get_one_random_line(std::ifstream& f, std::string& res_line, bool header)
{
    res_line.clear();
    f.seekg(0, f.beg);
    size_t get_this_line = 0;
    int start_at = n_metadata_lines;
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
        return get_this_line;
    }
    return -1;
}

bool USAMMv3_parameters::temporal_name_exists(std::string temporal_name)
{
    bool exists = false;
    for(auto& this_name : config->USAMM_temporal_order)
    {
        if(temporal_name.compare(this_name) == 0)
        {
            exists = true;
            break;
        }
    }
    return exists;
}

int USAMMv3_parameters::sample_post(std::string fpath)
{
    std::vector<std::string> metadata_vector;
    std::vector<std::string> header_vector;
    std::vector<std::string> data_vector;
    std::ifstream f(fpath);
    std::stringstream generation_data_ss;
    partstr_ptrs.reserve(100);
    if(f.is_open())
    {
        skipBOM(f);
        //Count the lines containing metadata.
        n_metadata_lines = 0;
        while(1)
        {
            std::string metadata_line;
            std::getline(f, metadata_line);
            if((*metadata_line.begin()) == '#')
            {
                ++n_metadata_lines;
            }
            else
            {
                break;
            }
        }
        f.seekg(0, f.beg); //Go back to beginning of file.

        //Skip the metadata.
        std::string temp_line;
        for(int i=0; i<n_metadata_lines; ++i)
        {
            std::getline(f, temp_line);
        }

        //Now ready to read header.
        std::string header;
        std::getline(f, header);
        generation_data_ss << header << std::endl;
        header_vector = split(header, '\t');
        std::string data;
        do
        {
            last_sampled_line = get_one_random_line(f, data, true);
        } while(data.compare("") == 0);

        generation_data_ss << data;
        generation_string = generation_data_ss.str();
        data_vector = split(data, '\t');
        f.close();
        for(size_t col_i = 0; col_i < header_vector.size(); col_i++)
        {
            double this_val = std::stod(data_vector.at(col_i));
            std::vector<std::string> colname_vector = split(header_vector.at(col_i), '_');

            if(colname_vector.size() == 3 and colname_vector.at(0) != "prior" and colname_vector.at(0) != "ll")
            {
                std::string parname = colname_vector[0];
                std::string period = colname_vector[2];
                if(!temporal_name_exists(period))
                {
                    std::cout << "The temporal identifier \"" << period << "\" in"
                              << config->USAMM_parameter_files.at(fty_idx) << " was not found "
                              << "among the temporal ordering in the config file (option 45). "
                              << "Exiting..." << std::endl;
                    exit(EXIT_FAILURE);
                }

                USAMM_time_periods.insert(period);

                if(parstruct_ptrs_by_period.find(period) == parstruct_ptrs_by_period.end())
                {
                    int period_idx = -1;
                    for(size_t i=0; i<config->USAMM_temporal_order.size(); ++i)
                    {
                        if(period == config->USAMM_temporal_order[i])
                        {
                            period_idx = int(i);
                            break;
                        }
                    }
                    if(period_idx < 0)
                    {
                        std::cout << "Failed to find period identifier " << period << " among the ones given in "
                                  << "the config file (option 45)." << std::endl;
                        exit(EXIT_FAILURE);

                    }
//                    USAMMv3_parstruct* parstr_ptr = new USAMMv3_parstruct(fty, period, period_idx);
                    parstruct_ptrs_by_period.emplace(period, USAMMv3_parstruct(fty, period, period_idx));
                    partstr_ptrs.push_back(&parstruct_ptrs_by_period.at(period));
                }
                USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(period);

                if(parname == "c" or parname == "kMu" or parname == "kNu" or parname == "kShape")
                {
                    std::vector<std::string> pcl_tags = split(colname_vector[1], '-');
                    int pcl_1_idx = pclass_tag_to_idx_map.at(pcl_tags[0]);
                    int pcl_2_idx = pclass_tag_to_idx_map.at(pcl_tags[1]);
                    if(parname == "c")
                    {
                        parstr.c_map[std::pair<int, int>(pcl_1_idx, pcl_2_idx)] = this_val / config->USAMM_temporal_n_timesteps[parstr.period_idx];
                    }
                    else if(parname == "kMu")
                    {
                        parstr.kMu_map[std::pair<int, int>(pcl_1_idx, pcl_2_idx)] = this_val;
                    }
                    else if(parname == "kNu" or parname == "kShape")
                    {
                        parstr.kNu_map[std::pair<int, int>(pcl_1_idx, pcl_2_idx)] = this_val;
                    }
                }
                else if(parname == "phiO")
                {
                    int pcl_idx = pclass_tag_to_idx_map.at(colname_vector[1]);
                    parstr.phiO_map[pcl_idx] = this_val;
                }
                else if(parname == "phiD")
                {
                    int pcl_idx = pclass_tag_to_idx_map.at(colname_vector[1]);
                    parstr.phiD_map[pcl_idx] = this_val;
                }
                else if(parname == "dhalf" or parname == "dhalfCVI")
                {
                    parstr.dhalf_map[colname_vector[1]] = this_val;
                }
                else if(parname == "hratio" or parname == "hratioCVI")
                {
                    parstr.hratio_map[colname_vector[1]] = this_val;
                }
                else if(parname == "dhalfCOM")
                {
                    parstr.dhalfCOM_map[colname_vector[1]] = this_val;
                }
                else if(parname == "hratioCOM")
                {
                    parstr.hratioCOM_map[colname_vector[1]] = this_val;
                }
                else if(parname == "inflow" or parname == "inflowCVI")
                {
                    parstr.inflow_map[colname_vector[1]] = this_val;
                }
                else if(parname == "outflow" or parname == "outflowCVI")
                {
                    parstr.outflow_map[colname_vector[1]] = this_val;
                }
                else if(parname == "inflowCOM")
                {
                    parstr.inflowCOM_map[colname_vector[1]] = this_val;
                }
                else if(parname == "outflowCOM")
                {
                    parstr.outflowCOM_map[colname_vector[1]] = this_val;
                }
                else if(parname == "ocov")
                {
                    parstr.ocov_map[colname_vector[1]] = this_val;
                }
                else if(parname == "dcov")
                {
                    parstr.dcov_map[colname_vector[1]] = this_val;
                }
                else if(parname == "z" or parname == "r")
                {
                    continue;
                }
                else
                {
                    //Unknown parameter name.
                    std::cout << "Unknown parameter name " << header_vector.at(col_i) <<
                                 " in " << fpath << ". Exiting..." << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    else
    {
        std::cout << "Failed to open " << fpath << ". Exting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    if(USAMM_time_periods != config_time_periods)
    {
        std::cout << "The temporal components of the parameters in "
                  << config->USAMM_parameter_files.at(fty_idx) << ":" << std::endl;
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
    n_time_periods = USAMM_time_periods.size();

    //Calculate and save kernel parameters a and b from dhalf and hratio.
    for(std::string period : USAMM_time_periods)
    {
        USAMMv3_parstruct& parstr = parstruct_ptrs_by_period.at(period);
        for(auto& state_dhalf_pair : parstr.dhalf_map)
        {
            std::string state = state_dhalf_pair.first;
            double dhalf = state_dhalf_pair.second;
            double hratio = parstr.hratio_map.at(state);
            std::pair<double, double> ab_pair = ab_fun(dhalf, hratio);
            parstr.a_map[state] = ab_pair.first;
            parstr.b_map[state] = ab_pair.second;
        }
    }
    return 0;
}

std::pair<double, double> USAMMv3_parameters::ab_fun(double dhalf, double hratio)
{
    static double x1 = 0.5;
    static double x2 = 0.05;
    static double one_over_x1_m1 = (1.0/x1) - 1.0;
    static double log_one_over_x1_m1 = std::log(one_over_x1_m1);
    static double log_one_over_x2_m1 = std::log((1.0/x2) - 1.0);

    double b = (log_one_over_x2_m1 - log_one_over_x1_m1) * (1.0 / std::log(hratio));
    double a = dhalf / std::pow(one_over_x1_m1, 1.0 / b);
    return std::make_pair(a, b);
}

int USAMMv3_parameters::farmSizeBinningFun(int p_original_size, Farm_type* ft, Prem_class* pcl)
{
    //Bins according to NASS categories for farm and feedlot.
    //Market bins are arbitrary.

    if(p_original_size < 1)
    {
        std::cout << "Error when binning the volume " << p_original_size << ". Exiting." << std::endl;
        exit(EXIT_FAILURE);
    }

    int binned_size = 0;
    if(pcl->tag == "Frm") //Farm
    {
        if(p_original_size < 10)
            { binned_size = 5; }
        else if(p_original_size >= 10 and p_original_size < 20)
            { binned_size = 15; }
        else if(p_original_size >= 20 and p_original_size < 50)
            { binned_size = 35; }
        else if(p_original_size >= 50 and p_original_size < 100)
            { binned_size = 75; }
        else if(p_original_size >= 100 and p_original_size < 200)
            { binned_size = 150; }
        else if(p_original_size >= 200 and p_original_size < 500)
            { binned_size = 350; }
        else if(p_original_size >= 500)
            { binned_size = 1200; } //This is the average of beef & dairy farms larger than 500 in flaps.
        else
        {
            std::cout << "Size binning failed with size = " <<
                         p_original_size << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(pcl->tag == "Fdl") //Feedlot
    {
        if(p_original_size < 20)
            { binned_size = 10; }
        else if(p_original_size >= 20 and p_original_size < 50)
            { binned_size = 35; }
        else if(p_original_size >= 50 and p_original_size < 100)
            { binned_size = 75; }
        else if(p_original_size >= 100 and p_original_size < 200)
            { binned_size = 150; }
        else if(p_original_size >= 200 and p_original_size < 500)
            { binned_size = 350; }
        else if(p_original_size >= 500 and p_original_size < 1000) //From the upper limit of this one I had to make something up since NASS is just 500+.
            { binned_size = 750; }
        else if(p_original_size >= 1000)
            { binned_size = 6500; } //This is the average size of feedlots larger than 1000 in flaps.
        else
        {
            std::cout << "Size binning failed with size = " <<
                         p_original_size << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(pcl->tag == "Mkt") //Market
    {
        p_original_size = p_original_size * 52; //Market size is weekly volume here.
        if(p_original_size < 10000)
            { binned_size = 5000; }
        else if(p_original_size >= 10000 and p_original_size < 20000)
            { binned_size = 15000; }
        else if(p_original_size >= 20000 and p_original_size < 40000)
            { binned_size = 30000; }
        else if(p_original_size >= 40000 and p_original_size < 60000)
            { binned_size = 50000; }
        else if(p_original_size >= 60000 and p_original_size < 80000)
            { binned_size = 70000; }
        else if(p_original_size >= 80000 and p_original_size < 100000)
            { binned_size = 90000; }
        else if(p_original_size >= 100000)
            { binned_size = 160000; } //This is the average of markets larger than 100000 in flaps.
        else
        {
            std::cout << "Size binning failed with size = " <<
                         p_original_size << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        std::cout << "Unknown prem type encountered in size binnning function: " <<
                     pcl->tag << std::endl;
        exit(EXIT_FAILURE);
    }

    if(USAMMv3_parameters::size_bins[ft][pcl].find(binned_size) ==
       USAMMv3_parameters::size_bins[ft][pcl].end())
    {
        USAMMv3_parameters::size_bins[ft][pcl].insert(binned_size);
    }

    return binned_size;
}

USAMMv3_parstruct::USAMMv3_parstruct(Farm_type* fty, std::string period, int period_idx) :
    fty(fty), period(period), period_idx(period_idx)
{

}
