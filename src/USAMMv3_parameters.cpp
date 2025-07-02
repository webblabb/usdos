//© 2024 Colorado State University
#include "USAMMv3_parameters.h"

#include "shared_functions.h"
#include "County.h"
#include "Farm.h"
#include "File_manager.h"
#include "State.h"

#include <math.h>
#include <string>
#include <vector>

std::map<Farm_type*, std::map<Prem_class*, std::set<int>>> USAMMv3_parameters::size_bins;

USAMMv3_parameters::USAMMv3_parameters(const Parameters* parameters, Farm_type* fty,
                                       std::set<Prem_class*> pclasses,
                                       std::map<Prem_class*, double> avg_prem_sizes,
                                       std::vector<County*>& counties,
                                       const Vec_d_2d& county_distance_matrix) :
    config(parameters), fty(fty), avg_prem_sizes(avg_prem_sizes),
    counties(counties), county_distance_matrix(county_distance_matrix)
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
        std::cout << pcl->tag << std::endl;
        pclass_tag_to_idx_map[pcl->tag] = pcl->idx;
    }
    fty_idx = fty->get_index();
    fty_name = fty->get_species();

    common_parstruct = USAMMv3_common_parstruct(fty);

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
        Prem_class* fdl_pcl = prem_classes[pclass_tag_to_idx_map["Fdl"]];
        USAMMv3_parameters::size_bins[fty][mkt_pcl] = USAMMv3_parameters::size_bins[beef_ft][mkt_pcl];
        USAMMv3_parameters::size_bins[fty][fdl_pcl] = USAMMv3_parameters::size_bins[beef_ft][fdl_pcl];
    }

    //Read parameters.
    post_fname = config->USAMM_parameter_files.at(fty_idx);
    n_lines_in_posterior_file = config->USAMM_post_lengths.at(fty_idx);
    sample_post(post_fname);
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
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
    return parstr->c_vec.at(pcl_comb_code) * config->shipment_rate_factor;
}

double USAMMv3_parameters::get_rho(std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->rho_vec.at(0);
}

double USAMMv3_parameters::get_kMu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->kMu_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_kNu(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->kNu_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_oShipSizeScalingCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeOScalingCVI_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_dShipSizeScalingCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeDScalingCVI_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_shipSizeShapeCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeShapeCVI_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_shipSizeInterceptCVI(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeInterceptCVI_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_oShipSizeScalingCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeOScalingCOM_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_dShipSizeScalingCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeDScalingCOM_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_shipSizeShapeCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeShapeCOM_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_shipSizeInterceptCOM(Prem_class* pcl_O, Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    int pcl_comb_code = pairing_function(pcl_O->idx, pcl_D->idx);
//    std::pair<int, int> pcl_pair = std::pair<int, int>(pcl_O->idx, pcl_D->idx);
    return parstr->shipSizeInterceptCOM_vec.at(pcl_comb_code);
}

double USAMMv3_parameters::get_phi_O(Prem_class* pcl_O, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->phiO_vec.at(pcl_O->idx);
}

double USAMMv3_parameters::get_phi_D(Prem_class* pcl_D, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->phiD_vec.at(pcl_D->idx);
}

std::pair<double, double> USAMMv3_parameters::get_abCVI(int state, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return std::make_pair(parstr->aCVI_vec.at(state), parstr->bCVI_vec.at(state));
}

std::pair<double, double> USAMMv3_parameters::get_abCOM(int state, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return std::make_pair(parstr->aCOM_vec.at(state), parstr->bCOM_vec.at(state));
}

double USAMMv3_parameters::get_outflowCVI(int state, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->outflowCVI_vec.at(state);
}

double USAMMv3_parameters::get_inflowCVI(int state, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->inflowCVI_vec.at(state);
}

double USAMMv3_parameters::get_outflowCOM(int state, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->outflowCOM_vec.at(state);
}

double USAMMv3_parameters::get_inflowCOM(int state, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->inflowCOM_vec.at(state);
}

double USAMMv3_parameters::get_ocov(std::string ocov_name, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->ocov_map.at(ocov_name);
}

double USAMMv3_parameters::get_dcov(std::string dcov_name, std::string time_period)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.at(time_period);
    return parstr->dcov_map.at(dcov_name);
}

double USAMMv3_parameters::get_producer_alpha(int producer_id)
{
    USAMMv3_period_parstruct* parstr = parstructs_by_period.begin()->second;
    return parstr->common_pars->alpha_map.at(producer_id);
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
    origin_size_weights_lookup.clear();
    dest_size_weights_lookup.clear();
    for(int period_idx=0; period_idx<int(USAMM_time_period_vec.size()); ++period_idx)
    {
        std::string period = USAMM_time_period_vec[period_idx];
        origin_size_weights[period_idx].clear();
        origin_size_weights[period_idx].resize(prem_classes.size(), std::vector<double>());
        destination_size_weights[period_idx].clear();
        destination_size_weights[period_idx].resize(prem_classes.size(), std::vector<double>());
        for(Prem_class* pcl : prem_classes)
        {
            if(USAMMv3_parameters::size_bins.at(fty).find(pcl) !=
               USAMMv3_parameters::size_bins.at(fty).end())
            {
                std::vector<int> bins(USAMMv3_parameters::size_bins.at(fty).at(pcl).begin(),
                                      USAMMv3_parameters::size_bins.at(fty).at(pcl).end());
                double avg_prem_size = avg_prem_sizes.at(pcl);
                for(size_t i=0; i<bins.size(); ++i)
                {
                    double binned_size = double(bins[i]);
                    double origin_size_weight = std::pow(binned_size / avg_prem_size, get_phi_O(pcl, period));
                    double destination_size_weight = std::pow(binned_size / avg_prem_size, get_phi_D(pcl, period));

                    origin_size_weights[period_idx][pcl->idx].push_back(origin_size_weight);
                    destination_size_weights[period_idx][pcl->idx].push_back(destination_size_weight);

                    origin_size_weights_lookup[period][pcl][bins[i]] = origin_size_weights[period_idx][pcl->idx][i];
                    dest_size_weights_lookup[period][pcl][bins[i]] = destination_size_weights[period_idx][pcl->idx][i];
                }
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

void USAMMv3_parameters::split_pclasses(const std::string& parname, int& pcl_idx_1, int& pcl_idx_2)
{
    std::vector<std::string> pcl_tags = split(parname, '-');
    pcl_idx_1 = pclass_tag_to_idx_map.at(pcl_tags[0]);
    pcl_idx_2 = pclass_tag_to_idx_map.at(pcl_tags[1]);
}

int USAMMv3_parameters::sample_post(std::string fpath)
{
    std::vector<std::string> metadata_vector;
    std::vector<std::string> header_vector;
    std::vector<std::string> data_vector;
    std::ifstream f(fpath);
    std::stringstream generation_data_ss;
    partstr_ptrs.reserve(12); //I don't expect higher resolution than months.
    if(f.is_open())
    {
        parstruct_vec.resize(config->USAMM_temporal_order.size());
        USAMM_time_period_vec.clear();
        USAMM_time_period_vec.resize(config->USAMM_temporal_order.size());
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
            metadata_vector.push_back(temp_line);
        }

        //Extract the species code from USAMM if it hasn't been done. I.e. b, d or s.
        if(USAMMv3_species_code.empty())
        {
            std::string varname = "<shipment type>";
            std::string metadata_first_line = metadata_vector[0];
            size_t found_at = metadata_first_line.find(varname);
            size_t first_pos = metadata_first_line.find("<", found_at + varname.size()); //Find first occurence of '<' after the name of the metadata variable.
            size_t second_pos = metadata_first_line.find(">", found_at + varname.size()); //And the first '>'. Between those is the value we're after.
            std::string val = metadata_first_line.substr(first_pos+1, second_pos - first_pos - 1);
            USAMMv3_species_code = trim(val);

            //If using a posterior from USAMM Swine, check to see if there is also a matching
            //file with COM ids etc.
            if(USAMMv3_species_code == "s")
            {
                size_t ending_pos = fpath.find(".");
                COM_id_fpath = fpath.substr(0, ending_pos) + ".ids";
            }
        }

        //If simulating swine, also get info from the COM id file.
        if(USAMMv3_species_code == "s" and !COM_id_fpath.empty())
        {
            readSwineExclusiveData(COM_id_fpath, last_sampled_line);
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

            if(colname_vector.at(0) != "prior" and colname_vector.at(0) != "ll")
            {
                std::string parname = colname_vector[0];
                std::string period;
                if(colname_vector.size() == 3)
                {
                    period = colname_vector[2];
                }
                else if(colname_vector.size() == 2 and colname_vector[0] == "alpha")
                {
                    //Alpha has no temporal dimension so the same value is shared between all parstructs.
                    period = "Q1";
                }

                if(!temporal_name_exists(period))
                {
                    std::cout << "The temporal identifier \"" << period << "\" in"
                              << config->USAMM_parameter_files.at(fty_idx) << " was not found "
                              << "among the temporal ordering in the config file (option 45). "
                              << "Exiting..." << std::endl;
                    exit(EXIT_FAILURE);
                }

                if(parstructs_by_period.find(period) == parstructs_by_period.end())
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
                    USAMM_time_periods.insert(period);
                    USAMM_time_period_vec[period_idx] = period;
                    parstruct_vec[period_idx] = USAMMv3_period_parstruct(fty, period, period_idx, &common_parstruct);
                    parstructs_by_period[period] = &parstruct_vec[period_idx];
                    partstr_ptrs.push_back(&parstruct_vec[period_idx]);
                }
                USAMMv3_period_parstruct& parstr = *parstructs_by_period.at(period);

                int pcl_idx_1;
                int pcl_idx_2;

                if(parname == "c")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.c_vec[pcl_comb_code] = this_val / config->USAMM_temporal_n_timesteps[parstr.period_idx];
                }
                else if(parname == "kMu")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.kMu_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "kNu" or parname == "kShape")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.kNu_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "rho")
                {
                    parstr.rho_vec[0] = this_val / config->USAMM_temporal_n_timesteps[parstr.period_idx];
                }
                else if(parname == "swineShipSizeOScalingCVI")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeOScalingCVI_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeDScalingCVI")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeDScalingCVI_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeShapeCVI")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeShapeCVI_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeInterceptCVI")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeInterceptCVI_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeOScalingCOM")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeOScalingCOM_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeDScalingCOM")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeDScalingCOM_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeShapeCOM")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeShapeCOM_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "swineShipSizeInterceptCOM")
                {
                    split_pclasses(colname_vector[1], pcl_idx_1, pcl_idx_2);
                    int pcl_comb_code = pairing_function(pcl_idx_1, pcl_idx_2);
                    parstr.shipSizeInterceptCOM_vec[pcl_comb_code] = this_val;
                }
                else if(parname == "phiO")
                {
                    int pcl_idx = pclass_tag_to_idx_map.at(colname_vector[1]);
                    parstr.phiO_vec[pcl_idx] = this_val;
                }
                else if(parname == "phiD")
                {
                    int pcl_idx = pclass_tag_to_idx_map.at(colname_vector[1]);
                    parstr.phiD_vec[pcl_idx] = this_val;
                }
                else if(parname == "dhalf" or parname == "dhalfCVI")
                {
                    parstr.dhalfCVI_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "hratio" or parname == "hratioCVI")
                {
                    parstr.hratioCVI_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "dhalfCOM")
                {
                    parstr.dhalfCOM_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "hratioCOM")
                {
                    parstr.hratioCOM_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "inflow") //With no specification, it is assumed that inflow is common for CVI and COM.
                {
                    parstr.inflowCVI_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                    parstr.inflowCOM_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "inflowCVI")
                {
                    parstr.inflowCVI_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "inflowCOM")
                {
                    parstr.inflowCOM_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "outflow") //With no specification, it is assumed that outflow is common for CVI and COM.
                {
                    parstr.outflowCVI_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                    parstr.outflowCOM_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "outflowCVI")
                {
                    parstr.outflowCVI_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "outflowCOM")
                {
                    parstr.outflowCOM_vec[abbr_to_fips_state_map.at(colname_vector[1])] = this_val;
                }
                else if(parname == "ocov")
                {
                    parstr.ocov_map[colname_vector[1]] = this_val;
                }
                else if(parname == "dcov")
                {
                    parstr.dcov_map[colname_vector[1]] = this_val;
                }
                else if(parname == "alpha")
                {
                    int prod_id = COM_producers_name_to_id.at(colname_vector[1]);
                    parstr.common_pars->alpha_map[prod_id] = this_val;
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
        USAMMv3_period_parstruct* parstr = parstructs_by_period.at(period);
        //First CVI
        for(auto& state_abbr_id_pair : abbr_to_fips_state_map)
        {
            int state_idx = state_abbr_id_pair.second;
            double dhalf = parstr->dhalfCVI_vec.at(state_idx);
            double hratio = parstr->hratioCVI_vec.at(state_idx);
            std::pair<double, double> ab_pair = ab_fun(dhalf, hratio);
            parstr->aCVI_vec.at(state_idx) = ab_pair.first;
            parstr->bCVI_vec.at(state_idx) = ab_pair.second;

            //Then COM if active.
            if(!COM_ids.empty()) //As a way to check if commuter shipments are active.
            {
                double dhalf_COM = parstr->dhalfCOM_vec.at(state_idx);
                double hratio_COM = parstr->hratioCOM_vec.at(state_idx);
                std::pair<double, double> ab_pair_COM = ab_fun(dhalf_COM, hratio_COM);
                parstr->aCOM_vec.at(state_idx) = ab_pair_COM.first;
                parstr->bCOM_vec.at(state_idx) = ab_pair_COM.second;
            }
        }

        //Also clear out any saved kernel values or county-county rates.
        size_t n_c = counties.size();
        size_t n_pcl = prem_classes.size();
        parstr->county_eval_kernel_CVI.clear();
        parstr->county_eval_kernel_CVI.resize(n_c, Vec_d_1d(n_c, -std::numeric_limits<double>::max()));
        parstr->county_eval_kernel_COM.clear();
        parstr->county_eval_kernel_COM.resize(n_c, Vec_d_1d(n_c, -std::numeric_limits<double>::max()));
        parstr->eval_cc_rates_CVI.clear();
        parstr->eval_cc_rates_CVI.resize(n_pcl, Vec_d_3d(n_pcl, Vec_d_2d(n_c, Vec_d_1d(n_c, -std::numeric_limits<double>::max()))));
        parstr->eval_c_o_rates_CVI.clear();
        parstr->eval_c_o_rates_CVI.resize(n_pcl, Vec_d_2d(n_pcl, Vec_d_1d(n_c, -std::numeric_limits<double>::max())));
    }

    //Update the evaluated size bins.
    evaluate_size_bins();
    return 0;
}

int USAMMv3_parameters::readSwineExclusiveData(std::string fpath, size_t posterior_line)
{
    std::ifstream f(fpath, std::ifstream::in);
    if(f.is_open())
    {
        skipBOM(f);
        if(COM_producers_name_to_id.empty())
        {
            std::string line;
            //First line is commuter agreement producer companies.
            std::getline(f, line);
            if(line != "NO_PRODUCERS")
            {
                std::map<int, std::string> temp_map;
                std::vector<std::string> line_vector = split(line, '\t');
                for(std::string& prod_str : line_vector)
                {
                    std::vector<std::string> prod_vec = split(prod_str, '_');
                    COM_producers_name_to_id[prod_vec[1]] = std::stoi(prod_vec[0]);
                }
            }

            //Second is extra unobserved premises added.
            std::getline(f, line);
            if(line != "NO UNOBSERVED PREMISES")
            {
                std::vector<std::string> line_vector = split(line, '\t');
                for(std::string& prem_str : line_vector)
                {
                    std::vector<std::string> prem_vec = split(prem_str, '_');
                    int prem_id = std::stoi(prem_vec[0]);
                    int prem_fips = std::stoi(prem_vec[1]);
                    unobserved_prem_county_association[prem_id] = prem_fips;
                }
            }

            //Read third and save for later when the rest of the file has been processed.
            std::string third_line;
            std::getline(f, third_line);


            //Fourth is the header for the rest of the
            //file and contains all COM_id ids which is useful when reading line three.
            size_t unobs_starts_at = std::numeric_limits<size_t>::max();
            if(third_line != "NO COMMUTER SHIPMENT IDS")
            {
                std::getline(f, line);
                std::vector<int> COM_id_vals;
                if(!line.empty())
                {
                    std::vector<std::string> line_vector = split(line, '\t');
                    COM_ids.reserve(line_vector.size());
                    unobserved_prem_ids.reserve(line_vector.size());
                    for(size_t id_idx=0; id_idx<line_vector.size(); ++id_idx)
                    {
                        std::vector<std::string> id_vector = split(line_vector[id_idx], '_');
                        int id = std::stoi(id_vector[0]);
                        if(id >= 0)
                        {
                            int producer_idx = std::stoi(id_vector[1]);
                            COM_id_vals.push_back(id);
                            COM_ids.emplace_back(id);
                            COM_ids.back().set_idx(id_idx);
                            COM_ids.back().set_producer(producer_idx);
                        }
                        else
                        {
                            //Negative id means it's an unobserved premises.
                            unobserved_prem_ids.push_back(id);
                            if(id_idx < unobs_starts_at)
                            {
                                unobs_starts_at = id_idx;
                            }
                        }
                    }

                    //Since COM_ids might reallocate the contents during building of the vector.
                    COM_ids_by_producer_idx.resize(COM_producers_name_to_id.size());
                    for(COM_prem_id& com_id : COM_ids)
                    {
                        COM_ids_by_id[com_id.get_id()] = &com_id;
//                        com_id.initialize_random_eff_vectors(COM_ids.size(), 4);
                        COM_ids_by_producer_idx[com_id.get_producer()].push_back(&com_id);
                    }
                }

                //Read the actual posterior data at the same line number as was
                //sampled from the posterior (beginning from the first line of
                //DATA, i.e. excluding headers & metadata).
                f.seekg(0, f.beg);
                //Skip ahead 4 lines for the metadata/header + the line number of the posterior sample.
                f.unsetf(std::ios_base::skipws);
                for(size_t i=0; i<4+posterior_line; ++i)
                {
                    f.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); //Discard characters until '\n' or EOF.
                }
                line.clear();
                std::getline(f, line); //Get the line.

                //Extract info from the relevant line.
                std::vector<std::string> line_vector = split(line, '\t');
                for(size_t i=0; i<line_vector.size(); ++i)
                {
                    if(i<unobs_starts_at)
                    {
                        //It's commuter id/prem info
                        int com_prem_id_val = COM_id_vals[i];
                        COM_prem_id* COM_id = COM_ids_by_id.at(com_prem_id_val);
                        COM_id->set_sampled_prem_id(std::stoi(line_vector[i]));
                    }
                    else
                    {
                        //It's unobserved prem size info.
                        int colname = unobserved_prem_ids[i-unobs_starts_at];
                        unobs_prem_sizes[colname] = std::stoi(line_vector[i]);
                    }
                }

                //Third is info about the number of shipments between all COM ids.
                if(!third_line.empty())
                {
                    std::vector<std::string> line_vector = split(third_line, '\t');
                    for(std::string& s : line_vector)
                    {
                        std::vector<std::string> str_vec = split(s, '_'); //oId_dId_quarter_n.
                        int oId = std::stoi(str_vec[0]);
                        int dId = std::stoi(str_vec[1]);
                        int period_idx = std::stoi(str_vec[2]);
                        int n = std::stoi(str_vec[3]);
                        COM_ids_by_id.at(oId)->set_n_shipments_to(COM_ids_by_id.at(dId), period_idx, n);
                    }
                }
            }
        } //End if COM_prodcers.empty()
    }
    else
    {
        std::cout << "Failed to read commuter id file from " << fpath << "." << std::endl;
        exit(EXIT_FAILURE);
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

double USAMMv3_parameters::kernel_fun(double d, double a, double b)
{
    return 1.0 / (1.0 + std::pow(d/a, b));
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
        if(ft->get_species() == "swine")
        {
            //Swine farm, separate binning.
            if(p_original_size < 25)
                { binned_size = 15; }
            else if(p_original_size >= 25 and p_original_size < 50)
                { binned_size = 35; }
            else if(p_original_size >= 50 and p_original_size < 100)
                { binned_size = 75; }
            else if(p_original_size >= 100 and p_original_size < 200)
                { binned_size = 150; }
            else if(p_original_size >= 200 and p_original_size < 500)
                { binned_size = 350; }
            else if(p_original_size >= 500 and p_original_size < 1000) //From the upper limit of this one I had to make something up since NASS is just 1000+.
                { binned_size = 750; }
            else if(p_original_size >= 1000 and p_original_size < 5000)
                { binned_size = 3000; }
            else if(p_original_size >= 5000 and p_original_size < 10000)
                { binned_size = 7500; }
            else if(p_original_size >= 10000)
                { binned_size = 19000; } //This is the average size of swine prems larger than 10000 in flaps.
            else
            {
                std::cout << "Size binning failed with size = " <<
                             p_original_size << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        else
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

///////////////////////////////////////////

USAMMv3_period_parstruct::USAMMv3_period_parstruct(Farm_type* fty, std::string period, int period_idx,
                                                   USAMMv3_common_parstruct* common_pars) :
    fty(fty), period(period), period_idx(period_idx), common_pars(common_pars)
{
    size_t max_n_ptypes = 3;
    size_t max_state_id = 56;

    size_t pty_max_size = pairing_function(max_n_ptypes, max_n_ptypes);
    c_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    kMu_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    kNu_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeOScalingCVI_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeDScalingCVI_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeShapeCVI_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeInterceptCVI_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeOScalingCOM_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeDScalingCOM_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeShapeCOM_vec.resize(pty_max_size, -std::numeric_limits<double>::max());
    shipSizeInterceptCOM_vec.resize(pty_max_size, -std::numeric_limits<double>::max());

    phiO_vec.resize(max_n_ptypes, -std::numeric_limits<double>::max());
    phiD_vec.resize(max_n_ptypes, -std::numeric_limits<double>::max());

    size_t state_max_size = pairing_function(max_state_id, max_state_id);
    dhalfCVI_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    hratioCVI_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    dhalfCOM_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    hratioCOM_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    outflowCVI_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    inflowCVI_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    outflowCOM_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    inflowCOM_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    aCVI_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    bCVI_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    aCOM_vec.resize(state_max_size, -std::numeric_limits<double>::max());
    bCOM_vec.resize(state_max_size, -std::numeric_limits<double>::max());

    //COM related only associated with one farm type.
    rho_vec.resize(1, -std::numeric_limits<double>::max());
    COM_random_effects.resize(1);
}

USAMMv3_common_parstruct::USAMMv3_common_parstruct(Farm_type* fty) :
    fty(fty)
{

}

////////////////////////////////////////////

COM_prem_id::COM_prem_id(int id) :
    id(id)
{
    id_str = std::to_string(id);
    n_shipments_to_others_by_period.resize(4);
}

COM_prem_id::~COM_prem_id()
{
    this->update_farm(nullptr); //Disconnect the currently associated farm.
}

//void COM_prem_id::initialize_random_eff_vectors(size_t n_com_ids, size_t n_periods)
//{
//    random_effects.resize(n_periods, std::vector<double>(n_com_ids, -std::numeric_limits<double>::max()));
//}

void COM_prem_id::update_farm(Farm* in_farm)
{
    //Disconnect the old farm
    if(farm != nullptr)
    {
        farm->setCOMIdentifier(nullptr);
    }
    farm = in_farm;
    if(in_farm != nullptr)
    {
        farm->setCOMIdentifier(this);
    }
}

int COM_prem_id::get_n_shipments_to(COM_prem_id* other, size_t period_idx)
{
    if(n_shipments_to_others_by_period.at(period_idx).find(other) !=
       n_shipments_to_others_by_period.at(period_idx).end())
    {
        return n_shipments_to_others_by_period[period_idx][other];
    }
    return 0; //If it's not in the map, there are no shipment to that COM id.
}

//double COM_prem_id::get_random_effect(size_t other_COM_id_idx, size_t period_idx)
//{
//    return random_effects[period_idx][other_COM_id_idx];
//}
//
//void COM_prem_id::set_random_effect(size_t other_COM_id_idx, size_t period_idx, double val)
//{
//    random_effects[period_idx][other_COM_id_idx] = val;
//}

void COM_prem_id::set_producer(int p)
{
    if(com_producer != -1)
    {
        std::cout << "Attempting to set producer on a COM id that already have a "
                  << "producer assigned to it." << std::endl;
        exit(EXIT_FAILURE);
    }
    com_producer = p;
}

void COM_prem_id::disconnect_farm()
{
    //Disconnect the old farm
    if(farm != nullptr)
    {
        farm->setCOMIdentifier(nullptr);
    }
    farm = nullptr;
}
