/*
 * DiscreteDirectionalWaveSpectrum.cpp
 *
 *  Created on: Jul 31, 2014
 *      Author: cady
 */

#include "DiscreteDirectionalWaveSpectrum.hpp"
#include "xdyn/exceptions/InvalidInputException.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>

#include <cmath>
#define PI M_PI

std::vector<double> FlatDiscreteDirectionalWaveSpectrum::get_periods() const
{
    std::vector<double> periods;
    periods.reserve(omega.size());
    for (size_t i=0;i<omega.size();++i)
    {
        periods.push_back(2.0 * PI / omega.at(i));
    }
    return periods;
}

void DiscreteDirectionalWaveSpectrum::check_sizes() const
{
    const size_t n_si = Si.size();
    const size_t n_dj = Dj.size();
    const size_t n_omega = omega.size();
    const size_t n_psi = psi.size();
    const size_t n_k = k.size();
    const size_t n_phase = phase.size();
    if (n_si==0)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "No data for Si");
    }
    if (n_omega==0)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "No data for omega");
    }
    if (n_dj==0)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "No data for Dj");
    }
    if (n_psi==0)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "No data for psi");
    }
    if (n_k==0)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "No data for k");
    }
    if (n_phase==0)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "No data for phase");
    }
    if (n_si!=n_omega)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException,
            "Si and omega should have the same size"+std::to_string(n_si)+"!="+std::to_string(n_omega));
    }
    if (n_si!=n_k)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException,
            "Si and k should have the same size"+std::to_string(n_si)+"!="+std::to_string(n_k));
    }
    if (n_dj!=n_psi)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException,
            "Dj and psi should have the same size"+std::to_string(n_dj)+"!="+std::to_string(n_psi));
    }
    if (n_si!=n_phase)
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException,
            "Si and phase should have the same size"+std::to_string(n_si)+"!="+std::to_string(n_phase));
    }
    for (size_t i_omega=0;i_omega<n_omega;++i_omega)
    {
        if (phase.at(i_omega).empty())
        {
            THROW(__PRETTY_FUNCTION__, InvalidInputException, "phase["+std::to_string(i_omega)+"] is empty");
        }
        if (phase.at(i_omega).size()!=n_dj)
        {
            THROW(__PRETTY_FUNCTION__, InvalidInputException,
                "phase["+std::to_string(i_omega)+"] and Dj should have the same size "+
                std::to_string(phase.at(i_omega).size())+"!="+std::to_string(n_dj));
        }
    }
}

double FlatDiscreteDirectionalWaveSpectrum::gain(const size_t i, const double t) const
{
    if (i >= fade_in_start.size()) return 1;
    const auto ramp = [t](const double start, const double duration)
    {
        if (t < start) return 0.;
        if (duration <= 0) return 1.;
        return std::min(1., (t - start) / duration);
    };
    const double u_in = ramp(fade_in_start[i], fade_in_duration[i]);
    const double u_out = ramp(fade_out_start[i], fade_out_duration[i]);
    return std::sqrt(u_in * (1 - u_out));
}

std::vector<double> FlatDiscreteDirectionalWaveSpectrum::amplitudes(const double t) const
{
    if (fade_in_start.empty()) return a;
    std::vector<double> ret(a.size());
    for (size_t i = 0 ; i < a.size() ; ++i) ret[i] = a[i] * gain(i, t);
    return ret;
}

bool read_rays_csv(const std::string& path, FlatDiscreteDirectionalWaveSpectrum& f, std::string& error)
{
    std::ifstream in(path);
    if (!in)
    {
        error = "cannot open '" + path + "'";
        return false;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        error = "'" + path + "' is empty";
        return false;
    }
    std::vector<std::string> names;
    {
        std::stringstream ss(line);
        std::string name;
        while (std::getline(ss, name, ','))
        {
            name.erase(std::remove_if(name.begin(), name.end(), [](const char c){return c == ' ' || c == '\r';}), name.end());
            names.push_back(name);
        }
    }
    std::map<std::string, std::vector<double> > columns;
    for (const auto& name: names) columns[name];
    size_t row = 1;
    while (std::getline(in, line))
    {
        ++row;
        if (line.empty() || line == "\r") continue;
        std::stringstream ss(line);
        std::string cell;
        size_t c = 0;
        for (; std::getline(ss, cell, ',') && c < names.size(); ++c)
        {
            try
            {
                columns[names[c]].push_back(std::stod(cell));
            }
            catch (const std::exception&)
            {
                error = "'";
                error += path;
                error += "', line ";
                error += std::to_string(row);
                error += ": '";
                error += cell;
                error += "' is not a number";
                return false;
            }
        }
        if (c != names.size())
        {
            error = "'" + path + "', line " + std::to_string(row) + ": " + std::to_string(c) + " values for " + std::to_string(names.size()) + " columns";
            return false;
        }
    }
    for (const char* required: {"a", "psi", "omega", "k", "phase"})
    {
        if (columns.find(required) == columns.end())
        {
            error = "'" + path + "' has no column '" + required + "'";
            return false;
        }
    }
    const char* fades[] = {"fade_in_start", "fade_in_duration", "fade_out_start", "fade_out_duration"};
    size_t nb_fades = 0;
    for (const char* name: fades) nb_fades += columns.count(name);
    if (nb_fades != 0 && nb_fades != 4)
    {
        error = "'" + path + "' needs all four fade columns or none";
        return false;
    }
    f.a = columns["a"];
    f.psi = columns["psi"];
    f.omega = columns["omega"];
    f.k = columns["k"];
    f.phase = columns["phase"];
    f.cos_psi.clear();
    f.sin_psi.clear();
    for (const auto& psi: f.psi)
    {
        f.cos_psi.push_back(std::cos(psi));
        f.sin_psi.push_back(std::sin(psi));
    }
    if (nb_fades == 4)
    {
        f.fade_in_start = columns["fade_in_start"];
        f.fade_in_duration = columns["fade_in_duration"];
        f.fade_out_start = columns["fade_out_start"];
        f.fade_out_duration = columns["fade_out_duration"];
    }
    else
    {
        f.fade_in_start.clear(); f.fade_in_duration.clear(); f.fade_out_start.clear(); f.fade_out_duration.clear();
    }
    return true;
}
