/*
 * EkmanUWCurrentModel.cpp
 */

#include "EkmanUWCurrentModel.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include "xdyn/yaml_parser/yaml_compat.h"
#include "xdyn/yaml_parser/parse_unit_value.hpp"
#include "xdyn/yaml_parser/external_data_structures_parsers.hpp"

EkmanUWCurrentModel::EkmanUWCurrentModel(const Input& input) :
    seabed(input.seabed),
    top_ekman_depth(input.top_ekman_depth),
    bottom_ekman_depth(input.bottom_ekman_depth),
    middle_velocity(input.middle_velocity),
    middle_orientation(input.middle_orientation),
    f_times_rho(input.f_times_rho),
    wind_angle(input.wind_angle),
    wind_stress(input.wind_stress),
    rho_air(input.rho_air),
    latitude(input.latitude),
    depth_from_wind(input.depth_from_wind),
    weather_file(input.weather_file),
    weather_from(layer_of(input.U10, input.wind_angle)),
    weather_to(weather_from)
{
    if (depth_from_wind)
    {
        top_ekman_depth = weather_from.top_ekman_depth;
    }
    reload_weather_if_changed();
}

double EkmanUWCurrentModel::layer_depth_of(const double U10, const double latitude_)
{
    return 7.6*U10/std::sqrt(std::sin(std::max(10.*M_PI/180., std::abs(latitude_))));
}

EkmanUWCurrentModel::WindLayer EkmanUWCurrentModel::layer_of(const double U10, const double angle) const
{
    return WindLayer{wind_stress_of(U10, rho_air), angle, depth_from_wind ? layer_depth_of(U10, latitude) : top_ekman_depth};
}

double EkmanUWCurrentModel::wind_stress_of(const double U10, const double rho_air_)
{
    const double drag_coefficient = (U10 < 20.2) ? 0.79e-3 + 0.08e-3*U10 : 0.002423;
    return drag_coefficient * rho_air_ * U10 * U10;
}

void EkmanUWCurrentModel::reload_weather_if_changed() const
{
    if (weather_file.empty()) return;
    const auto now = std::chrono::steady_clock::now();
    if (now < next_weather_check) return;
    next_weather_check = now + std::chrono::milliseconds(200);
    std::error_code ec;
    const auto time = std::filesystem::last_write_time(weather_file, ec);
    if (ec || (has_weather && time == weather_file_time)) return;
    // A file that cannot be read keeps the current weather: the writer renames whole files.
    std::ifstream in(weather_file);
    std::string line;
    std::vector<double> v;
    while (std::getline(in, line))
    {
        if (line.empty() || !(std::isdigit(static_cast<unsigned char>(line[0])) || line[0] == '-' || line[0] == '.')) continue;
        std::stringstream ss(line);
        std::string cell;
        v.clear();
        try
        {
            while (std::getline(ss, cell, ',')) v.push_back(std::stod(cell));
        }
        catch (const std::exception&)
        {
            return;
        }
    }
    if (v.size() != 6) return;
    weather_start = v[0];
    weather_duration = v[1];
    weather_from = layer_of(v[2], v[3]*M_PI/180.);
    weather_to = layer_of(v[4], v[5]*M_PI/180.);
    weather_file_time = time;
    has_weather = true;
}

std::string EkmanUWCurrentModel::model_name() {return "ekman current";}

EkmanUWCurrentModel::Input EkmanUWCurrentModel::parse(const std::string& yaml_input)
{
    YAML::Node node = YAML::Load(yaml_input);
    Input input;
    double rho = 1026.0;        // kg/m3
    try_to_parse(node, "rho", rho);
    double omega = 7.2921e-5;   // rad/s, Earth rotation rate
    try_to_parse(node, "omega", omega);
    double air_density = 1.225; // kg/m3
    try_to_parse(node, "rho_air", air_density);
    if (node["seabed file"])
    {
        std::string seabed_file;
        node["seabed file"] >> seabed_file;
        input.seabed = Seabed(seabed_file);
    }
    else if (node["seabed depth"])
    {
        double seabed_depth;
        node["seabed depth"] >> seabed_depth;
        input.seabed = Seabed(seabed_depth);
    }
    else
    {
        THROW(__PRETTY_FUNCTION__, InvalidInputException, "Please provide either a 'seabed file' or a 'seabed depth' for the Ekman current model.");
    }
    double U10;
    xdyn::yaml_parser::parse_uv(node["wind orientation"], input.wind_angle);
    xdyn::yaml_parser::parse_uv(node["U10"], U10);
    input.rho_air = air_density;
    input.U10 = U10;
    input.wind_stress = wind_stress_of(U10, air_density);
    input.depth_from_wind = false;
    if (node["layer thickness follows wind"])
    {
        node["layer thickness follows wind"] >> input.depth_from_wind;
    }
    if (node["weather file"])
    {
        node["weather file"] >> input.weather_file;
    }
    xdyn::yaml_parser::parse_uv(node["current velocity"], input.middle_velocity);
    xdyn::yaml_parser::parse_uv(node["current orientation"], input.middle_orientation);
    double lat;
    xdyn::yaml_parser::parse_uv(node["latitude"], lat);
    input.latitude = lat;
    xdyn::yaml_parser::parse_uv(node["top layer thickness"], input.top_ekman_depth);
    xdyn::yaml_parser::parse_uv(node["bottom layer thickness"], input.bottom_ekman_depth);
    // Coriolis parameter times water density. Classical Ekman theory gives the
    // surface velocity as V0 = tau/(rho*sqrt(Az*f)) with D_s = pi*sqrt(2*Az/f),
    // hence V0 = sqrt(2)*pi*tau/(rho*f*D_s): rho, not its square root. The
    // sqrt(rho) form is also dimensionally inconsistent -- it resolves to
    // kg^0.5 m^0.5 / s rather than a velocity -- and overestimates V0 by
    // sqrt(rho) = 32x at rho = 1026 kg/m3.
    input.f_times_rho = 2*omega*std::sin(lat)*rho;
    return input;
}

Eigen::Vector3d EkmanUWCurrentModel::get_UWCurrent(const Eigen::Vector3d& position, const double t, const double wave_height) const
{
    reload_weather_if_changed();
    if (!has_weather)
    {
        return current_for(WindLayer{wind_stress, wind_angle, top_ekman_depth}, position, wave_height);
    }
    const double u = (weather_duration <= 0) ? (t >= weather_start ? 1. : 0.)
                                             : std::min(1., std::max(0., (t - weather_start)/weather_duration));
    if (u <= 0) return current_for(weather_from, position, wave_height);
    if (u >= 1) return current_for(weather_to, position, wave_height);
    // Two currents add as velocities: the change moves linearly from one to the other.
    return (1 - u)*current_for(weather_from, position, wave_height) + u*current_for(weather_to, position, wave_height);
}

Eigen::Vector3d EkmanUWCurrentModel::current_for(const WindLayer& layer, const Eigen::Vector3d& position, const double wave_height) const
{
    const double top_depth = layer.top_ekman_depth;
    const double z = position(2);
    const double seabed_height = seabed.get_height(position(0), position(1));
    const Eigen::Vector3d middle_current = get_middle_layer_current();
    // Below two Ekman depths the spiral contributes less than 0.2%, so each layer stops there
    if (z > wave_height && z < (wave_height + 2*top_depth))
    {
        return get_top_layer_current(layer, position, middle_current, wave_height);
    }
    else if (z >= (wave_height + 2*top_depth) && z <= (seabed_height - 2*bottom_ekman_depth))
    {
        return middle_current;
    }
    else if (z > (seabed_height - 2*bottom_ekman_depth) && z < seabed_height)
    {
        return get_bottom_layer_current(position, middle_current, seabed_height);
    }
    return Eigen::Vector3d::Zero();
}

Eigen::Vector3d EkmanUWCurrentModel::get_top_layer_current(const WindLayer& layer,
                                                           const Eigen::Vector3d& position,
                                                           const Eigen::Vector3d& middle_current,
                                                           const double) const
{
    const double z = position(2);
    const double sgn_f = (f_times_rho >= 0) ? 1.0 : -1.0;
    const double V0 = std::sqrt(2)*M_PI*layer.wind_stress/(layer.top_ekman_depth*f_times_rho);
    const double decay_factor = layer.top_ekman_depth/M_PI;
    // The wind orientation must enter both components with the SAME sign. The
    // horizontal plane is isotropic (f and Az are scalars), so rotating the wind
    // rotates the whole solution rigidly and leaves its magnitude unchanged.
    // The previous form used cos(... - wind_angle) against sin(... + wind_angle),
    // which is not a rotation: the perturbation's magnitude then varied with wind
    // direction AND with depth (+28% at the surface to -40% at z/D ~ 1.5 for a
    // 20 deg wind), distorting the vertical shear this model exists to provide.
    const double phase = M_PI/4 - z/decay_factor + sgn_f*layer.wind_angle;
    const double u = middle_current(0) + sgn_f*V0*std::exp(-z/decay_factor)*std::cos(phase);
    const double v = middle_current(1) +       V0*std::exp(-z/decay_factor)*std::sin(phase);
    // Vertical velocity would have to be composed with the orbital velocities, which is not done yet
    return Eigen::Vector3d(u, v, 0.);
}

Eigen::Vector3d EkmanUWCurrentModel::get_middle_layer_current() const
{
    return Eigen::Vector3d(middle_velocity*std::cos(middle_orientation), middle_velocity*std::sin(middle_orientation), 0.);
}

Eigen::Vector3d EkmanUWCurrentModel::get_bottom_layer_current(const Eigen::Vector3d& position,
                                                              const Eigen::Vector3d& middle_current,
                                                              const double seabed_height) const
{
    const double depth_factor = M_PI*(seabed_height - position(2))/bottom_ekman_depth;
    const double u = middle_current(0)*(1 - std::exp(-depth_factor)*std::cos(depth_factor)) - middle_current(1)*std::exp(-depth_factor)*std::sin(depth_factor);
    const double v = middle_current(0)*std::exp(-depth_factor)*std::sin(depth_factor) + middle_current(1)*(1 - std::exp(-depth_factor)*std::cos(depth_factor));
    return Eigen::Vector3d(u, v, 0.);
}
