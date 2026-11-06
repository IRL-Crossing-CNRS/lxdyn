/*
 * EkmanUWCurrentModel.hpp
 *
 * Three-layer current: a wind-driven Ekman spiral at the top, a uniform background current in the
 * middle, and a seabed-influenced layer at the bottom.
 */

#ifndef ENVIRONMENT_MODELS_INC_EKMANUWCURRENTMODEL_HPP_
#define ENVIRONMENT_MODELS_INC_EKMANUWCURRENTMODEL_HPP_

#include <chrono>
#include <filesystem>
#include <string>

#include "xdyn/environment_models/UWCurrentModel.hpp"
#include "xdyn/environment_models/Seabed.hpp"

class EkmanUWCurrentModel : public UWCurrentModel
{
public:
    struct Input
    {
        Input();
        Input(Seabed);
        Seabed seabed;
        double top_ekman_depth;     ///< Top layer thickness (m)
        double bottom_ekman_depth;  ///< Bottom layer thickness (m)
        double middle_velocity;     ///< Background current velocity (m/s)
        double middle_orientation;  ///< Background current direction (rad)
        double f_times_rho;      ///< Coriolis parameter times water density
        double wind_angle;          ///< Wind direction (rad)
        double wind_stress;         ///< Wind stress (N/m2)
        double rho_air;             ///< Air density (kg/m3)
        double latitude;            ///< rad
        double U10;                 ///< 10 m wind (m/s)
        bool depth_from_wind;       ///< "layer thickness follows wind": D = 7.6 U10 / sqrt(sin|latitude|)
        std::string weather_file;   ///< Optional: a change of wind during the run (see below)
    };
    EkmanUWCurrentModel(const Input& input);
    virtual ~EkmanUWCurrentModel();

    virtual Eigen::Vector3d get_UWCurrent(const Eigen::Vector3d& position, const double time, const double wave_height) const override;

    static std::string model_name();
    static Input parse(const std::string&);

    /// Wind stress (N/m2) of the 10 m wind U10 (m/s): C_D rho_air U10^2.
    static double wind_stress_of(const double U10, const double rho_air);

    /// Depth (m) of the wind-driven layer that a 10 m wind U10 (m/s) sets at a latitude (rad):
    /// Ekman's empirical D = 7.6 U10 / sqrt(sin|latitude|) (Stewart, Introduction to Physical
    /// Oceanography, ch. 9), the latitude floored at 10 deg. With a fixed depth the surface current
    /// would grow as U10^2.
    static double layer_depth_of(const double U10, const double latitude);

    /// The wind-driven top layer: the part of the model a change of weather moves.
    struct WindLayer
    {
        double wind_stress;      ///< N/m2
        double wind_angle;       ///< rad, direction the wind blows towards
        double top_ekman_depth;  ///< m
    };

private:
    Seabed seabed;
    double top_ekman_depth;
    double bottom_ekman_depth;
    double middle_velocity;
    double middle_orientation;
    double f_times_rho;
    double wind_angle;
    double wind_stress;
    double rho_air;
    double latitude;
    bool depth_from_wind;
    WindLayer layer_of(const double U10, const double angle) const;

    // "weather file": a CSV with one row, start_s, duration_s, U10_from, wind_from_deg, U10_to,
    // wind_to_deg (U10 in m/s, wind towards, from north to east). Between the two simulation times the current is
    // (1 - u) current(from) + u current(to), u going linearly from 0 to 1; before, the "from" wind;
    // after, the "to" wind. The file is re-read when it changes (checked at most every 200 ms), so a
    // live weather can rewrite it during the run, as the rays file of the waves.
    std::string weather_file;
    mutable WindLayer weather_from;
    mutable WindLayer weather_to;
    mutable double weather_start;
    mutable double weather_duration;
    mutable bool has_weather;
    mutable std::chrono::steady_clock::time_point next_weather_check;
    mutable std::filesystem::file_time_type weather_file_time;
    void reload_weather_if_changed() const;

    Eigen::Vector3d current_for(const WindLayer& layer, const Eigen::Vector3d& position, const double wave_height) const;
    Eigen::Vector3d get_top_layer_current(const WindLayer&, const Eigen::Vector3d&, const Eigen::Vector3d&, const double) const;
    Eigen::Vector3d get_middle_layer_current() const;
    Eigen::Vector3d get_bottom_layer_current(const Eigen::Vector3d&, const Eigen::Vector3d&, const double) const;
};

#endif /* ENVIRONMENT_MODELS_INC_EKMANUWCURRENTMODEL_HPP_ */
