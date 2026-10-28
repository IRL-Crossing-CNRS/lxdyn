/*
 * DiscreteDirectionalWaveSpectrum.hpp
 *
 *  Created on: Jul 31, 2014
 *      Author: cady
 */

#ifndef DISCRETEDIRECTIONALWAVESPECTRUM_HPP_
#define DISCRETEDIRECTIONALWAVESPECTRUM_HPP_

#include <functional>
#include <string>
#include <vector>

/** \author cady
 *  \date Jul 31, 2014, 1:08:15 PM
 *  \brief Used by 'discretize'
 *  \ingroup wave_models
 *  \section ex1 Example
 *  \snippet environment_models/unit_tests/DirectionalWaveSpectrumTest.cpp DirectionalWaveSpectrumTest example
 *  \section ex2 Expected output
 *  \snippet environment_models/unit_tests/DirectionalWaveSpectrumTest.cpp DirectionalWaveSpectrumTest expected output
 */
struct DiscreteDirectionalWaveSpectrum
{
    DiscreteDirectionalWaveSpectrum();

    std::vector<double> Si;                                     //!< Discretized spectral density (in s m^2/rad)
    std::vector<double> Dj;                                     //!< Spatial spreading (in 1/rad)
    std::vector<double> omega;                                  //!< Angular frequencies the spectrum was discretized at (in rad/s)
    std::vector<double> psi;                                    //!< Directions between 0 & 2pi the spatial spreading was discretized at (in rad)
    std::vector<double> k;                                      //!< Discretized wave number (for each frequency) in rad/m
    std::vector<std::vector<double> > phase;                    //!< Random phases, for each (frequency, direction) couple (but time invariant) in radian phases *phase[i_freq][i_dir]*


    std::function<double(double,double,double)> pdyn_factor;    //!< Factor used when computing the dynamic pressure (no unit)
    std::function<double(double,double,double)> pdyn_factor_sh; //!< Factor used when computing the orbital velocity (no unit)
    void check_sizes() const;
};

/** \brief Used by the wave models (eg. Airy, Stokes, etc.)
 *  \details Avoids having to double loop on frequency and direction. All vectors should have the same size as each element corresponds to a (direction, omega) pair: that way we only need one loop (over, say, i) to go over all possible (omega[i], psi[i]).
 *  \ingroup wave_models
 *  \section ex1 Example
 *  \snippet environment_models/unit_tests/DirectionalWaveSpectrumTest.cpp DirectionalWaveSpectrumTest example
 *  \section ex2 Expected output
 *  \snippet environment_models/unit_tests/DirectionalWaveSpectrumTest.cpp DirectionalWaveSpectrumTest expected output
 */
struct FlatDiscreteDirectionalWaveSpectrum
{
    FlatDiscreteDirectionalWaveSpectrum();
    std::vector<double> a;       //!< Amplitude (in m), for each angular frequency omega, and direction
    std::vector<double> omega;   //!< Angular frequencies the spectrum was discretized at (in rad/s), for each angular frequency omega, and direction
    std::vector<double> psi;     //!< Directions between 0 & 2pi the spatial spreading was discretized at (in rad)
    std::vector<double> cos_psi; //!< Cosinus directions between 0 & 2pi the spatial spreading was discretized at (so we do not compute it each time), for each angular frequency omega, and direction
    std::vector<double> sin_psi; //!< Sinus directions between 0 & 2pi the spatial spreading was discretized at (so we do not compute it each time), for each angular frequency omega, and direction
    std::vector<double> k;       //!< Discretized wave number (for each frequency) (in 1/m), for each angular frequency omega, i.e. same size as omega
    std::vector<double> phase;   //!< Random phases, for each (frequency, direction) couple (but time invariant) in radian, for each angular frequency omega, and direction
    // Optional fade of each ray while the sea state changes, in simulation time (s). Empty:
    // every ray at full amplitude. gain(t) = sqrt(u_in (1 - u_out)), u the clamped ramps.
    std::vector<double> fade_in_start, fade_in_duration, fade_out_start, fade_out_duration;
    std::string rays_file;       //!< CSV re-read when it changes ("reload rays file: true"), empty otherwise
    double gain(const size_t i, const double t) const; //!< Fade of ray i at time t, 1 without fades
    std::vector<double> amplitudes(const double t) const; //!< a times the fade of each ray at time t
    std::function<double(double,double,double)> pdyn_factor;    //!< Factor used when computing the dynamic pressure (no unit)
    std::function<double(double,double,double)> pdyn_factor_sh; //!< Factor used when computing the orbital velocity (no unit)
    std::vector<double> get_periods() const; //< Get the ray periods as a vector, from omega attribute (in s)
};

/** \brief Reads the rays of a CSV with columns a,psi,omega,k,phase and, optionally,
 *  fade_in_start,fade_in_duration,fade_out_start,fade_out_duration (any order, header line first).
 *  \returns false, with the reason in error, if the file cannot be read: f is then unchanged.
 */
bool read_rays_csv(const std::string& path, FlatDiscreteDirectionalWaveSpectrum& f, std::string& error);

#endif /* DISCRETEDIRECTIONALWAVESPECTRUM_HPP_ */
