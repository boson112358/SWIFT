/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2022 Doug Rennehan (douglas.rennehan@gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_RENNEHAN_BLACK_HOLES_H
#define SWIFT_RENNEHAN_BLACK_HOLES_H

/* Local includes */
#include "black_holes_properties.h"
#include "black_holes_struct.h"
#include "cooling.h"
#include "cosmology.h"
#include "dimension.h"
#include "exp10.h"
#include "gravity.h"
#include "kernel_hydro.h"
#include "minmax.h"
#include "physical_constants.h"
#include "random.h"

/* Standard includes */
#include <float.h>
#include <math.h>
#include <gsl/gsl_poly.h>


/**
 * @brief How much of the feedback actually couples to the medium?
 *
 * @param props The properties of the black hole scheme.
 * @param BH_state The current state of the black hole.
 */
__attribute__((always_inline)) INLINE static float get_black_hole_coupling(
    const struct black_holes_props* props, const int BH_state) {
  switch (BH_state) {
    case BH_states_adaf:
        return props->adaf_coupling;
        break;
    case BH_states_quasar:
        return props->quasar_coupling;
        break;
    case BH_states_slim_disk:
        return props->slim_disk_coupling;
        break;
    default:
        error("Invalid black hole state.");
        return 0.f;
        break;
  }
}

/**
 * @brief Computes the radiative efficiency in the slim disk mode.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd M_dot,BH / M_dot,Edd
 */
__attribute__((always_inline)) INLINE static 
float get_black_hole_slim_disk_efficiency(
    const struct black_holes_props* props, const double f_Edd) {
  if (f_Edd <= 0.f) return 0.f;
  const float R = 1.f / f_Edd;
  /* Efficiency from Lupi et al. (2014), 
   * super eddington accretion and feedback */
  return (R / 16.f) * props->A_lupi * 
         (0.985f / (R + (5.f / 8.f) * props->B_lupi) + 0.015f / 
            (R + (5.f / 8.f) * props->C_lupi));
}

/**
 * @brief Computes the radiative efficiency in the ADAF mode.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd M_dot,BH / M_dot,Edd
 */
__attribute__((always_inline)) INLINE static 
float get_black_hole_adaf_efficiency(
    const struct black_holes_props* props, const double f_Edd) {
  return props->epsilon_r * f_Edd;  /* scales with M_dot,BH */
}

/**
 * @brief Chooses and calls the proper radiative efficiency function for the state.
 *
 * @param props The properties of the black hole scheme.
 * @param f_Edd The accretion rate over the Eddington rate.
 * @param BH_state The current state of the BH.
 */
__attribute__((always_inline)) INLINE static 
float get_black_hole_radiative_efficiency(
    const struct black_holes_props* props, 
    const double f_Edd, const int BH_state) {
  switch(BH_state) {
    case BH_states_adaf:
        return get_black_hole_adaf_efficiency(props, f_Edd);
    case BH_states_quasar:
        return props->epsilon_r;
    case BH_states_slim_disk:
        return get_black_hole_slim_disk_efficiency(props, f_Edd);
    default:
        error("Invalid black hole state.");
        break;
  }

  return 0.f;
}

/**
 * @brief Computes the fraction of M_dot,inflow that should go into the BH.
 *
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param m_dot_inflow_m_dot_edd M_dot,inflow scaled to M_dot,Edd for the BH.
 */
__attribute__((always_inline)) INLINE static float get_black_hole_upper_mdot_medd(
    const struct black_holes_props* props, 
    const struct phys_const* constants,
    const float m_dot_inflow_m_dot_edd) {
  if (m_dot_inflow_m_dot_edd <= 0.f) return 0.f;

  double x1, x2, x3;
  double a3, a2, a1, a0;
  const double phi = props->slim_disk_wind_mass_loading;
                      
  /* Kinetic constrained models don't work at these resolutions 
    * const double phi = 2.0 * props->slim_disk_coupling * 
    * pow(c / props->slim_disk_wind_speed, 2);
    */

  int num_roots;

  a3 = ((5. * 5.) / (8. * 8.)) * props->B_lupi * props->C_lupi;
  a2 = (5. / 8.) * ((props->B_lupi + props->C_lupi) + (phi / 16.) * 
       props->A_lupi * (0.015 * props->B_lupi + 0.985 * props->C_lupi) - 
       (5. / 8.) * props->B_lupi * props->C_lupi * m_dot_inflow_m_dot_edd);
  a1 = 1. + (phi / 16.) * props->A_lupi - (5. / 8.) * 
       (props->B_lupi + props->C_lupi) * m_dot_inflow_m_dot_edd;
  a0 = -m_dot_inflow_m_dot_edd;

  a2 /= a3;
  a1 /= a3;
  a0 /= a3;

  num_roots = gsl_poly_solve_cubic(a2, a1, a0, &x1, &x2, &x3);
  if (num_roots == 1) {
    if (x1 >= 0.) {
      return x1;
    } else {
      error("num_roots=1 m_dot_inflow_m_dot_edd=%g phi=%g a3=%g a2=%g "
            "a1=%g a0=%g",
            m_dot_inflow_m_dot_edd, phi, a3, a2, a1, a0);
    }
  }
  if (x3 >= 0.) {
    return x3;
  } else {
    error("m_dot_inflow_m_dot_edd=%g phi=%g a3=%g a2=%g a1=%g a0=%g",
          m_dot_inflow_m_dot_edd, phi, a3, a2, a1, a0);
  }

  return 0.f;
}

/**
 * @brief Computes the fraction of M_dot,inflow that should go into the BH.
 *
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param m_dot_inflow M_dot,inflow in internal units.
 * @param BH_mass The subgrid mass of the BH in internal units.
 * @param BH_state The current state of the BH.
 * @param Eddington_rate M_dot,Edd in internal units.
 */
__attribute__((always_inline)) INLINE static 
float get_black_hole_accretion_factor(
    const struct black_holes_props* props, 
    const struct phys_const* constants,
    const float m_dot_inflow, const float BH_mass, 
    const int BH_state, 
    const float Eddington_rate) {
  
  if (m_dot_inflow <= 0.f || BH_mass <= 0.f) return 0.f;

  switch (BH_state) {
    case BH_states_adaf:
      return props->adaf_f_accretion;
    case BH_states_quasar:
      return props->quasar_f_accretion;
    case BH_states_slim_disk:
      /* This is the FRACTION of the total so divide by M_dot,inflow */
      return get_black_hole_upper_mdot_medd(
               props, constants, m_dot_inflow / Eddington_rate
             ) * 
             Eddington_rate / m_dot_inflow;
    default:
      error("Invalid black hole state.");
      return 0.f;
      break;
  }
}

/**
 * @brief Computes the time-step of a given black hole particle.
 *
 * @param bp Pointer to the s-particle data.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 */
__attribute__((always_inline)) INLINE static float black_holes_compute_timestep(
    const struct bpart* const bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo) {

  /* Allow for finer timestepping if necessary! */
  float dt_accr = FLT_MAX;
  if (bp->accretion_rate > 0.f) {
    dt_accr = props->dt_accretion_factor * bp->mass / bp->accretion_rate;
  }

  if (dt_accr < props->time_step_min) {
    message(
        "Warning! BH_TIMESTEP_LOW: id=%lld (%g Myr) is below time_step_min (%g "
        "Myr).",
        bp->id, dt_accr * props->time_to_Myr,
        props->time_step_min * props->time_to_Myr);
  }

  return max(dt_accr, props->time_step_min);
}

/**
 * @brief Initialises the b-particles for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions.
 *
 * @param bp The particle to act upon
 * @param props The properties of the black holes model.
 */
__attribute__((always_inline)) INLINE static void black_holes_first_init_bpart(
    struct bpart* bp, const struct black_holes_props* props) {

  bp->time_bin = 0;
  if (props->use_subgrid_mass_from_ics == 0) {
    bp->subgrid_mass = bp->mass;
  } else if (props->with_subgrid_mass_check && bp->subgrid_mass <= 0) {
    error(
        "Black hole %lld has a subgrid mass of %f (internal units).\n"
        "If this is because the ICs do not contain a 'SubgridMass' data "
        "set, you should set the parameter "
        "'RennehanAGN:use_subgrid_mass_from_ics' to 0 to initialize the "
        "black hole subgrid masses to the corresponding dynamical masses.\n"
        "If the subgrid mass is intentionally set to this value, you can "
        "disable this error by setting 'RennehanAGN:with_subgrid_mass_check' "
        "to 0.",
        bp->id, bp->subgrid_mass);
  }
  bp->total_accreted_mass = 0.f;
  bp->accretion_rate = 0.f;
  bp->formation_time = -1.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->number_of_time_steps = 0;
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;
  bp->accreted_angular_momentum[0] = 0.f;
  bp->accreted_angular_momentum[1] = 0.f;
  bp->accreted_angular_momentum[2] = 0.f;
  bp->last_repos_vel = 0.f;
  bp->dt_accr = FLT_MAX;
  bp->radiative_luminosity = 0.f;
  bp->delta_energy_this_timestep = 0.f;
  bp->AGN_number_of_AGN_events = 0;
  bp->AGN_number_of_energy_injections = 0;
  bp->state = BH_states_slim_disk;
  bp->radiative_efficiency = 0.f;
  bp->f_accretion = 0.f;
  bp->m_dot_inflow = 0.f;
  bp->jet_mass_reservoir = 0.f;
  bp->adaf_energy_to_dump = 0.f;
  bp->dm_mass = 0.f;
  bp->dm_mass_low_vel = 0.f;
  bp->relative_velocity_to_dm_com2 = 0.f;
  bp->dm_com_velocity[0] = 0.f;
  bp->dm_com_velocity[1] = 0.f;
  bp->dm_com_velocity[2] = 0.f;
  bp->relative_velocity_to_dm_com[0] = 0.f;
  bp->relative_velocity_to_dm_com[1] = 0.f;
  bp->relative_velocity_to_dm_com[2] = 0.f;

#ifdef WITH_FOF_GALAXIES
  bp->group_data.mass = 0.f;
  bp->group_data.stellar_mass = 0.f;
  bp->group_data.ssfr = 0.f;
#endif
}

/**
 * @brief Prepares a b-particle for its interactions
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_init_bpart(
    struct bpart* bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_density[i] = -1;
  bp->num_ngb_density = 0;
#endif

  bp->density.wcount = 0.f;
  bp->density.wcount_dh = 0.f;
  bp->rho_gas = 0.f;
  bp->sound_speed_gas = 0.f;
  bp->internal_energy_gas = 0.f;
  bp->hot_gas_mass = 0.f;
  bp->cold_gas_mass = 0.f;
  bp->hot_gas_internal_energy = 0.f;
  bp->rho_subgrid_gas = -1.f;
  bp->sound_speed_subgrid_gas = -1.f;
  bp->velocity_gas[0] = 0.f;
  bp->velocity_gas[1] = 0.f;
  bp->velocity_gas[2] = 0.f;
  bp->circular_velocity_gas[0] = 0.f;
  bp->circular_velocity_gas[1] = 0.f;
  bp->circular_velocity_gas[2] = 0.f;
  bp->angular_momentum_gas[0] = 0.f;
  bp->angular_momentum_gas[1] = 0.f;
  bp->angular_momentum_gas[2] = 0.f;
  bp->specific_angular_momentum_stars[0] = 0.f;
  bp->specific_angular_momentum_stars[1] = 0.f;
  bp->specific_angular_momentum_stars[2] = 0.f;
  bp->stellar_mass = 0.f;
  bp->stellar_bulge_mass = 0.f;
  bp->radiative_luminosity = 0.f;
  bp->ngb_mass = 0.f;
  bp->num_ngbs = 0;
  bp->reposition.delta_x[0] = -FLT_MAX;
  bp->reposition.delta_x[1] = -FLT_MAX;
  bp->reposition.delta_x[2] = -FLT_MAX;
  bp->reposition.min_potential = FLT_MAX;
  bp->reposition.potential = FLT_MAX;
  bp->accretion_rate = 0.f; /* Optionally accumulated ngb-by-ngb */
  bp->cold_disk_mass = 0.f;
  bp->mass_at_start_of_step = bp->mass; /* bp->mass may grow in nibbling mode */
  bp->m_dot_inflow = 0.f; /* reset accretion rate */
  bp->adaf_energy_to_dump = 0.f;
  bp->dm_mass = 0.f;
  bp->dm_mass_low_vel = 0.f;
  bp->relative_velocity_to_dm_com2 = 0.f;
  bp->dm_com_velocity[0] = 0.f;
  bp->dm_com_velocity[1] = 0.f;
  bp->dm_com_velocity[2] = 0.f;
  bp->relative_velocity_to_dm_com[0] = 0.f;
  bp->relative_velocity_to_dm_com[1] = 0.f;
  bp->relative_velocity_to_dm_com[2] = 0.f;
}

/**
 * @brief Predict additional particle fields forward in time when drifting
 *
 * The fields do not get predicted but we move the BH to its new position
 * if a new one was calculated in the repositioning loop.
 *
 * @param bp The particle
 * @param dt_drift The drift time-step for positions.
 */
__attribute__((always_inline)) INLINE static void black_holes_predict_extra(
    struct bpart* restrict bp, float dt_drift) {

  /* Are we doing some repositioning? */
  if (bp->reposition.min_potential != FLT_MAX) {

#ifdef SWIFT_DEBUG_CHECKS
    if (bp->reposition.delta_x[0] == -FLT_MAX ||
        bp->reposition.delta_x[1] == -FLT_MAX ||
        bp->reposition.delta_x[2] == -FLT_MAX) {
      error("Something went wrong with the new repositioning position");
    }

    const double dx = bp->reposition.delta_x[0];
    const double dy = bp->reposition.delta_x[1];
    const double dz = bp->reposition.delta_x[2];
    const double d = sqrt(dx * dx + dy * dy + dz * dz);
    if (d > 1.01 * kernel_gamma * bp->h)
      error("Repositioning BH beyond the kernel support!");
#endif

    /* Move the black hole */
    bp->x[0] += bp->reposition.delta_x[0];
    bp->x[1] += bp->reposition.delta_x[1];
    bp->x[2] += bp->reposition.delta_x[2];

    /* Move its gravity properties as well */
    bp->gpart->x[0] += bp->reposition.delta_x[0];
    bp->gpart->x[1] += bp->reposition.delta_x[1];
    bp->gpart->x[2] += bp->reposition.delta_x[2];

    /* Store the delta position */
    bp->x_diff[0] -= bp->reposition.delta_x[0];
    bp->x_diff[1] -= bp->reposition.delta_x[1];
    bp->x_diff[2] -= bp->reposition.delta_x[2];

    /* Reset the reposition variables */
    bp->reposition.delta_x[0] = -FLT_MAX;
    bp->reposition.delta_x[1] = -FLT_MAX;
    bp->reposition.delta_x[2] = -FLT_MAX;
    bp->reposition.min_potential = FLT_MAX;

    /* Count the jump */
    bp->number_of_repositions++;
  }
}

/**
 * @brief Sets the values to be predicted in the drifts to their values at a
 * kick time
 *
 * @param bp The particle.
 */
__attribute__((always_inline)) INLINE static void
black_holes_reset_predicted_values(struct bpart* bp) {}

/**
 * @brief Kick the additional variables
 *
 * @param bp The particle to act upon
 * @param dt The time-step for this kick
 */
__attribute__((always_inline)) INLINE static void black_holes_kick_extra(
    struct bpart* bp, float dt) {}

/**
 * @brief Finishes the calculation of density on black holes
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_end_density(
    struct bpart* bp, const struct cosmology* cosmo) {

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                       /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv);       /* 1/h^d */
  const float h_inv_dim_plus_one = h_inv_dim * h_inv; /* 1/h^(d+1) */

  /* --- Finish the calculation by inserting the missing h factors --- */
  bp->density.wcount *= h_inv_dim;
  bp->density.wcount_dh *= h_inv_dim_plus_one;
  bp->rho_gas *= h_inv_dim;
  float rho_inv = 1.f;
  if (bp->rho_gas > 0.f) rho_inv = 1.f / bp->rho_gas;
  /* All mass-weighted quantities are for the hot & cold gas */
  float m_hot_inv = 1.f;
  if (bp->hot_gas_mass > 0.f) m_hot_inv /= bp->hot_gas_mass;

  /* For the following, we also have to undo the mass smoothing
   * (N.B.: bp->velocity_gas is in BH frame, in internal units). */
  bp->sound_speed_gas *= h_inv_dim * rho_inv;
  bp->internal_energy_gas *= h_inv_dim * rho_inv;
  bp->hot_gas_internal_energy *= m_hot_inv;
  bp->velocity_gas[0] *= h_inv_dim * rho_inv; /* h_inv_dim * rho_inv;*/
  bp->velocity_gas[1] *= h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */
  bp->velocity_gas[2] *= h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */
  bp->circular_velocity_gas[0] *= h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */
  bp->circular_velocity_gas[1] *= h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */
  bp->circular_velocity_gas[2] *= h_inv_dim * rho_inv; /* h_inv_dim * rho_inv; */

  /* Calculate circular velocity at the smoothing radius from specific
   * angular momentum (extra h_inv) */
  bp->circular_velocity_gas[0] *= h_inv;
  bp->circular_velocity_gas[1] *= h_inv;
  bp->circular_velocity_gas[2] *= h_inv;
}

/**
 * @brief Sets all particle fields to sensible values when the #spart has 0
 * ngbs.
 *
 * @param bp The particle to act upon
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void
black_holes_bpart_has_no_neighbours(struct bpart* bp,
                                    const struct cosmology* cosmo) {

  warning(
      "BH particle with ID %lld treated as having no neighbours (h: %g, "
      "wcount: %g).",
      bp->id, bp->h, bp->density.wcount);

  /* Some smoothing length multiples. */
  const float h = bp->h;
  const float h_inv = 1.0f / h;                 /* 1/h */
  const float h_inv_dim = pow_dimension(h_inv); /* 1/h^d */

  /* Re-set problematic values */
  bp->density.wcount = kernel_root * h_inv_dim;
  bp->density.wcount_dh = 0.f;

  bp->velocity_gas[0] = FLT_MAX;
  bp->velocity_gas[1] = FLT_MAX;
  bp->velocity_gas[2] = FLT_MAX;

  bp->internal_energy_gas = -FLT_MAX;
  bp->hot_gas_internal_energy = -FLT_MAX;
}

/**
 * @brief Return the current instantaneous accretion rate of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accretion_rate(const struct bpart* bp) {
  return bp->accretion_rate;
}

/**
 * @brief Return the total accreted gas mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_accreted_mass(const struct bpart* bp) {
  return bp->total_accreted_mass;
}

/**
 * @brief Return the subgrid mass of this BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_subgrid_mass(const struct bpart* bp) {
  return bp->subgrid_mass;
}

/**
 * @brief Return the current bolometric luminosity of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_bolometric_luminosity(const struct bpart* bp,
                                      const struct phys_const* constants) {
  const double c = constants->const_speed_light_c;
  return bp->accretion_rate * bp->radiative_efficiency * c * c;
}

/**
 * @brief Return the current kinetic jet power of the BH.
 *
 * @param bp the #bpart.
 */
__attribute__((always_inline)) INLINE static double black_holes_get_jet_power(
    const struct bpart* bp, const struct phys_const* constants,
    const struct black_holes_props* props) {
  const double c = constants->const_speed_light_c;
  const double eta_jet = props->jet_efficiency;
  const double psi_jet = props->jet_mass_loading;

  /* accretion_rate is M_dot,acc from the paper */
  return ((eta_jet * eta_jet) / psi_jet) * bp->accretion_rate * c * c;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a gas particle.
 *
 * @param bp The #bpart to update.
 * @param p The #part that is swallowed.
 * @param xp The #xpart that is swallowed.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_part(
    struct bpart* bp, const struct part* p, const struct xpart* xp,
    const struct cosmology* cosmo) {

  /* Get the current dynamical masses */
  const float gas_mass = hydro_get_mass(p);
  const float BH_mass = bp->mass;

  /* Increase the dynamical mass of the BH. */
  bp->mass += gas_mass;
  bp->gpart->mass += gas_mass;

  /* Physical velocity difference between the particles */
  const float dv[3] = {(bp->v[0] - p->v[0]) * cosmo->a_inv,
                       (bp->v[1] - p->v[1]) * cosmo->a_inv,
                       (bp->v[2] - p->v[2]) * cosmo->a_inv};

  /* Physical distance between the particles */
  const float dx[3] = {(bp->x[0] - p->x[0]) * cosmo->a,
                       (bp->x[1] - p->x[1]) * cosmo->a,
                       (bp->x[2] - p->x[2]) * cosmo->a};

  /* Collect the swallowed angular momentum */
  bp->swallowed_angular_momentum[0] +=
      gas_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
  bp->swallowed_angular_momentum[1] +=
      gas_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
  bp->swallowed_angular_momentum[2] +=
      gas_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

  /* Update the BH momentum */
  const float BH_mom[3] = {BH_mass * bp->v[0] + gas_mass * p->v[0],
                           BH_mass * bp->v[1] + gas_mass * p->v[1],
                           BH_mass * bp->v[2] + gas_mass * p->v[2]};

  bp->v[0] = BH_mom[0] / bp->mass;
  bp->v[1] = BH_mom[1] / bp->mass;
  bp->v[2] = BH_mom[2] / bp->mass;
  bp->gpart->v_full[0] = bp->v[0];
  bp->gpart->v_full[1] = bp->v[1];
  bp->gpart->v_full[2] = bp->v[2];

  const float dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
  message(
      "BH %lld swallowing gas particle %lld "
      "(Delta_v = [%f, %f, %f] U_V, "
      "Delta_x = [%f, %f, %f] U_L, "
      "Delta_v_rad = %f)",
      bp->id, p->id, -dv[0], -dv[1], -dv[2], -dx[0], -dx[1], -dx[2],
      (dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2]) / dr);

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_add_part_to_bpart(bp_chem, p_chem, gas_mass);

  /* This BH swallowed a gas particle */
  bp->number_of_gas_swallows++;
  bp->number_of_direct_gas_swallows++;

  /* This BH lost a neighbour */
  bp->num_ngbs--;
  bp->ngb_mass -= gas_mass;
}

/**
 * @brief Update the properties of a black hole particles by swallowing
 * a BH particle.
 *
 * @param bpi The #bpart to update.
 * @param bpj The #bpart that is swallowed.
 * @param cosmo The current cosmological model.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param props The properties of the black hole scheme.
 */
__attribute__((always_inline)) INLINE static void black_holes_swallow_bpart(
    struct bpart* bpi, const struct bpart* bpj, const struct cosmology* cosmo,
    const double time, const int with_cosmology,
    const struct black_holes_props* props, const struct phys_const* constants) {

  /* Get the current dynamical masses */
  const float bpi_dyn_mass = bpi->mass;
  const float bpj_dyn_mass = bpj->mass;

  /* Is this merger ratio above the threshold for recording? */
  const double merger_ratio = bpj->subgrid_mass / bpi->subgrid_mass;
  if (merger_ratio > props->major_merger_threshold) {
    if (with_cosmology) {
      bpi->last_major_merger_scale_factor = cosmo->a;
    } else {
      bpi->last_major_merger_time = time;
    }
  } else if (merger_ratio > props->minor_merger_threshold) {
    if (with_cosmology) {
      bpi->last_minor_merger_scale_factor = cosmo->a;
    } else {
      bpi->last_minor_merger_time = time;
    }
  }

  /* Increase the masses of the BH. */
  bpi->mass += bpj->mass;
  bpi->gpart->mass += bpj->mass;
  bpi->subgrid_mass += bpj->subgrid_mass;

  /* Collect the swallowed angular momentum */
  bpi->swallowed_angular_momentum[0] += bpj->swallowed_angular_momentum[0];
  bpi->swallowed_angular_momentum[1] += bpj->swallowed_angular_momentum[1];
  bpi->swallowed_angular_momentum[2] += bpj->swallowed_angular_momentum[2];

  /* Update the BH momentum */
  const float BH_mom[3] = {bpi_dyn_mass * bpi->v[0] + bpj_dyn_mass * bpj->v[0],
                           bpi_dyn_mass * bpi->v[1] + bpj_dyn_mass * bpj->v[1],
                           bpi_dyn_mass * bpi->v[2] + bpj_dyn_mass * bpj->v[2]};

  bpi->v[0] = BH_mom[0] / bpi->mass;
  bpi->v[1] = BH_mom[1] / bpi->mass;
  bpi->v[2] = BH_mom[2] / bpi->mass;
  bpi->gpart->v_full[0] = bpi->v[0];
  bpi->gpart->v_full[1] = bpi->v[1];
  bpi->gpart->v_full[2] = bpi->v[2];

  /* Update the BH metal masses */
  struct chemistry_bpart_data* bpi_chem = &bpi->chemistry_data;
  const struct chemistry_bpart_data* bpj_chem = &bpj->chemistry_data;
  chemistry_add_bpart_to_bpart(bpi_chem, bpj_chem);

  /* Update the energy reservoir */
  bpi->jet_mass_reservoir += bpj->jet_mass_reservoir;

  /* Add up all the BH seeds */
  bpi->cumulative_number_seeds += bpj->cumulative_number_seeds;

  /* Add up all the gas particles we swallowed */
  bpi->number_of_gas_swallows += bpj->number_of_gas_swallows;

  /* Add the subgrid angular momentum that we swallowed */
  bpi->accreted_angular_momentum[0] += bpj->accreted_angular_momentum[0];
  bpi->accreted_angular_momentum[1] += bpj->accreted_angular_momentum[1];
  bpi->accreted_angular_momentum[2] += bpj->accreted_angular_momentum[2];

  /* We had another merger */
  bpi->number_of_mergers++;
}

/**
 * @brief Compute the wind launch speed for this feedback step.
 *
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param m_dot_bh The accretion rate onto the black hole.
 * @param m_dot_inflow The accretion rate near the black hole.
 * @param Eddington_rate M_dot,Edd for the black hole.
 * @param BH_state The current state of the black hole.
 */
__attribute__((always_inline)) INLINE static float get_black_hole_wind_speed(
    const struct black_holes_props* props,
    const struct phys_const* constants,
    const float m_dot_bh, const float m_dot_inflow, 
    const float Eddington_rate, const int BH_state) {
  if (m_dot_bh < 0.f || m_dot_inflow < 0.f) return 0.f;
  switch (BH_state) {   
    case BH_states_adaf:
      return props->adaf_wind_speed;
      break;
    case BH_states_quasar:
      return props->quasar_wind_speed;
      break;
    case BH_states_slim_disk:
      return props->slim_disk_wind_speed;
      break;
    default:
      error("Invalid black hole state.");
      return 0.f;
      break;
  }
}

/**
 * @brief Compute the accretion rate of the black hole and all the quantities
 * required for the feedback loop.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param cooling Properties of the cooling model.
 * @param floor_props Properties of the entropy fllor.
 * @param time Time since the start of the simulation (non-cosmo mode).
 * @param with_cosmology Are we running with cosmology?
 * @param dt The time-step size (in physical internal units).
 * @param ti_begin The time at which the step begun (ti_current).
 */
__attribute__((always_inline)) INLINE static void black_holes_prepare_feedback(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct cooling_function_data* cooling,
    const struct entropy_floor_properties* floor_props, const double time,
    const int with_cosmology, const double dt, const integertime_t ti_begin) {

  /* Record that the black hole has another active time step */
  bp->number_of_time_steps++;

  if (dt == 0. || bp->rho_gas == 0.) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bp->group_data.mass <= 0.f) return;

  /* Gather some physical constants (all in internal units) */
  const double G = constants->const_newton_G;
  const double c = constants->const_speed_light_c;
  const double proton_mass = constants->const_proton_mass;
  const double sigma_Thomson = constants->const_thomson_cross_section;

  /* Gather the parameters of the model */
  const double f_Edd_maximum = props->f_Edd_maximum;

  /* (Subgrid) mass of the BH (internal units) */
  const double BH_mass = bp->subgrid_mass;

  /* Convert the quantities we gathered to physical frame (all internal units).
   * Note: for the velocities this means peculiar velocities */
  const double gas_v_norm2 = 0.;

  const float h_inv = 1.f / bp->h;
  double gas_rho = 0.;
  if (props->bondi_use_all_gas) {
    gas_rho = bp->rho_gas;
  } else {
    gas_rho =
        bp->hot_gas_mass * (3. / (4. * M_PI)) * pow_dimension(h_inv);
  }

  const double gas_rho_phys = gas_rho * cosmo->a3_inv;

  double gas_c_phys = 0.;

  if (props->bondi_use_all_gas) {
    if (bp->internal_energy_gas > 0.) {
      gas_c_phys = 
            gas_soundspeed_from_internal_energy(gas_rho, bp->internal_energy_gas)
              * cosmo->a_factor_sound_speed;
    }
  } else {
    if (bp->hot_gas_internal_energy > 0.) {
      gas_c_phys = 
            gas_soundspeed_from_internal_energy(
                gas_rho, bp->hot_gas_internal_energy) *
                    cosmo->a_factor_sound_speed;
    }
  }

  double gas_c_phys2 = gas_c_phys * gas_c_phys;

  /* We can now compute the Bondi accretion rate (internal units) 
   * D. Rennehan: In Simba, we only consider the hot gas within
   * the kernel for the Bondi rate, and the cold gas using the
   * torque accretion estimator.
   */
  double Bondi_rate = 0.;

  const double denominator2 = gas_v_norm2 + gas_c_phys2;
#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure that the denominator is strictly positive */
  if (denominator2 < 0)
    error(
        "Invalid denominator for black hole particle %lld in Bondi rate "
        "calculation. den2=%g v_norm2-%g c_phys2=%g rho=%g",
        bp->id, denominator2, gas_v_norm2, gas_c_phys2, gas_rho);
#endif

  if (denominator2 > 0.) {
    const double denominator_inv = 1. / sqrt(denominator2);
    Bondi_rate = 4. * M_PI * G * G * BH_mass * BH_mass * gas_rho_phys *
                  denominator_inv * denominator_inv * denominator_inv;
  }


  /* Compute the Eddington rate (internal units).
   * IMPORTANT: epsilon_r = 0.1 is the SET value for the Eddington rate.
   * It is assumed in all of the derivations for the state model.
   */
  const double Eddington_rate =
      4. * M_PI * G * BH_mass * proton_mass / (0.1 * c * sigma_Thomson);

  /* In the case of standard Bondi, the accretion rate can runaway
   * significantly, and we MUST limit it. */
  if (props->bondi_use_all_gas) {
    /* Now we can Eddington limit. */
    Bondi_rate = min(Bondi_rate, f_Edd_maximum * Eddington_rate);
  }

  /* The accretion rate estimators give Mdot,inflow  
   * (Mdot,BH = f_acc * Mdot,inflow) */
  double accr_rate = props->bondi_alpha * Bondi_rate;

#ifdef RENNEHAN_DEBUG_CHECKS
  message("BH_ACCRETION: bondi accretion rate id=%lld, %g Msun/yr", 
      bp->id, accr_rate * props->mass_to_solar_mass / props->time_to_yr);
#endif

  /* Here the accretion rate is only based on Mgas / tdyn.
   * We do not use the DM mass to compute tdyn since it probably
   * doesn't contribute much near the core of the system. We also
   * assume that the gas fraction is constant in the galaxy in
   * order to compute Mstar within the kernel of the black hole.
   * Therefore, Mdot = Mgas / tdyn = Mgas / sqrt(3pi/(32 G rho))
   * and rho = (Mgas + Mstar + Mdm) / (4pi h^3 / 3) where
   * Mstar = Mgas / fgas, Mdm = 0. Therefore,
   * rho = 3 * ((1 + fgas) / fgas) * Mgas / (4 * pi * h^3)
   * and
   * Mdot = Mgas * sqrt(32 * G * 3 * ((1 + fgas) / fgas) * Mgas)) /
   *    sqrt(3 * pi * 4 * pi * h^3)
   *      = sqrt(96 * G * ((1 + fgas) / fgas) * Mgas^3) /
   *    sqrt(12 * pi^2 * h^3)
   *      = (1 / pi) * sqrt(8 * G * ((1 + fgas) / fgas) * (Mgas / h)^3))
   */
  float torque_accr_rate = 0.f;
  double f_corr_stellar = 10.;
  if (props->torque_accretion_norm > 0.f) {
    const float galaxy_gas_mass
          = bp->group_data.mass - bp->group_data.stellar_mass;
    if (galaxy_gas_mass > 0) {
      f_corr_stellar = min(
        1. + bp->group_data.stellar_mass / galaxy_gas_mass, 
        f_corr_stellar
      );
    }

    const double rho_bh = bp->subgrid_mass / (4.18879 * bp->h * bp->h * bp->h);
    const double tdyn_inv = sqrt(
      G * f_corr_stellar * (bp->rho_gas + rho_bh) * cosmo->a3_inv
    );

    torque_accr_rate 
          = props->torque_accretion_norm * bp->cold_disk_mass * tdyn_inv;
    accr_rate += torque_accr_rate;
  }

#ifdef RENNEHAN_DEBUG_CHECKS
  message("BH_TORQUE: f_corr_stellar=%g, cold_disk_mass=%g, "
          "rho_gas=%g, torque_accretion_norm=%g",
          f_corr_stellar, bp->cold_disk_mass * props->mass_to_solar_mass, 
          gas_rho_phys * props->rho_to_n_cgs, 
          props->torque_accretion_norm);

  message("BH_ACCRETION: torque accretion rate id=%lld, %g Msun/yr", 
          bp->id, 
          torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr);
#endif

  if (props->bh_characteristic_suppression_mass > 0.) {
     accr_rate *= 1. - exp(-bp->subgrid_mass * props->mass_to_solar_mass /
                                 props->bh_characteristic_suppression_mass);
  }

  /* Right now this is M_dot,inflow. We will multiply by 
   * f_accretion later to make it M_dot,acc */
  bp->accretion_rate = accr_rate;

#ifdef RENNEHAN_DEBUG_CHECKS
  message("BH_STATES: id=%lld, old_state=%d",
          bp->id, bp->state);
#endif

  /* We will use eddington_fraction_lower_boundary and 
   * eddington_fraction_upper_boundary to divide up the accretion rate 
   * in three regimes.
   * 
   * In order to switch out of a regime (i.e. a state), it is necessary 
   * for the true accretion rate (compared to Eddington rate) to switch 
   * over the boundary. Therefore, before we switch a state we must calculate 
   * what the previous state predicts the true accretion rate onto the SMBH is, 
   * and then update the state if it crosses a boundary.
   */
  float predicted_mdot_medd = 0.f;
  switch (bp->state) {
    case BH_states_adaf:
      predicted_mdot_medd 
            = accr_rate * props->adaf_f_accretion / Eddington_rate;

      if (predicted_mdot_medd > props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_slim_disk;
        break;
      }
      if (predicted_mdot_medd > props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_quasar;
      }

      break; /* end case ADAF */
    case BH_states_quasar:
      predicted_mdot_medd 
            = accr_rate * props->quasar_f_accretion / Eddington_rate;

      if (BH_mass > props->adaf_mass_limit &&
          predicted_mdot_medd < props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_adaf;
        break;
      }

      if (predicted_mdot_medd > props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_slim_disk;
      }
  
      break; /* end case quasar */
    case BH_states_slim_disk:
      predicted_mdot_medd = 
          get_black_hole_upper_mdot_medd(props, constants, 
                                         accr_rate / Eddington_rate);

      if (BH_mass > props->adaf_mass_limit &&
          predicted_mdot_medd < props->eddington_fraction_lower_boundary) {
        bp->state = BH_states_adaf;
        break;
      }

      if (predicted_mdot_medd < props->eddington_fraction_upper_boundary) {
        bp->state = BH_states_quasar;
      }

      break; /* end case slim disk */
    default:
      error("Invalid black hole state.");
      break;
  }

  /* We need to store the full M_dot,inflow rate to calculate the 
   * fraction at high accretion rate */
  bp->m_dot_inflow = accr_rate;

  /* This depends on the new state */
  bp->f_accretion = 
      get_black_hole_accretion_factor(props, constants, bp->m_dot_inflow,
                                      BH_mass, bp->state, Eddington_rate);

  /* This accretion rate is M_dot,acc. */
  bp->accretion_rate *= bp->f_accretion;
  accr_rate *= bp->f_accretion;

  if (!props->bondi_use_all_gas) {
    /* Now we can Eddington limit. */
    bp->accretion_rate 
          = min(bp->accretion_rate, f_Edd_maximum * Eddington_rate);
  }

  bp->eddington_fraction = accr_rate / Eddington_rate;

  /* Get the new radiative efficiency based on the new state */
  bp->radiative_efficiency = 
      get_black_hole_radiative_efficiency(props, bp->eddington_fraction, 
                                          bp->state);
      
  double mass_rate = 0.;
  const double luminosity = bp->radiative_efficiency * accr_rate * c * c;

  /* Factor in the radiative efficiency, don't subtract 
   * jet BZ efficiency (spin is fixed) */
  mass_rate = (1. - bp->radiative_efficiency) * accr_rate;

  /* This is used for X-ray feedback later */
  bp->radiative_luminosity = luminosity;

  /* Integrate forward in time */
  bp->subgrid_mass += mass_rate * dt;
  bp->total_accreted_mass += mass_rate * dt;

  if (bp->state == BH_states_adaf) {
    /* ergs to dump in a kernel-weighted fashion */
    bp->adaf_energy_to_dump = get_black_hole_coupling(props, bp->state) *
                              props->adaf_disk_efficiency *
                              bp->accretion_rate * c * c;
  }

  if (bp->state == BH_states_adaf || 
        (props->slim_disk_jet_active && bp->state == BH_states_slim_disk)) { 
    /* Psi_jet*M_dot,acc*dt is the total mass expected in the jet this step */
    bp->jet_mass_reservoir += props->jet_mass_loading * bp->accretion_rate * dt;
  }

  if (bp->subgrid_mass < bp->mass) {
    /* In this case, the BH is still accreting from its (assumed) subgrid gas
     * mass reservoir left over when it was formed. There is some loss in this
     * due to radiative losses, so we must decrease the particle mass
     * in proportion to its current accretion rate. We do not account for this
     * in the swallowing approach, however. */
    bp->mass -= bp->radiative_efficiency * accr_rate * dt;

    if (bp->mass < 0)
      error("Black hole %lld reached negative mass (%g). Trouble ahead...",
            bp->id, bp->mass);
  }

  /* Increase the subgrid angular momentum according to what we accreted
   * Note that this is already in physical units, a factors from velocity and
   * radius cancel each other. Also, the circular velocity contains an extra
   * smoothing length factor that we undo here. */
  bp->accreted_angular_momentum[0] +=
      bp->circular_velocity_gas[0] * mass_rate * dt / bp->h;
  bp->accreted_angular_momentum[1] +=
      bp->circular_velocity_gas[1] * mass_rate * dt / bp->h;
  bp->accreted_angular_momentum[2] +=
      bp->circular_velocity_gas[2] * mass_rate * dt / bp->h;

  /* Keep v_kick physical, there are a lot of comparisons */
  bp->v_kick = get_black_hole_wind_speed(
                  props, constants, bp->accretion_rate, bp->m_dot_inflow,
                  Eddington_rate, bp->state);

#ifdef RENNEHAN_DEBUG_CHECKS
  message("BH_STATES: id=%lld, new_state=%d, predicted_mdot_medd=%g, "
          "eps_r=%g, f_Edd=%g, f_acc=%g, "
          "luminosity=%g, accr_rate=%g Msun/yr, coupling=%g, v_kick=%g km/s, "
          "jet_mass_reservoir=%g Msun",
          bp->id,
          bp->state, 
          predicted_mdot_medd, 
          bp->radiative_efficiency,
          bp->eddington_fraction,
          bp->f_accretion, 
          bp->radiative_luminosity * props->conv_factor_energy_rate_to_cgs, 
          accr_rate * props->mass_to_solar_mass / props->time_to_yr,  
          get_black_hole_coupling(props, bp->state), 
          bp->v_kick / props->kms_to_internal,
          bp->jet_mass_reservoir * props->mass_to_solar_mass);
#endif

  printf("BH_DETAILS "
         "%2.12f %lld "
         " %g %g %g %g %g %g %g "
         " %g %g %g %g " 
         " %g %g %g %g %g "
         " %2.10f %2.10f %2.10f "
         " %2.7f %2.7f %2.7f "
         " %g %g %g  %g %g %g"
         " %g %d %g %g"
         " %g %g\n",
         cosmo->a,
         bp->id,
         bp->mass * props->mass_to_solar_mass, 
         bp->subgrid_mass * props->mass_to_solar_mass, 
         0.f, 
         bp->accretion_rate * props->mass_to_solar_mass / props->time_to_yr, 
         Bondi_rate * props->mass_to_solar_mass / props->time_to_yr, 
         torque_accr_rate * props->mass_to_solar_mass / props->time_to_yr, 
         dt * props->time_to_Myr,
         (bp->rho_gas * cosmo->a3_inv) * props->rho_to_n_cgs, 
         bp->hot_gas_internal_energy * cosmo->a_factor_internal_energy * 
             props->conv_factor_specific_energy_to_cgs, 
         0.f /* SFR */, 
         (bp->hot_gas_mass + bp->cold_gas_mass) * props->mass_to_solar_mass, 
         bp->hot_gas_mass * props->mass_to_solar_mass, 
         bp->stellar_mass * props->mass_to_solar_mass, 
         0.f /* Mgas,bulge */, 
         bp->stellar_bulge_mass * props->mass_to_solar_mass, 
         0.f,
         bp->x[0] * cosmo->a * props->length_to_parsec / 1.0e3f, 
         bp->x[1] * cosmo->a * props->length_to_parsec / 1.0e3f, 
         bp->x[2] * cosmo->a * props->length_to_parsec / 1.0e3f, 
         bp->v[0] * cosmo->a_inv / props->kms_to_internal, 
         bp->v[1] * cosmo->a_inv / props->kms_to_internal, 
         bp->v[2] * cosmo->a_inv / props->kms_to_internal,
         bp->angular_momentum_gas[0], 
         bp->angular_momentum_gas[1], 
         bp->angular_momentum_gas[2],
         bp->specific_angular_momentum_stars[0], 
         bp->specific_angular_momentum_stars[1], 
         bp->specific_angular_momentum_stars[2],
         bp->radiative_luminosity * props->conv_factor_energy_rate_to_cgs,
         bp->state,
         bp->f_accretion,
         bp->radiative_efficiency,
         bp->eddington_fraction,
         (gas_rho * cosmo->a3_inv) * props->rho_to_n_cgs);
}

/**
 * @brief Computes the (maximal) repositioning speed for a black hole.
 *
 * Calculated as upsilon * (m_BH / m_ref) ^ beta_m * (n_H_BH / n_ref) ^ beta_n
 * where m_BH = BH subgrid mass, n_H_BH = physical gas density around BH
 * and upsilon, m_ref, beta_m, n_ref, and beta_n are parameters.
 *
 * @param bp The #bpart.
 * @param props The properties of the black hole model.
 * @param cosmo The current cosmological model.
 */
__attribute__((always_inline)) INLINE static double
black_holes_get_repositioning_speed(const struct bpart* restrict bp,
                                    const struct black_holes_props* props,
                                    const struct cosmology* cosmo) {

  const double n_gas_phys = bp->rho_gas * cosmo->a3_inv * props->rho_to_n_cgs;
  const double v_repos =
      props->reposition_coefficient_upsilon *
      pow(bp->subgrid_mass / props->reposition_reference_mass,
          props->reposition_exponent_mass) *
      pow(n_gas_phys / props->reposition_reference_n_H,
          props->reposition_exponent_n_H);

  /* Make sure the repositioning is not back-firing... */
  if (v_repos < 0)
    error(
        "BH %lld wants to reposition at negative speed (%g U_V). Do you "
        "think you are being funny? No-one is laughing.",
        bp->id, v_repos);

  return v_repos;
}

/**
 * @brief Finish the calculation of the new BH position.
 *
 * Here, we check that the BH should indeed be moved in the next drift.
 *
 * @param bp The black hole particle.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants (in internal units).
 * @param cosmo The cosmological model.
 * @param dt The black hole particle's time step.
 * @param ti_begin The time at the start of the temp
 */
__attribute__((always_inline)) INLINE static void black_holes_end_reposition(
    struct bpart* restrict bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const double dt, const integertime_t ti_begin) {

  /* First check: did we find any eligible neighbour particle to jump to? */
  if (bp->reposition.min_potential != FLT_MAX) {

    /* Record that we have a (possible) repositioning situation */
    bp->number_of_reposition_attempts++;

    /* Is the potential lower (i.e. the BH is at the bottom already)
     * OR is the BH massive enough that we don't reposition? */
    const float potential = gravity_get_comoving_potential(bp->gpart);
    if (potential < bp->reposition.min_potential ||
        bp->subgrid_mass > props->max_reposition_mass) {

      /* No need to reposition */
      bp->reposition.min_potential = FLT_MAX;
      bp->reposition.delta_x[0] = -FLT_MAX;
      bp->reposition.delta_x[1] = -FLT_MAX;
      bp->reposition.delta_x[2] = -FLT_MAX;

    } else if (props->set_reposition_speed) {

      /* If we are re-positioning, move the BH a fraction of delta_x, so
       * that we have a well-defined re-positioning velocity (repos_vel
       * cannot be negative). */
      double repos_vel = black_holes_get_repositioning_speed(bp, props, cosmo);

      /* Convert target reposition velocity to a fractional reposition
       * along reposition.delta_x */
      const double dx = bp->reposition.delta_x[0];
      const double dy = bp->reposition.delta_x[1];
      const double dz = bp->reposition.delta_x[2];
      const double d = sqrt(dx * dx + dy * dy + dz * dz);

      /* Exclude the pathological case of repositioning by zero distance */
      if (d > 0) {
        double repos_frac = repos_vel * dt / d;

        /* We should never get negative repositioning fractions... */
        if (repos_frac < 0)
          error("Wanting to reposition by negative fraction (%g)?", repos_frac);

        /* ... but fractions > 1 can occur if the target velocity is high.
         * We do not want this, because it could lead to overshooting the
         * actual potential minimum. */
        if (repos_frac > 1) {
          repos_frac = 1.;
          repos_vel = repos_frac * d / dt;
        }

        bp->last_repos_vel = (float)repos_vel;
        bp->reposition.delta_x[0] *= repos_frac;
        bp->reposition.delta_x[1] *= repos_frac;
        bp->reposition.delta_x[2] *= repos_frac;
      }

      /* ends section for fractional repositioning */
    } else {

      /* We _should_ reposition, but not fractionally. Here, we will
       * reposition exactly on top of another gas particle - which
       * could cause issues, so we add on a small fractional offset
       * of magnitude 0.001 h in the reposition delta. */

      /* Generate three random numbers in the interval [-0.5, 0.5[; id,
       * id**2, and id**3 are required to give unique random numbers (as
       * random_unit_interval is completely reproducible). */
      const float offset_dx =
          random_unit_interval(bp->id, ti_begin, random_number_BH_reposition) -
          0.5f;
      const float offset_dy =
          random_unit_interval(bp->id * bp->id, ti_begin,
                               random_number_BH_reposition) -
          0.5f;
      const float offset_dz =
          random_unit_interval(bp->id * bp->id * bp->id, ti_begin,
                               random_number_BH_reposition) -
          0.5f;

      const float length_inv =
          1.0f / sqrtf(offset_dx * offset_dx + offset_dy * offset_dy +
                       offset_dz * offset_dz);

      const float norm = 0.001f * bp->h * length_inv;

      bp->reposition.delta_x[0] += offset_dx * norm;
      bp->reposition.delta_x[1] += offset_dy * norm;
      bp->reposition.delta_x[2] += offset_dz * norm;
    }
  } /* ends section if we found eligible repositioning target(s) */
}

/**
 * @brief Reset acceleration fields of a particle
 *
 * This is the equivalent of hydro_reset_acceleration.
 * We do not compute the acceleration on black hole, therefore no need to use
 * it.
 *
 * @param bp The particle to act upon
 */
__attribute__((always_inline)) INLINE static void black_holes_reset_feedback(
    struct bpart* restrict bp) {

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
  for (int i = 0; i < MAX_NUM_OF_NEIGHBOURS_STARS; ++i)
    bp->ids_ngbs_force[i] = -1;
  bp->num_ngb_force = 0;
#endif
}

/**
 * @brief Store the gravitational potential of a black hole by copying it from
 * its #gpart friend.
 *
 * @param bp The black hole particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_bpart(struct bpart* bp, const struct gpart* gp) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->gpart != gp) error("Copying potential to the wrong black hole!");
#endif

  bp->reposition.potential = gp->potential;
}

/**
 * @brief Store the gravitational potential of a particle by copying it from
 * its #gpart friend.
 *
 * @param p_data The black hole data of a gas particle.
 * @param gp The black hole's #gpart.
 */
__attribute__((always_inline)) INLINE static void
black_holes_store_potential_in_part(struct black_holes_part_data* p_data,
                                    const struct gpart* gp) {
  p_data->potential = gp->potential;
}

/**
 * @brief Initialise a BH particle that has just been seeded.
 *
 * @param bp The #bpart to initialise.
 * @param props The properties of the black hole scheme.
 * @param constants The physical constants in internal units.
 * @param cosmo The current cosmological model.
 * @param p The #part that became a black hole.
 * @param xp The #xpart that became a black hole.
 */
INLINE static void black_holes_create_from_gas(
    struct bpart* bp, const struct black_holes_props* props,
    const struct phys_const* constants, const struct cosmology* cosmo,
    const struct part* p, const struct xpart* xp,
    const integertime_t ti_current) {

  /* All the non-basic properties of the black hole have been zeroed
   * in the FOF code. We update them here.
   * (i.e. position, velocity, mass, time-step have been set) */

  /* Birth time and density */
  bp->formation_scale_factor = cosmo->a;
  bp->formation_gas_density = hydro_get_physical_density(p, cosmo);

  /* Initial seed mass */
  bp->subgrid_mass = props->subgrid_seed_mass;

  /* We haven't accreted anything yet */
  bp->total_accreted_mass = 0.f;
  bp->cumulative_number_seeds = 1;
  bp->number_of_mergers = 0;
  bp->number_of_gas_swallows = 0;
  bp->number_of_direct_gas_swallows = 0;
  bp->number_of_time_steps = 0;

  /* We haven't repositioned yet, nor attempted it */
  bp->number_of_repositions = 0;
  bp->number_of_reposition_attempts = 0;
  bp->last_repos_vel = 0.f;

  /* Copy over the splitting struct */
  bp->split_data = xp->split_data;

  /* Initial metal masses */
  const float gas_mass = hydro_get_mass(p);
  struct chemistry_bpart_data* bp_chem = &bp->chemistry_data;
  const struct chemistry_part_data* p_chem = &p->chemistry_data;
  chemistry_bpart_from_part(bp_chem, p_chem, gas_mass);

  /* No swallowed angular momentum */
  bp->swallowed_angular_momentum[0] = 0.f;
  bp->swallowed_angular_momentum[1] = 0.f;
  bp->swallowed_angular_momentum[2] = 0.f;

  /* Last time of mergers */
  bp->last_minor_merger_time = -1.;
  bp->last_major_merger_time = -1.;

  /* First initialisation */
  black_holes_init_bpart(bp);

  bp->state = BH_states_slim_disk;
  
  black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);
}

/**
 * @brief Should this bh particle be doing any stars looping?
 *
 * @param bp The #bpart.
 * @param e The #engine.
 */
__attribute__((always_inline)) INLINE static int bh_stars_loop_is_active(
    const struct bpart* bp, const struct engine* e) {
  /* Active bhs never do the stars loop for the Rennehan model */
  return 0;
}

#endif /* SWIFT_RENNEHAN_BLACK_HOLES_H */