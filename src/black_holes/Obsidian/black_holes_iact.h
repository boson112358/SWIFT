/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2021 Edo Altamura (edoardo.altamura@manchester.ac.uk)
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
#ifndef SWIFT_OBSIDIAN_BH_IACT_H
#define SWIFT_OBSIDIAN_BH_IACT_H

/* Local includes */
#include "black_holes_parameters.h"
#include "black_holes_properties.h"
#include "entropy_floor.h"
#include "equation_of_state.h"
#include "gravity.h"
#include "gravity_iact.h"
#include "hydro.h"
#include "random.h"
#include "space.h"
#include "timestep_sync_part.h"
#include "tracers.h"

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param bi First particle (black hole).
 * @param sj Second particle (stars, not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_stars_density(
    const float r2, const float dx[3],
    struct bpart *bi, const struct spart *sj) {

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the black hole) */
  const float dv[3] = {sj->v[0] - bi->v[0], sj->v[1] - bi->v[1],
                       sj->v[2] - bi->v[2]};

  bi->specific_angular_momentum_stars[0] += dx[1] * dv[2] - dx[2] * dv[1];
  bi->specific_angular_momentum_stars[1] += dx[2] * dv[0] - dx[0] * dv[2];
  bi->specific_angular_momentum_stars[2] += dx[0] * dv[1] - dx[1] * dv[0];  
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param bi First particle (black hole).
 * @param sj Second particle (stars, not updated).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_stars_bulge(
    const float r2, const float dx[3],
    struct bpart *bi, const struct spart *sj) {

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the black hole) */
  const float dv[3] = {sj->v[0] - bi->v[0], sj->v[1] - bi->v[1],
                       sj->v[2] - bi->v[2]};

  const float star_angular_momentum[3] = {dx[1] * dv[2] - dx[2] * dv[1],
                                          dx[2] * dv[0] - dx[0] * dv[2],
                                          dx[0] * dv[1] - dx[1] * dv[0]};

  const float proj =
      star_angular_momentum[0] * bi->specific_angular_momentum_stars[0] +
      star_angular_momentum[1] * bi->specific_angular_momentum_stars[1] +
      star_angular_momentum[2] * bi->specific_angular_momentum_stars[2];

  bi->stellar_mass += sj->mass;
  /* The bulge mass is twice the counter rotating mass */
  if (proj < 0.f) bi->stellar_bulge_mass += 2.f * sj->mass;
}

/**
 * @brief Density interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas, not updated).
 * @param xpj The extended data of the second particle (not updated).
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time Current physical time in the simulation.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_density(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time,
    const double time_base) {

  /* Ignore decoupled winds in density computation */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  float wi, wi_dx;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the number of neighbours */
  bi->density.wcount += wi;
  bi->density.wcount_dh -= (hydro_dimension * wi + ui * wi_dx);

  /* Contribution to the number of neighbours */
  bi->num_ngbs += 1;

  /* Neighbour gas mass */
  const float mj = hydro_get_mass(pj);

  /* Contribution to the BH gas density */
  bi->rho_gas += mj * wi;

  /* Contribution to the total neighbour mass */
  bi->ngb_mass += mj;

  /* Contribution to the smoothed sound speed */
  const float cj = hydro_get_comoving_soundspeed(pj);
  bi->sound_speed_gas += mj * wi * cj;

  /* Neighbour internal energy */
  const float uj = hydro_get_drifted_comoving_internal_energy(pj);

  /* Contribution to the smoothed internal energy */
  bi->internal_energy_gas += mj * uj * wi;

  /* Neighbour's (drifted) velocity in the frame of the black hole
   * (we don't include a Hubble term since we are interested in the
   * velocity contribution at the location of the black hole) */
  const float dv[3] = {pj->v[0] - bi->v[0], pj->v[1] - bi->v[1],
                       pj->v[2] - bi->v[2]};

  /* Account for hot and cold gas surrounding the SMBH */
  const float Tj =
      uj * cosmo->a_factor_internal_energy / bh_props->temp_to_u_factor;
  int is_hot_gas = 0;
  const float rho_crit_0 = cosmo->critical_density_0;
  const float rho_crit_baryon = cosmo->Omega_b * rho_crit_0;
  const double rho_com = hydro_get_comoving_density(pj);
  const double rho_phys = hydro_get_physical_density(pj, cosmo);
  /* Are we in the regime of the Jeans equation of state? */
  if ((rho_com >= rho_crit_baryon * floor_props->Jeans_over_density_threshold) &&
      (rho_phys >= floor_props->Jeans_density_threshold)) {
    const float T_EoS = entropy_floor_temperature(pj, cosmo, floor_props);
    /* Only hot if above a small region of the EoS */
    if (Tj > T_EoS * bh_props->fixed_T_above_EoS_factor) {
      is_hot_gas = 1;
    }
  }
  else {
    /* If off of the EoS always hot above the temperature cut */
    if (Tj > bh_props->environment_temperature_cut) {
      is_hot_gas = 1;
    }
  }

  /* Star forming gas is never considered "hot" */
  if (pj->sf_data.SFR > 0.f) is_hot_gas = 0;

  if (is_hot_gas) {
    bi->hot_gas_mass += mj;
    bi->hot_gas_internal_energy += mj * uj; /* Not kernel weighted */
  } else {
    bi->cold_gas_mass += mj;
  }

  /* Gas angular momentum in kernel */
  bi->angular_momentum_gas[0] -= mj * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->angular_momentum_gas[1] -= mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->angular_momentum_gas[2] -= mj * (dx[2] * dv[0] - dx[0] * dv[2]);

  /* Contribution to the smoothed velocity (gas w.r.t. black hole) */
  bi->velocity_gas[0] += mj * wi * dv[0];
  bi->velocity_gas[1] += mj * wi * dv[1];
  bi->velocity_gas[2] += mj * wi * dv[2];

  /* Contribution to the specific angular momentum of gas, which is later
   * converted to the circular velocity at the smoothing length */
  bi->circular_velocity_gas[0] -= mj * wi * (dx[1] * dv[2] - dx[2] * dv[1]);
  bi->circular_velocity_gas[1] -= mj * wi * (dx[2] * dv[0] - dx[0] * dv[2]);
  bi->circular_velocity_gas[2] -= mj * wi * (dx[0] * dv[1] - dx[1] * dv[0]);

#ifdef DEBUG_INTERACTIONS_BH
  /* Update ngb counters */
  if (si->num_ngb_density < MAX_NUM_OF_NEIGHBOURS_BH)
    bi->ids_ngbs_density[si->num_ngb_density] = pj->id;

  /* Update ngb counters */
  ++si->num_ngb_density;
#endif

}

/**
 * @brief Repositioning interaction between two particles (non-symmetric).
 *
 * Function used to identify the gas particle that this BH may move towards.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time Current physical time in the simulation.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_repos(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, const struct part *pj, const struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time,
    const double time_base) {

  /* Ignore cooling shut off particles */
  if (pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  /* Ignore decoupled wind particles */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  float wi;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
  const float hi_inv = 1.0f / hi;
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* (Square of) Max repositioning distance allowed based on the softening */
  const float max_dist_repos2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      bh_props->max_reposition_distance_ratio *
      bh_props->max_reposition_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* Is this gas neighbour close enough that we can consider its potential
     for repositioning? */
  if (r2 < max_dist_repos2) {

    /* Flag to check whether neighbour is slow enough to be considered
     * as repositioning target. Always true if velocity cut is switched off. */
    int neighbour_is_slow_enough = 1;
    if (bh_props->with_reposition_velocity_threshold) {

      /* Compute relative peculiar velocity between the two BHs
       * Recall that in SWIFT v is (v_pec * a) */
      const float delta_v[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                                bi->v[2] - pj->v[2]};
      const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                       delta_v[2] * delta_v[2];
      const float v2_pec = v2 * cosmo->a2_inv;

      /* Compute the maximum allowed velocity */
      float v2_max = bh_props->max_reposition_velocity_ratio *
                     bh_props->max_reposition_velocity_ratio *
                     bi->sound_speed_gas * bi->sound_speed_gas *
                     cosmo->a_factor_sound_speed * cosmo->a_factor_sound_speed;

      /* If desired, limit the value of the threshold (v2_max) to be no
       * smaller than a user-defined value */
      if (bh_props->min_reposition_velocity_threshold > 0) {
        const float v2_min_thresh =
            bh_props->min_reposition_velocity_threshold *
            bh_props->min_reposition_velocity_threshold;
        v2_max = max(v2_max, v2_min_thresh);
      }

      /* Is the neighbour too fast to jump to? */
      if (v2_pec >= v2_max) neighbour_is_slow_enough = 0;
    }

    if (neighbour_is_slow_enough) {
      float potential = pj->black_holes_data.potential;

      if (bh_props->correct_bh_potential_for_repositioning) {

        /* Let's not include the contribution of the BH
         * itself to the potential of the gas particle */

        /* Note: This assumes the BH and gas have the same
         * softening, which is currently true */
        const float eps = gravity_get_softening(bi->gpart, grav_props);
        const float eps2 = eps * eps;
        const float eps_inv = 1.f / eps;
        const float eps_inv3 = eps_inv * eps_inv * eps_inv;
        const float BH_mass = bi->mass;

        /* Compute the Newtonian or truncated potential the BH
         * exherts onto the gas particle */
        float dummy, pot_ij;
        runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, BH_mass, &dummy,
                                 &pot_ij);

        /* Deduct the BH contribution */
        potential -= pot_ij * grav_props->G_Newton;
      }

      /* Is the potential lower? */
      if (potential < bi->reposition.min_potential) {

        /* Store this as our new best */
        bi->reposition.min_potential = potential;
        bi->reposition.delta_x[0] = -dx[0];
        bi->reposition.delta_x[1] = -dx[1];
        bi->reposition.delta_x[2] = -dx[2];
      }
    }
  }
}

/**
 * @brief Swallowing interaction between two particles (non-symmetric).
 *
 * Function used to flag the gas particles that will be swallowed
 * by the black hole particle.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time Current physical time in the simulation.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_swallow(
    const float r2, const float dx[3], const float hi, const float hj,
    struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time,
    const double time_base) {

  /* Skip cooling delayed particles */
  if (pj->feedback_data.cooling_shutoff_delay_time > 0.f) return;

  /* Do not even consider wind particles for accretion/feedback */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bi->group_data.mass <= 0.f) return;

  /* If there is no gas, skip */
  if (bi->rho_gas <= 0.f) return;

  float wi;

  /* Compute the kernel function; note that r cannot be optimised
   * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
  const float r = sqrtf(r2);
  const float hi_inv = 1.0f / hi;
  const float hi_inv_dim = pow_dimension(hi_inv);
  const float ui = r * hi_inv;
  kernel_eval(ui, &wi);

  /* Sum up cold disk mass corotating relative to total angular momentum. */
  /* Neighbour internal energy */
  const float uj = hydro_get_drifted_comoving_internal_energy(pj);
  const float mj = hydro_get_mass(pj);

  /* Account for hot and cold gas surrounding the SMBH */
  const float Tj =
      uj * cosmo->a_factor_internal_energy / bh_props->temp_to_u_factor;
  int is_hot_gas = 0;
  const float rho_crit_0 = cosmo->critical_density_0;
  const float rho_crit_baryon = cosmo->Omega_b * rho_crit_0;
  const double rho_com = hydro_get_comoving_density(pj);
  const double rho_phys = hydro_get_physical_density(pj, cosmo);

  /* Are we in the regime of the Jeans equation of state? */
  if ((rho_com >= rho_crit_baryon * floor_props->Jeans_over_density_threshold) &&
      (rho_phys >= floor_props->Jeans_density_threshold)) {
    const float T_EoS = entropy_floor_temperature(pj, cosmo, floor_props);
    /* Only hot if above a small region of the EoS */
    if (Tj > T_EoS * bh_props->fixed_T_above_EoS_factor) {
      is_hot_gas = 1;
    }
  }
  else {
    /* If off of the EoS always hot above the temperature cut */
    if (Tj > bh_props->environment_temperature_cut) {
      is_hot_gas = 1;
    }
  }

  /* Star forming gas is never considered "hot" */
  if (pj->sf_data.SFR > 0.f) is_hot_gas = 0;

  const float dv[3] = {bi->v[0] - pj->v[0], bi->v[1] - pj->v[1],
                       bi->v[2] - pj->v[2]};
  const float Lx = mj * (dx[1] * dv[2] - dx[2] * dv[1]);
  const float Ly = mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  const float Lz = mj * (dx[2] * dv[0] - dx[0] * dv[2]);
  const float proj = Lx * bi->angular_momentum_gas[0] 
                    + Ly * bi->angular_momentum_gas[1] 
                    + Lz * bi->angular_momentum_gas[2];
  if ((proj > 0.f) && (is_hot_gas == 0)) {
    bi->cold_disk_mass += mj;
  }

  /* Probability to swallow this particle */
  float prob = -1.f;
  float f_accretion = bi->f_accretion;
  if (f_accretion <= 0.f) return;

  /* Radiation was already accounted for in bi->subgrid_mass
   * so if is is bigger than bi->mass we can simply
   * flag particles to eat and satisfy the mass constraint.
   *
   * If bi->subgrid_mass < bi->mass then there is a problem,
   * but we use the full accretion rate to select particles
   * and then don't actually take anything away from them.
   * The bi->mass variable is decreased previously to account
   * for the radiative losses.
   * 
   * The probability only needs to be set if we are not in the ADAF mode,
   * since the ADAF mode only heats gas surrounding the black hole. The jet 
   * identification happens below.
   */
  if (bi->state != BH_states_adaf) {
    /* Get particle time-step */
    double dt;
    if (with_cosmology) {
      const integertime_t ti_step = get_integer_timestep(bi->time_bin);
      const integertime_t ti_begin =
          get_integer_time_begin(ti_current - 1, bi->time_bin);

      dt = cosmology_get_delta_time(cosmo, ti_begin,
                                    ti_begin + ti_step);
    } else {
      dt = get_timestep(bi->time_bin, time_base);
    }

    const float mass_deficit = bi->subgrid_mass - bi->mass_at_start_of_step;
    if (mass_deficit >= 0.f) {
      /* Don't nibble from particles that are too small already */
      if (mj < bh_props->min_gas_mass_for_nibbling) return;

      prob = (mass_deficit / bi->f_accretion) * hi_inv_dim * wi / bi->rho_gas;
    } else {
      prob = ((1.f - bi->f_accretion) / bi->f_accretion) * bi->accretion_rate *
            dt * (hi_inv_dim * wi / bi->rho_gas);
      /* We do NOT accrete when subgrid_mass < physical_mass
      * but we still kick.
      */
      f_accretion = 0.f;
    }
  }

  /* Draw a random number (Note mixing both IDs) */
  const float rand = random_unit_interval(bi->id + pj->id, ti_current,
                                          random_number_BH_swallow);
  const float pj_mass_orig = mj;
  float new_gas_mass = pj_mass_orig;
  /* Are we lucky? */
  if (rand < prob) {

    if (f_accretion > 0.f) {
      const float bi_mass_orig = bi->mass;
      const float nibbled_mass = f_accretion * pj_mass_orig;
      new_gas_mass = pj_mass_orig - nibbled_mass;
      /* Don't go below the minimum for stability */
      if (new_gas_mass < bh_props->min_gas_mass_for_nibbling) return;

      bi->mass += nibbled_mass;
      hydro_set_mass(pj, new_gas_mass);

      /* Add the angular momentum of the accreted gas to the BH total.
       * Note no change to gas here. The cosmological conversion factors for
       * velocity (a^-1) and distance (a) cancel out, so the angular momentum
       * is already in physical units. */
      bi->swallowed_angular_momentum[0] +=
          nibbled_mass * (dx[1] * dv[2] - dx[2] * dv[1]);
      bi->swallowed_angular_momentum[1] +=
          nibbled_mass * (dx[2] * dv[0] - dx[0] * dv[2]);
      bi->swallowed_angular_momentum[2] +=
          nibbled_mass * (dx[0] * dv[1] - dx[1] * dv[0]);

      /* Update the BH momentum and velocity. Again, no change to gas here. */
      const float bi_mom[3] = {
          bi_mass_orig * bi->v[0] + nibbled_mass * pj->v[0],
          bi_mass_orig * bi->v[1] + nibbled_mass * pj->v[1],
          bi_mass_orig * bi->v[2] + nibbled_mass * pj->v[2]};

      bi->v[0] = bi_mom[0] / bi->mass;
      bi->v[1] = bi_mom[1] / bi->mass;
      bi->v[2] = bi_mom[2] / bi->mass;

      /* Update the BH and also gas metal masses */
      struct chemistry_bpart_data *bi_chem = &bi->chemistry_data;
      struct chemistry_part_data *pj_chem = &pj->chemistry_data;
      chemistry_transfer_part_to_bpart(bi_chem, pj_chem, nibbled_mass,
                                       nibbled_mass / pj_mass_orig);

    }

    /* This particle is swallowed by the BH with the largest ID of all the
     * candidates wanting to swallow it */
    if (pj->black_holes_data.swallow_id < bi->id) {
      pj->black_holes_data.swallow_id = bi->id;
    } else {
      message(
          "BH %lld wants to swallow gas particle %lld BUT CANNOT (old "
          "swallow id=%lld)",
          bi->id, pj->id, pj->black_holes_data.swallow_id);
    }
  }

  if (bi->jet_mass_reservoir >= bh_props->jet_minimum_reservoir_mass) {

#ifdef OBSIDIAN_DEBUG_CHECKS
    message("BH_JET: bid=%lld, jet_mass_reservoir=%g",
            bi->id, bi->jet_mass_reservoir);
#endif

    const float jet_prob = bi->jet_mass_reservoir * 
                           (hi_inv_dim * wi / bi->rho_gas);
    const float rand_jet = random_unit_interval(bi->id + pj->id, ti_current,
                                                  random_number_BH_kick);

#ifdef OBSIDIAN_DEBUG_CHECKS
    message("BH_JET: bid=%lld, pid=%lld, rand_jet=%g, jet_prob=%g",
            bi->id, pj->id, rand_jet, jet_prob);
#endif

    /* Here the particle is also identified to be kicked out as a jet */
    if (rand_jet < jet_prob) {

      /* If we also are accreting above, the mass loss is already taken
       * into account */

      if (pj->black_holes_data.jet_id < bi->id) {
        bi->jet_mass_kicked_this_step += new_gas_mass;
        pj->black_holes_data.jet_id = bi->id;
      } else {
        message(
            "BH %lld wants to kick jet particle %lld BUT CANNOT (old "
            "swallow id=%lld)",
            bi->id, pj->id, pj->black_holes_data.jet_id);
      }
    }
  }
}

/**
 * @brief Swallowing interaction between two BH particles (non-symmetric).
 *
 * Function used to identify the BH particle that this BH may move towards.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param bj Second particle (black hole)
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_bh_repos(const float r2, const float dx[3],
                               const float hi, const float hj, struct bpart *bi,
                               const struct bpart *bj,
                               const struct cosmology *cosmo,
                               const struct gravity_props *grav_props,
                               const struct black_holes_props *bh_props,
                               const integertime_t ti_current) {

  /* Compute relative peculiar velocity between the two BHs
   * Recall that in SWIFT v is (v_pec * a) */
  const float delta_v[3] = {bi->v[0] - bj->v[0], bi->v[1] - bj->v[1],
                            bi->v[2] - bj->v[2]};
  const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                   delta_v[2] * delta_v[2];

  const float v2_pec = v2 * cosmo->a2_inv;

  /* (Square of) Max repositioning distance allowed based on the softening */
  const float max_dist_repos2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      bh_props->max_reposition_distance_ratio *
      bh_props->max_reposition_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  /* Is this BH neighbour close enough that we can consider its potential
     for repositioning? */
  if (r2 < max_dist_repos2) {

    /* Flag to check whether neighbour is slow enough to be considered
     * as repositioning target. Always true if velocity cut switched off */
    int neighbour_is_slow_enough = 1;
    if (bh_props->with_reposition_velocity_threshold) {

      /* Compute the maximum allowed velocity */
      float v2_max = bh_props->max_reposition_velocity_ratio *
                     bh_props->max_reposition_velocity_ratio *
                     bi->sound_speed_gas * bi->sound_speed_gas *
                     cosmo->a_factor_sound_speed * cosmo->a_factor_sound_speed;

      /* If desired, limit the value of the threshold (v2_max) to be no
       * smaller than a user-defined value */
      if (bh_props->min_reposition_velocity_threshold > 0) {
        const float v2_min_thresh =
            bh_props->min_reposition_velocity_threshold *
            bh_props->min_reposition_velocity_threshold;
        v2_max = max(v2_max, v2_min_thresh);
      }

      /* Is the neighbour too fast to jump to? */
      if (v2_pec >= v2_max) neighbour_is_slow_enough = 0;
    }

    if (neighbour_is_slow_enough) {
      float potential = bj->reposition.potential;

      if (bh_props->correct_bh_potential_for_repositioning) {

        /* Let's not include the contribution of the BH i
         * to the potential of the BH j */
        const float eps = gravity_get_softening(bi->gpart, grav_props);
        const float eps2 = eps * eps;
        const float eps_inv = 1.f / eps;
        const float eps_inv3 = eps_inv * eps_inv * eps_inv;
        const float BH_mass = bi->mass;

        /* Compute the Newtonian or truncated potential the BH
         * exherts onto the gas particle */
        float dummy, pot_ij;
        runner_iact_grav_pp_full(r2, eps2, eps_inv, eps_inv3, BH_mass, &dummy,
                                 &pot_ij);

        /* Deduct the BH contribution */
        potential -= pot_ij * grav_props->G_Newton;
      }

      /* Is the potential lower? */
      if (potential < bi->reposition.min_potential) {

        /* Store this as our new best */
        bi->reposition.min_potential = potential;
        bi->reposition.delta_x[0] = -dx[0];
        bi->reposition.delta_x[1] = -dx[1];
        bi->reposition.delta_x[2] = -dx[2];
      }
    }
  }
}

/**
 * @brief Swallowing interaction between two BH particles (non-symmetric).
 *
 * Function used to flag the BH particles that will be swallowed
 * by the black hole particle.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param bj Second particle (black hole)
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_bh_swallow(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct bpart *bi, struct bpart *bj,
                                 const struct cosmology *cosmo,
                                 const struct gravity_props *grav_props,
                                 const struct black_holes_props *bh_props,
                                 const integertime_t ti_current) {

  /* Compute relative peculiar velocity between the two BHs
   * Recall that in SWIFT v is (v_pec * a) */
  const float delta_v[3] = {bi->v[0] - bj->v[0], bi->v[1] - bj->v[1],
                            bi->v[2] - bj->v[2]};
  const float v2 = delta_v[0] * delta_v[0] + delta_v[1] * delta_v[1] +
                   delta_v[2] * delta_v[2];

  const float v2_pec = v2 * cosmo->a2_inv;

  /* Find the most massive of the two BHs */
  float M = bi->subgrid_mass;
  float h = hi;
  if (bj->subgrid_mass > M) {
    M = bj->subgrid_mass;
    h = hj;
  }

  /* (Square of) max swallowing distance allowed based on the softening */
  const float max_dist_merge2 =
      kernel_gravity_softening_plummer_equivalent_inv *
      kernel_gravity_softening_plummer_equivalent_inv *
      bh_props->max_merging_distance_ratio *
      bh_props->max_merging_distance_ratio * grav_props->epsilon_baryon_cur *
      grav_props->epsilon_baryon_cur;

  const float G_Newton = grav_props->G_Newton;

  /* The BH with the smaller mass will be merged onto the one with the
   * larger mass.
   * To avoid rounding issues, we additionally check for IDs if the BHs
   * have the exact same mass. */
  if ((bj->subgrid_mass < bi->subgrid_mass) ||
      (bj->subgrid_mass == bi->subgrid_mass && bj->id < bi->id)) {

    /* Merge if gravitationally bound AND if within max distance
     * Note that we use the kernel support here as the size and not just the
     * smoothing length */

    /* Maximum velocity difference between BHs allowed to merge */
    float v2_threshold;

    if (bh_props->merger_threshold_type == BH_mergers_circular_velocity) {

      /* 'Old-style' merger threshold using circular velocity at the
       * edge of the more massive BH's kernel */
      v2_threshold = G_Newton * M / (kernel_gamma * h);
    } else {

      /* Arguably better merger threshold using the escape velocity at
       * the distance between the BHs */

      if (bh_props->merger_threshold_type == BH_mergers_escape_velocity) {

        /* Standard formula (not softening BH interactions) */
        v2_threshold = 2.f * G_Newton * M / sqrt(r2);
      } else if (bh_props->merger_threshold_type ==
                 BH_mergers_dynamical_escape_velocity) {

        /* General two-body escape velocity based on dynamical masses */
        v2_threshold = 2.f * G_Newton * (bi->mass + bj->mass) / sqrt(r2);
      } else {
        /* Cannot happen! */
#ifdef SWIFT_DEBUG_CHECKS
        error("Invalid choice of BH merger threshold type");
#endif
        v2_threshold = 0.f;
      }
    } /* Ends sections for different merger thresholds */

    if ((v2_pec < v2_threshold) && (r2 < max_dist_merge2)) {

      /* This particle is swallowed by the BH with the largest mass of all the
       * candidates wanting to swallow it (we use IDs to break ties)*/
      if ((bj->merger_data.swallow_mass < bi->subgrid_mass) ||
          (bj->merger_data.swallow_mass == bi->subgrid_mass &&
           bj->merger_data.swallow_id < bi->id)) {

#ifdef SWIFT_DEBUG_CHECKS
        message("BH %lld wants to swallow BH particle %lld", bi->id, bj->id);
#endif

        bj->merger_data.swallow_id = bi->id;
        bj->merger_data.swallow_mass = bi->subgrid_mass;

      } else {

        message(
            "BH %lld wants to swallow bh particle %lld BUT CANNOT (old "
            "swallow id=%lld)",
            bi->id, bj->id, bj->merger_data.swallow_id);
      }
    }
  }
}

/**
 * @brief Feedback interaction between two particles (non-symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param bi First particle (black hole).
 * @param pj Second particle (gas)
 * @param xpj The extended data of the second particle.
 * @param with_cosmology Are we doing a cosmological run?
 * @param cosmo The cosmological model.
 * @param grav_props The properties of the gravity scheme (softening, G, ...).
 * @param bh_props The properties of the BH scheme
 * @param ti_current Current integer time value (for random numbers).
 * @param time current physical time in the simulation
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_bh_gas_feedback(
    const float r2, const float dx[3], const float hi, const float hj,
    const struct bpart *bi, struct part *pj, struct xpart *xpj,
    const int with_cosmology, const struct cosmology *cosmo,
    const struct gravity_props *grav_props,
    const struct black_holes_props *bh_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_current, const double time,
    const double time_base) {

  /* This shouldn't happen, but just be sure anyway */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return;

  /* A black hole should never accrete/feedback if it is not in a galaxy */
  if (bi->group_data.mass <= 0.f) return;

  /* A black hole should have gas surrounding it. */
  if (bi->rho_gas <= 0.f) return;

  /* Save gas density and entropy before feedback */
  tracers_before_black_holes_feedback(pj, xpj, cosmo->a);

  float v_kick = bi->v_kick;  /* PHYSICAL */
  int jet_flag = 0;
  float dir[3];
  float norm;

  /* In the swallow loop the particle was marked as a jet particle */
  if (pj->black_holes_data.jet_id == bi->id) {
    v_kick = bh_props->jet_velocity; 
    jet_flag = 1;
  }

  /* Only heat this particle if it is NOT a jet particle */
  if (!jet_flag
      && (bi->state == BH_states_adaf && bi->adaf_energy_to_dump > 0.f)) {
    float wi;

    /* Compute the kernel function; note that r cannot be optimised
     * to r2 / sqrtf(r2) because of 1 / 0 behaviour. */
    const float r = sqrtf(r2);
    const float hi_inv = 1.0f / hi;
    const float hi_inv_dim = pow_dimension(hi_inv);
    const float ui = r * hi_inv;
    kernel_eval(ui, &wi);

    /* Compute new energy per unit mass of this particle */
    const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
    /* Below is equivalent to 
     * E_inject_i = E_ADAF * (w_ij * m_i) / Sum(w_ij * mj) */
    const double E_inject 
        = bi->adaf_energy_to_dump * 
          (hydro_get_mass(pj) * hi_inv_dim * wi / bi->rho_gas);

    double u_new = u_init + E_inject / hydro_get_mass(pj); 
    double T_new = u_new / bh_props->temp_to_u_factor;

    /* There can sometimes be VERY large amounts of energy to deposit */
    if (T_new > bh_props->adaf_maximum_temperature) {
      u_new = bh_props->adaf_maximum_temperature * bh_props->temp_to_u_factor;
    }

#ifdef OBSIDIAN_DEBUG_CHECKS
    message("BH_ADAF_HEAT: bid=%lld heating pid=%lld to T=%g K.",
            bi->id, pj->id, 
            u_new / (bh_props->T_K_to_int * bh_props->temp_to_u_factor));
#endif

    /* If no heating, just mark as not swallowed and continue */
    if (u_new > u_init) {

      /* Check whether we are close to the entropy floor. If we are, we
      * classify the gas as cold regardless of temperature.
      * All star forming gas is considered cold.
      * Only shut off cooling for cold ISM gas.
      */

      /* Hydrogen number density (X_H * rho / m_p) [cm^-3] */
      const float n_H_cgs =
          hydro_get_physical_density(pj, cosmo) * bh_props->rho_to_n_cgs;
      const float T_gas_cgs =
          u_init / (bh_props->temp_to_u_factor * bh_props->T_K_to_int);
      const float T_EoS_cgs = entropy_floor_temperature(pj, cosmo, floor_props)
                                  / bh_props->T_K_to_int;
      if ((n_H_cgs > bh_props->adaf_heating_n_H_threshold_cgs &&
            (T_gas_cgs < bh_props->adaf_heating_T_threshold_cgs ||
                T_gas_cgs < T_EoS_cgs * bh_props->fixed_T_above_EoS_factor)) ||
            pj->sf_data.SFR > 0.f) {

        /* Kick with half the energy when not shutting off cooling */
        if (bh_props->adaf_cooling_shutoff_factor <= 0.f) {
          const double random_for_theta = 
              random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
          const double random_for_phi = 
              random_unit_interval(bi->id, ti_current, random_number_BH_feedback);

          const float theta = acosf(2.f * random_for_theta - 1.f);
          const float phi = 2.f * M_PI * random_for_phi;

          dir[0] = sinf(theta) * cosf(phi);
          dir[1] = sinf(theta) * sinf(phi);
          dir[2] = cosf(theta);

          norm = sqrtf(dir[0] * dir[0] +
                       dir[1] * dir[1] +
                       dir[2] * dir[2]);
          const double E_kick = bh_props->adaf_kick_factor * E_inject;
          v_kick = sqrtf(2.f * E_kick / hydro_get_mass(pj));
          if (v_kick > bh_props->adaf_wind_speed) {
            v_kick = bh_props->adaf_wind_speed;
          }

          double E_heat = E_inject - E_kick;
          if (E_heat <= 0.f) E_heat = 0.f;

          u_new = u_init + E_heat / hydro_get_mass(pj); 
          T_new = u_new / bh_props->temp_to_u_factor;

          /* There can sometimes be VERY large amounts of energy to deposit */
          if (T_new > bh_props->adaf_maximum_temperature) {
            u_new = bh_props->adaf_maximum_temperature * bh_props->temp_to_u_factor;
          }

          /* Never cool */
          if (u_new < u_init) u_new = u_init;

  #ifdef OBSIDIAN_DEBUG_CHECKS
          const float pj_vel_norm = sqrtf(
              pj->gpart->v_full[0] * pj->gpart->v_full[0] + 
              pj->gpart->v_full[1] * pj->gpart->v_full[1] + 
              pj->gpart->v_full[2] * pj->gpart->v_full[2]
          );
  #endif

          const float prefactor = v_kick * cosmo->a / norm;

          pj->v[0] += prefactor * dir[0];
          pj->v[1] += prefactor * dir[1];
          pj->v[2] += prefactor * dir[2];

  #ifdef OBSIDIAN_DEBUG_CHECKS
      message("BH_ADAF_KICK: bid=%lld kicking id=%lld, v_kick=%g km/s, v_kick/v_part=%g",
          bi->id, pj->id, v_kick / bh_props->kms_to_internal, v_kick * cosmo->a / pj_vel_norm);
  #endif

          /* Update the signal velocity of the particle based on the velocity kick. */
          hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, v_kick);
        }
        else {
          /* Get particle time-step */
          double dt;
          if (with_cosmology) { 
            const integertime_t ti_step = get_integer_timestep(bi->time_bin);
            const integertime_t ti_begin =
              get_integer_time_begin(ti_current - 1, bi->time_bin);

            dt = cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
          } else {
            dt = get_timestep(bi->time_bin, time_base);
          }

          /* u_init is physical so cs_physical is physical */
          const double cs_physical 
              = gas_soundspeed_from_internal_energy(pj->rho, u_new);

          /* a_factor_sound_speed converts cs_physical to comoving units */
          pj->feedback_data.cooling_shutoff_delay_time = 
              bh_props->adaf_cooling_shutoff_factor 
                * max(cosmo->a_factor_sound_speed * (pj->h / cs_physical),
                      dt); /* BH timestep as a lower limit */
        }

        /* Wind cannot be star forming */
        if (pj->sf_data.SFR > 0.f) {

          /* Record the current time as an indicator of when this particle was last
            star-forming. */
          if (with_cosmology) {
            pj->sf_data.SFR = -cosmo->a;
          } else {
            pj->sf_data.SFR = -time;
          }
        }
      }

      /* We are overwriting the internal energy of the particle */
      hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
      hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

      /* Impose maximal viscosity */
      hydro_diffusive_feedback_reset(pj);

      /* Synchronize the particle on the timeline */
      timestep_sync_part(pj);
    }

    /* IMPORTANT: The particle MUST NOT be swallowed. 
      * We are taking a f_accretion from each particle, and then
      * kicking the rest. We used the swallow marker as a temporary
      * passer in order to remember which particles have been "nibbled"
      * so that we can kick them out.
      */
    black_holes_mark_part_as_not_swallowed(&pj->black_holes_data);

    return;
  }

  if (pj->black_holes_data.swallow_id == bi->id || jet_flag) {
    dir[0] = bi->angular_momentum_gas[0];
    dir[1] = bi->angular_momentum_gas[1];
    dir[2] = bi->angular_momentum_gas[2];

    /* Kick along the angular momentum axis of gas in the kernel */
    norm = sqrtf(dir[0] * dir[0] +
                 dir[1] * dir[1] +
                 dir[2] * dir[2]);

    /* TODO: random_uniform() won't work here?? */
    /*const float dirsign = (random_uniform(-1.0, 1.0) > 0. ? 1.f : -1.f);*/
    const double random_number = 
        random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
    float dirsign = (random_number > 0.5) ? 1.f : -1.f;

    /* Heat the particle and set kinetic kick information if jet particle */
    if (jet_flag) {
      float new_Tj = 0.f;
      /* Use the halo Tvir? */
      if (bh_props->jet_temperature < 0.f) {
        /* TODO: Get the halo Tvir for pj */
        new_Tj = 1.0e8f * bh_props->T_K_to_int; /* internal */
      } else {
        new_Tj = bh_props->jet_temperature; /* internal */
      }

      /* Compute new energy per unit mass of this particle */
      const double u_init = hydro_get_physical_internal_energy(pj, xpj, cosmo);
      double u_new = bh_props->temp_to_u_factor * new_Tj;

#ifdef OBSIDIAN_DEBUG_CHECKS
      message("BH_JET: bid=%lld heating pid=%lld to T=%g K and kicking to v=%g km/s",
        bi->id, pj->id, u_new / (bh_props->T_K_to_int * bh_props->temp_to_u_factor),
        v_kick / bh_props->kms_to_internal);
#endif

      /* Don't decrease the gas temperature if it's already hotter */
      if (u_new > u_init) {
        /* We are overwriting the internal energy of the particle */
        hydro_set_physical_internal_energy(pj, xpj, cosmo, u_new);
        hydro_set_drifted_physical_internal_energy(pj, cosmo, NULL, u_new);

        const double delta_energy = (u_new - u_init) * hydro_get_mass(pj);

        tracers_after_black_holes_feedback(pj, xpj, with_cosmology, cosmo->a,
                                            time, delta_energy);
      }

      if (bh_props->jet_is_isotropic) {
        const double random_for_theta = 
            random_unit_interval(bi->id, ti_current, random_number_BH_feedback);
        const double random_for_phi = 
            random_unit_interval(bi->id, ti_current, random_number_BH_feedback);

        const float theta = acosf(2.f * random_for_theta - 1.f);
        const float phi = 2.f * M_PI * random_for_phi;

        dir[0] = sinf(theta) * cosf(phi);
        dir[1] = sinf(theta) * sinf(phi);
        dir[2] = cosf(theta);

        norm = sqrtf(dir[0] * dir[0] +
                      dir[1] * dir[1] +
                      dir[2] * dir[2]);
        dirsign = 1.f;
      }
    }

#ifdef OBSIDIAN_DEBUG_CHECKS
    const float pj_vel_norm = sqrtf(
        pj->gpart->v_full[0] * pj->gpart->v_full[0] + 
        pj->gpart->v_full[1] * pj->gpart->v_full[1] + 
        pj->gpart->v_full[2] * pj->gpart->v_full[2]
    );
#endif

    const float prefactor = v_kick * cosmo->a * dirsign / norm;

    pj->v[0] += prefactor * dir[0];
    pj->v[1] += prefactor * dir[1];
    pj->v[2] += prefactor * dir[2];

#ifdef OBSIDIAN_DEBUG_CHECKS
    message("BH_KICK: bid=%lld kicking id=%lld, v_kick=%g km/s, v_kick/v_part=%g",
        bi->id, pj->id, v_kick / bh_props->kms_to_internal, v_kick * cosmo->a / pj_vel_norm);
#endif

    /* Set delay time */
    pj->feedback_data.decoupling_delay_time = bh_props->wind_decouple_time_factor * 
        cosmology_get_time_since_big_bang(cosmo, cosmo->a);

    if (jet_flag) {
      pj->feedback_data.number_of_times_decoupled += 100000;
    } else {
      pj->feedback_data.number_of_times_decoupled += 1000;
    }

    /* Update the signal velocity of the particle based on the velocity kick. */
    hydro_set_v_sig_based_on_velocity_kick(pj, cosmo, v_kick);

    /* Wind cannot be star forming */
    if (pj->sf_data.SFR > 0.f) {

      /* Record the current time as an indicator of when this particle was last
        star-forming. */
      if (with_cosmology) {
        pj->sf_data.SFR = -cosmo->a;
      } else {
        pj->sf_data.SFR = -time;
      }

    }

    /* Impose maximal viscosity */
    hydro_diffusive_feedback_reset(pj);
    
    /* Synchronize the particle on the timeline */
    timestep_sync_part(pj);

    /* IMPORTANT: The particle MUST NOT be swallowed. 
      * We are taking a f_accretion from each particle, and then
      * kicking the rest. We used the swallow marker as a temporary
      * passer in order to remember which particles have been "nibbled"
      * so that we can kick them out.
      */
    black_holes_mark_part_as_not_swallowed(&pj->black_holes_data);
  }
}

#endif /* SWIFT_OBSIDIAN_BH_IACT_H */