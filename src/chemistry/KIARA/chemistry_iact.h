/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_KIARA_CHEMISTRY_IACT_H
#define SWIFT_KIARA_CHEMISTRY_IACT_H

/**
 * @brief Sums ambient quantities for the firehose wind model
 *
 * This is called from runner_iact_chemistry, which is called during the density loop
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void firehose_compute_ambient_sym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj) {

  struct chemistry_part_data* chi = &pi->chemistry_data;
  struct chemistry_part_data* chj = &pj->chemistry_data;
  if (chi->rho_ambient < 0.f || chj->rho_ambient < 0.f) return;  // this signifies that firehose model is not on

  /* Do accumulation of ambient quantities */
  if (r2 > 0.f) {
    /* Compute the kernel function for pi */
    const float r = sqrtf(r2);
    const float hi_inv = 1. / hi;
    const float hj_inv = 1. / hj;
    const float ui = r * hi_inv;
    const float uj = r * hj_inv;
    float wi, wj;
    kernel_eval(ui, &wi);
    kernel_eval(uj, &wj);
    /* Accumulate ambient neighbour quantities with an SPH gather operation */
    if (wi > 0.f) {
      chi->u_ambient += pj->u * wi * pj->chemistry_data.volume_factor;
      chi->rho_ambient += pj->mass * wi * hi_inv * hi_inv * hi_inv;
      chj->u_ambient += pi->u * wj * pi->chemistry_data.volume_factor;
      chj->rho_ambient += pi->mass * wj * hj_inv * hj_inv * hj_inv;
      //message("FIREHOSE sum: %lld %lld %g %g %g %g %g %g\n",pi->id, pj->id, sqrt(r2), wi, pj->mass, pj->rho, chi->u_ambient, chi->rho_ambient);
    }
  }
}

/**
 * @brief Sums ambient quantities for the firehose wind model
 *
 * This is called from runner_iact_chemistry, which is called during the density loop
 *
 * @param p The particle to act upon
 * @param cd #chemistry_global_data containing chemistry informations.
 */
__attribute__((always_inline)) INLINE static void firehose_compute_ambient_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj) {


  struct chemistry_part_data* chi = &pi->chemistry_data;
  if (chi->rho_ambient < 0.f) return;  // this signified that firehose model is not on

  /* Do accumulation of ambient quantities */
  if (r2 > 0.f) {
    /* Compute the kernel function for pi */
    const float r = sqrtf(r2);
    const float h_inv = 1. / hi;
    const float ui = r * h_inv;
    float wi;
    kernel_eval(ui, &wi);
    /* Accumulate ambient neighbour quantities with an SPH gather operation */
    if (wi > 0.f) {
      chi->u_ambient += pj->u * wi * pj->chemistry_data.volume_factor;
      chi->rho_ambient += pj->mass * wi * h_inv * h_inv * h_inv;
    }
  }
}

/**
 * @file KIARA/chemistry_iact.h
 * @brief Smooth metal interaction functions following the KIARA model.
 */

/**
 * @brief do chemistry computation after the runner_iact_density (symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, compute ambient quantities for firehose wind diffusion */
    firehose_compute_ambient_sym(r2, dx, hi, hj, pi, pj);
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  float wi, wi_dx;
  float wj, wj_dx;

  /* Get the masses. */
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute the kernel function for pj */
  const float uj = r / hj;
  kernel_deval(uj, &wj, &wj_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  const float wj_dr = wj_dx * r_inv;
  const float mi_wj_dr = mi * wj_dr;

    /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    const float dxi_mi_wj_dr = dx[i] * mi_wj_dr;

    chi->shear_tensor[i][0] += (pj->v[0] - pi->v[0]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][1] += (pj->v[1] - pi->v[1]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][2] += (pj->v[2] - pi->v[2]) * dxi_mj_wi_dr;

    chj->shear_tensor[i][0] -= (pj->v[0] - pi->v[0]) * dxi_mi_wj_dr;
    chj->shear_tensor[i][1] -= (pj->v[1] - pi->v[1]) * dxi_mi_wj_dr;
    chj->shear_tensor[i][2] -= (pj->v[2] - pi->v[2]) * dxi_mi_wj_dr;
  }
}

/**
 * @brief do chemistry computation after the runner_iact_density (non symmetric
 * version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_chemistry(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, const struct part *restrict pj, const float a,
    const float H) {

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, compute ambient quantities for firehose wind diffusion */
    firehose_compute_ambient_nonsym(r2, dx, hi, hj, pi, pj);
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;

  float wi, wi_dx;

  /* Get the masses. */
  const float mj = hydro_get_mass(pj);

  /* Get r */
  const float r = sqrtf(r2);
  const float r_inv = 1.f / r;

  /* Compute the kernel function for pi */
  const float ui = r / hi;
  kernel_deval(ui, &wi, &wi_dx);

  const float wi_dr = wi_dx * r_inv;
  const float mj_wi_dr = mj * wi_dr;

  /* Compute the shear tensor */
  for (int i = 0; i < 3; i++) {
    const float dxi_mj_wi_dr = dx[i] * mj_wi_dr;
    chi->shear_tensor[i][0] += (pj->v[0] - pi->v[0]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][1] += (pj->v[1] - pi->v[1]) * dxi_mj_wi_dr;
    chi->shear_tensor[i][2] += (pj->v[2] - pi->v[2]) * dxi_mj_wi_dr;
  }
}


/**
 * @brief Computes the mass exchanged between the firehose stream and the ambient medium
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle (not updated).
 * @param pj Gas particle.
 * @param xpi Extra particle data (wind)
 * @param xpj Extra particle data (gas)
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static float firehose_compute_mass_exchange(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd, 
    float *v2) {

  /* Are we using firehose model? If not, return */
  if (pi->chemistry_data.radius_stream <= 0.f) return 0.;
  /* Ignore COUPLED particles */
  if (pi->feedback_data.decoupling_delay_time <= 0.f) return 0.;
  /* No wind-wind interaction */
  if (pj->feedback_data.decoupling_delay_time > 0.f) return 0.;

  /* Compute the kernel function */
  const float r = sqrtf(r2);
  const float h_inv = 1. / hi;
  const float ui = r * h_inv;
  float wi;
  kernel_eval(ui, &wi);
  if (wi <= 0.) return 0.;  // Too far out to interact

  /* Get the timestep */ 
  const float dt = get_timestep(pi->time_bin, time_base);
  if (dt <= 0.f) return 0.;  // zero timestep so no mixing

  /* Compute the velocity of the stream relative to ambient gas */
  *v2 = 0.f;
  for (int i=0; i<3; i++){
    *v2 += (pi->v_full[i] - pj->v_full[i]) * (pi->v_full[i] - pj->v_full[i]);
  }

  /* Don't apply above some velocity to avoid jets */
  if (*v2 > cd->firehose_max_velocity) return 0.;

  /* Compute thermal energy ratio for stream and ambient */
  const float chi = pi->chemistry_data.u_ambient / pi->u;
  const float v_stream = sqrtf(*v2);
  const float c_stream = sqrtf(pi->u * hydro_gamma * hydro_gamma_minus_one);
  const float c_amb = sqrtf(pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one);
  const float Mach = v_stream / (c_stream+c_amb);
  const float alpha = 0.21 * (0.8 * exp(-3 * Mach * Mach) + 0.2);

  /* Define timescales for streams*/
  /* Get mixing layer cooling time, which is negative if cooling */
  float tcoolmix = 1.e10 * dt;   // some large number that will result in no mixing
  if (pi->cooling_data.mixing_layer_cool_time < 0.f) tcoolmix = fabs(pi->cooling_data.mixing_layer_cool_time);
  const double tshear = pi->chemistry_data.radius_stream / (alpha * v_stream);
  const float tsc = 2.f * pi->chemistry_data.radius_stream / c_stream;  // sound-crossing time

  //message("FIREHOSE cool: %lld tc=%g ts=%g r=%g td=%g", pi->id, tcoolmix, tshear, pi->chemistry_data.radius_stream, pi->chemistry_data.destruction_time/dt);

  /* Mass change is growth due to cooling minus loss due to shearing, kernel-weighted */
  float dm = 0.f, delta_shear = 0.f, delta_growth = 0.f;
  if (tshear < tcoolmix) delta_shear = (1.f - exp(-dt / tshear));
  if (pi->cooling_data.mixing_layer_cool_time < 0.f) delta_growth = 4. / chi / tsc * pow(tcoolmix / tsc, -0.25) * dt;
  dm = wi * pi->mass * (delta_growth - delta_shear);

  /* Limit amount of mixing per neighbor */
  float fmix_max = 0.02;  // with ~50 neighbours, this limits total loss/gain to a particle's mass in a single step.
  if (dm > fmix_max * pi->mass) dm = fmix_max * pi->mass;
  if (dm < -fmix_max * pi->mass) dm = -fmix_max * pi->mass;

  /* Update stream radius */
  float stream_growth_factor = (pi->mass + dm) / pi->mass;
  if (stream_growth_factor > 0.f) pi->chemistry_data.radius_stream *= sqrtf(stream_growth_factor);
  else pi->chemistry_data.radius_stream = 0.;

  /* If stream is growing, don't mix */
  if (dm > 0.f) dm = 0.f;

  if (dm < 0.f) message("FIREHOSE: %lld %lld m=%g rhoi=%g rhoamb=%g rhoj=%g ui=%g uamb=%g uj=%g Ti/Tamb=%g grow=%g shear=%g tshear=%g tcmix=%g tcool=%g fexch=%g", pi->id, pj->id, pi->mass, pi->rho, pi->chemistry_data.rho_ambient, pj->rho, pi->u, pi->chemistry_data.u_ambient, pj->u, pi->u/pi->chemistry_data.u_ambient, delta_growth, delta_shear, tshear, tcoolmix/tshear, pi->cooling_data.mixing_layer_cool_time, dm/pi->mass);

  return dm;
}


/**
 * @brief Check recoupling criterion for firehose stream particle 
 *
 * @param pi Wind particle (not updated).
 * @param Mach Stream Mach number vs ambient
 * @param r_stream Current radius of stream 
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static float firehose_recoupling_criterion(
	struct part *pi, const float Mach, const float r_stream, const struct chemistry_global_data* cd) {

  // negative value indicates it should recouple
  float u_diff = fabs(pi->u - pi->chemistry_data.u_ambient) / max(pi->u, pi->chemistry_data.u_ambient);
  //message("FIREHOSE recouple: %lld M=%g u=%g r=%g dm=%g tdel=%g", pi->id, Mach, u_diff, r_stream, pi->chemistry_data.exchanged_mass/pi->mass, pi->feedback_data.decoupling_delay_time);
  if (Mach < cd->firehose_recoupling_mach && u_diff < cd->firehose_recoupling_u_factor ) return -1.;  
  if (pi->chemistry_data.exchanged_mass / pi->mass > cd->firehose_recoupling_fmix) return -1.;  
  if (r_stream == 0.f) return -1.;
  return pi->chemistry_data.radius_stream;
}


/**
 * @brief Computes the particle interaction via the firehose stream model
 *
 * This is called from runner_iact_diffusion, which is called during the force loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle (not updated).
 * @param pj Gas particle.
 * @param xpi Extra particle data (wind)
 * @param xpj Extra particle data (gas)
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_particle_sym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd) {

  if (r2 <= 0.f) return;

  /* Compute the amount of mass mixed between stream particle and ambient gas */
  float v2=0.f;
  const float dm = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, time_base, ti_current, phys_const, cd, &v2);
  float delta_m = fabs(dm);
  /* This is the fraction of each quantity which particle j will exchange, from among all ambient gas */
  if (delta_m <= 0.) return; // no interaction, nothing to do
  pi->chemistry_data.exchanged_mass += delta_m;  // track amount of gas mixed

  /* Constants */ 
  const float Mach = sqrtf(v2 / (pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* set weights for averaging i and j */
  float pii_weight = (pi->mass - delta_m) / pi->mass;
  float pij_weight = delta_m / pi->mass;
  float pji_weight = delta_m / pj->mass;
  float pjj_weight = (pj->mass - delta_m) / pj->mass;
  if (pij_weight < 1.e-10 && pji_weight < 1.e-10) return;  // mixing is negligibly small, avoid underflows

  /* 1) Update chemistry */
  pi->chemistry_data.metal_mass_fraction_total = 0.f;
  pj->chemistry_data.metal_mass_fraction_total = 0.f;
  float pi_metal_frac;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    pi_metal_frac = pi->chemistry_data.metal_mass_fraction[elem];
    pi->chemistry_data.metal_mass_fraction[elem] = (pii_weight * pi->mass * pi->chemistry_data.metal_mass_fraction[elem] + pij_weight * pj->mass * pj->chemistry_data.metal_mass_fraction[elem]) / pi->mass;
    pj->chemistry_data.metal_mass_fraction[elem] = (pjj_weight * pj->mass * pj->chemistry_data.metal_mass_fraction[elem] + pji_weight * pi->mass * pi_metal_frac) / pj->mass;
    if (elem > chemistry_element_He) {
      pi->chemistry_data.metal_mass_fraction_total += pi->chemistry_data.metal_mass_fraction[elem];
      pj->chemistry_data.metal_mass_fraction_total += pj->chemistry_data.metal_mass_fraction[elem];
    }
  }

  int spread_dust=1; // flag (for testing) whether to spread dust in addition to metals
//  message("FIREHOSE_1dust: %lld %g %g %g %g %g %g", pi->id, pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2], delta_m, dm);
  if (spread_dust && (pi->cooling_data.dust_mass+pj->cooling_data.dust_mass) > 0.f) {
    /* Spread dust mass between particles */
    float pi_dust_mass = pi->cooling_data.dust_mass;
    float pj_dust_mass = pj->cooling_data.dust_mass;
    pi->cooling_data.dust_mass = pii_weight * pi_dust_mass + pij_weight * pj_dust_mass;
    pj->cooling_data.dust_mass = pji_weight * pi_dust_mass + pjj_weight * pj_dust_mass;
    /* Spread individual dust elements */
    float pi_dust_frac;
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      pi_dust_frac = pi->cooling_data.dust_mass_fraction[elem];
      if (pi->cooling_data.dust_mass > 0.f) pi->cooling_data.dust_mass_fraction[elem] = (pii_weight * pi_dust_mass * pi->cooling_data.dust_mass_fraction[elem] + pij_weight * pj_dust_mass * pj->cooling_data.dust_mass_fraction[elem]) / pi->cooling_data.dust_mass;
      else pi->cooling_data.dust_mass_fraction[elem] = 0.f;
      if (pj->cooling_data.dust_mass > 0.f) pj->cooling_data.dust_mass_fraction[elem] = (pjj_weight * pj_dust_mass * pj->cooling_data.dust_mass_fraction[elem] + pji_weight * pi_dust_mass * pi_dust_frac) / pj->cooling_data.dust_mass;
      else pj->cooling_data.dust_mass_fraction[elem] = 0.f;
    }
  }
//  message("FIREHOSE_2dust: %lld %g %g %g %g %g %g", pi->id, pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2], delta_m, dm);

  /* 2) Update particles' internal energy per unit mass */
  float pi_u = pi->u;
  //float pj_u = pj->u;
  pi->u = (pii_weight * pi->mass * pi_u + pij_weight * pj->mass * pj->u) / pi->mass;
  pj->u = (pji_weight * pi->mass * pi_u + pjj_weight * pj->mass * pj->u) / pj->mass;
  //message("FIREHOSE_u: %lld %g %g %g %g %g", pi->id, wi, pi_u, pi->u, pj_u, pj->u);


    /* 3) Update particles' velocities, conserving momentum */
  float pi_vfull, new_v2 = 0.f;
  for (int i=0; i<3; i++){
    pi_vfull = pi->v_full[i];
    pi->v_full[i] = (pii_weight * pi->mass * pi_vfull + pij_weight * pj->mass * pj->v_full[i]) / pi->mass;
    pj->v_full[i] = (pji_weight * pi->mass * pi_vfull + pjj_weight * pj->mass * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i]-pj->v_full[i]) * (pi->v_full[i]-pj->v_full[i]);
  }

   /* 4) Deposit excess energy onto stream */
  float delE = 0.5f * delta_m * (v2 - new_v2);
  pi->u += delE / pi->mass;

  /* Check if particle should recouple */
  pi->chemistry_data.radius_stream = firehose_recoupling_criterion(pi, Mach, pi->chemistry_data.radius_stream, cd);  // Negative value signifies particle should be recoupled (which happens in feedback.h)

  return;
}

/**
 * @brief Computes non-symmetric (single) particle interaction via the firehose stream model
 *
 * This is called from runner_iact_nonsym_diffusion, which is called during the force loop
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi Wind particle (not updated).
 * @param pj Gas particle.
 * @param xpi Extra particle data (wind)
 * @param xpj Extra particle data (gas)
 * @param ti_current Current integer time used value for seeding random number generator   
 * 
 * @param phys_const Physical constants
 * @param us Unit system
 *
 */
__attribute__((always_inline)) INLINE static void firehose_evolve_particle_nonsym(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *pi, struct part *pj,
    const float time_base, const integertime_t ti_current,
    const struct phys_const* phys_const, const struct chemistry_global_data* cd) {

  if (r2 <= 0.f) return;

  /* Compute the interaction terms between stream particle and ambient gas */
  float v2=0.f;
  float dm = firehose_compute_mass_exchange(r2, dx, hi, hj, pi, pj, time_base, ti_current, phys_const, cd, &v2);
  float delta_m = fabs(dm);
  if (delta_m <= 0.) return; // no interaction, nothing to do
  pi->chemistry_data.exchanged_mass += delta_m;  // track amount of gas mixed

  /* Constants */ 
  const float Mach = sqrtf(v2/(pi->chemistry_data.u_ambient * hydro_gamma * hydro_gamma_minus_one));  

  /* set weights for averaging i and j */
  float pii_weight = (pi->mass - delta_m) / pi->mass;
  float pij_weight = delta_m / pi->mass;
  float pji_weight = delta_m / pj->mass;
  float pjj_weight = (pj->mass - delta_m) / pj->mass;
  if (pij_weight < 1.e-10 && pji_weight < 1.e-10) return;  // mixing is negligibly small

  /* 1) Update chemistry */
  pi->chemistry_data.metal_mass_fraction_total = 0.f;
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    pi->chemistry_data.metal_mass_fraction[elem] = (pii_weight * pi->mass * pi->chemistry_data.metal_mass_fraction[elem] + pij_weight * pj->mass * pj->chemistry_data.metal_mass_fraction[elem]) / pi->mass;
    //message("FIREHOSE met: %d %g %g %g %g %g",elem,pi->chemistry_data.metal_mass_fraction[elem],pii_weight,pij_weight,pj->chemistry_data.metal_mass_fraction[elem],pi->mass);
    if (elem > chemistry_element_He) {
      pi->chemistry_data.metal_mass_fraction_total += pi->chemistry_data.metal_mass_fraction[elem];
    }
  }

  int spread_dust=1;
//  message("FIREHOSE_1dust: %lld %g %g %g %g %g %g", pi->id, pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2], delta_m, dm);
  if (spread_dust && (pi->cooling_data.dust_mass+pj->cooling_data.dust_mass) > 0.f) {
    /* Spread dust mass between particles */
    float pi_dust_mass = pi->cooling_data.dust_mass;
    float pj_dust_mass = pj->cooling_data.dust_mass;
    pi->cooling_data.dust_mass = pii_weight * pi_dust_mass + pij_weight * pj_dust_mass;
    /* Spread individual dust elements */
    for (int elem = 0; elem < chemistry_element_count; ++elem) {
      pi->cooling_data.dust_mass_fraction[elem] = (pii_weight * pi_dust_mass * pi->cooling_data.dust_mass_fraction[elem] + pij_weight * pj_dust_mass * pj->cooling_data.dust_mass_fraction[elem]) / pi->cooling_data.dust_mass;
    }
  }
//  message("FIREHOSE_2dust: %lld %g %g %g %g %g %g", pi->id, pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2], delta_m, dm);
  //message("FIREHOSE after: %g %g %g %g", pi->cooling_data.dust_mass, pj->cooling_data.dust_mass, pi->cooling_data.dust_mass_fraction[2], pj->cooling_data.dust_mass_fraction[2]);

  /* 2) Update particles' internal energy per unit mass */
  float pi_u = pi->u;
  pi->u = (pii_weight * pi->mass * pi_u + pij_weight * pj->mass * pj->u) / pi->mass;


    /* 3) Update particles' velocities, conserving momentum */
  float pi_vfull, pj_vfull, new_v2 = 0.f;
  for (int i=0; i<3; i++){
    pi_vfull = pi->v_full[i];
    pj_vfull = pj->v_full[i];
    pi->v_full[i] = (pii_weight * pi->mass * pi_vfull + pij_weight * pj->mass * pj->v_full[i]) / pi->mass;
    pj->v_full[i] = (pji_weight * pi->mass * pi_vfull + pjj_weight * pj->mass * pj->v_full[i]) / pj->mass;
    new_v2 += (pi->v_full[i]-pj->v_full[i]) * (pi->v_full[i]-pj->v_full[i]);
    //message("FIREHOSE v2: %d %g %g %g %g %g %g %g %g %g %g %g %g %g",i,pi_vfull,pj_vfull,pi->v_full[i],pj->v_full[i],v2,new_v2,pii_weight,pij_weight,pji_weight,pjj_weight,delta_m/pi->mass,dm,pow(pi->h,3.));
    pj->v_full[i] = pj_vfull;
  }

   /* 4) Deposit excess energy into ambient */
  float delE = 0.5f * delta_m * (v2 - new_v2);
  //if (delE<0.f) error("delE<0! %g %g %g %g\n",delta_m,v2,new_v2,delE);
  pi->u += delE / pi->mass;

  /* Check if particle should recouple */
  pi->chemistry_data.radius_stream = firehose_recoupling_criterion(pi, Mach, pi->chemistry_data.radius_stream, cd);  // Negative value signifies particle should be recoupled (which happens in feedback.h)

  return;
}


/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (symmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology, 
    const struct phys_const* phys_const, const struct chemistry_global_data *cd) {

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, do firehose wind diffusion */
    firehose_evolve_particle_sym(r2, dx, hi, hj, pi, pj, time_base, t_current, phys_const, cd);
    return;
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  /* No need to diffuse if both particles are not diffusing. */
  if (chj->diffusion_coefficient > 0 && chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float mi = hydro_get_mass(pi);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, wj, dwi_dx, dwj_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part j */
    /* Get the kernel for hj */
    const float hj_inv = 1.0f / hj;

    /* Compute the kernel function for pj */
    const float xj = r * hj_inv;
    kernel_deval(xj, &wj, &dwj_dx);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);

    const float wi_dr = dwi_dx * r_inv;
    const float wj_dr = dwj_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;
    const float mi_dw_r = mi * wj_dr;

    const float rhoij_inv = 1.f / (rhoi * rhoj);

    /**
     * Compute the diffusion following Eq. 2.14
     * from Monaghan, Huppert, & Worster (2006).
     */
    float coef = 4.f * chi->diffusion_coefficient * chj->diffusion_coefficient;
    coef /= chi->diffusion_coefficient + chj->diffusion_coefficient;

    const float coef_i = coef * mj_dw_r * rhoij_inv;
    const float coef_j = coef * mi_dw_r * rhoij_inv;

    /* Compute the time derivative of metals due to diffusion */
    const float dZ_ij_tot = chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;
    chj->dZ_dt_total -= coef_j * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
      chj->dZ_dt[elem] -= coef_j * dZ_ij;
    }
  }
}

/**
 * @brief do metal diffusion computation in the <FORCE LOOP>
 * (nonsymmetric version)
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param time_base The time base used in order to convert integer to float
 * time.
 * @param ti_current The current time (in integer)
 * @param cosmo The #cosmology.
 * @param with_cosmology Are we running with cosmology?
 *
 */
__attribute__((always_inline)) INLINE static void runner_iact_nonsym_diffusion(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const float time_base, const integertime_t t_current,
    const struct cosmology *cosmo, const int with_cosmology,
    const struct phys_const* phys_const, const struct chemistry_global_data *cd) {

  if (pi->feedback_data.decoupling_delay_time > 0.f) {
    /* If in wind mode, do firehose wind diffusion */
    firehose_evolve_particle_nonsym(r2, dx, hi, hj, pi, pj, time_base, t_current, phys_const, cd);
    return;
  }

  struct chemistry_part_data *chi = &pi->chemistry_data;
  struct chemistry_part_data *chj = &pj->chemistry_data;

  if (chj->diffusion_coefficient > 0 && chi->diffusion_coefficient > 0) {

    /* Get mass */
    const float mj = hydro_get_mass(pj);
    const float rhoj = hydro_get_physical_density(pj, cosmo);
    const float rhoi = hydro_get_physical_density(pi, cosmo);

    float wi, dwi_dx;

    /* Get r */
    const float r = sqrtf(r2);

    /* part i */
    /* Get the kernel for hi */
    const float hi_inv = 1.0f / hi;

    /* Compute the kernel function for pi */
    const float xi = r * hi_inv;
    kernel_deval(xi, &wi, &dwi_dx);

    /* Get 1/r */
    const float r_inv = 1.f / sqrtf(r2);
    const float wi_dr = dwi_dx * r_inv;

    const float mj_dw_r = mj * wi_dr;

    const float rhoij_inv = 1.f / (rhoi * rhoj);

    /**
     * Compute the diffusion following Eq. 2.14
     * from Monaghan, Huppert, & Worster (2006).
     */
    float coef = 4.f * chi->diffusion_coefficient * chj->diffusion_coefficient;
    coef /= chi->diffusion_coefficient + chj->diffusion_coefficient;

    const float coef_i = coef * mj_dw_r * rhoij_inv;

    /* Compute the time derivative */
    const float dZ_ij_tot = chi->metal_mass_fraction_total - chj->metal_mass_fraction_total;
    chi->dZ_dt_total += coef_i * dZ_ij_tot;

    for (int elem = 0; elem < chemistry_element_count; elem++) {
      const float dZ_ij = chi->metal_mass_fraction[elem] - chj->metal_mass_fraction[elem];
      chi->dZ_dt[elem] += coef_i * dZ_ij;
    }
  }
}


#endif /* SWIFT_KIARA_CHEMISTRY_IACT_H */