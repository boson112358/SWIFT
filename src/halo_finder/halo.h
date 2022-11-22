/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2022 Will Roper (w.roper@sussex.ac.uk)
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
#ifndef SWIFT_PHASE_CELL_H
#define SWIFT_PHASE_CELL_H

/* Config parameters. */
#include <config.h>

/* Includes. */
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/* Local includes. */
#include "cell.h"
#include "fof.h"

/**
 * @brief What kind of halo search are we doing?
 *
 * 0 = A 3D FOF group.
 * 1 = A 6D Host halo.
 * 2 = A 6D Subhalo
 */
enum halo_types {fof_group, host_halo, sub_halo};


/* Prototypes */
void halo_finder_search_self_cell_gpart(const struct fof_props *props,
                                        const double l_x2,
                                        const int halo_level,
                                        const struct gpart *const space_gparts,
                                        struct cell *c);
void halo_finder_search_pair_cells_gpart(const struct fof_props *props,
                                         const double dim[3],
                                         const double l_x2,
                                         const int halo_level,
                                         const int periodic,
                                         const struct gpart *const space_gparts,
                                         struct cell *restrict ci,
                                         struct cell *restrict cj);

#endif /* SWIFT_FOF_STRUCT_H */