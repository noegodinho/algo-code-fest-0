/* init/default.c
 *
 * (C) 2021 Andreia P. Guerreiro <andreia.guerreiro@tecnico.ulisboa.pt>
 *          Carlos M. Fonseca <cmfonsec@dei.uc.pt>
 *          Samuel B. Outeiro <souteiro@student.dei.uc.pt>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License, version 3, as
 * published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include "init.h"
#include "problem.h"
#include <stdio.h>
#include <stdlib.h>

/*
 * Allocate a problem structure and initialise it with problem-specific data
 * contained in a specified input file.
 */
struct problem *initProblem(int argc, char **argv)
{
    struct problem *p = NULL;

    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file name>\n", argv[0]);
        return NULL;
    }

    p = newProblem(argv[1]);

    return p;
}
