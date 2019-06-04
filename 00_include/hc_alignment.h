/*
This source code is part of specs, an application for
protein residue specialization and conservation scoring.
Written by Ivana Mihalek.
Copyright (C) 2010-2013 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/
# define ALMT_NAME_LENGTH 50

typedef struct{
    char group_name[MEDSTRING];
    int  no_members;
    int  * group_member;
    int  *gaps; /* number of gaps at each position, in this group */
    int  *unks; /* number of unknowns at each position, in this group */
} Group;

typedef struct{
    int number_of_seqs;
    int length;
    char ** sequence;
    char ** name;
    int  *  seq_gaps;
    int  *  column_gaps;
    int  *  seq_unks;
    int  *  column_unks;
    int  * sunk; /* positions "sunk" to the bottom bcs of  too many gaps */
    int number_of_sunk_positions;
    double **seq_dist;
    int  **aligned_sites;
    int  **identical_sites;
    int  **similar_sites;
    char * struct_seq;
    int  * tp_flag; /* positions in the lmt tagged for some reason as "true positive" */

    int no_groups;
    Group * group;
    int   * belongs_to_group;

    double *** trust_score;

    double ** per_node_spec_score;
    int    ** per_node_pool_size;

}  Alignment;
