#pragma once

/* Functions to execute in A* search algorithm */

/* Include external modules */
#include <iostream>
#include <vector>

/* Include library API */
#include "SL_Comptime_Interface.h"
#include "SL_Runtime_Interface.h"

/* Include any internal APIs*/
#include "Cell_ds.h"

/* Macro for debugging
- DEBUG 0 = Debug -> Real time printing of pathfinding, boundary condition cell output,
- DEBUG 1 = No Debugging. 
*/
#define DEBUG 1

/* Grid generation function*/
void grid_generation(RMF::DYN_C2D<Cell_ds>& grid, size_t row_PT, size_t col_PT, float spacing); 

/* Obstacle generation */
void obstacle_generation(RMF::DYN_C2D<Cell_ds>& grid, int row_PT, int col_PT, int no_of_objects);

/* Choosing correct parent cell based on F and G costs */
std::vector<Cell_ds> min_FHcost(std::vector<Cell_ds>& priority_list, Cell_ds*& parent_pt, Cell_ds& pt_A);

/* A* algorithm implementation */
void A_star_algorithm(size_t(&loop_index)[4], RMF::DYN_C2D<Cell_ds>& grid, std::vector<Cell_ds>& priority_list, size_t& priolist_index, Cell_ds*& parent_pt, Cell_ds& pt_B, const size_t& p_col, const size_t& p_row, const int& diag_movement, const int& adj_movement);
