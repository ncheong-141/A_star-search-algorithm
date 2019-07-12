/* Implementation for A* functions */

#include "A_star_functions.h"

/* Grid generation function (Cartesian uniform grid) */
void grid_generation(RMF::DYN_C2D<Cell_ds>& grid, size_t row_PT, size_t col_PT, float spacing) {
	
	
	// Pre-allocate memory and construct Cell_ds entities. Note, DYN_C2D::emplace_back(object, rows) reserves memory in function, however, with
	//															 DYN_C2D::emplace_back(objects(parameters)) this has to be done manually. 

	grid.reserve(row_PT, col_PT);	// Reserve memory to avoid unneccesary copying.
									// Furthermore, emplace_back also constructs Cell_ds entity in the vector to further reduce copying

	// Predefine counters 
	int ID = 0;
	float ypos = 0.0f;
	float xpos = 0.0f;

	// Predefine collumn and row end points (so only need to calculate once 
	size_t row_end_pt = row_PT - 1; 
	size_t col_end_pt = col_PT - 1; 

	for (size_t i = 0; i < row_PT; i++) {
		for (size_t j = 0; j < col_PT; j++) {

			grid.element_emplace_back(Cell_ds(ID, i, j, col_PT, xpos, ypos), j);

			ID += 1;				// Update ID; 
			xpos += spacing;		// Update X position

			/* Setting boundary points*/
			/* B_IDs :	BL corner	= 1
						B edge		= 2
						BR corner	= 3
						L edge		= 4 
						R edge		= 5 
						TL corner	= 6
						T edge		= 7 
						TR corner	= 8
			*/
			if (i == 0 && j == 0) {										// Bottom left corner point
				grid(i, j).set_Boundary_pt(); 
				grid(i, j).set_Boundary_pt_val(1);
			}
			else if (i == 0 && j < col_end_pt) {						// Bottom edge without BR corner
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(2);
			}
			else if (i == 0 && j == col_end_pt) {						// Bottom right corner
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(3);
			}
			else if (i > 0 && i < row_end_pt && j == 0) {				// Left edge
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(4);
			}
			else if (i > 0 && i < row_end_pt && j == col_end_pt) {		// Right edge
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(5);
			}
			else if (i == row_end_pt && j == 0) {						// Top left corner
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(6);
			}
			else if (i == row_end_pt && j > 0 && j < col_end_pt) {		// Top edge
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(7);
			}
			else if (i == row_end_pt && j == col_end_pt) {				// Top right edge
				grid(i, j).set_Boundary_pt();
				grid(i, j).set_Boundary_pt_val(8);
			}


		}

		xpos = 0.0f;				// Re-zero xpos as starting from next row 
		ypos += spacing;			// Update ypos 
	}

	/* Debug print grid */
	#if DEBUG==0
	for (int i_plot = 0; i_plot < row_PT; i_plot++) {
		std::cout << "| ";
		for (int j_plot = 0; j_plot < col_PT; j_plot++) {
			std::cout << grid(i_plot, j_plot).get_Boundary_pt_val() << " ";
		}
		std::cout << "|\n";
	}
	#endif
}

/* Obstacle generation */
void obstacle_generation(RMF::DYN_C2D<Cell_ds>& grid, int row_PT, int col_PT, int obstacle_setting) {

	switch (obstacle_setting) {
		case (1):			// Uniform cells as objects 
		{
			for (int i = 0; i < row_PT; i++) {
				for (int j = 0; j < col_PT; j++) {
					if (RMF::isEven(i)) {
						if (RMF::isOdd(j)) {
							// Do not make the end point an obstacle!! 
							if (grid(i, j).get_isEnd() == false) {
								grid(i, j).set_Obstacle();
							}
						}
					}
				}
			}
			break;
		}
		case (2):			// Diagonal line through middle 
		{
			int diag_counter = 0; 
			for (int i = 0; i < row_PT; i++) {
				for (int j = 2; j < col_PT; j++) {
					if (j == diag_counter || j == diag_counter + 1) {
						grid(i, j).set_Obstacle();
					}
				}
				diag_counter += 1;
			}
			break;
		}
	}
	// Loop over each grid cell and set 
}
/* Choosing correct parent cell based on F and H costs */
std::vector<Cell_ds> min_FHcost(std::vector<Cell_ds>& priority_list, Cell_ds*& parent_pt, Cell_ds& pt_A) {

	// This is initial conditions and depending on this, the algorithm can get stuck
	// (i.e. if you use the parent point it gets stuck if it cant find a better g and h cost)
	int min_fcost = 2*priority_list[0].get_Fcost();	// Set min f_cost
	//int min_fcost = pt_A.get_Fcost();

	// Initialize a cell in a vector with  minimum fcost for processing (vector as I also may want to return multiple cells in future for optomizing the algortithm)
	std::vector<Cell_ds> fcost_min_cell;			
	fcost_min_cell.reserve(1);						// Reserve memory for avoiding initial de-alloc and alloc (copying) 
	fcost_min_cell.push_back(priority_list[0]);		// Initialize cell with any point

	// Predefine size so not constantly calling function when checking loop condition
	size_t priolist_size = priority_list.size();

	// Loop over priority list (Could use iterator instead as its slightly faster)
	for (size_t j = 0; j < priolist_size; j++) {

		// Since current parent is on priority list, DO NOT, consider it again 
		if (priority_list[j].get_Parent_eval() == false) {

			if (priority_list[j].get_Fcost() < min_fcost) {

				min_fcost = priority_list[j].get_Fcost();

				// Reset the fcost_min list as a new minimum is found 
				if (fcost_min_cell.empty() == false) {		// if the list is not empty
					fcost_min_cell.clear();					// Clears all elements and returns the vector to zero size. Note, capacity is NOT zeroed! 
				}
				fcost_min_cell.push_back(priority_list[j]);	// Insert the cell into the fcost min vector

			}
			// If there is more than one cell on priority list which has the minimum fcost value
			else if (priority_list[j].get_Fcost() == min_fcost) {

				// Since fcosts are the same; compare Hcosts and input the cell with the smallest gcost! 
				if (priority_list[j].get_Hcost() < fcost_min_cell[0].get_Hcost()) {
					fcost_min_cell.clear();						// Clears all elements and returns the vector to zero size. Note, capacity is NOT zeroed! 
					fcost_min_cell.push_back(priority_list[j]);	// Insert the cell into list of cells with minimum g and f cost. 
				}
				// If it occurs that the Hcosts are also the same. Just use any of them... 
				else if (priority_list[j].get_Hcost() == fcost_min_cell[0].get_Hcost()) {
					continue;
				}
				else {
					std::cout << "\nNo better cells Fcost and Gcost-wise\n";
				}
			}
		}

	}


	// Return the minimum f and hcost cell to become the next parent!
	return fcost_min_cell; 
}

/* A* algorithm implementation */
void A_star_algorithm(size_t(&loop_index)[4], RMF::DYN_C2D<Cell_ds>& grid, std::vector<Cell_ds>& priority_list, size_t& priolist_index, Cell_ds*& parent_pt, Cell_ds& pt_B, const size_t& p_col, const size_t& p_row, const int& diag_movement, const int& adj_movement) {

	// Loop for upper and lower limits 
	for (size_t i = loop_index[0]; i <= loop_index[1]; i = i + 1) {
		
		for (size_t j = loop_index[2]; j <= loop_index[3]; j++) {

			// Condition if the cell has been a parent before (don't evaluate) and don't evaulate the actual current parent cell 
			if (grid(i, j).get_Parent_eval() == true || (j == p_col && i == p_row) || grid(i, j).get_Obstacle() == true) {
				continue; // Dont evaluate statement which goes to next iteration of for loop 
						  // In terms of obstacles, not evaluating means it will not be put on the priolist and therefore not considered in the path
			}
			// Check if it has been a child cell before --> compare the movement costs (gcosts)
			else if (grid(i, j).get_On_priolist() == true) {

				// Predefine the index which the cell is allocated on the priority list
				size_t priolist_cell_index = grid(i, j).get_Priolist_index();

				// Conditions for adjacent movement (Combinations of indices)
				if ((i == p_row - 1 && j == p_col) || (i == p_row + 1 && j == p_col) || (i == p_row && j == p_col - 1) || (i == p_row && j == p_col + 1)) {

					// Select adjacent movement
					int movement = adj_movement;

					// Calculate the potential new gcost 
					int potential_new_gcost = grid(i, j).generate_gcost(parent_pt, movement);

					// Compare the movement costs and select the lowest one 
					if (potential_new_gcost < priority_list[priolist_cell_index].get_Gcost()) {

						// Update the cell data 
						grid(i, j).set_Gcost(potential_new_gcost);

						// Let the currently evaluated cell know its new parent
						grid(i, j).set_From_Cell_ID(parent_pt->get_ID());

						// Update the priority list 
						priority_list[priolist_cell_index].set_Gcost(potential_new_gcost);
						priority_list[priolist_cell_index].set_From_Cell_ID(parent_pt->get_ID());
					}
				}
				else {	// else, diagonal movement
					int movement = diag_movement;
					int potential_new_gcost = grid(i, j).generate_gcost(parent_pt, movement);

					if (potential_new_gcost < priority_list[priolist_cell_index].get_Gcost()) {
						grid(i, j).set_Gcost(potential_new_gcost);
						grid(i, j).set_From_Cell_ID(parent_pt->get_ID());
						priority_list[priolist_cell_index].set_Gcost(potential_new_gcost);
						priority_list[priolist_cell_index].set_From_Cell_ID(parent_pt->get_ID());
					}
				}
			}
			// The cell is new and not been considered before
			else {

				 // Generate all the costs 
				 // Condition adjacent movement 
				if ((i == p_row - 1 && j == p_col) || (i == p_row + 1 && j == p_col) || (i == p_row && j == p_col - 1) || (i == p_row && j == p_col + 1)) {
					// Select adjacent movement
					int movement = adj_movement;

					// Calculate the potential new gcost 
					grid(i, j).set_Gcost(grid(i, j).generate_gcost(parent_pt, movement));
					grid(i, j).generate_heurcost(pt_B, diag_movement, adj_movement);
					grid(i, j).generate_fcost();

					// Let the currently evaluated cell know which parent it is from
					grid(i, j).set_From_Cell_ID(parent_pt->get_ID());
				}
				else {	// else, diagonal movement
					int movement = diag_movement;
					grid(i, j).set_Gcost(grid(i, j).generate_gcost(parent_pt, movement));
					grid(i, j).generate_heurcost(pt_B, diag_movement, adj_movement);
					grid(i, j).generate_fcost();
					grid(i, j).set_From_Cell_ID(parent_pt->get_ID());
				}

				// Put cell on priority list, assign it an index on the priority list for future reference and increment the priority list index
				grid(i, j).set_On_priolist();
				grid(i, j).set_Priolist_index(priolist_index);
				priolist_index += 1;
				priority_list.push_back(grid(i, j));
			}
		} // End of jth loop of neighbouring cells
	} // End of ith loop of neighbouring cellsuring cells
}