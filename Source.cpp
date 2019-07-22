/* A* search algorithm for pathfinding */
/* Finds the shortest paths from point A to B past different obstacles. */

/* How it works: 

- At a evals node (point A for start), the algorithm looks at all nodes surrounding it and calculations
	-> gcost () =	distance from starting node (point A)
					(typical for abstract A* algorithms is to have distance between nodes as 1, and multiply by 10 for irradicating FP numbers)
	-> hcost () =  heuristic function / distance from end node (point B). 
	-> fcost()  = total cost -> gcost() + hcost() 
- The algorithm finds the lowest fcosts of the surrounding nodes and then choses that to be the next parent node and repeat the above. 
- If there is obstacles points, they are not evaluated. 
- If there is two points with the same fcost, the smallest hcost takes preference. 
- All fcosts are stored and the losest cost is stored each time
- Do not consider neighbouring nodes which have already been evaluated which have a high fcost() 
- A ndoe that has been a parent node is not evaluated again as a parent node
- The path it took to get to the parent node has to be stored as this will hold the optimaoal path. 
 
 - Cells and lsits to determine best path
 - Each Cell has a h, g and f value and a parent node
 - Lists: Open and closed list. 
 - Open lists -> a list which holds the cells which need to be checked, as we check nodes we put them on the closed list. 
 
 h cost: 
 - h value is the distance from the cell to the end point (B)
 - Pre-compute all the h values if all the cells are not moving (static). 
 - You can use the manhatten formula which is the summation of the adjacent and opposite sides

 g cost: 
 - g value is the movement cost, which is the distance between the cell to another cell (diag is 1.4, adjacent is 1 if the spacing is 1) 
 - The total g value of a cell is the value of the parent cell plus the movement cost to the cell. 
 - Need a condition if the gcost is greater or less than previous parent
*/

/* Include external modules */
#include <iostream>
#include <vector>

/* Include library API */
#include "SL_Comptime_Interface.h"
#include "SL_Runtime_Interface.h"

/* Include any internal APIs*/
#include "Cell_ds.h"
#include "A_star_functions.h"

/* Notes:
Steps:
1.) Discretize a 2D cartesian grid (origin bottom left corner) 
2.) Establish obstacles, point A and B.
3.) create two lists, open and closed. (or model each node with a bool (true for eval it, false for dont) 
	- Open - Nodes which can be evaluated as a parent node
	- Closed - Nodes which have already been evalled. 
4.) will need a priority list to choose the lowest fcost -> lowest hcost to be the next parent node	
5.) Set current node to start node
In main loop: 
5.) If current not eval == true 
	-> if current node is poitn B, stop and generate path. 
	-> Loop over each neighbouring cell
	In neighbour cell loop 
		-> if neighbour cell != obstacle && eval != false  
			-> Evalute, else dont
		-> Calculate fcost
		-> if new path to neighbour is shorter than old OR(and?) neighbour is not in OPEN 
			-> set fcosts
	
6.) Find smallest fcost which has not been a parent node -> if same fcosts check hcosts -> set as new current and loop again 
*/


/* Start of A* search algorithm code */
int main() {

	/* Grid settings */
	const int	row_PT = 10;			// Number of row grid points
	const int	col_PT = 10;			// Number of col  grid points 
	const float	spacing = 1.0f;			// Spacing between cell centres (or nodes basically) 
	const int diag_movement = 14; 
	const int adj_movement  = 10; 

	/* Start and end point. Note, these are indices so size -1 is max size */
	size_t start_pt[2]	{ 9,0 };
	size_t end_pt[2]	{ 0,4 };
	
	/* ------------------------- Discretize the a cartesian grid with cells (or nodes) -------------------------------------
	Each cell has all the information about itself such as ID, the path to the cell, etc.
	Note, test this with a stack array variant as interested in difference of speed */ 

	// DYN_C2D is a contigious memory 2D dynamic array (std::vector< std::vector <primitive type> > is NOT contigious)
	RMF::DYN_C2D<Cell_ds> grid(col_PT);		// Construct 2D contiguous array from MF library with collumns known
											// (rows defined when constructing Cell_ds entites with parameters when defining grid)
	
	
	/* -------------------------- Function to generate grid with cell data -------------------------- */
	grid_generation(grid, row_PT, col_PT, spacing);


	/* -------------------------- Set start and end points (point A and B) -------------------------- */
	
	// This is so the point knows it is the start and end  
	grid(start_pt[0], start_pt[1]).set_Start();	
	grid(end_pt[0]	, end_pt[1]).set_End();

	// Let the start point know its a parent (For plot checks that the algorithm ended and started at the correct points (A and B),
	// really only relevant for pt B as its guaranteed for A but just for plot tidiness)
	grid(start_pt[0], start_pt[1]).set_Parent_eval(true);

	// Containers for start and end point values to reference
	Cell_ds pt_A = grid(start_pt[0], start_pt[1]);		
	Cell_ds pt_B = grid(end_pt[0], end_pt[1]); 

	/* -------------------------- Set first parent cell as start cell -------------------------- */
	Cell_ds* parent_pt = &pt_A;				// Parent_pt is a pointer which switches from difernt cells (should do that, currently pointing to container for parent cell after the initial condition)..
	std::vector<Cell_ds> container_for_parent_cell;		// Need to optomize this out if using a vector of possible parent cells is not good. 

	// Calculate costs from initial condition parent and start end points for referencing
	pt_A.set_Gcost(0); 
	pt_A.generate_heurcost(pt_B, diag_movement, adj_movement);
	pt_A.generate_fcost();
	
	parent_pt->set_Gcost(0);
	parent_pt->generate_heurcost(pt_B, diag_movement, adj_movement); 
	parent_pt->generate_fcost();

	/* -------------------------------  Obstacle generation function ------------------------------------------- */
	// Note, this is after the start and end points are defined as you do not, do not, want to make them an obstacle.. (I did)
	obstacle_generation(grid, row_PT, col_PT, 1);

	/*------------------------------- Initalize a priority list ----------------------------------*/
	/*Choose which cell to evaluate next based on smallest fcost/hcosts
	This has elements added to it (cells) as the algorithm progresses */

	std::vector<Cell_ds> priority_list;	
	size_t				 priolist_index = 0;
	
	// Reserve memory of the maximum size the priority list can be
	priority_list.reserve(col_PT * row_PT);		
	
	// Initialize iteration counter
	int iteration = 0; 

	/* ------------------------------ Start of A* search algorithm ------------------------------ */

	while (parent_pt->get_isEnd() == false) {

		// Predefine index of parent cell 
		size_t p_row = parent_pt->get_Row_pos(); 
		size_t p_col = parent_pt->get_Col_pos();
		
		/* Get points surrounding current_pt (will be 8 points including diagonals if parent is not a boundary cell) */
		/* Boundary cell filtering -> boundary cells must be treated differently as the available neighbouring cells are different*/

		// Switch for different types of boundaries
		/* B_IDs :	Normal cell = 0
					BL corner	= 1
					B edge		= 2
					BR corner	= 3
					L edge		= 4
					R edge		= 5
					TL corner	= 6
					T edge		= 7
					TR corner	= 8
		*/

		switch (parent_pt->get_Boundary_pt_val()) {
			case 0: {	// Normal cell

				// Indexed lrow_lim[0], urow_lim[1], lcol_lim[2], lcol_lim[3]
				size_t neighbour_limits[4]{ p_row - 1, p_row + 1, p_col - 1, p_col + 1 };

				// void A_star_algorithm(size_t(&loop_index)[4], RMF::DYN_C2D<Cell_ds> & grid, std::vector<Cell_ds> & priority_list, size_t & priolist_index, Cell_ds & parent_pt, Cell_ds & pt_B, const size_t & p_col, const size_t & p_row, const int& diag_movement, const int& adj_movement)				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 1: {	// Bottom left corner cell (do not consider points at row - 1, col - 1 area) 

				size_t neighbour_limits[4]{ p_row, p_row + 1, p_col , p_col + 1 };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 2: {	// Bottom edge 

				size_t neighbour_limits[4]{ p_row, p_row + 1, p_col - 1, p_col + 1 };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 3: {	// Bottom right corner 

				size_t neighbour_limits[4]{ p_row, p_row + 1, p_col - 1, p_col };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 4: {	// Left edge 

				size_t neighbour_limits[4]{ p_row - 1, p_row + 1, p_col, p_col + 1 };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 5: {	// Right edge

				size_t neighbour_limits[4]{ p_row - 1, p_row + 1, p_col -1, p_col };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 6: {	// Top left corner 

				size_t neighbour_limits[4]{ p_row - 1, p_row, p_col, p_col + 1 };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 7: {	// Top edge

				size_t neighbour_limits[4]{ p_row - 1, p_row, p_col - 1, p_col + 1 };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
			case 8: {	// Top right corner

				size_t neighbour_limits[4]{ p_row - 1, p_row, p_col - 1, p_col };
				A_star_algorithm(neighbour_limits, grid, priority_list, priolist_index, parent_pt, pt_B, p_col, p_row, diag_movement, adj_movement);
				break;
			}
		}

		/* Select new parent cell for the minimum fcost -> minimum hcost if fcosts are the same */
		container_for_parent_cell = min_FHcost(priority_list, parent_pt, pt_A);	// Not ideal, but using vector for possible optimizations of code structure later
		parent_pt = &container_for_parent_cell[0];				// Set parent cell pointer to the new parent cell. 

#if DEBUG==0 
		// Print fcost values of priority list 
		for (size_t ip = 0; ip < priority_list.size(); ip++) {
			std::cout << "Fcost of cell " << priority_list[ip].get_ID() << " : " << priority_list[ip].get_Fcost() << "\n"; 
			std::cout << "Hcost of cell " << priority_list[ip].get_ID() << " : " << priority_list[ip].get_Hcost() << "\n";
		}
#endif
		/* Set new parent cell to acknoledge it is a parent cell.*/
		// Note, cell ID is the same as contiguous index for DYN_C2D object "grid" since ID is 0 -> N cells and contig ind is 0 -> N cells.
		grid(parent_pt->get_ID()).set_Parent_eval(true);
		priority_list[parent_pt->get_Priolist_index()].set_Parent_eval(true);

		/* NOTES: 
		-	When the end cell (target or whatever) is met, the while condition is checked. Therefore, the end cell is actually not calculated.
			But the new parent (which will be the end cell) knows the parent it comes from, i.e. the path, so is fine. Just weird for plotting i guess
		-   NEED A FUNCTION TO CLEAR THE PRIORITY LIST OF REDUNDENT CELLS to avoid excess looping over list*/

		// Iteration logger 
		iteration += 1;
		std::cout << "Iteration: " << iteration << "\n\n";
		if (iteration > 10000) {
			std::cin.get(); 
		}

		/* "Real time" plot grid using ASCII characters for console.. need a better plotter */
#if DEBUG==0
		for (int i_plot = 0; i_plot < row_PT; i_plot++) {

			std::cout << "| ";

			for (int j_plot = 0; j_plot < col_PT; j_plot++) {

				// Plot start and end poitns with them beign a parent at one point (i.e. the algorithm works, this is just to check for bugs) 
				if ((grid(i_plot, j_plot).get_isStart() == true && grid(i_plot, j_plot).get_Parent_eval() == true) || (grid(i_plot, j_plot).get_isEnd() == true && grid(i_plot, j_plot).get_Parent_eval() == true)) {
					std::cout << "@ "; 
				}
				// Plot start and end points (If not the above then there is a problem)
				else if (grid(i_plot, j_plot).get_isStart() == true || grid(i_plot, j_plot).get_isEnd() == true) {
					std::cout << "% ";
				}
				// If cell is a parent
				else if (grid(i_plot, j_plot).get_Parent_eval() == true) {
					std::cout << "P ";
				}
				// If cell has been evaluated 
				else if (grid(i_plot, j_plot).get_On_priolist() == true) {
					std::cout << "A ";
				}
				else if (grid(i_plot, j_plot).get_Obstacle() == true) {
					std::cout << (char)254u << " ";
				}

				else {
					std::cout << ". ";
				}
			}
			// Next row 
			std::cout << "|\n";
		}
#endif

	} // End of while loop for finding target/end point of algorithm. 


	/* ----------------------------- POST PROCESSING ------------------------------------------------ */
	/* -------------------------- Establish path from cell IDs -------------------------------------- */

	// Pointer to current end grid cell (which is, hopefully, in this container, as all parents were)
	Cell_ds* backtrack_cell_pointer = &container_for_parent_cell[0]; 
	
	// Initialize vector to contain all path cell IDs and reserve memory 
	std::vector<int> path_cell_IDs; 
	path_cell_IDs.reserve(row_PT* col_PT);

	// Store cell ID of target cell for plotting
	path_cell_IDs.push_back(backtrack_cell_pointer->get_ID());

	// While loop where the condition is if it is the start cell (which the grid cells are aware of) 
	while (backtrack_cell_pointer->get_isStart() == false) {

		// Get the ID of the parent cell of current cell backtrack_cell_pointer is pointing too. 
		int from_parent_cell_ID = backtrack_cell_pointer->get_From_Cell_ID();
		path_cell_IDs.push_back(from_parent_cell_ID);				// Store path cell ID 

		// Point to the parent cell of the current cell and repeat until the start point is reached. 
		backtrack_cell_pointer = &grid(from_parent_cell_ID);
		backtrack_cell_pointer->set_Plot_path_activator();			// Set the cell to plot when called by the plotter
	}

	/* ----------------------------------------------- Plot the path --------------------------------------------- */

	// Print display output verbose 
	std::cout << "\n\n Path the A* algorithm found: \n";

	/* Plot grid using ASCII characters for console.. need a better plotter */
	for (int i_plot = 0; i_plot < row_PT; i_plot++) {

		std::cout << "| ";

		for (int j_plot = 0; j_plot < col_PT; j_plot++) {
		
			// Plot start and end poitns with them beign a parent at one point (i.e. the algorithm works, this is just to check for bugs) 
			if ((grid(i_plot, j_plot).get_isStart() == true && grid(i_plot, j_plot).get_Parent_eval() == true) || (grid(i_plot, j_plot).get_isEnd() == true && grid(i_plot, j_plot).get_Parent_eval() == true)) {
				std::cout << "@ ";
			}
			// Plot start and end points  (If not the above then there is a problem)
			else if (grid(i_plot, j_plot).get_isStart() == true || grid(i_plot, j_plot).get_isEnd() == true) {
				std::cout << "% ";
			}
			// Plot the path points
			else if (grid(i_plot, j_plot).get_Plot_path_activator() == true) {
				std::cout << "& "; 
			}
			// If cell was a parent
			else if (grid(i_plot, j_plot).get_Parent_eval() == true) {
				std::cout << "P ";
			}
			// If cell has been evaluated 
			else if (grid(i_plot, j_plot).get_On_priolist() == true) {
				std::cout << "A ";
			}
			else if (grid(i_plot, j_plot).get_Obstacle() == true) {
				std::cout << (char)254u << " ";
			}
			else {
				std::cout << ". ";
			}
		}

		// Next row 
		std::cout << "|\n";
	}
}


