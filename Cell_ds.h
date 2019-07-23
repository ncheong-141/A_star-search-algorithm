#pragma once
/* Interface for Node_ds*/
/* Class for cell data structure containing all information about a node*/


struct Cell_ds {

private:
	const int		cell_ID;				// cell ID to identify each cell for path calculations (same as the contigious index!)
	int			from_cell_ID;				// From where the cell has been defined from (the path it took to reach the cell)
	const size_t		row_pos, col_pos, contig_pos;		// Cell row and col position for index manipulation
	const float		x_pos, y_pos;				// The spatial position of the cell centres 
	int			gcost, heurcost, fcost;			// Cost variables for A* algorithm 

	/* Cell and boundary cell settings */
	bool			parent_eval, start_pt, end_pt;		// Logic variables for defining cell types (e.g. if it is the start cell or been evalled)
	bool			boundary_pt;				// Logic variable for defining if the cell is a boundary_pt (only needs one boolean condition then) 
	int			boundary_pt_val;			// Each boundary cell will get a value corresponding to left, right edge etc etc. so calculations will not include cells which dont exist 
	bool			obstacle_pt;				// Logic variable to determine if the cell is an obstacle or not. 
	
	bool			on_priolist;				// Logic variable to determine if the cell is on the priority list
	size_t			priolist_index;				// Integer variable for where the cell is on the priority list (so can update it)

	// Plot variables (Doesnt need these to run) 
	bool			plot_path_activator;			// Just for plotting using a simple loop. 

public:
	
	/* Constructors and destructors */
	
	// Defuault
	Cell_ds() :

		cell_ID(0), from_cell_ID(0),
		row_pos(0), col_pos(0), contig_pos(0),
		x_pos(0), y_pos(0),
		gcost(0), heurcost(0), fcost(0),
		parent_eval(false), start_pt(false), end_pt(false),
		boundary_pt(false), boundary_pt_val(0), obstacle_pt(false), on_priolist(false), priolist_index(0), 
		plot_path_activator(false) 
	{
		std::cout << "Error. Do not use the default constructor for this data structure atm.\n";
		exit(-1);

	}

	// Parameters input for std::vector parameter instansation of grid 
	Cell_ds(int ID, size_t row_ind, size_t col_ind, size_t col_size, float x_position, float y_position) : 
		
		cell_ID(ID), from_cell_ID(0),
		row_pos(row_ind), col_pos(col_ind), contig_pos((row_ind*col_size) + col_ind),
		x_pos(x_position), y_pos(y_position),
		gcost(0), heurcost(0), fcost(0), 
		parent_eval(false), start_pt(false), end_pt(false),
		boundary_pt(false), boundary_pt_val(0), obstacle_pt(false), on_priolist(false), priolist_index(0), 
		plot_path_activator(false) {}

	// Copy constructor
	//Cell_ds operator=(const Cell_ds& new_cell) {
	//	return new_cell; 
	//}
	

	/* Getters and setters (encapsulating class for code safety) */
	// Start and end point getters and setter
	void	set_Start() { start_pt = true;	}
	void	set_End()	{ end_pt = true;	}
	bool	get_isStart() { return start_pt;   }
	bool	get_isEnd() { return end_pt; }
	
	// Getting and setting cell condition (if it was a parent or not, boundary conditions etcs) 
	void	set_Parent_eval(bool new_cell_cond) { parent_eval = new_cell_cond; }			
	bool	get_Parent_eval() { return parent_eval; }

	void	set_Boundary_pt() { boundary_pt = true;  }
	bool	get_Boundary_pt() { return boundary_pt; }
	void	set_Boundary_pt_val(int val) { boundary_pt_val = val;  }
	int		get_Boundary_pt_val() { return boundary_pt_val;  }
	void	set_Obstacle() { obstacle_pt = true; }
	bool	get_Obstacle() { return obstacle_pt; }

	void	set_On_priolist() { on_priolist = true;  }
	bool	get_On_priolist() { return on_priolist; }
	void	set_Priolist_index(size_t ind) { priolist_index = ind; }
	size_t	get_Priolist_index() { return priolist_index; }

	// From_node_ID to establish paths 
	void	set_From_Cell_ID(int setCellIDval)	{ from_cell_ID = setCellIDval;	}		
	int	get_From_Cell_ID()			{ return from_cell_ID;		}

	// Get cost variables of cell
	int		get_Gcost() { return gcost; }			
	int		get_Hcost() { return heurcost; }
	int		get_Fcost() { return fcost; }

	// Only one setter fucntion for costs (for gcost) as the lower gcost of newer paths can overwrite older paths (updating movement paths)
	void	set_Gcost(int new_cost) { gcost = new_cost; }

	// Get cell charcteristic data, Note, xpos and ypos does not change
	float	get_Xpos()		{ return x_pos;	  }
	float	get_Ypos()		{ return y_pos;	  }
	int		get_ID()	{ return cell_ID; }
	size_t	get_Row_pos()		{ return row_pos; }
	size_t	get_Col_pos()		{ return col_pos; }
	size_t	get_Contig_pos()	{ return contig_pos; }

	// Text plotter setter (condition to mark cell to plot as a path) 
	bool get_Plot_path_activator() { return plot_path_activator; }
	void set_Plot_path_activator() { plot_path_activator = true;  }

	/* Operator overload to access cells based on their ID */

	/* Operator overload to allow copy?? */
	//void operator=(Cell_ds& overwriting_values) {
	//	this = overwriting_values;
	//}


	/* Member functions for manipulating node data*/

	/* gcost -> distance of parent to start + movement to cell */
	int generate_gcost(Cell_ds*& parent, int movement) {
		return  parent->gcost + movement;
	}

	/* Heuristic cost (decided abitrarily) */
	int generate_heurcost(Cell_ds& pt_end, int diag_movement, int adj_movement) {

		// Casting to int as size_t is unsigned and negative values can occur in difference (in which if size_t is negative, it becomes its max number - 1..)
		int pt_end_row_pos = (int)pt_end.get_Row_pos(); 
		int pt_end_col_pos = (int)pt_end.get_Col_pos(); 

		int difference_row_ind = RMF::abs(pt_end_row_pos - (int)row_pos); 
		int difference_col_ind = RMF::abs(pt_end_col_pos - (int)col_pos);

		/* Resolve row-wise movement as diagonal and adjacent 
		- Dealing in absolutes (not a star wars 3 reference) so doesnt matter if current point is left right above below etc. 
		- Only allowed to move diagonally till difference_col_ind value and differnece_row_ind
		- Index variables since row_pos and col_pos do not actually change.  */

		int row_indexer = 0;
		int col_indexer = 0;
		while ((row_indexer < difference_row_ind) && (col_indexer < difference_col_ind)) {

			heurcost += diag_movement;	// Add diagonal movement 
			col_indexer += 1;			// Diagonal movement has incrememented the row and collumn position 
			row_indexer += 1;
		}

		// Remainder of movement is adjacent; conditions for vertical/horizontal adjacent movement
		if (col_indexer == difference_col_ind) {		// for adjecent vertical movement
			for (int adj_m = 0; adj_m < (difference_row_ind - row_indexer); adj_m++) {
				heurcost += adj_movement;
			}
		}
		else if (row_indexer == difference_row_ind) {	// for adjacent horizonal movement
			// difference_col_ind - col_indexer is the remainder of collumn-wise movements
			for (int adj_m = 0; adj_m < (difference_col_ind - col_indexer); adj_m++) {
				heurcost += adj_movement;
			}
		}
		return heurcost; 
	}

	int generate_fcost() {
		fcost = gcost + heurcost;
		return fcost;
	}
};

