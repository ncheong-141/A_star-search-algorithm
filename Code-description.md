# A_start-search-algorithm
A star search algorithm to find the optimal path from A (start point) to B (target point) around various objects. 


Program objectives:

      - A star search algorithm to find the optimal path from A (start point) to B (target point) around various objects. 
        
	How it works: 
      - At a cell(point A for start), the algorithm looks at all cells surrounding it and calculates:
	          -> gcost () =	distance from starting node (point A)
		    (typical for abstract A* algorithms is to have distance between nodes as 1, and multiply by 10 for 
		    irradicating FP numbers)
	          -> hcost () =  heuristic function / distance from end node (point B). 
	          -> fcost()  = total cost -> gcost() + hcost() 
      - The algorithm finds the lowest fcosts of the surrounding nodes and then choses that to be the next parent node and repeat the 		above until the end point is found. 
      - If there is obstacles points, they are not evaluated. 
      - If there is two points with the same fcost, the smallest hcost takes preference. 
      - All fcosts are stored and the lowest cost is stored each time
      - Do not consider neighbouring nodes which have already been evaluated.
      - A node that has been a parent node is not evaluated again as a parent node
      - The path it took to get to the parent node has to be stored as this will hold the optimal path.

Notes on code structure: 

        - Algorithm linked to numerical function static library. 
        
        - The spatial discretization (mesh) is created using an 2D array of cell data structures (Cell_ds).
        
        - The array is the contiguous dynamic std::vector(DYN_C2D) which is implemented on the static library and 
	  generated in a function for modularity. 
        
        - Each cell holds data which holds its gcost, hcost, fcosts, cell ID, from_cell_ID (for path determination), 
          x_position, y_position etc.
          
        - Additionally, the cell also holds algorithm flow control parameters (bools) such as:
          -> Boundary point value, this is an integer value to establish if the cell is on the left, right, top etc wall. 
             since the calculation at these cells are DIFFERENT from other cells. (e.g. you cannot evaluate cells to the
             left on the left wall since they do not exist. 
          -> Obstacle point, if it is an obstacle do not evaluate
          -> parent_eval, if true then don't re-evaluate as its been a parent before
          -> on_priolist and priolist_index, if true then its on the priority list 
          
        - Calculates for g,h,f costs are done using the member functions of Cell_ds. 
	
	- Cell_ds class is encapsulated for code safety. 
        
        - A priority list is used for containing data of cells which have been evaluated. The algorithm selects the next 
          best cell on the priority list with every iteration until the end point is found. 
          
        - For every iteration, the current cell is pointed to by a pointer, which is de-referenced to obtain the current
          cell values. 

Learning objectives: 

	- Practisting code and data structures 
