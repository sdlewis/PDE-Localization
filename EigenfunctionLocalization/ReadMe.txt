Read Me for Eigenvalue Localization Package

The following is a package in development for the localization of eigenfunctions of elliptic partial differential equations.

The package includes the following files:

	twoDPartition.py --- Provides the partition object, an object which implements a 2 dimensional subregion of a square in the first quadrant with lower left corner (0,0), and a partition of it into subregions Omega. Each region has certain data attached to it, including its interior (self.Omega[key]), its boundary (self.bOmega[key]), and its interior points which are adjacent to the boundary (self.abOmega[key]).

Its principle pieces of data are:
	--n the number of grid ticks in each direction including the boundary (e.g., 102)
	--meshSize the thickness of the meshgrid (e.g., 1./102)
	--Domain the underlying domain points
	--bDomain the boundary of the domain (does not intersect with Domain)
	--Omega a dictionary whose ith index is the ith set in the partition
	--bOmega a dictionary whose ith index is the boundary of the ith set
	--abOmega a dict whose ith index is the set of points in Omega adjacent to bOmega

Its principle methods are:
	--evolve which updates the regions probabilistically (by using the deterministic grow and shrink methods)
	--merge which merges two sets into one (and handles the necessary updates to the boundary points self.bOmega and the adjacent boundary points self.abOmega)
	--updateBoundary which updates the boundary after a change to the underyling sets is made
	--globalOperator returns the matrix representing the elliptic partial differential operator imported in setup.py for the whole Domain
	--localOperator returns the matrix representing the elliptic partial differential operator imported in setup.py for the specified subregion

	DFJM.py --- Responsible for handling the partition object and running running a search for the optimal partition. Methods for search include simulatedAnnealing, in which steps are taken probabilistically, and lessThan, in which steps are taken if and only if they improve the current energy.

	utilities.py and graphingUtilities.py --- Responsible for simple data manipulation and graphing features. All necessary methods have a handler in the partition object.

	harmonicSquare.py  and quantumSquare.py --- Two different options for partial differential operators. Responsible for setting up the matrix representing the global operator and provides the parameters to used by the rest of the package. Imported through setup.py

	
