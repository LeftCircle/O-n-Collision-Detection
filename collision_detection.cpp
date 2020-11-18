#include "ContestDetector.h"
#include "TwoDScene.h"
#include <iostream>
#include <cmath>
#include <list>
#include <set>

// Given particle positions, computes lists of *potentially* overlapping object
// pairs. How exactly to do this is up to you.
// Inputs: 
//   scene:  The scene object. Get edge information, radii, etc. from here. If 
//           for some reason you'd also like to use particle velocities in your
//           algorithm, you can get them from here too.
//   x:      The positions of the particle.
// Outputs:
//   pppairs: A list of (particle index, particle index) pairs of potentially
//            overlapping particles. IMPORTANT: Each pair should only appear
//            in the list at most once. (1, 2) and (2, 1) count as the same 
//            pair.
//   pepairs: A list of (particle index, edge index) pairs of potential
//            particle-edge overlaps.
//   phpairs: A list of (particle index, halfplane index) pairs of potential
//            particle-halfplane overlaps.
void ContestDetector::findCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs)
{
    // Your code goes here!
    // If this is the first frame, then initialize grid. 
    const std::vector<double> r = scene.getRadii();
    const std::vector<std::pair<int, int>>& edges = scene.getEdges();
    int n_edges = scene.getNumEdges();
    int n_particles = scene.getNumParticles();
    int n_halfplanes = scene.getNumHalfplanes();
    // zero out the sets 
    pppairs.clear();
    pepairs.clear();
    // note, might not have to do this since we are always checking everything against halfplanes. 
    phpairs.clear();
    
    if (!(firstListConstruction)){
        // Just throwing everything used once here. 
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
        // Begin by determining the min/max X
        double maxX = x[0];
        double minX = x[0];
        double totalSize = maxX - minX;
        for (int i = 0; i < n_particles; i++){
            Vector2s aabb;
            boundingBoxX(scene, x, i, true, r, aabb);
            if (aabb[0] < minX){
                minX = aabb[0] - 0.1;
            }
            if (aabb[0] > maxX){
                maxX = aabb[0] + 0.1;
            }
            if (aabb[1] < minX){
                minX = aabb[1] - 0.1;
            }
            if (aabb[1] > maxX){
                maxX = aabb[1] + 0.1; 
            }    
        }

        for (int i = 0; i < n_edges; i++){
            Vector2s aabb;
            boundingBoxX(scene, x, i, false, r, aabb);
            if (aabb[0] < minX){
                minX = aabb[0] - 0.1;
            }
            if (aabb[0] > maxX){
                maxX = aabb[0] + 0.1;
            }
            if (aabb[1] < minX){
                minX = aabb[1] - 0.1;
            }
            if (aabb[1] > maxX){
                maxX = aabb[0] + 0.1;
            }
        } 
        // Max and min X have been found. Box size can be initialized. 
        double columnSize = (maxX - minX) / NCELLS;

        // Now loop through all particles and edges to find the column they fit into 
        for (int i = 0; i < n_particles; i++){
            // build the bounding box
            Vector2s aabb;
            boundingBoxX(scene, x, i, true, r, aabb);
            for (int j = 0; j < 2; j++){
                ListComponent particleListPiece;
                // setting max first then min 
                particleListPiece.setIsMin(j);
                particleListPiece.setIsParticle(true);
                particleListPiece.setIdx(i);
                if (j == 0){
                    // set max 
                    particleListPiece.setPaddedPos(aabb[1]);
                }else{
                    particleListPiece.setPaddedPos(aabb[0]);
                }
                // Add list piece to appropriate colomn
                int columnN = floor((aabb[0] - minX) / columnSize);
                backgroundMatrix[columnN].push_back(particleListPiece);
            }
        } // end looping through particles
        // Now initialize edges 
        for (int i = 0; i < n_edges; i++){
            Vector2s aabb;
            boundingBoxX(scene, x, i, false, r, aabb);
            for (int j = 0; j < 2; j++){
                ListComponent edgeListPiece;
                edgeListPiece.setIsMin(j);
                edgeListPiece.setIsParticle(false);
                edgeListPiece.setIdx(i);
                if (j == 0){
                    edgeListPiece.setPaddedPos(aabb[1]);
                }else{
                    edgeListPiece.setPaddedPos(aabb[0]);
                }
                int columN = floor((aabb[0] - minX) / columnSize) ;
                backgroundMatrix[columN].push_back(edgeListPiece);
            }
        }


        // Now intialize the bigList;
        for (int i = 0; i < NCELLS; i++){
            // add each list to the big list then sort. 
            if(backgroundMatrix[i].size() > 0){
                bubbleSort(backgroundMatrix[i]);
                theBigListX.merge(backgroundMatrix[i]);
            }
        }
        bubbleSort(theBigListX);
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        possibleCollisions(theBigListX, pppairs, pepairs);
        
        checkCollisions(scene, x, pppairs, pepairs);
        
    }else{
        updateList(theBigListX, x, scene);
        possibleCollisions(theBigListX, pppairs, pepairs);
        checkCollisions(scene, x, pppairs, pepairs); 
    }  
    
    // Now add every particle and edge to potential halfplane collisions 
    for (int i = 0; i < n_halfplanes; i++){
        for (int j = 0; j < n_particles; j++){
            std::pair<int, int> phpair;
            phpair.first = j;
            phpair.second = i;
            phpairs.insert(phpair);
        }
    }
    
    firstListConstruction = true;
    
    
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     Function for initializing the grid
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  
void initializeColumns(BackgroundGrid &grid, const VectorXs &x, const TwoDScene &scene, const std::vector<double> &r, const std::vector<std::pair<int, int>>& edges){
    int n_edges = scene.getNumEdges();
    int n_particles = scene.getNumParticles();
    int n_halfplanes = scene.getNumHalfplanes();
    
    // Begin by determining the min/max X
    double maxX = x[0];
    double minX = x[0];
    double totalSize = maxX - minX;
    for (int i = 0; i < n_particles; i++){
        Vector2s aabb;
        boundingBoxX(scene, x, i, true, r, aabb);
        if (aabb[0] < minX){
            minX = aabb[0] - 0.1;
        }
        if (aabb[0] > maxX){
            maxX = aabb[0] + 0.1;
        }
        if (aabb[1] < minX){
            minX = aabb[1] - 0.1;
        }
        if (aabb[1] > maxX){
            maxX = aabb[1] + 0.1; 
        }    
    }
    
    for (int i = 0; i < n_edges; i++){
        Vector2s aabb;
        boundingBoxX(scene, x, i, false, r, aabb);
        if (aabb[0] < minX){
            minX = aabb[0] - 0.1;
        }
        if (aabb[0] > maxX){
            maxX = aabb[0] + 0.1;
        }
        if (aabb[1] < minX){
            minX = aabb[1] - 0.1;
        }
        if (aabb[1] > maxX){
            maxX = aabb[0] + 0.1;
        }
    } 
    // Max and min X have been found. Box size can be initialized. 
    double columnSize = (maxX - minX) / NCELLS;
    
    // Now loop through all particles and edges to find the column they fit into 
    for (int i = 0; i < n_particles; i++){
        // build the bounding box
        Vector2s aabb;
        boundingBoxX(scene, x, i, true, r, aabb);
        for (int j = 0; j < 2; j++){
            ListComponent particleListPiece;
            // setting max first then min 
            particleListPiece.setIsMin(j);
            particleListPiece.setIsParticle(true);
            particleListPiece.setIdx(i);
            if (j == 0){
                // set max 
                particleListPiece.setPaddedPos(aabb[1]);
            }else{
                particleListPiece.setPaddedPos(aabb[0]);
            }
            // Add list piece to appropriate colomn
            int columnN = floor((aabb[0] - minX) / columnSize);
            backgroundMatrix[columnN].push_back(particleListPiece);
        }
    } // end looping through particles
    // Now initialize edges 
    for (int i = 0; i < n_edges; i++){
        Vector2s aabb;
        boundingBoxX(scene, x, i, false, r, aabb);
        for (int j = 0; j < 2; j++){
            ListComponent edgeListPiece;
            edgeListPiece.setIsMin(j);
            edgeListPiece.setIsParticle(false);
            edgeListPiece.setIdx(i);
            if (j == 0){
                edgeListPiece.setPaddedPos(aabb[1]);
            }else{
                edgeListPiece.setPaddedPos(aabb[0]);
            }
            int columN = floor((aabb[0] - minX) / columnSize) ;
            backgroundMatrix[columN].push_back(edgeListPiece);
        }
    }
    
    // Now intialize the bigList;
    for (int i = 0; i < NCELLS; i++){
        // add each list to the big list then sort. 
        if(backgroundMatrix[i].size() > 0){
            theBigListX.merge(backgroundMatrix[i]);
        }
    }
}    

// Generates the bounding box for a particle or an edge 
void boundingBoxX(const TwoDScene &scene, const VectorXs &x, int idx, bool isParticle, const std::vector<double>& r, Vector2s& aabb){
    const std::vector<std::pair<int,int>>& edges = scene.getEdges();
    if (isParticle){
        double radius = r[idx];
        double pad = radius / 4;
        double position = x[2 * idx];
        double minX = position - radius - pad;
        double maxX = position + radius + pad;
        aabb[0] = minX;
        aabb[1] = maxX;
    }else{
        // Generate bounding box for edge
        int pidx1 = edges[idx].first;
        int pidx2 = edges[idx].second;
        
        double x1 = x[2 * pidx1];
        double x2 = x[2 * pidx2];
        
        double radius = r[pidx1];
        double pad = radius / 4;
        
        double minX = x1 < x2 ? x1 - radius - pad : x2 - radius - pad;
        double maxX = x1 < x2 ? x2 + radius + pad : x1 + radius + pad;
        
        aabb[0] = minX;
        aabb[1] = maxX;
    }
}

void boundingBoxY(const TwoDScene &scene, const VectorXs &x, int idx, bool isParticle, const std::vector<double>& r, Vector2s& aabb){
    const std::vector<std::pair<int,int>>& edges = scene.getEdges();
    if (isParticle){
        double radius = r[idx];
        double pad = radius / 4;
        double position = x[2 * idx + 1];
        double minX = position - radius - pad;
        double maxX = position + radius + pad;
        aabb[0] = minX;
        aabb[1] = maxX;
        
    }else{
        // Generate bounding box for edge
        int pidx1 = edges[idx].first;
        int pidx2 = edges[idx].second;
        
        double x1 = x[2 * pidx1 + 1];
        double x2 = x[2 * pidx2 + 1];
        
        double radius = r[pidx1];
        double pad = radius / 4;
        
        double minX = x1 < x2 ? x1 - radius - pad : x2 - radius - pad;
        double maxX = x1 < x2 ? x2 + radius + pad : x1 + radius + pad;
        
        aabb[0] = minX;
        aabb[1] = maxX;
    }
}

// Loops through all components within the list, and checks the 
// position value. It then starts moving back within the list from the current position 
// until the insertion position is found. 
// All other positions are moved up one
void bubbleSort(std::list<ListComponent>& listToSort) {
	int nOperations = 0;
	int sizeMinusOne = listToSort.size() - 1;

	for (int i = sizeMinusOne; i > 0; --i) {
		bool swapOccured = false;

		std::list<ListComponent>::iterator it = listToSort.begin();
		std::list<ListComponent>::iterator it2 = listToSort.begin();
		++it2;

		ListComponent currentLC = *it;
		ListComponent compareLC = *it2;
		double currentPos = currentLC.getPaddedPos();
		double comparePos = compareLC.getPaddedPos();
		for (int j = 0; j < i; j++) {
			currentLC = *it;
			compareLC = *it2;
			currentPos = currentLC.getPaddedPos();
			comparePos = compareLC.getPaddedPos();

			if (comparePos < currentPos) {
				std::swap(*it2, *it);
				swapOccured = true;
			}
			it++;
			it2++;
			nOperations += 1;
		}
		if (swapOccured == false) {
			break;
		}
	}
}

//prints the list
void printList(std::list<ListComponent>& listToPrint) {
	std::list<ListComponent>::iterator it2 = listToPrint.begin();
	for (std::list<ListComponent>::iterator it = listToPrint.begin(); it != listToPrint.end(); ++it) {
		ListComponent currentListComponent = *it;
		it2++;
		double pos = currentListComponent.getPaddedPos();
		int idx = currentListComponent.getIdx();
		int type = currentListComponent.getIsParticle();
		printf("pos = ");
		std::cout << pos;
		printf(" ");
		printf("idx = ");
		std::cout << idx;
		printf(" type = ");
		std::cout << type;
		printf(" Min or max = ");
		std::cout << currentListComponent.getIsMin() << std::endl;

	}
}

void updateList(std::list<ListComponent>& listToUpdate, const VectorXs& x, const TwoDScene& scene){
    
    const std::vector<std::pair<int, int> >& edges = scene.getEdges();
	const std::vector<double> r = scene.getRadii();
	// Go through the list and update positions of the particles and edges
	for (std::list<ListComponent>::iterator it = listToUpdate.begin(); it != listToUpdate.end(); ++it) {
		// look at the current list component, and check if it is a particle or an edge 
		ListComponent currentListComponent = *it;

		if (currentListComponent.getIsParticle()){
			// Particle has been found. Find particle with this index. 
			int pIndex = currentListComponent.getIdx();
			double radius = r[pIndex];
			Vector2s aabb;
            boundingBoxX(scene, x, pIndex, true, r, aabb);
			// Check to see if this is a max or a min;
			if (currentListComponent.getIsMin()) {
				it->setPaddedPos(aabb[0]);
			}else{
				// Max found
				it->setPaddedPos(aabb[1]);	
			}        
		}// end particle found
		else {
            // edge has been found. Find edge index
			int eIdx = currentListComponent.getIdx();
			Vector2s aabb;
            boundingBoxX(scene, x, eIdx, false, r, aabb);
			if (currentListComponent.getIsMin()){
				it->setPaddedPos(aabb[0]);
			}else {
				it->setPaddedPos(aabb[1]);
			}
		} // End if particle or edge 
	} // End loop through list 

	// Now sort the list 
	bubbleSort(listToUpdate);
}

void possibleCollisions(const std::list<ListComponent>& list, PPList& pppairs, PEList& pepairs) {
	// Loop through the big list. If the component is a min, set it to active. 
	// Create a vector of active indeces. Each time a component is added to the vector, loop through 
	// the vector and say that the new component can collide with all already active components. 
	// 
	// Loop through the big list. If a min is found, set the second iterator to the next position.
	// Iterate the second iterator until the current particle is no longer active, then continue looping 
	// through the list with the first iterator. For every other min we find, state a possible collision. 
	
	std::list<ListComponent>::const_iterator activeIterator = list.begin();
	for (std::list<ListComponent>::const_iterator it = list.begin(); it != list.end(); ++it) {

		ListComponent clc = *it;
		if (clc.getIsMin()) {
			// Min found. Set the activeIterator equal to it 
			activeIterator = it;
			++activeIterator;

			while (true) {
				// First check if the next component has the same index and type as the current list component. 
				// if so, set the tag to inactive and break 
				ListComponent nlc = *activeIterator;
				if (nlc.getIdx() == clc.getIdx() && nlc.getIsParticle() == clc.getIsParticle()) {
					break;
				}

				// if we find a min, then we found a possible collision 
				if (nlc.getIsMin()) {
					// Check the type of clc and nlc and add possible collision to appropriate list
					std::pair<int, int> collisionPair;

					if (clc.getIsParticle()) {
						// current component is a particle
						if (nlc.getIsParticle()) {
							// pp collision 
							collisionPair.first = clc.getIdx();
							collisionPair.second = nlc.getIdx();
							pppairs.insert(collisionPair);
						}
						else {
							// pe collision. 
							collisionPair.first = clc.getIdx();
							collisionPair.second = nlc.getIdx();
							pepairs.insert(collisionPair);
						}
					}
					else {
						// current is an edge
						if (nlc.getIsParticle()) {
							// pe collision 
							collisionPair.first = nlc.getIdx();
							collisionPair.second = clc.getIdx();
							pepairs.insert(collisionPair);
						}
						else {
							// edge edge collision  -> this should never occur.
						}
					}
				} // end finding of another min (active) object 
				++activeIterator;
			} // Particle is no longer active. Continue iterating starting with the next component. 

		} // End check for if the clc is min. 
	}
}

void checkCollisions(const TwoDScene& scene, const VectorXs& x, PPList& pppairs, PEList& pepairs) {
	// Loop through the PP list and check to see if the Y bounding boxes overlap.
	// We will recalculate the bounding box for each potential overlap. 
	// If they do not overlap, remove the item from the list. 
	const std::vector<std::pair<int, int> >& edges = scene.getEdges();
	const std::vector<double> r = scene.getRadii();
    // Start with the PPList 
	std::set<std::pair<int, int> >::iterator it = pppairs.begin();
	while (it != pppairs.end()) {
		// Check the y components of the two particles bounding boxes.
		int idx1 = it->first;
		int idx2 = it->second;
        
		double y1min = x[2 * idx1 + 1] - r[idx1];
		double y1max = y1min + 2 * r[idx1];

		double y2min = x[2 * idx2 + 1] - r[idx2];
		double y2max = y2min + 2 * r[idx2];
        
		if ((y2min <= y1max && y2max >= y1min) || (y1min <= y2max && y1max >= y2min)) {
			// Collision found! move on to the next one
            ++it;
			continue;
		}
		else {
			//remove this element from the set
			pppairs.erase(it++);
		}
	}
	// Check PE list 
	std::set<std::pair<int, int>>::iterator eit = pepairs.begin();
	while (eit != pepairs.end()) {
		int pidx = eit->first;
		int eidx = eit->second;

		double pmin = x[2 * pidx + 1] - r[pidx];
		double pmax = pmin + 2 * r[pidx];

		// setting up edge bounding box 
		int eidx1 = edges[eidx].first;
		int eidx2 = edges[eidx].second;
		//double radius = scene.getEdgeRadii()[i];
		Vector2s p1 = x.segment<2>(2 * eidx1);
		Vector2s p2 = x.segment<2>(2 * eidx2);

		// Extending AABB by particle radius instead of edge radius (particle radius is larger and at the ends anyway)
		double p_radius = r[eidx1];
		double emin = p1[1] < p2[1] ? p1[1] - p_radius : p2[1] - p_radius;
		double emax = p1[1] < p2[1] ? p2[1] + p_radius : p1[1] + p_radius;

		if ((pmin <= emax && pmax >= emin) || (emin <= pmax && emax >= pmin)) {
			// Collision found! move on to the next one 
            ++eit;
            continue;
		}
		else {
			//remove this element from the set 
			pepairs.erase(eit++);
		}
	}
}