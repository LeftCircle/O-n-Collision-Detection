#ifndef CONTEST_DETECTOR_H
#define CONTEST_DETECTOR_H

#include "CollisionDetector.h"
#include <iostream>
#include <vector>
#include <list>
#include <set>

typedef std::set<std::pair<int, int> > PPList;
typedef std::set<std::pair<int, int> > PEList;
typedef std::set<std::pair<int, int> > PHList;

class ContestDetector : public CollisionDetector
{
 public:
  ContestDetector() {}

  virtual void performCollisionDetection(const TwoDScene &scene, const VectorXs &qs, const VectorXs &qe, DetectionCallback &dc);

 private:
  void findCollidingPairs(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs, PHList &phpairs);
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//               
//                Personal Classes
//             
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// List Components
class ListComponent
{
public:
    void setIsMin(bool minMax){
        isMin = minMax;
    }
    
    void setIsParticle(bool isP){
        isParticle = isP;
    }
        
    void setIdx(int index){
        idx = index;
    }
    
    void setPaddedPos(double position){
        paddedPos = position;
    }
    // Getters
    bool getIsMin(){
        return isMin;
    }
    
    bool getIsParticle(){
        return isParticle;
    }
    
    int getIdx(){
        return idx;
    }
    
    double getPaddedPos(){
        return paddedPos;
    }
    
    // Overloaded operators 
    bool operator< (ListComponent& rhs) {
		return this->getPaddedPos() < rhs.getPaddedPos();
	}
    
private:
    bool isMin;
    bool isParticle;
    int idx;
    double paddedPos;
};

std::list<ListComponent> theBigList;
std::list<ListComponent>& theBigListX = theBigList;

// Functions for sorting the list 

void bubbleSort(std::list<ListComponent>& listToSort);
void updateList(std::list<ListComponent>& listToUpdate, const VectorXs& x, const TwoDScene& scene);
void printList(std::list<ListComponent>& listToPrint);

// Constructing initial conditions

const double NCELLS = 100;

typedef Eigen::Matrix<std::list<ListComponent>, (int)NCELLS, 1> BackgroundGrid;
BackgroundGrid backgroundMatrixFinal;
BackgroundGrid &backgroundMatrix = backgroundMatrixFinal;
bool firstListConstruction = false;

// Functions for generating the grid
void initializeColumns(BackgroundGrid &grid, const VectorXs &x, const TwoDScene &scene, const std::vector<double> &r, const std::vector<std::pair<int, int>>& edges);


// Other functions

// Goes into a cell and returns the x or y list from that cell
void possibleCollisions(const std::list<ListComponent>& list, PPList &pppairs, PEList &pepairs);
void checkCollisions(const TwoDScene &scene, const VectorXs &x, PPList &pppairs, PEList &pepairs);
void boundingBoxX(const TwoDScene &scene, const VectorXs &x, int idx, bool isParticle, const std::vector<double>& r, Vector2s& aabb);
void boundingBoxY(const TwoDScene &scene, const VectorXs &x, int idx, bool isParticle, const std::vector<double>& r, Vector2s& aabb);


#endif
