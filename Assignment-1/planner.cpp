/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include <random>
#include <vector>
#include <array>
#include <algorithm>
#include <queue>
using std::queue;
#include <stdlib.h>
using namespace std;
#include <tuple>
#include <string>
#include <stdexcept>
#include <regex> // For regex and split logic
#include <iostream> // cout, endl
#include <fstream> // For reading/writing files
#include <assert.h> 
//#include "../../../../../usr/include/c++/9/bits/algorithmfwd.h"

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTSTAR     1
#define PRM         2
#define ALL         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654
#define TIMELIMIT 60
#define min(a,b) ((a)<(b)?(a):(b))

//the length of each link in the arm
#define LINKLENGTH_CELLS 10

#ifdef _BLAS64_
#define int_F ptrdiff_t
#else
#define int_F int
#endif

// Some potentially helpful imports
using std::vector;
using std::array;
using std::string;
using std::runtime_error;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                GIVEN FUNCTIONS                                                    //
//                                                                                                                   //
//*******************************************************************************************************************//

/// @brief 
/// @param filepath 
/// @return map, x_size, y_size
tuple<double*, int, int> loadMap(string filepath) {
	std::FILE *f = fopen(filepath.c_str(), "r");
	if (f) {
	}
	else {
		printf("Opening file failed! \n");
		throw runtime_error("Opening map file failed!");
	}
	int height, width;
	if (fscanf(f, "height %d\nwidth %d\n", &height, &width) != 2) {
		throw runtime_error("Invalid loadMap parsing map metadata");
	}
	
	////// Go through file and add to m_occupancy
	double* map = new double[height*width];

	double cx, cy, cz;
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			char c;
			do {
				if (fscanf(f, "%c", &c) != 1) {
					throw runtime_error("Invalid parsing individual map data");
				}
			} while (isspace(c));
			if (!(c == '0')) { 
				map[y+x*width] = 1; // Note transposed from visual
			} else {
				map[y+x*width] = 0;
			}
		}
	}
	fclose(f);
	return make_tuple(map, width, height);
}

// Splits string based on deliminator
vector<string> split(const string& str, const string& delim) {   
		// https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c/64886763#64886763
		const std::regex ws_re(delim);
		return { std::sregex_token_iterator(str.begin(), str.end(), ws_re, -1), std::sregex_token_iterator() };
}


double* doubleArrayFromString(string str) {
	vector<string> vals = split(str, ",");
	double* ans = new double[vals.size()];
	for (int i = 0; i < vals.size(); ++i) {
		ans[i] = std::stod(vals[i]);
	}
	return ans;
}

bool equalDoubleArrays(double* v1, double *v2, int size) {
    for (int i = 0; i < size; ++i) {
        if (abs(v1[i]-v2[i]) > 1e-3) {
            cout << endl;
            return false;
        }
    }
    return true;
}

typedef struct {
	int X1, Y1;
	int X2, Y2;
	int Increment;
	int UsingYIndex;
	int DeltaX, DeltaY;
	int DTerm;
	int IncrE, IncrNE;
	int XIndex, YIndex;
	int Flipped;
} bresenham_param_t;


void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size) {
	double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params) {
	params->UsingYIndex = 0;

	if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
		(params->UsingYIndex)++;

	if (params->UsingYIndex)
		{
			params->Y1=p1x;
			params->X1=p1y;
			params->Y2=p2x;
			params->X2=p2y;
		}
	else
		{
			params->X1=p1x;
			params->Y1=p1y;
			params->X2=p2x;
			params->Y2=p2y;
		}

	 if ((p2x - p1x) * (p2y - p1y) < 0)
		{
			params->Flipped = 1;
			params->Y1 = -params->Y1;
			params->Y2 = -params->Y2;
		}
	else
		params->Flipped = 0;

	if (params->X2 > params->X1)
		params->Increment = 1;
	else
		params->Increment = -1;

	params->DeltaX=params->X2-params->X1;
	params->DeltaY=params->Y2-params->Y1;

	params->IncrE=2*params->DeltaY*params->Increment;
	params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
	params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

	params->XIndex = params->X1;
	params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y) {
	if (params->UsingYIndex) {
        *y = params->XIndex;
        *x = params->YIndex;
        if (params->Flipped)
            *x = -*x;
    }
	else {
        *x = params->XIndex;
        *y = params->YIndex;
        if (params->Flipped)
            *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params) {
	if (params->XIndex == params->X2) {
        return 0;
    }
	params->XIndex += params->Increment;
	if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
		params->DTerm += params->IncrE;
	else {
        params->DTerm += params->IncrNE;
        params->YIndex += params->Increment;
	}
	return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
			 int x_size, int y_size) {
	bresenham_param_t params;
	int nX, nY; 
	short unsigned int nX0, nY0, nX1, nY1;

	//printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
		
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

	//printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
			return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
			 int x_size, int y_size) {
    double x0,y0,x1,y1;
    int i;
		
	 //iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
	y1 = 0;
	for(i = 0; i < numofDOFs; i++){
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}    
	return 1;
}
struct Node {
    double* joint;
    Node* parent;
    int nodeNum;
    double cost;
};

struct ExperimentResult {
    double planningTime;
    int numNodes;
    int planLength;
    double planQuality;
};

static int getAngleDiscretizationFactor(int numofDOFs) {
    return round((2 * PI) / (2 * asin(sqrt(2)/(2 * LINKLENGTH_CELLS * numofDOFs))));
}

// Given currJoint and joints (all the existing joints),
// sets closestNeighbor to the joint within joints that's closest to currJoint
// and also returns the distance between closestNeighbor and currJoint
double getClosestNeighborFromTree(double* currJoint, vector<Node*>* tree, int numofDOFs, Node** closestNeighbor) {

    double closestNeighborDistance = (pow(2 * PI, 2) * numofDOFs);

    for (int i = 0; i < tree->size(); i++) {
        Node* neighbor = (*tree)[i];
        double* neighborJoint = neighbor->joint;
        double currNeighborDistance = 0;
        for (int j = 0; j < numofDOFs; j++) {
            currNeighborDistance += pow(fabs(neighborJoint[j] - currJoint[j]), 2);
        }
        //printf("***neighborJoint is , [%f, %f, %f, %f, %f]\n",
        //            neighborJoint[0], neighborJoint[1], neighborJoint[2], neighborJoint[3], neighborJoint[4]);
        if(currNeighborDistance < closestNeighborDistance) {
            *closestNeighbor = neighbor;
            closestNeighborDistance = currNeighborDistance;
        }
    }
    return sqrt(closestNeighborDistance);
}

// Given currJoint and joints (all the existing joints),
// sets closestNeighbor to the joint within joints that's closest to currJoint
// and also returns the distance between closestNeighbor and currJoint
double getClosestNeighborFromTreeAndNearNodes(double* currJoint, vector<Node*>* tree, int numofDOFs, Node** closestNeighbor,
        vector<Node*>* nearNodes, vector<double>* nearNodeDistances, double radius) {

    radius = pow(radius, 2);
    double closestNeighborDistance = (pow(2 * PI, 2) * numofDOFs);

    for (int i = 0; i < tree->size(); i++) {
        Node* neighbor = (*tree)[i];
        double* neighborJoint = neighbor->joint;
        double currNeighborDistance = 0;
        for (int j = 0; j < numofDOFs; j++) {
            currNeighborDistance += pow(fabs(neighborJoint[j] - currJoint[j]), 2);
        }
        //printf("***neighborJoint is , [%f, %f, %f, %f, %f]\n",
        //            neighborJoint[0], neighborJoint[1], neighborJoint[2], neighborJoint[3], neighborJoint[4]);
        if(currNeighborDistance < closestNeighborDistance) {
            *closestNeighbor = neighbor;
            closestNeighborDistance = currNeighborDistance;
        }
        if (currNeighborDistance <= radius) {
            nearNodes->push_back(neighbor);
            nearNodeDistances->push_back(sqrt(currNeighborDistance));
        }
    }
    return sqrt(closestNeighborDistance);
}

static int isJointTransitionValid(double distance, double discretizationStep, int numofDOFs, double* currJoint, double* closestNeighbor,
        double* worldMap, int x_size, int y_size) {
    double* tempJoint = (double*) malloc(numofDOFs * sizeof(double));
    int numSteps = ((int) (distance/discretizationStep));
    for (int i = 1; i <= numSteps; i++) {
        for (int j = 0; j < numofDOFs; j++) {
            tempJoint[j] = closestNeighbor[j] + (i * discretizationStep) * ((currJoint[j] - closestNeighbor[j])/distance);
        }
        if (!IsValidArmConfiguration(tempJoint, numofDOFs, worldMap, x_size, y_size)) {
            return 0;
        }
    }
    return 1;
}

static void generateRandomJoint(double** joint, int numofDOFs) {
    for (int i = 0; i < numofDOFs; i++) {
        (*joint)[i] = (rand() / (RAND_MAX/(2 * PI )));
    }
}

static double getPlanQuality(double*** plan, int* planlength, int numofDOFs) {
    double distance = 0;
    for (int i = 0; i < *planlength - 1; i++) {
        double* currPlan = (*plan)[i];
        double* nextPlan = (*plan)[i+1];
        double currDistance = 0;
        for (int j = 0; j < numofDOFs; j++) {
            currDistance += pow(fabs(currPlan[j] - nextPlan[j]), 2);
        }
        distance += sqrt(currDistance);
    }
    return distance;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                          DEFAULT PLANNER FUNCTION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//

static void planner(
			double* map,
			int x_size,
			int y_size,
			double* armstart_anglesV_rad,
			double* armgoal_anglesV_rad,
            int numofDOFs,
            double*** plan,
            int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
		
    //for now just do straight interpolation between start and goal checking for the validity of samples

    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf) {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;
    
    return;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              RRT IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//

static ExperimentResult plannerRRT(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
	clock_t start = clock();
	//no plan by default
	*plan = NULL;
	*planlength = 0;

    int discretizationFactor = getAngleDiscretizationFactor(numofDOFs);
    double discretizationStep = (2 * PI)/discretizationFactor;
    double epsilon = PI/4;
    //printf("Discretization factor is %d and epsilon is %f\n", discretizationFactor, epsilon);

	Node* startNode = (Node*) malloc(sizeof(Node));
    double* startJoint = (double*) malloc(numofDOFs * sizeof(double));
    for (int i = 0; i < numofDOFs; i++) {
        startJoint[i] = armstart_anglesV_rad[i];
    }
    startNode->joint = startJoint;
    startNode->parent = 0;
    startNode->nodeNum = 1;
    vector<Node*>* nodes = new vector<Node*>();
    nodes->push_back(startNode);
    //printf("Created startTree and added startNode to it.\n");

    double* currJoint;
    Node* closestNeighbor;
    int isGoalJoint = 0;
    while (1) {
        // if (((clock() - start ) / (double) CLOCKS_PER_SEC) > TIMELIMIT) {
        //     ExperimentResult result;
        //     result.planningTime = -1;
        //     return result;
        // }
        if (rand() % 2 == 1) {
            currJoint = (double*) malloc(numofDOFs * sizeof(double));
            for (int i = 0; i < numofDOFs; i++) {
                currJoint[i] = armgoal_anglesV_rad[i];
            }
            
            isGoalJoint = 1;
            
        } else {
            currJoint = (double*) malloc(numofDOFs * sizeof(double));
            generateRandomJoint(&currJoint, numofDOFs);
            isGoalJoint = 0;
        }
        if(!IsValidArmConfiguration(currJoint, numofDOFs, map, x_size, y_size))
            continue;

        // Calculate closest neighbor
        double closestNeighborDistance = getClosestNeighborFromTree(currJoint, nodes, numofDOFs, &closestNeighbor);
        //printf("closestNeighbor to currNode is, [%f, %f, %f, %f, %f]\n",
        //            closestNeighbor[0], closestNeighbor[1], closestNeighbor[2], closestNeighbor[3], closestNeighbor[4]);
        //printf("Distance between closestNeighbor and currNode is %f\n", closestNeighborDistance);
        //printf("(%lf)/n",closestNeighborDistance);

        if (closestNeighborDistance > epsilon) {
            isGoalJoint = 0;
            for (int j = 0; j < numofDOFs; j++) {
                currJoint[j] = closestNeighbor->joint[j] + epsilon * ((currJoint[j] - closestNeighbor->joint[j])/closestNeighborDistance);
            }
            closestNeighborDistance = epsilon;
            //printf("Distance was greater than epsilon(%f), so updated currJoint to [%f, %f, %f, %f, %f]\n",
            //        epsilon, currJoint[0], currJoint[1], currJoint[2], currJoint[3], currJoint[4]);
        }
        
        int jointTransitionValid = isJointTransitionValid(closestNeighborDistance, discretizationStep, numofDOFs,
                currJoint, closestNeighbor->joint, map, x_size, y_size);
        
        // cout << jointTransitionValid << endl;
        if (jointTransitionValid) {
            Node* currNode = (Node*) malloc(sizeof(Node));
            currNode->joint = currJoint;
            currNode->parent = closestNeighbor;
            currNode->nodeNum = closestNeighbor->nodeNum + 1;
            nodes->push_back(currNode);
            
            if (isGoalJoint) {
                
                // printf("Reached goalJoint -- building plan of length %d.\n", 1);
              
                *plan = (double**) malloc(currNode->nodeNum * sizeof(double*));
                *planlength = currNode->nodeNum;

                for (int i = *planlength - 1; i >= 0; i--) {
                    (*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
                    for(int j = 0; j < numofDOFs; j++){
                        (*plan)[i][j] = currNode->joint[j];
                    }
                    currNode = currNode->parent;
                }
                ExperimentResult result;
                result.planningTime = (clock() - start ) / (double) CLOCKS_PER_SEC;
                result.numNodes = nodes->size();
                result.planLength = *planlength;
                result.planQuality = getPlanQuality(plan, planlength, numofDOFs);
                for(int i = 0; i < nodes->size(); i++) {
                    free((*nodes)[i]->joint);
                    free((*nodes)[i]);
                }
                
                //printf("Path has %d nodes and %d total nodes were generated in %f seconds with planQuality %f.\n", result.planLength, result.numNodes, result.planningTime, result.planQuality);
                return result;
            }
        }
    }
	
}


//*******************************************************************************************************************//
//                                                                                                                   //
//                                           RRT STAR IMPLEMENTATION                                                 //
//                                                                                                                   //
//*******************************************************************************************************************//
static double getRRTStarRadius(int numVertices, int numofDOFs, double epsilon) {
    double calcRad = pow((1000 * log(numVertices)/numVertices), (1.0/numofDOFs));
    return min(calcRad, epsilon);
}

static ExperimentResult plannerRRTStar(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
	clock_t start = clock();

	//no plan by default
	*plan = NULL;
	*planlength = 0;

    int discretizationFactor = getAngleDiscretizationFactor(numofDOFs);
    double discretizationStep = (2 * PI)/discretizationFactor;
    double epsilon = PI/4;

	Node* startNode = (Node*) malloc(sizeof(Node));
    double* startJoint = (double*) malloc(numofDOFs * sizeof(double));
    for (int i = 0; i < numofDOFs; i++) {
        startJoint[i] = armstart_anglesV_rad[i];
    }
    startNode->joint = startJoint;
    startNode->parent = 0;
    startNode->nodeNum = 1;
    startNode->cost = 0;
    vector<Node*>* nodes = new vector<Node*>();
    nodes->push_back(startNode);

    double* currJoint;
    Node* closestNeighbor;
    Node* goalNode;
    int isGoalJoint = 0;
    int numAfterGoal = -1;
    while (1) {
        // if (((clock() - start ) / (double) CLOCKS_PER_SEC) > TIMELIMIT) {
        //     ExperimentResult result;
        //     result.planningTime = -1;
        //     return result;
        // }
        if (rand() % 2 == 1 && numAfterGoal == -1) {
            currJoint = (double*) malloc(numofDOFs * sizeof(double));
            for (int i = 0; i < numofDOFs; i++) {
                currJoint[i] = armgoal_anglesV_rad[i];
            }
            isGoalJoint = 1;
        } else {
            currJoint = (double*) malloc(numofDOFs * sizeof(double));
            generateRandomJoint(&currJoint, numofDOFs);
            isGoalJoint = 0;
        }
        if(!IsValidArmConfiguration(currJoint, numofDOFs, map, x_size, y_size)) {
            free(currJoint);
            continue;
        }

        // Calculate closest neighbor
        double radius = getRRTStarRadius(nodes->size(), numofDOFs, epsilon);
        vector<Node*>* nearNodes = new vector<Node*>();
        vector<double>* nearNodeDistances = new vector<double>();

        double closestNeighborDistance = getClosestNeighborFromTreeAndNearNodes(
                currJoint, nodes, numofDOFs, &closestNeighbor, nearNodes, nearNodeDistances, radius);

        // printf("Radius = %f, Num nearest nodes = %d, total num nodes = %d\n", radius, nearNodes->size(), nodes->size());

        if (closestNeighborDistance > epsilon) {
            isGoalJoint = 0;
            for (int j = 0; j < numofDOFs; j++) {
                currJoint[j] = closestNeighbor->joint[j] + epsilon * ((currJoint[j] - closestNeighbor->joint[j])/closestNeighborDistance);
            }
            closestNeighborDistance = epsilon;
            //printf("Distance was greater than epsilon(%f), so updated currJoint to [%f, %f, %f, %f, %f]\n",
            //        epsilon, currJoint[0], currJoint[1], currJoint[2], currJoint[3], currJoint[4]);
        }

        int jointTransitionValid = isJointTransitionValid(closestNeighborDistance, discretizationStep, numofDOFs,
                currJoint, closestNeighbor->joint, map, x_size, y_size);

        if (jointTransitionValid) {

            vector<int>* nearNodeObstacleFree = new vector<int>();

            Node* minNode = closestNeighbor;
            double minCost = closestNeighbor->cost + closestNeighborDistance;
            for (int i = 0; i < nearNodes->size(); i++) {
                int jointTransitionValid = isJointTransitionValid((*nearNodeDistances)[i], discretizationStep, numofDOFs,
                    currJoint, (*nearNodes)[i]->joint, map, x_size, y_size);
                
                nearNodeObstacleFree->push_back(jointTransitionValid);
                if (jointTransitionValid) {
                    double currCost = (*nearNodes)[i]->cost + (*nearNodeDistances)[i];
                    if (currCost < minCost) {
                        minNode = (*nearNodes)[i];
                        minCost = currCost;
                    }
                }
            }
            
            Node* currNode = (Node*) malloc(sizeof(Node));
            currNode->joint = currJoint;
            currNode->parent = minNode;
            currNode->nodeNum = minNode->nodeNum + 1;
            currNode->cost = minCost;
            nodes->push_back(currNode);

            for (int i = 0; i < nearNodes->size(); i++) {
                if ((*nearNodes)[i] == minNode)
                    continue;
                double currCost = currNode->cost + (*nearNodeDistances)[i];
                if ((*nearNodeObstacleFree)[i] && currCost < (*nearNodes)[i]->cost) {
                    (*nearNodes)[i]->cost = currCost;
                    (*nearNodes)[i]->parent = currNode;
                    (*nearNodes)[i]->nodeNum = currNode->nodeNum + 1;
                }
            }
            
            if (numAfterGoal > 0){
                
                numAfterGoal--;
            }
            if (numAfterGoal == 0){
                
                break;
            }
            if (isGoalJoint) {
                numAfterGoal = 1000; // Start the countdown!
                
                //printf("Reached goalJoint -- expanding %d more nodes to improve path quality.\n", numAfterGoal);
                goalNode = currNode;
            }
            
        }
        
    }
    Node* tempNode = goalNode;
    int planLength = 1;
    while(tempNode->parent != 0) {
        planLength++;
        tempNode = tempNode->parent;
    }
    
    //printf("Reached goalJoint -- building plan of length %d.\n", goalNode->nodeNum);
    *plan = (double**) malloc(planLength * sizeof(double*));
    *planlength = planLength;

    //printf("startNode = [%f, %f, %f, %f, %f]\n",
    //        startNode->joint[0], startNode->joint[1], startNode->joint[2],
    //        startNode->joint[3], startNode->joint[4]);
    for (int i = *planlength - 1; i >= 0; i--) {
        //printf("plan[%d] = [%f, %f, %f, %f, %f]\n",
        //    i, goalNode->joint[0], goalNode->joint[1], goalNode->joint[2],
        //        goalNode->joint[3], goalNode->joint[4]);
        (*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
        for(int j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = goalNode->joint[j];
        }
        goalNode = goalNode->parent;
    }
    ExperimentResult result;
    result.planningTime = (clock() - start ) / (double) CLOCKS_PER_SEC;
    result.numNodes = nodes->size();
    result.planLength = *planlength;
    result.planQuality = getPlanQuality(plan, planlength, numofDOFs);
    for(int i = 0; i < nodes->size(); i++) {
        free((*nodes)[i]->joint);
        free((*nodes)[i]);
    }
    //printf("Path has %d nodes and %d total nodes were generated in %f seconds with planQuality %f.\n", result.planLength, result.numNodes, result.planningTime, result.planQuality);
    return result;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                              PRM IMPLEMENTATION                                                   //
//                                                                                                                   //
//*******************************************************************************************************************//
struct PRMNode {
    double* joint;
    vector<PRMNode*>* neighbors;
    int connectedToStart;
    int connectedToGoal;

    PRMNode* neighborToStart;
    int nodeNum;
};

static void getNearPRMNodes(double* joint, vector<PRMNode*>* nodes, vector<PRMNode*>* nearNodes, vector<double>* nearNodeDistances, double radius, int numofDOFs) {
    radius = pow(radius, 2);
    for (int i = 0; i < nodes->size(); i++) {
        PRMNode* neighbor = (*nodes)[i];
        double* neighborJoint = neighbor->joint;
        double currNeighborDistance = 0;
        for (int j = 0; j < numofDOFs; j++) {
            currNeighborDistance += pow(fabs(neighborJoint[j] - joint[j]), 2);
        }
        
        if (currNeighborDistance <= radius) {
            nearNodes->push_back(neighbor);
            //printf("Added neighbor [%f, %f, %f, %f, %f]\n",
            //    neighbor->joint[0], neighbor->joint[1], neighbor->joint[2], neighbor->joint[3], neighbor->joint[4]);
            nearNodeDistances->push_back(sqrt(currNeighborDistance));
        }
    }
}

static int propagateStartGoalConnected(PRMNode* node) {
    int startGoalConnected = (node->connectedToStart && node->connectedToGoal);
    PRMNode* currNeighbor;
    int currNeighborConnectedToStart;
    int currNeighborConnectedToGoal;
    for (int i = 0; i < node->neighbors->size(); i++) {
        currNeighbor = (*(node->neighbors))[i];
        currNeighborConnectedToStart = currNeighbor->connectedToStart;
        currNeighborConnectedToGoal = currNeighbor->connectedToGoal;
        currNeighbor->connectedToStart = (currNeighbor->connectedToStart || node->connectedToStart);
        currNeighbor->connectedToGoal = (currNeighbor->connectedToGoal || node->connectedToGoal);
        if (currNeighborConnectedToStart != node->connectedToStart || currNeighborConnectedToGoal != node->connectedToGoal) {
            if(propagateStartGoalConnected(currNeighbor)) {
                startGoalConnected = 1;
            }
        }
    }
    return startGoalConnected;
}

static ExperimentResult plannerPRM(
    double *map,
    int x_size,
    int y_size,
    double *armstart_anglesV_rad,
    double *armgoal_anglesV_rad,
    int numofDOFs,
    double ***plan,
    int *planlength)
{
    /* TODO: Replace with your implementation */
    clock_t start = clock();

	//no plan by default
	*plan = NULL;
	*planlength = 0;

    int discretizationFactor = getAngleDiscretizationFactor(numofDOFs);
    double discretizationStep = (2 * PI)/discretizationFactor;
    double epsilon = PI/4;

    vector<PRMNode*>* nodes = new vector<PRMNode*>();

    
    PRMNode* startNode = (PRMNode*) malloc(sizeof(PRMNode));
    double* startJoint = (double*) malloc(numofDOFs * sizeof(double));
    for (int i = 0; i < numofDOFs; i++) {
        startJoint[i] = armstart_anglesV_rad[i];
    }
    startNode->joint = startJoint;
    startNode->connectedToStart = 1;
    startNode->connectedToGoal = 0;
    startNode->neighbors = new vector<PRMNode*>();
    startNode->nodeNum = 1;
    nodes->push_back(startNode);

    PRMNode* goalNode = (PRMNode*) malloc(sizeof(PRMNode));
    double* goalJoint = (double*) malloc(numofDOFs * sizeof(double));
    for (int i = 0; i < numofDOFs; i++) {
        goalJoint[i] = armgoal_anglesV_rad[i];
    }
    goalNode->joint = goalJoint;
    goalNode->connectedToStart = 0;
    goalNode->connectedToGoal = 1;
    goalNode->neighbors = new vector<PRMNode*>();
    goalNode->nodeNum = -1;
    nodes->push_back(goalNode);

    double* currJoint;
    while(1) {
        if (((clock() - start ) / (double) CLOCKS_PER_SEC) > TIMELIMIT) {
            ExperimentResult result;
            result.planningTime = -1;
            return result;
        }
        currJoint = (double*) malloc(numofDOFs * sizeof(double));
        generateRandomJoint(&currJoint, numofDOFs);
        if(!IsValidArmConfiguration(currJoint, numofDOFs, map, x_size, y_size))
            continue;
        //printf("currJoint is , [%f, %f, %f, %f, %f]\n",
        //    currJoint[0], currJoint[1], currJoint[2], currJoint[3], currJoint[4]);

        vector<PRMNode*>* nearNodes = new vector<PRMNode*>();
        vector<double>* nearNodeDistances = new vector<double>();
        double radius = getRRTStarRadius(nodes->size(), numofDOFs, epsilon);
        getNearPRMNodes(currJoint, nodes, nearNodes, nearNodeDistances, radius, numofDOFs);
        // printf("Radius = %f, Num nearest nodes = %d, total num nodes = %d\n", radius, nearNodes->size(), nodes->size());

        PRMNode* currNode = (PRMNode*) malloc(sizeof(PRMNode));
        currNode->joint = currJoint;
        currNode->connectedToStart = 0;
        currNode->connectedToGoal = 0;
        currNode->neighbors = new vector<PRMNode*>();
        currNode->nodeNum = -1;
        nodes->push_back(currNode);
        // printf("currNode added to nodes\n");

        PRMNode* neighbor;
        double neighborDistance;
        for(int i = 0; i < nearNodes->size(); i++) {
            neighbor = (*nearNodes)[i];
            neighborDistance = (*nearNodeDistances)[i];
            //printf("neighbor is , [%f, %f, %f, %f, %f]\n",
            //    neighbor->joint[0], neighbor->joint[1], neighbor->joint[2], neighbor->joint[3], neighbor->joint[4]);
            if(!isJointTransitionValid(neighborDistance, discretizationStep, numofDOFs, currNode->joint, neighbor->joint, map, x_size, y_size))
                continue;
            neighbor->neighbors->push_back(currNode);
            currNode->neighbors->push_back(neighbor);
        }
        if (propagateStartGoalConnected(currNode))
            break;
    }
    // printf("Start goal connected!  %d nodes expanded!", nodes->size());
    
    queue<PRMNode*> prmQueue;
    prmQueue.push(startNode);
            
    PRMNode* currNode;
    while(prmQueue.size() != 0) {
        currNode = prmQueue.front();
        prmQueue.pop();
        if (currNode == goalNode) {
            //printf("Found path to goalNode!\n");
            break;
        }
        PRMNode* neighbor;
        for(int i = 0; i < currNode->neighbors->size(); i++) {
            neighbor = (*(currNode->neighbors))[i];
            if (neighbor->nodeNum == -1) {
                neighbor->neighborToStart = currNode;
                neighbor->nodeNum = currNode->nodeNum + 1;
                prmQueue.push(neighbor);
            }
        }
    }

    *plan = (double**) malloc(currNode->nodeNum * sizeof(double*));
    *planlength = currNode->nodeNum;

    for (int i = *planlength - 1; i >= 0; i--) {
        (*plan)[i] = (double*) malloc(numofDOFs * sizeof(double));
        for(int j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = currNode->joint[j];
        }
        currNode = currNode->neighborToStart;
    }
    
    ExperimentResult result;
    result.planningTime = (clock() - start ) / (double) CLOCKS_PER_SEC;
    result.numNodes = nodes->size();
    result.planLength = *planlength;
    result.planQuality = getPlanQuality(plan, planlength, numofDOFs);
    //printf("Path has %d nodes and %d total nodes were generated in %f seconds with planQuality %f.\n", result.planLength, result.numNodes, result.planningTime, result.planQuality);
    for(int i = 0; i < nodes->size(); i++) {
        free((*nodes)[i]->joint);
        free((*nodes)[i]);
    }
    return result;
}

//*******************************************************************************************************************//
//                                                                                                                   //
//                                                MAIN FUNCTION                                                      //
//                                                                                                                   //
//*******************************************************************************************************************//

/** Your final solution will be graded by an grading script which will
 * send the default 6 arguments:
 *    map, numOfDOFs, commaSeparatedStartPos, commaSeparatedGoalPos, 
 *    whichPlanner, outputFilePath
 * An example run after compiling and getting the planner.out executable
 * >> ./planner.out map1.txt 5 1.57,0.78,1.57,0.78,1.57 0.392,2.35,3.14,2.82,4.71 0 output.txt
 * See the hw handout for full information.
 * If you modify this for testing (e.g. to try out different hyper-parameters),
 * make sure it can run with the original 6 commands.
 * Programs that do not will automatically get a 0.
 * */
int main(int argc, char** argv) {
	double* map;
	int x_size, y_size;

	tie(map, x_size, y_size) = loadMap(argv[1]);
	const int numOfDOFs = std::stoi(argv[2]);
	double* startPos = doubleArrayFromString(argv[3]);
	double* goalPos = doubleArrayFromString(argv[4]);
	int whichPlanner = std::stoi(argv[5]);
	string outputFile = argv[6];

	if(!IsValidArmConfiguration(startPos, numOfDOFs, map, x_size, y_size)||
			!IsValidArmConfiguration(goalPos, numOfDOFs, map, x_size, y_size)) {
		throw runtime_error("Invalid start or goal configuration!\n");
	}

	///////////////////////////////////////
	//// Feel free to modify anything below. Be careful modifying anything above.

	double** plan = NULL;
	int planlength = 0;

    // Call the corresponding planner function
    if (whichPlanner == PRM)
    {
        plannerPRM(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRT)
    {
        plannerRRT(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == RRTSTAR)
    {
        plannerRRTStar(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }
    else if (whichPlanner == ALL)
    {
        printf("Running All Planners\n");
        int numIterations = 20;
        
        double rrtPlanningTime = 0;
        int rrtNumNodes = 0;
        double rrtPlanQuality = 0;

        double rrtStarPlanningTime = 0;
        int rrtStarNumNodes = 0;
        double rrtStarPlanQuality = 0;

        double prmPlanningTime = 0;
        int prmNumNodes = 0;
        double prmPlanQuality = 0;
        
        srand(time(NULL));
        int i = 1;
        while (i <= numIterations) {
            printf("Iteration %d\n", i);
            double* start = (double*) malloc(numOfDOFs * sizeof(double));
            while(1) {
                generateRandomJoint(&start, numOfDOFs);
                if(IsValidArmConfiguration(start, numOfDOFs, map, x_size, y_size))
                    break;
            }
            double* goal = (double*) malloc(numOfDOFs * sizeof(double));
            while(1) {
                generateRandomJoint(&goal, numOfDOFs);
                if(IsValidArmConfiguration(goal, numOfDOFs, map, x_size, y_size))
                    break;
            }
            printf("start is  [");
            for (int i = 0; i < numOfDOFs-1; i++) {
                printf("%f, ", start[i]);
            }
            printf("%f]\n", start[numOfDOFs-1]);

            printf("goal is  [");
            for (int i = 0; i < numOfDOFs-1; i++) {
                printf("%f, ", goal[i]);
            }
            printf("%f]\n", goal[numOfDOFs-1]);

            printf("Running RRT\n");
            ExperimentResult rrtResult = plannerRRT(map,x_size,y_size, start, goal, numOfDOFs, &plan, &planlength);
            // if (rrtResult.planningTime == -1) {
            //     printf("RRT took more than %d seconds, retrying iteration...\n\n", TIMELIMIT);
            //     continue;
            // }
            printf("Running RRTStar\n");
            ExperimentResult rrtStarResult = plannerRRTStar(map,x_size,y_size, start, goal, numOfDOFs, &plan, &planlength);
            // if (rrtStarResult.planningTime == -1) {
            //     printf("RRTStar took more than %d seconds, retrying iteration...\n\n", TIMELIMIT);
            //     continue;
            // }
            printf("Running PRM\n");
            ExperimentResult prmResult = plannerPRM(map,x_size,y_size, start, goal, numOfDOFs, &plan, &planlength);
            // if (prmResult.planningTime == -1) {
            //     printf("PRM took more than %d seconds, retrying iteration...\n\n", TIMELIMIT);
            //     continue;
            // }

            printf("Algorithm | planningTime | numNodes | planLength | planQuality\n");
            rrtPlanningTime += rrtResult.planningTime;
            rrtNumNodes += rrtResult.numNodes;
            rrtPlanQuality += rrtResult.planQuality;
            printf("RRT | %f | %d | %d | %f\n", rrtResult.planningTime, rrtResult.numNodes, rrtResult.planLength, rrtResult.planQuality);
            
            rrtStarPlanningTime += rrtStarResult.planningTime;
            rrtStarNumNodes += rrtStarResult.numNodes;
            rrtStarPlanQuality += rrtStarResult.planQuality;
            printf("RRTStar | %f | %d | %d | %f\n", rrtStarResult.planningTime, rrtStarResult.numNodes, rrtStarResult.planLength, rrtStarResult.planQuality);

            prmPlanningTime += prmResult.planningTime;
            prmNumNodes += prmResult.numNodes;
            prmPlanQuality += prmResult.planQuality;
            printf("PRM | %f | %d | %d | %f\n", prmResult.planningTime, prmResult.numNodes, prmResult.planLength, prmResult.planQuality);
            printf("-----------------------------------\n\n");

            i++;
        }
        rrtPlanningTime /= numIterations;
        rrtNumNodes /= numIterations;
        rrtPlanQuality /= numIterations;
        rrtStarPlanningTime /= numIterations;
        rrtStarNumNodes /= numIterations;
        rrtStarPlanQuality /= numIterations;
        prmPlanningTime /= numIterations;
        prmNumNodes /= numIterations;
        prmPlanQuality /= numIterations;
        printf("Final Results!\n");
        printf("Algorithm | avgPlanningTime | avgNumNodes |avgPlanQuality\n");
        printf("RRT | %f | %d | %f\n", rrtPlanningTime, rrtNumNodes, rrtPlanQuality);
        printf("RRTStar | %f | %d | %f\n", rrtStarPlanningTime, rrtStarNumNodes, rrtStarPlanQuality);
        printf("PRM | %f | %d | %f\n", prmPlanningTime, prmNumNodes, prmPlanQuality);
        printf("-----------------------------------\n\n");
    }
    else
    {
        planner(map, x_size, y_size, startPos, goalPos, numOfDOFs, &plan, &planlength);
    }

	//// Feel free to modify anything above.
	//// If you modify something below, please change it back afterwards as my 
	//// grading script will not work and you will recieve a 0.
	///////////////////////////////////////

    // Your solution's path should start with startPos and end with goalPos
    if (!equalDoubleArrays(plan[0], startPos, numOfDOFs) || 
    	!equalDoubleArrays(plan[planlength-1], goalPos, numOfDOFs)) {
		throw std::runtime_error("Start or goal position not matching");
	}

	/** Saves the solution to output file
	 * Do not modify the output log file output format as it is required for visualization
	 * and for grading.
	 */
	std::ofstream m_log_fstream;
	m_log_fstream.open(outputFile, std::ios::trunc); // Creates new or replaces existing file
	if (!m_log_fstream.is_open()) {
		throw std::runtime_error("Cannot open file");
	}
	m_log_fstream << argv[1] << endl; // Write out map name first
	/// Then write out all the joint angles in the plan sequentially
	for (int i = 0; i < planlength; ++i) {
		for (int k = 0; k < numOfDOFs; ++k) {
			m_log_fstream << plan[i][k] << ",";
		}
		m_log_fstream << endl;
	}
}
