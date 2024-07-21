/*=================================================================
 *
 * planner.cpp
 *
 *=================================================================*/
#include "planner.h"
#include <cmath>
#include <iostream>
using namespace std;
#include <queue>
#include <vector>
#include <unordered_map>
#include <chrono>
using namespace std::chrono;

#define GET_MAP_INDEX(X, Y, XSIZE, YSIZE) ((Y - 1) * XSIZE + (X - 1))

#if !defined(MAX_VAL)
#define MAX_VAL(A, B) ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN_VAL)
#define MIN_VAL(A, B) ((A) < (B) ? (A) : (B))
#endif

#define NUM_DIRS 8

class NodeInfo {
public:
    int posX, posY;
    double costFromStart;
    int costValue;
    NodeInfo* parentNode;
    int timeStep;
    double heuristicCost;

    double totalCost() const {
        return costFromStart + heuristicCost;
    }
};

struct CompareNodeInfo {
    bool operator()(const NodeInfo* lhs, const NodeInfo* rhs) const {
        return lhs->totalCost() > rhs->totalCost();
    }
};

void planner(
    int* map,
    int collision_thresh,
    int x_size,
    int y_size,
    int robotposeX,
    int robotposeY,
    int target_steps,
    int* target_traj,
    int targetposeX,
    int targetposeY,
    int curr_time,
    int* action_ptr
) {
    static const int directionX[NUM_DIRS] = {-1, -1, -1, 0, 0, 1, 1, 1};
    static const int directionY[NUM_DIRS] = {-1, 0, 1, -1, 1, -1, 0, 1};

    priority_queue<NodeInfo*, vector<NodeInfo*>, CompareNodeInfo> openList;
    unordered_map<int, NodeInfo*> openListMap;
    unordered_map<int, bool> closedListMap;
    int neighborX, neighborY, neighborKey, currentKey, key;
    static unordered_map<int, int> actionMapX;
    static unordered_map<int, int> actionMapY;

    NodeInfo* startNode = new NodeInfo();
    startNode->posX = robotposeX;
    startNode->posY = robotposeY;
    startNode->costFromStart = 0;
    startNode->heuristicCost = abs(robotposeX - targetposeX) + abs(robotposeY - targetposeY); // Manhattan distance heuristic
    openList.push(startNode);
    openListMap[GET_MAP_INDEX(startNode->posX, startNode->posY, x_size, y_size)] = startNode;

    if (curr_time == 0) {
        while (!openList.empty()) {
            NodeInfo* currentNode = openList.top();
            openList.pop();
            currentKey = GET_MAP_INDEX(currentNode->posX, currentNode->posY, x_size, y_size);

            if (closedListMap[currentKey]) {
                continue;
            }

            closedListMap[currentKey] = true;

            for (int dir = 0; dir < NUM_DIRS; dir++) {
                neighborX = currentNode->posX + directionX[dir];
                neighborY = currentNode->posY + directionY[dir];

                if ((neighborX >= 1 && neighborX <= x_size && neighborY >= 1 && neighborY <= y_size) &&
                    ((int)map[GET_MAP_INDEX(neighborX, neighborY, x_size, y_size)] < collision_thresh)) {

                    neighborKey = GET_MAP_INDEX(neighborX, neighborY, x_size, y_size);

                    if (closedListMap[neighborKey]) {
                        continue;
                    }

                    if (openListMap[neighborKey]) {
                        NodeInfo* neighborNode = openListMap[neighborKey];
                        double newCostFromStart = currentNode->costFromStart + neighborNode->costValue;
                        if (newCostFromStart < neighborNode->costFromStart) {
                            neighborNode->costFromStart = newCostFromStart;
                            neighborNode->parentNode = currentNode;
                            neighborNode->timeStep = currentNode->timeStep + 1;
                            neighborNode->heuristicCost = abs(neighborX - targetposeX) + abs(neighborY - targetposeY); // Manhattan distance heuristic
                            openList.push(neighborNode);
                        }
                    } else {
                        NodeInfo* neighborNode = new NodeInfo();
                        neighborNode->posX = neighborX;
                        neighborNode->posY = neighborY;
                        neighborNode->costValue = (int)map[GET_MAP_INDEX(neighborX, neighborY, x_size, y_size)];
                        neighborNode->costFromStart = currentNode->costFromStart + neighborNode->costValue;
                        neighborNode->parentNode = currentNode;
                        neighborNode->timeStep = currentNode->timeStep + 1;
                        neighborNode->heuristicCost = abs(neighborX - targetposeX) + abs(neighborY - targetposeY); // Manhattan distance heuristic
                        openListMap[neighborKey] = neighborNode;
                        openList.push(neighborNode);
                    }
                }
            }
        }

        int minCost = 0;
        NodeInfo* goalNode = nullptr;

        for (int i = 0; i < target_steps; i++) {
            currentKey = GET_MAP_INDEX(target_traj[i], target_traj[i + target_steps], x_size, y_size);
            NodeInfo* currentNode = openListMap[currentKey];

            if (i > currentNode->timeStep) {
                int newCost = currentNode->costFromStart + (i - currentNode->timeStep) * (currentNode->costValue);
                if (minCost == 0 || newCost < minCost) {
                    minCost = newCost;
                    goalNode = currentNode;
                }
            }
        }

        currentKey = GET_MAP_INDEX(goalNode->posX, goalNode->posY, x_size, y_size);
        NodeInfo* currentNode = goalNode;
        actionMapX[currentKey] = currentNode->posX;
        actionMapY[currentKey] = currentNode->posY;
        NodeInfo* previousNode = currentNode->parentNode;

        while (previousNode != nullptr) {
            key = GET_MAP_INDEX(previousNode->posX, previousNode->posY, x_size, y_size);
            actionMapX[key] = currentNode->posX;
            actionMapY[key] = currentNode->posY;
            currentNode = previousNode;
            previousNode = currentNode->parentNode;
        }
    }

    key = GET_MAP_INDEX(robotposeX, robotposeY, x_size, y_size);

    action_ptr[0] = actionMapX[key];
    action_ptr[1] = actionMapY[key];

    return;
}