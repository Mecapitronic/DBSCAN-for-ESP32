/*
    DBSCAN unsupervised classification algorithm
    inspired by https://penseeartificielle.fr/clustering-avec-lalgorithme-dbscan/
    https://openclassrooms.com/fr/courses/4379436-explorez-vos-donnees-avec-des-algorithmes-non-supervises/4379571-partitionnez-vos-donnees-avec-dbscan

    (c) 2021 Lesept
    contact: lesept777@gmail.com

*/
#ifndef dbscan_h
#define dbscan_h

#include <Arduino.h>
#include <vector>
#include "Debugger.h"
#include "A010.h"
using namespace std;

#define CLOUD_SIZE PICTURE_SIZE
#define CLOUD_WiDTH PICTURE_RES

enum DISTANCE_TYPE
{
    NONE,
    EUCLIDEAN,   // Euclidean distance
    MINKOWSKI,   // Minkowski distance (is Euclidean if param = 2)
    MANHATTAN,   // Manhattan distance
    TCHEBYCHEV,  // Tchebychev distance
    CANBERRA,    // Canberra distance
    PROJECTION   // Projection on view axe
};

enum POINT_STATUS
{
    NOT_VISITED,
    VISITED,
    NOISE
};

class Dbscan
{
   public:
    struct ClusterPoint4D
    {
        Point4D point;
        POINT_STATUS type;
        uint8_t cluster;
    };

    struct ConfigDbscan
    {
        float epsilon;
        uint16_t minPts;
        DISTANCE_TYPE distanceType;
        float mink;
    };

   private:
    ConfigDbscan dbscanConfig = {0.0, 0, NONE, 0.0};

    uint16_t _nData = 0;
    uint16_t _resolution = 0;
    ClusterPoint4D _dataset[CLOUD_SIZE];

    uint16_t _nNoise = 0;
    uint16_t _nClusters = 0;
    vector<vector<uint16_t>> _clusters;

   private:
    Point4D computeCentroid(vector<uint16_t> const &);
    // vector<uint16_t> findNeighbours(uint16_t);
    vector<uint16_t> findClosestNeighbours(uint16_t);
    float computeTightness(vector<uint16_t> const &, Point4D const &);
    float getDistance(Point4D point1, Point4D point2, DISTANCE_TYPE = NONE);
    // int countNeighbours(ClusterPoint4D point1);
    bool isNeighbour(Point4D point1, Point4D point2);
    void enlargeCluster(vector<uint16_t>, vector<uint16_t> &);

   public:
    Dbscan(void);
    void Config(float, int, DISTANCE_TYPE, float = 1.0f);
    vector<vector<uint16_t>> Process(Point4D *);
    void displayStats();
};

#endif
