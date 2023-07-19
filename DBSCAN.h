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
    EUCLIDEAN,   // Euclidean distance
    MINKOWSKI,   // Minkowski distance (is Euclidean if param = 2)
    MANHATTAN,   // Manhattan distance
    TCHEBYCHEV,  // Tchebychev distance
    CANBERRA     // Canberra distance
};

enum TYPE
{
    NOT_VISITED,
    VISITED,
    NOISE
};

class Dbscan
{
   public:
    struct Point3D
    {
        int16_t x;
        int16_t y;
        int16_t z;
    };

    struct Point3DCluster
    {
        Point3D point;
        uint8_t type;
        uint8_t cluster;
    };

    struct ConfigDbscan
    {
        float epsilon;
        uint16_t minPts;
        uint8_t distanceType;
        float mink;
    };

   private:
    ConfigDbscan dbscanConfig = {0.0, 0, 0, 0.0};

    uint16_t _nNoise = 0;
    uint8_t _distanceType = 0;
    uint16_t _nData = 0;
    uint16_t _nClusters = 0;
    Point3DCluster *_dataset;
    vector<vector<uint16_t>> _clusters;

    Dbscan::Point3D computeCentroid(vector<uint16_t> const &);
    vector<uint16_t> findNeighbours(uint16_t);
    float computeTightness(vector<uint16_t> const &, Point3D const &);
    float distance(Point3D point1, Point3D point2);
    int countNeighbours(Point3DCluster point1);
    bool isNeighbour(Point3DCluster point1, Point3DCluster point2);
    void enlargeCluster(vector<uint16_t>, vector<uint16_t> &);

   public:
    Dbscan(void);
    void Config(uint16_t, float, int, DISTANCE_TYPE, float = 1.0f);
    vector<vector<uint16_t>> Process(Point3D *);
    void displayStats();
};

#endif
