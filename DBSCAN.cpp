#include "DBSCAN.h"

/* Constructor
 */
Dbscan::Dbscan()
{
    _epsilon = 0;
    _minPts = 0;
    _distanceType = 0;
    _mink = 0;
}

/* Configuration
        Arguments :
                epsilon : maximum distance between neighbours
                minPts  : minimum number of points in a cluster
                type    : the type of distance
                mink    : the exponent in case of the Minkovski distance (optionnal)
*/
void Dbscan::Config(uint16_t size, float epsilon, int minPts, uint8_t distance, float mink)
{
    _nData = size;

    // Delete the array
    delete[] _dataset;

    // Resize the array
    _dataset = new Point3DCluster[_nData];

    Serial.printf(" %i ", _nData);

    _epsilon = epsilon;
    _minPts = minPts;
    _distanceType = distance;
    _mink = mink;

    Serial.printf(" %f ", _epsilon);
    Serial.printf(" %i ", _minPts);
    Serial.printf(" %u ", _distanceType);
    Serial.printf(" %f ", _mink);
    Serial.println();
}

/* Process the dataset */
vector<vector<uint16_t>> Dbscan::Process(Point3D *dataset)
{
    // Fill the array
    for (uint16_t i = 0; i < _nData; ++i)
    {
        _dataset[i].point.x = dataset[i].x;
        _dataset[i].point.y = dataset[i].y;
        _dataset[i].point.z = dataset[i].z;
        //_dataset[i].point = dataset[i];
        _dataset[i].type = NOT_VISITED;
        _dataset[i].cluster = 0;
    }

    // Average distance
    float averageDistance = 0.0f;
    for (uint16_t i = 0; i < _nData; ++i)
        for (uint16_t j = i + 1; j < _nData; ++j) averageDistance += distance(_dataset[i].point, _dataset[j].point);
    averageDistance /= ((_nData - 1) * _nData / 2);
    Serial.printf("Average distance : %f\n", averageDistance);

    // Process dataset
    vector<uint16_t> noise;
    for (uint16_t i = 0; i < _nData; ++i)
    {
        if (_dataset[i].type == NOT_VISITED)
        {
            _dataset[i].type = VISITED;
            vector<uint16_t> currentCluster;
            vector<uint16_t> neighbours = findNeighbours(i);
            // If the point has too few neighbours : set to noise
            if (neighbours.size() < _minPts)
            {
                _dataset[i].type = NOISE;
                noise.push_back(i);
                ++_nNoise;
                // Serial.println ("Noise!");
            }
            else
            {
                // create	 a cluster with this point
                currentCluster.push_back(i);
                enlargeCluster(neighbours, currentCluster);
                // Mark all points in the cluster as VISITED
                for (uint16_t j = 0; j < currentCluster.size(); ++j) _dataset[currentCluster[j]].type = VISITED;
                // Add current cluster to clusters list
                _clusters.push_back(currentCluster);
                ++_nClusters;
            }
        }
    }
    // Noise cluster is inserted at position 0 (first cluster) even if empty
    _clusters.insert(_clusters.begin(), noise);
    displayStats();
    return _clusters;
}

void Dbscan::displayStats()
{
    // Print statistics about the clusters
    vector<Point3D> centroid;
    vector<float> tightness;
    Serial.printf("Created %d clusters.\n", _nClusters);
    for (uint16_t i = 0; i < _nClusters; ++i)
    {
        Serial.printf("Cluster %d : %d points\n", i, _clusters[i + 1].size() - 1);

        // Centroid
        Point3D c;
        c = computeCentroid(_clusters[i + 1]);
        Serial.print("\tCentroid: ");
        Serial.printf("%f ", c.x);
        Serial.printf("%f ", c.y);
        Serial.printf("%f ", c.z);

        Serial.println();
        centroid.push_back(c);

        // Tightness (mean distance to centroid)
        float t = computeTightness(_clusters[i + 1], c);
        Serial.printf("\tTightness = %.3f\n", t);
        tightness.push_back(t);
    }

    // Separation of clusters (mean distance of centroids)
    float separation = 0.0f;
    float indexDB = 0.0;
    for (uint16_t i = 0; i < _nClusters; ++i)
    {
        for (uint16_t j = i + 1; j < _nClusters; ++j)
        {
            separation += distance(centroid[i], centroid[j]);
            float index = (tightness[i] + tightness[j]) / separation;
            indexDB = max(indexDB, index);
        }
    }
    separation = separation * 2.0f / _nClusters / (_nClusters - 1.0);
    Serial.printf("\nSeparation = %.3f", separation);

    // Davies-Bouldin index (max ratio tightness over separation)
    Serial.printf("\nDavies-Bouldin index = %.3f\n", indexDB);

    Serial.printf("%d noise points\n", _nNoise);
}

/* Compute the coordinates of the centroid of a cluster */
Dbscan::Point3D Dbscan::computeCentroid(vector<uint16_t> const &cluster)
{
    Point3D centroid;
    for (uint16_t j = 0; j < cluster.size(); ++j)
    {
        centroid.x += _dataset[cluster[j]].point.x / cluster.size();
        centroid.y += _dataset[cluster[j]].point.y / cluster.size();
        centroid.z += _dataset[cluster[j]].point.z / cluster.size();
    }
    return centroid;
}

/* Compute the tightness of a cluster */
float Dbscan::computeTightness(vector<uint16_t> const &cluster, Point3D const &centroid)
{
    float tightness = 0.0f;
    for (uint16_t j = 0; j < cluster.size(); ++j) tightness += distance(_dataset[cluster[j]].point, centroid) / cluster.size();
    return tightness;
}

/* Enlarge an existing cluster */
void Dbscan::enlargeCluster(vector<uint16_t> neighbours, vector<uint16_t> &currentCluster)
{
    uint16_t i = 0;
    while (i < neighbours.size())
    {
        uint16_t index = neighbours[i++];
        if (_dataset[index].type == NOT_VISITED)
        {
            vector<uint16_t> neighbours2 = findNeighbours(index);
            if (neighbours2.size() > _minPts)
            {
                // make union of both neighbourhoods
                for (uint16_t j = 0; j < neighbours2.size(); ++j)
                {
                    bool isInNeighbours = false;
                    for (uint16_t k = 0; k < neighbours.size(); ++k)
                    {
                        if (neighbours2[j] == neighbours[k])
                        {
                            isInNeighbours = true;
                            break;
                        }
                    }
                    if (!isInNeighbours)
                        neighbours.push_back(neighbours2[j]);
                }
            }
        }
        // add current point to current cluster is not already part of a cluster
        bool isInCluster = false;
        for (uint16_t j = 1; j < _nClusters; ++j)
            for (uint16_t k = 0; k < _clusters[j].size(); ++k)
                if (_clusters[j][k] == index)
                {
                    isInCluster = true;
                    break;
                }
        if (!isInCluster)
            currentCluster.push_back(index);
    }
}

/* Find the neighbours of a point in the dataset */
vector<uint16_t> Dbscan::findNeighbours(uint16_t n)
{
    vector<uint16_t> neighbours;
    for (uint16_t i = 0; i < _nData; ++i)
        if (isNeighbour(_dataset[n], _dataset[i]))
            neighbours.push_back(i);
    return neighbours;
}

/*
        Compute the distance between 2 vectors
*/
float Dbscan::distance(Point3D point1, Point3D point2)
{
    float distance = 0.0f;
    switch (_distanceType)
    {
        case EUCLIDIAN:
            distance = pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2) + pow(point1.z - point2.z, 2);
            distance = sqrt(distance);
            break;
            /* TODO im lazy ...
         case MINKOVSKI:
             for (uint8_t i = 0; i < vector1.size(); ++i) distance += pow(abs(vector1[i] - vector2[i]), _mink);
             distance = pow(distance, 1. / _mink);
             break;
         case MANHATTAN:
             for (uint8_t i = 0; i < vector1.size(); ++i) distance += abs(vector1[i] - vector2[i]);
             break;
         case CHEBYCHEV:
             for (uint8_t i = 0; i < vector1.size(); ++i) distance += max(distance, abs(vector1[i] - vector2[i]));
             break;
         case CANBERRA:
             for (uint8_t i = 0; i < vector1.size(); ++i) distance += abs(vector1[i] - vector2[i]) / (abs(vector1[i]) + abs(vector2[i]));
             break;*/
        default:
            Serial.println("Distance type problem !");
            distance = 2.0e10;
    }
    return distance;
}

int Dbscan::countNeighbours(Point3DCluster point1)
{
    int neighbours = 0.0;
    for (uint8_t i = 0; i < _nData; ++i)
        if (isNeighbour(point1, _dataset[i]))
            ++neighbours;
    return neighbours;
}

bool Dbscan::isNeighbour(Point3DCluster point1, Point3DCluster point2) { return (distance(point1.point, point2.point) <= _epsilon); }
