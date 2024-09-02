#include <iostream>
#include <vector>
#include <cmath>
#include <map>

class Point
{
public:
    Point(double x, double y) : x(x), y(y), clusterID(0), isBorderPointFlg(0),nearestCorePointDistance(100){}
    Point(double x, double y, int id) : x(x), y(y), clusterID(id), isBorderPointFlg(0), nearestCorePointDistance(100){}

    double x, y;
    int clusterID;
    bool isBorderPointFlg;
    double nearestCorePointDistance;
};

class Grid
{
public:
    Grid(double size, double sp_lon, double sp_lat, int num_lon, int num_lat) : cellsize(size), startpoint_lon(sp_lon), startpoint_lat(sp_lat), cellnum_lon(num_lon), cellnum_lat(num_lat) {}
    double cellsize;
    double startpoint_lon;
    double startpoint_lat;
    int cellnum_lon;
    int cellnum_lat;
    std::map<int, std::vector<int>> GridMapTable;
    std::map<int, std::vector<int>> GridMapCoreTable;

    void RegisterPoint2Grid(std::vector<Point> &points);
    void DefineSearchingRange(int cellidx, std::vector<int> &neighborcells);
};

void Grid::RegisterPoint2Grid(std::vector<Point> &points)
{
    int lon_i;
    int lat_i;
    int cellID;

    std::vector<int> corevector;

    for (int i = 0; i < points.size(); i++)
    {
        lon_i = (int)floor((points[i].x - startpoint_lon) / Grid::cellsize);
        lat_i = (int)floor((points[i].y - startpoint_lat) / Grid::cellsize);

        if (lon_i >= 0 && lat_i >= 0 && lon_i < Grid::cellnum_lon && lat_i < Grid::cellnum_lat)
        {
            cellID = Grid::cellnum_lat * lon_i + lat_i;
            corevector.push_back(i);
        }
    }
    Grid::GridMapTable.emplace(cellID, corevector);
}

void Grid::DefineSearchingRange(int cellidx, std::vector<int> &neighborcells)
{
    neighborcells = {cellidx - Grid::cellnum_lat * 2 - 1, cellidx - Grid::cellnum_lat * 2, cellidx - Grid::cellnum_lat * 2 + 1,
                     cellidx - Grid::cellnum_lat - 2, cellidx - Grid::cellnum_lat - 1, cellidx - Grid::cellnum_lat, cellidx - Grid::cellnum_lat + 1, cellidx - Grid::cellnum_lat + 2,
                     cellidx - 2, cellidx - 1, cellidx + 1, cellidx + 2,
                     cellidx + Grid::cellnum_lat - 2, cellidx + Grid::cellnum_lat - 1, cellidx + Grid::cellnum_lat, cellidx + Grid::cellnum_lat + 1, cellidx + Grid::cellnum_lat + 2,
                     cellidx + Grid::cellnum_lat * 2 + 1, cellidx + Grid::cellnum_lat * 2, cellidx + Grid::cellnum_lat * 2 + 1};
};

double P2Pdistance(const Point &a, const Point &b)
{
    return std::sqrt((a.x - b.x) * (a.y - b.y));
}

std::vector<int> FindNeighbors(const std::vector<Point> &points, int pointIdx, double eps)
{
    std::vector<int> neighbors;
    for (int i = 0; i < points.size(); i++)
    {
        if (P2Pdistance(points[pointIdx], points[i]) <= eps && i != pointIdx)
        {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

void ExpandCluster(std::vector<Point> &points, int pointIdx, int &clusterID, double eps, int minPts)
{
    std::vector<int> seeds = FindNeighbors(points, pointIdx, eps);
    if (seeds.size() < minPts - 1)
    {
        points[pointIdx].clusterID = -1; // Noise
        return;
    }
    clusterID++;
    points[pointIdx].clusterID = clusterID; // Start Point of cluster

    for (int idx : seeds)
    {
        if (points[idx].clusterID == 0)
        {
            points[idx].clusterID = clusterID;
            points[idx].nearestCorePointDistance = P2Pdistance(points[idx],points[pointIdx]);
        }
        else if (points[idx].clusterID == -1 || points[idx].isBorderPointFlg)
        {
            points[idx].isBorderPointFlg = 1;
            double tmpdist = P2Pdistance(points[idx],points[pointIdx]);
            if (tmpdist < points[idx].nearestCorePointDistance){
                    points[idx].clusterID = clusterID;
                    points[idx].nearestCorePointDistance = tmpdist;
            }
        }
    }

    while (!seeds.empty())
    {
        int currentP = seeds.back();
        seeds.pop_back();
        std::vector<int> result = FindNeighbors(points, currentP, eps);
        if (result.size() >= minPts - 1)
        {
            for (int i : result)
            {
                if (points[i].clusterID == 0)
                {
                    seeds.push_back(i);
                    points[i].clusterID = clusterID;
                    points[i].nearestCorePointDistance = P2Pdistance(points[i],points[currentP]);;
                }
                else if (points[i].clusterID == -1 || points[i].isBorderPointFlg)
                {
                    points[i].isBorderPointFlg = 1;
                    double tmpdist = P2Pdistance(points[i],points[currentP]);
                    if (tmpdist < points[i].nearestCorePointDistance){
                        points[i].clusterID = clusterID;
                        points[i].nearestCorePointDistance = tmpdist;
                    }
                }
            }
        }
        else
        {
            points[currentP].isBorderPointFlg = 1;
        }
    }
}

void DBSCAN(std::vector<Point> &points, double eps, int minPts)
{
    int clusterID = 0;
    for (int i = 0; i < points.size(); i++)
    {
        if (points[i].clusterID == 0)
        {
            ExpandCluster(points, i, clusterID, eps, minPts);
        }
    }
}





void DeterminateCorePoints(std::vector<Point> &points, Grid &mygrid, double eps, int minPts)
{
    int cnt;
    bool fnd;
    std::vector<int> neighborcells;

    for (auto cell : mygrid.GridMapTable)
    {
        if (cell.second.size() > minPts)
        {
            std::vector<int> corevector;
            for (int ptidx : cell.second)
            {
                // points[ptidx].clusterID = 0;
                corevector.push_back(ptidx);
            }
            mygrid.GridMapCoreTable.emplace(cell.first, corevector);
        }
        else
        {
            std::vector<int> corevector;
            fnd = 0;
            for (int ptidx : cell.second)
            { // 20 iteration
                points[ptidx].clusterID = -1;
                cnt = cell.second.size() - 1;
                mygrid.DefineSearchingRange(cell.first, neighborcells);
                for (int cidx : neighborcells)
                {
                    for (int tgtptidx : mygrid.GridMapTable[cidx])
                    {
                        if (P2Pdistance(points[tgtptidx], points[ptidx]) < eps)
                        {
                            cnt++;
                            if (cnt >= minPts)
                            {
                                points[ptidx].clusterID = 0;
                                corevector.push_back(ptidx);
                                fnd = 1;
                                break;
                            }
                        }
                    }
                }
            }
            if (fnd)
            {
                mygrid.GridMapCoreTable.emplace(cell.first, corevector);
            }
        }
    }
}

void DeterminateNoiseBorderPoint(std::vector<Point> &points, Grid &mygrid, double eps, int minPts)
{
    int cnt;
    std::vector<int> neighborcells;
    for (auto cell : mygrid.GridMapTable)
        {
        mygrid.DefineSearchingRange(cell.first, neighborcells);
            for (int ptidx : cell.second){
                if (points[ptidx].clusterID == -1){
                for (int cellidx : neighborcells){
                    for (int coreptidx :mygrid.GridMapCoreTable[cellidx]){
                        double tmpdist = P2Pdistance(points[ptidx],points[coreptidx]);
                        if (tmpdist < eps && tmpdist < points[ptidx].nearestCorePointDistance ){
                            points[ptidx].isBorderPointFlg = 1;
                            points[ptidx].clusterID = points[coreptidx].clusterID;
                            points[ptidx].nearestCorePointDistance = tmpdist;
                        }
                    }
                }
            }
        }

    }
}

bool FindReachPoint(int cellidx, int targetcellidx, std::vector<Point> &points, Grid &mygrid, double eps)
{
    if (mygrid.GridMapTable[targetcellidx][0] >= 0)
    {
        return 0; // already have visted
    }
    for (int ptidx : mygrid.GridMapCoreTable[cellidx])
    {
        for (int tgtptidx : mygrid.GridMapCoreTable[targetcellidx])
        {
            if (P2Pdistance(points[ptidx], points[tgtptidx]) < eps)
            {
                return 1;
            }
        }
    }
    return 0;
}

void MergingCluster(std::vector<Point> &points, Grid &mygrid, double eps, int minPts)
{
    std::vector<int> seeds;
    int tmpcellidx;
    std::vector<int> neighborcells;

    int clusterID = 0;
    for (auto cell : mygrid.GridMapCoreTable)
    {
        tmpcellidx  = cell.first;
        if (mygrid.GridMapTable[tmpcellidx][0] <= 0)//未到達
        {
            seeds = {cell.first};
            clusterID++;
            for (int pntidx : mygrid.GridMapCoreTable[tmpcellidx])
            {
                points[pntidx].clusterID = clusterID;
            }

        while (!seeds.empty())
            tmpcellidx = seeds.back();
            seeds.pop_back();
            mygrid.DefineSearchingRange(tmpcellidx, neighborcells);
            for (int cidx : neighborcells)
            {
                if (FindReachPoint(tmpcellidx, cidx, points, mygrid, eps))
                {
                    seeds.push_back(cidx);
                    for (int pntidx : mygrid.GridMapCoreTable[cidx])
                    {
                        points[pntidx].clusterID = clusterID;
                    }
                }
            }
        }
    }
}

void GridDBSCAN(std::vector<Point> &points, double eps, int minPts, double gridsp_lon, double gridsp_lat, int gridnum_lon, int gridnum_lat)
{
    Grid mygrid = Grid(eps / 1.4141356, gridsp_lon, gridsp_lat, gridnum_lon, gridnum_lat);
    mygrid.RegisterPoint2Grid(points);
    DeterminateCorePoints(points, mygrid, eps, minPts);
    MergingCluster(points, mygrid, eps, minPts);
    DeterminateNoiseBorderPoint(points, mygrid, eps, minPts);
    int clusterID = 0;
}

int main()
{
    std::vector<Point> points = {{1, 2}, {2, 3}, {2, 2}, {8, 7}, {8, 8}, {25, 80}};
    double eps = 2.0;
    int minPts = 2;
    DBSCAN(points, eps, minPts);
    return 0;
}