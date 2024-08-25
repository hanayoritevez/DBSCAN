#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <set>

struct Point {
    double x, y;
};

struct Edge {
    Point start, end;
    Point direction;
    Point* neighbor;
};

struct Event {
    double y;
    Point point;
    bool isSiteEvent;
    bool operator<(const Event& other) const {
        return y > other.y; // Priority queue is max-heap by default, so we invert the comparison
    }
};

struct BeachLine {
    Point site;
    BeachLine* left;
    BeachLine* right;
    Edge* edge;
    Event* circleEvent;
};

class Voronoi {
public:
    Voronoi(const std::vector<Point>& points) : points(points) {}

    void generate() {
        for (const auto& point : points) {
            events.push({point.y, point, true});
        }

        while (!events.empty()) {
            Event event = events.top();
            events.pop();

            if (event.isSiteEvent) {
                handleSiteEvent(event.point);
            } else {
                handleCircleEvent(event.point);
            }
        }
    }

    void printEdges() const {
        for (const auto& edge : edges) {
            std::cout << "Edge from (" << edge.start.x << ", " << edge.start.y << ") to ("
                      << edge.end.x << ", " << edge.end.y << ")\n";
        }
    }

private:
    std::vector<Point> points;
    std::priority_queue<Event> events;
    std::vector<Edge> edges;
    BeachLine* beachLine = nullptr;

    void handleSiteEvent(const Point& point) {
        if (!beachLine) {
            beachLine = new BeachLine{point, nullptr, nullptr, nullptr, nullptr};
            return;
        }

        BeachLine* arc = findArcAbove(point);
        if (arc->circleEvent) {
            events.erase(std::remove(events.begin(), events.end(), *arc->circleEvent), events.end());
            delete arc->circleEvent;
            arc->circleEvent = nullptr;
        }

        BeachLine* newArc = new BeachLine{point, arc, arc->right, nullptr, nullptr};
        arc->right = newArc;
        if (newArc->right) {
            newArc->right->left = newArc;
        }

        Edge* newEdge = new Edge{arc->site, point, {0, 0}, nullptr};
        edges.push_back(*newEdge);
        arc->edge = newEdge;
        newArc->edge = newEdge;

        checkCircleEvent(arc);
        checkCircleEvent(newArc);
        checkCircleEvent(newArc->right);
    }

    void handleCircleEvent(const Point& point) {
        BeachLine* arc = findArcAbove(point);
        if (!arc || !arc->left || !arc->right) return;

        BeachLine* leftArc = arc->left;
        BeachLine* rightArc = arc->right;

        if (leftArc->circleEvent) {
            events.erase(std::remove(events.begin(), events.end(), *leftArc->circleEvent), events.end());
            delete leftArc->circleEvent;
            leftArc->circleEvent = nullptr;
        }
        if (rightArc->circleEvent) {
            events.erase(std::remove(events.begin(), events.end(), *rightArc->circleEvent), events.end());
            delete rightArc->circleEvent;
            rightArc->circleEvent = nullptr;
        }

        edges.push_back({leftArc->site, rightArc->site, {0, 0}, nullptr});
        leftArc->right = rightArc;
        rightArc->left = leftArc;
        delete arc;

        checkCircleEvent(leftArc);
        checkCircleEvent(rightArc);
    }

    BeachLine* findArcAbove(const Point& point) {
        BeachLine* arc = beachLine;
        while (arc) {
            if (point.x < arc->site.x) {
                if (arc->left) {
                    arc = arc->left;
                } else {
                    break;
                }
            } else {
                if (arc->right) {
                    arc = arc->right;
                } else {
                    break;
                }
            }
        }
        return arc;
    }

    void checkCircleEvent(BeachLine* arc) {
        if (!arc || !arc->left || !arc->right) return;

        Point a = arc->left->site;
        Point b = arc->site;
        Point c = arc->right->site;

        if ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) >= 0) return;

        double bx = b.x - a.x;
        double by = b.y - a.y;
        double cx = c.x - a.x;
        double cy = c.y - a.y;

        double d = 2 * (bx * cy - by * cx);
        double x = (cy * (bx * bx + by * by) - by * (cx * cx + cy * cy)) / d + a.x;
        double y = (bx * (cx * cx + cy * cy) - cx * (bx * bx + by * by)) / d + a.y;

        Point center = {x, y};
        double radius = sqrt((x - a.x) * (x - a.x) + (y - a.y) * (y - a.y));

        Point circleEventPoint = {center.x, center.y - radius};
        Event* circleEvent = new Event{circleEventPoint.y, circleEventPoint, false};
        arc->circleEvent = circleEvent;
        events.push(*circleEvent);
    }
};

int main() {
    std::vector<Point> points = {{2.0, 3.0}, {5.0, 4.0}, {1.0, 1.0}, {6.0, 2.0}};
    Voronoi voronoi(points);
    voronoi.generate();
    voronoi.printEdges();
    return 0;
}