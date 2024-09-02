#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>
#include <set>


class Point {
    public:
    double x, y;
    Point(){};
    Point(double x,double y):x(x),y(y){};
    auto operator<(const Point& other) const {
        return  y < other.y;
    }
    auto operator>(const Point& other) const {
        return  y > other.y;
    }
};


class BeachLine{
    public:
    std::vector<Point>* points;//参照点一覧
    std::pair<int,int> sites;//beach lineの左側の切れ端を生成している2点のindex
    double* swwpliney;
    
    // relation parameter in tree data structure
    BeachLine* left;
    BeachLine* right;
    BeachLine* parent;
    bool isleft_edge;
    bool isright_edge;
    bool isroot;
    bool isleaf;
    int height;

    //real world relation parameter 
    BeachLine* prev;
    BeachLine* next;
    bool ismostleft;
    bool ismostright;
    //managing circle event
    Event* circle_event;
    
    BeachLine(std::vector<Point>* points,std::pair<int,int> sites):points(points),sites(sites){};
    double xvalue();
    double balance();
    void update_height();
};


struct Event {
    enum { SITE = 0, CIRCLE = 1, SKIP = 2, };

    Point eventpoint;//eventが発生するラインのyとx
    int pointidx;// 参照点(site eventでのみ有効)
    int event_type;
    BeachLine* manager;// 自分を制御しているbeachline (circle eventでのみ有効

    bool operator<(const Event& other) const {
        return eventpoint.y < other.eventpoint.y;
    }
    Event(){};
    Event(Point eventpoint,int pointidx,int type):eventpoint(eventpoint),pointidx(pointidx),event_type(type){};//site event用
    Event(Point eventpoint,BeachLine* manager,int type):eventpoint(eventpoint),manager(manager),event_type(type){};//circle event用
};



double BeachLine::balance(){
    int leftheight = 0;
    int rightheight = 0;
    if(~isleft_edge){
        leftheight = left->height +1;
    }
    if(~isright_edge){
        rightheight = right->height + 1;
    }
    return   leftheight-rightheight ;
}

void BeachLine::update_height(){
    int leftheight = 0;
    int rightheight = 0;
    if(~isleft_edge){
        leftheight = left->height +1;
    }
    if(~isright_edge){
        rightheight = right->height + 1;
    }
    height =  std::max(leftheight,rightheight);
    return;
}


double BeachLine::xvalue(){
    Point p1 = (*points)[sites.first];
    Point p2 = (*points)[sites.second];
    std::vector<Point> bp = calcBreakPoints(p1,p2,*swwpliney);
    
    if (bp.size() == 2) {
                if (p1.y < p2.y) {
                    return bp[1].x;
                } else {
                    return bp[0].x;
                }
            } else {
                return bp[0].x;
    }
}



std::vector<Point> calcBreakPoints(Point &p1, Point &p2, double d) {
    std::vector<Point> ans;
    if (std::fabs(p1.x - p2.x) < POINT_EPSILON) {
        double y = 0.5 * (p1.y + p2.y);
        double D = std::sqrt(d * d - d * (p1.y + p2.y) + p1.y * p2.y);
        ans.push_back(Point(p1.x - D, y));
        ans.push_back(Point(p1.x + D, y));
    } else if (std::fabs(p1.y - p2.y) < POINT_EPSILON) { //殆ど同じ高さ　交点は1点のみ
        double x = 0.5 * (p1.x + p2.x);
        ans.push_back(Point(x, 0.5 * ((x - p1.x) * (x - p1.x) + p1.y * p1.y  - d * d) / (p1.y - d)));
    } else {
        double D = 2. * std::sqrt(std::pow(p1.x - p2.x, 2) * (d - p1.y) * (d - p2.y) * (std::pow(p1.x - p2.x, 2) + std::pow(p1.y - p2.y, 2)));
        double T = -2. * d * std::pow(p1.x - p2.x, 2) + (p1.y + p2.y) * (std::pow(p2.x - p1.x, 2) + std::pow(p2.y - p1.y, 2));
        double Q = 2. * std::pow(p1.y - p2.y, 2);
        
        double y1 = (T - D) / Q, y2 = (T + D) / Q;
        double x1 = 0.5 * (p1.x * p1.x - p2.x * p2.x + (2 * y1 - p2.y - p1.y) * (p2.y - p1.y)) / (p1.x - p2.x);
        double x2 = 0.5 * (p1.x * p1.x - p2.x * p2.x + (2 * y2 - p2.y - p1.y) * (p2.y - p1.y)) / (p1.x - p2.x);
        
        if (x1 > x2) {
            ans.push_back(Point(x2, y2));//上から順番
            ans.push_back(Point(x1, y1));
        }
        else{
            ans.push_back(Point(x1, y1));
            ans.push_back(Point(x2, y2));
        }
    }
    return ans;
}

bool findCircleCenter(const Point &p1, const Point &p2, const Point &p3, Point &center) {   
    // get normalized vectors
    Point u1 = (p1 - p2).normalized(), u2 = (p3 - p2).normalized();
    
    double cross = crossProduct(u1, u2);
    
    // check if vectors are collinear
    if (std::fabs(cross) < CIRCLE_CENTER_EPSILON) {
        return false;
    }
    
    // get cental points
    Point pc1 = 0.5 * (p1 + p2), pc2 = 0.5 * (p2 + p3);
    
    // get free components
    double b1 = dotProduct(u1, pc1), b2 = dotProduct(u2, pc2);
    // calculate the center of a circle
    center.x = (b1 * u2.y - b2 * u1.y) / cross;
    center.y = (u1.x * b2 - u2.x * b1) / cross;
    
    return true;
}


class BeachLineHandler{
    BeachLine* root;
    std::vector<Point>* points;
    BeachLine* findAboveArc(Point site);
    public:
    std::vector<BeachLine*> insertarc(int i);
    BeachLine* BeachLineHandler::removearc(BeachLine* arc);

    double horizonline;

    void remove(BeachLine* arc);
    void rotate_left(BeachLine* arc);
    void rotate_right(BeachLine* arc);
    void balancecheck(BeachLine* newbeachline);
    void insert(BeachLine* newbeachline, BeachLine* targetbeachline, bool left);

    BeachLineHandler(std::vector<Point>* points):points(points){};
};



BeachLine* BeachLineHandler::findAboveArc(Point site){
    BeachLine* arc = BeachLineHandler::root;
        while (arc != nullptr) {
            }
            if (site.x < arc->xvalue()) { //左にある
                if (!arc->isleft_edge) {
                    arc = arc->left;
                } else {
                    return arc;
                }
            } else {
                if (!arc->isright_edge) {
                    arc = arc->right;
                } else {
                    return arc;
                }
            }  
        return arc;
}


//制約: targetbeachlineは絶対に　leftかrightが0
void BeachLineHandler::insert(BeachLine* newbeachline,  BeachLine* targetbeachline, bool left){
    double newx = newbeachline->xvalue();

    if (left){ //左に挿入
        targetbeachline->isleft_edge = 0;
        targetbeachline->left = newbeachline;
        newbeachline->parent = targetbeachline; //親子の盃を交わす
        
        BeachLineHandler::balancecheck(newbeachline);
    }
    else{
        targetbeachline->isright_edge = 0;
        targetbeachline->right = newbeachline;
        newbeachline->parent = targetbeachline; //親子の盃を交わす
        targetbeachline->update_height(); 
        BeachLineHandler::balancecheck(newbeachline);
    }
}


void BeachLineHandler::remove(BeachLine* rmbeachline){

    BeachLine* replacebeachline;
    BeachLine*  tmparc;
    if(rmbeachline->isleft_edge == 0 && rmbeachline->isright_edge == 0){
        if (rmbeachline->isroot){
            BeachLineHandler::root = nullptr;
            return;//　何にもつながっていない（多分起こらない）
        }else{
        if (rmbeachline->parent->left == rmbeachline){
            rmbeachline->parent->left = nullptr;
            rmbeachline->parent->isleft_edge = 1;
        }else{
            rmbeachline->parent->right = nullptr;
            rmbeachline->parent->isright_edge = 1;
        }
           tmparc = rmbeachline->parent;
        }
    }else{
    if(rmbeachline->isleft_edge == 0){//左に存在している
            replacebeachline = rmbeachline->left;
            while(!replacebeachline->isright_edge){
                replacebeachline = replacebeachline->right; 
            }
            if(!replacebeachline->isleft_edge){
                replacebeachline->parent->right = replacebeachline->left;
                replacebeachline->left->parent = replacebeachline->parent;//親交代
            }else{
                replacebeachline->parent->right = nullptr;
                replacebeachline->parent->isright_edge = 1;
            }
    }
    else if(rmbeachline->isright_edge == 0){//右に存在している
        replacebeachline = rmbeachline->right;
            while(!replacebeachline->isleft_edge){
                replacebeachline = replacebeachline->left; 
            }
            if(!replacebeachline->isright_edge){
                replacebeachline->parent->left = replacebeachline->right;
                replacebeachline->right->parent = replacebeachline->parent;//親交代
            }else{
                replacebeachline->parent->left = nullptr;
                replacebeachline->parent->isleft_edge = 1;
            }
    }
    tmparc = replacebeachline->parent;//replace元の親から上のbalanceを図る
    
    replacebeachline->parent = rmbeachline->parent;
    replacebeachline->isroot = rmbeachline->isroot;
    if(!rmbeachline->isroot){
        if(rmbeachline->parent->left == rmbeachline){
            rmbeachline->parent->left = replacebeachline;
        }
        else{
            rmbeachline->parent->right = replacebeachline;
        }
        replacebeachline->parent = rmbeachline->parent;
    }else{
        replacebeachline->isroot = 1;
        replacebeachline->parent = nullptr;
        BeachLineHandler::root = replacebeachline;
    }

    // replace 
    if(!rmbeachline->isleft_edge){//元データに左が存在するなら
        if(rmbeachline->left != replacebeachline){//自分自身じゃないなら
            rmbeachline->left->parent = replacebeachline;
            replacebeachline->left = rmbeachline->left;
            replacebeachline->isleft_edge = 0;
        }
        else{
            replacebeachline->isleft_edge = 1;
            replacebeachline->left = nullptr;
            tmparc = replacebeachline;//replace元の親から上と言ったがそれが自分になってしまった
        }
    }else{
        replacebeachline->isleft_edge = 1;
        replacebeachline->left = nullptr;
    }
    if(!rmbeachline->isright_edge){//元データに右が存在するなら
        if(rmbeachline->right != replacebeachline){//自分自身じゃないなら
            rmbeachline->right->parent = replacebeachline;
            replacebeachline->right = rmbeachline->right;
            replacebeachline->isright_edge = 0;
        }else{
            replacebeachline->isright_edge = 1;
            replacebeachline->right = nullptr;
            tmparc = replacebeachline;
        }
    }else{
        replacebeachline->isright_edge = 1;
        replacebeachline->right = nullptr;
    }
    }
       BeachLineHandler::balancecheck(tmparc);
       return;
    }






void BeachLineHandler::balancecheck(BeachLine* newbeachline){
    BeachLine* arc = newbeachline;

    bool stopflg = 0;
    while(!stopflg){
        arc->update_height();
        if(BeachLineHandler::root = arc){
            stopflg = 1;
        }
        else{
            arc = arc->parent;
        }
        
        if (arc->balance() >1){ //左over load
           BeachLineHandler::rotate_right(arc);
        }
        else if(arc->balance()< -1){
            BeachLineHandler::rotate_left(arc);
        }
        
    }
    return;
}


void BeachLineHandler::rotate_right(BeachLine* arc){//left over loadのときのみ作動
    BeachLine* tmparc;
    BeachLine* tmparc2;
    if(arc->left->balance()==-1){//右over load 二重回転
     tmparc = arc->left;//値を保存
     tmparc2 = arc->left->right;     
       if (BeachLineHandler::root  == arc){
          tmparc2->isroot = 1;
          arc->isroot = 0;
          BeachLineHandler::root = tmparc2;
       }
       else{
        tmparc2->parent = arc->parent;
        arc->parent = tmparc2;
        tmparc->parent = tmparc2;
       }

       arc->left = tmparc2->right;
       arc->isleft_edge = tmparc2->isright_edge;
       
       tmparc->right = tmparc2->left;
       tmparc->isright_edge =tmparc2->isleft_edge;
       
       tmparc2->left = tmparc;
       tmparc2->right = arc;
       tmparc2->isleft_edge = 0;
       tmparc2->isright_edge = 0;


       arc->update_height();
       tmparc->update_height();
       tmparc2->update_height();
    }
    else{//一重回転
       tmparc = arc->left;//値を保存     
       if (BeachLineHandler::root  == arc){
          tmparc->isroot = 1;
          arc->isroot = 0;
          BeachLineHandler::root = tmparc;
       }
       else{
        tmparc->parent = arc->parent;
        arc->parent = tmparc;
       }
       arc->left = tmparc->right;
       arc->isleft_edge = tmparc->isright_edge;
       
       tmparc->right = arc;
       tmparc->isright_edge = 0;
       arc->update_height();
       tmparc->update_height();
    }
    return;
}

void BeachLineHandler::rotate_left(BeachLine* arc){//left over loadのときのみ作動
    BeachLine* tmparc;
    BeachLine* tmparc2;
    if(arc->right->balance()==1){//左over load 二重回転
     tmparc = arc->right;//値を保存
     tmparc2 = arc->right->left;     
       if (BeachLineHandler::root  == arc){
          tmparc2->isroot = 1;
          arc->isroot = 0;
          BeachLineHandler::root = tmparc2;
       }
       else{
        tmparc2->parent = arc->parent;
        arc->parent = tmparc2;
        tmparc->parent = tmparc2;
       }
       
       arc->right = tmparc2->left;
       arc->isright_edge = tmparc2->isleft_edge;
       
       tmparc->left = tmparc2->right;
       tmparc->isleft_edge =tmparc2->isright_edge;
       
       tmparc2->right = tmparc;
       tmparc2->left = arc;
       tmparc2->isleft_edge = 0;
       tmparc2->isright_edge = 0;

       arc->update_height();
       tmparc->update_height();
       tmparc2->update_height();
    }
    else{//一重回転
       tmparc = arc->right;//値を保存     
       if (BeachLineHandler::root  == arc){
          tmparc->isroot = 1;
          arc->isroot = 0;
          BeachLineHandler::root = tmparc;
       }
       else{
        tmparc->parent = arc->parent;
        arc->parent = tmparc;
       }
       arc->right = tmparc->left;
       arc->isright_edge = tmparc->isleft_edge;
       
       tmparc->left = arc;
       tmparc->isleft_edge = 0;
       arc->update_height();
       tmparc->update_height();
    }
    return;
}

std::vector<BeachLine*> BeachLineHandler::insertarc(int i){
    std::vector<BeachLine*> results;
    BeachLine* arc = BeachLineHandler::findAboveArc((*points)[i]);//挿入先のbeach line
    if ((*points)[i].x< arc->xvalue()){
        arc->isleft_edge = 0;
        std::pair<int,int> thiscite =  {i,arc->sites.first};
        BeachLine* newbeachline1 =  new BeachLine(points,thiscite);
        BeachLineHandler::insert(newbeachline1,arc,1);
        thiscite =  {arc->sites.first,i};
        BeachLine* newbeachline2 =  new BeachLine(points,thiscite);
        BeachLineHandler::insert(newbeachline2,newbeachline1,1);
        // pudate next prev relatioon  @@@prev edgeとかいらないのか？
        newbeachline2->prev = arc->prev;
        newbeachline2->next = newbeachline1;
        newbeachline1->prev = newbeachline2;
        newbeachline1->next = arc;
        arc->prev->next = newbeachline2;
        arc->prev = newbeachline1;
        results.push_back(newbeachline2);
        results.push_back(newbeachline1);
    }
    else{
        arc->isright_edge = 0;
        std::pair<int,int> thiscite =  {arc->sites.second,i};
        BeachLine* newbeachline1 =  new BeachLine(points,thiscite);
        BeachLineHandler::insert(newbeachline1,arc,0);
        thiscite =  {i,arc->sites.second};
        BeachLine* newbeachline2 =  new BeachLine(points,thiscite);
        BeachLineHandler::insert(newbeachline2,newbeachline1,0);
        // pudate next prev relatioon
        newbeachline2->next = arc->next;
        newbeachline2->prev = newbeachline1;
        newbeachline1->next = newbeachline2;
        newbeachline1->prev = arc;
        arc->next->prev = newbeachline2;
        arc->next = newbeachline1;
        results.push_back(newbeachline1);
        results.push_back(newbeachline2);
    }
    return results;
}


BeachLine* BeachLineHandler::removearc(BeachLine* arc){

    std::pair<int,int> thiscite =  {arc->sites.first,arc->next->sites.second};
    arc->next->sites = thiscite;// 一個目はreplace 
    arc->prev->next = arc->next;
    arc->next->prev = arc->prev;

    BeachLineHandler::remove(arc);
    return arc->next;
}



void execSiteEvent(Event* event,BeachLineHandler &bh, std::vector<Point> &points){
    std::vector<BeachLine*> beachlines = bh.insertarc(event->pointidx);//beachlineに点を挿入
    Event* circleevent;//eventは　2つのarcの交点にひも付き　自分から見て右側の点とのcircle eventを管理する
    if(checkCircleEvent(beachlines[0]->prev,beachlines[0],circleevent,points)){
       eventqueue.push(*circleevent);
    } //左方向
    
    if(checkCircleEvent(beachlines[1],beachlines[1]->next,circleevent,points)){
        eventqueue.push(*circleevent);        
        }//右方向
}


void execCircleEvent(Event* event,BeachLineHandler &bh, std::vector<Point> &points){
    //arcを削除 
    BeachLine* beachline = bh.removearc(event->manager);
    Event circleevent;
    if(checkCircleEvent(beachline->prev,beachline,circleevent,points)){// 左側のbeachlineのevent
        eventqueue.push(circleevent);
    } //左方向
    if (checkCircleEvent(beachline,beachline->next,circleevent,points)){//eventを更新
        eventqueue.push(circleevent);
    }
}


bool checkCircleEvent(BeachLine* leftarc, BeachLine* rightarc ,Event &circleevent,std::vector<Point> &points){

    Point p1 = points[leftarc->sites.first];
    Point p2 = points[leftarc->sites.second];
    Point p3 = points[rightarc->sites.second];

    Point centerpoint(0,0),eventpoint(0,0);

    if(leftarc->circle_event!=nullptr){
            leftarc->circle_event->event_type= Event::SKIP; //いらない子
    }

    if (p2.y < p1.y && p2.y < p3.y){ //  already passed threepoint あり得るのか？
        return false;
    }
    if(!findCircleCenter(p1, p2, p3, centerpoint)){//center pointが見つからない並びの可能性もある
        return false;
    }
    eventpoint.x = centerpoint.x;
    eventpoint.y = centerpoint.y - P2Pdistance(center,p2);
    circleevent = Event(eventpoint,leftarc,Event::CIRCLE);
    leftarc->circle_event = &circleevent;
    return true;
}

void build_voronoi(std::vector<Point> &points){
    std::sort(points.rbegin(),points.rend());//点群を降順sort
    std::priority_queue<Event> Event_Queue;
    BeachLineHandler bh = BeachLineHandler(&points);
    int cnt = 0;
    for(Point p: points){
        Event(p,cnt,Event::SITE);
        
        if(cnt==0){
            bh.horizonline =p.y;
        }
        cnt++;
    }

    while(!Event_Queue.empty()){
        Event event = Event_Queue.top();

        execCircleEvent();
        
        Event_Queue.pop();
    }

}
    




int main() {
    std::vector<Point> points = {{2.0, 3.0}, {5.0, 4.0}, {1.0, 1.0}, {6.0, 2.0}};
    Voronoi voronoi(points);
    voronoi.generate();
    voronoi.printEdges();
    return 0;
}