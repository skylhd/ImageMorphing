#pragma once
#ifndef MORPH
#define MORPH

#include "CImg.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace cimg_library;

class Point {
public:
    Point(double _x, double _y){
        x = _x;
        y = _y;
    }
    Point& operator = (const Point a) {
        x = a.x;
        y = a.y;
        return *this;
    }
    Point operator - (Point a) {
        Point res(x - a.x, y - a.y);
        return res;
    }
    Point operator + (Point a) {
        Point res(x + a.x, y + a.y);
        return res;
    }
    int operator * (Point a) {
        return x*a.x + y*a.y;
    }

    Point operator * (int a) {
        Point res(x*a, y*a);
        return res;
    }
    Point operator * (float a) {
        Point res(x*a, y*a);
        return res;
    }
    Point operator * (double a) {
        Point res(x*a, y*a);
        return res;
    }

    Point operator / (int a) {
        Point res(x/a, y/a);
        return res;
    }
    Point operator / (float a) {
        Point res(x/a, y/a);
        return res;
    }
    Point operator / (double a) {
        Point res(x/a, y/a);
        return res;
    }


    Point perpendicular() {
        Point res(-y, x);
        return res;
    }
    double abs() {
        return sqrt(x*x + y*y);
    }
    double Sabs() {
        return x*x + y*y;
    }
    double x, y;
};

class Line {
public:
    Line(Point s, Point e) : sp(s),ep(e) {}
    Point sp, ep;
};

template<typename T>
class Morph {
public:
    /*
    the file path of two images
    */
    Morph(const char * s, const char * s2) {
        _img1.load_bmp(s);
        _img2.load_bmp(s2);
        unifySize();
    }
    Morph(CImg<T> &img1, CImg<T> &img2) {
        _img1 = img1;
        _img2 = img2;
        unifySize();
    }
    /*
    draw lines on the Images
    You must press F after your any click to sure the operation
    You can press right mouse button to cancel your operation
    You should keep the direction and the order of the lines in two images are corresponding
    */
    void drawPQ() {
        CImgDisplay dis, dis2;
        CImg<T> copy1(_img1), copy2(_img2);
        dis.display(copy1);
        dis2.display(copy2);
        while (!dis.is_closed()) {
            dis.wait(dis, dis2);
            if ((dis.button()&1)&&dis.mouse_y() >= 0) {
                Point startPos(dis.mouse_x(), dis.mouse_y());
                FtoSure(dis);
                while (true) {
                    dis.wait();
                    if ((dis.button() & 1) && dis.mouse_y() >= 0) {
                        Point endPos(dis.mouse_x(), dis.mouse_y());
                        unsigned char white[3] = { 255,255,255 };
                        copy1.draw_line(startPos.x, startPos.y, endPos.x, endPos.y, white);
                        dis.display(copy1);
                        printf("img1 : (%lf, %lf) -> (%lf, %lf)\n", startPos.x, startPos.y, endPos.x, endPos.y);
                        P1.push_back(startPos);
                        Q1.push_back(endPos);
                        FtoSure(dis);
                        break;
                    }
                    if ((dis.button() & 2)) break;
                }
            }
            if ((dis2.button() & 1) && dis2.mouse_y() >= 0) {
                Point startPos(dis2.mouse_x(), dis2.mouse_y());
                FtoSure(dis2);
                while (true) {
                    dis2.wait();
                    if ((dis2.button() & 1) && dis2.mouse_y() >= 0) {
                        Point endPos(dis2.mouse_x(), dis2.mouse_y());
                        unsigned char white[3] = { 255,255,255 };
                        copy2.draw_line(startPos.x, startPos.y, endPos.x, endPos.y, white);
                        dis2.display(copy2);
                        printf("img2 : (%lf, %lf) -> (%lf, %lf)\n", startPos.x, startPos.y, endPos.x, endPos.y);
                        P2.push_back(startPos);
                        Q2.push_back(endPos);
                        FtoSure(dis2);
                        break;
                    }
                    if ((dis2.button() & 2)) break;
                }
            }
        }
        SavePointsToFile("Data/temp/p1.txt", P1);
        SavePointsToFile("Data/temp/q1.txt", Q1);
        SavePointsToFile("Data/temp/p2.txt", P2);
        SavePointsToFile("Data/temp/q2.txt", Q2);
    }
    /*
    Beier&Neely Algorithm(SIGGRAPH 1992)

    param:
    src : the source image
    p1 : the start Point of the line in the source image
    q1 : the end Point of the line in the source image
    p2 : the start Point of the line in the destination image
    q2 : the end Point of the line in the destination image
    a, p, b : the param to decide the weight
    */
    CImg<unsigned char> BNAlgorithm(const CImg<unsigned char> &src, vector<Point> p1,
        vector<Point> q1, vector<Point> p2, vector<Point> q2, int a = 1, int p = 1, int b = 1) {
        CImg<unsigned char> ans(src._width, src._height, src._depth, src._spectrum);
        cimg_forXY(ans, x, y) {
            Point DSUM(0, 0);
            double weightsum = 0;
            Point X(x, y);
            for (int k = 0; k < p2.size(); k++) {
                double u = (X - p2[k])*(q2[k] - p2[k])/((q2[k]-p2[k]).Sabs());
                double v = (X - p2[k])*((q2[k] - p2[k]).perpendicular()) / ((q2[k] - p2[k]).abs());
                Point sX = p1[k] + (q1[k] - p1[k])*u + ((q1[k] - p1[k]).perpendicular())*v / ((q1[k] - p1[k]).abs());
                Point D = sX - X;
                double length = (q2[k] - p2[k]).abs();
                double dist;
                if (u > 0 && u < 1) dist = abs(v);
                else if (u < 0) dist = (X - p2[k]).abs();
                else dist = (X - q2[k]).abs();
                double weight = pow(pow(length, p) / (a + dist), b);
                DSUM = DSUM + D*weight;
                weightsum += weight;
            }
            Point Xd = X + DSUM / weightsum;
            float *num = Bilinear(src, Xd.x, Xd.y);
            ans(x, y, 0) = num[0];
            ans(x, y, 1) = num[1];
            ans(x, y, 2) = num[2];
            if (num != NULL) delete num;
        }
        return ans;
    }
    /*
    show the Morphing animation
    frames means the number of images bewteen two images
    */
    void ShowMorph(int frames) {
        CImgList<unsigned char> list;
        list.push_back(_img2);
        for (int i = 1; i < frames; i++) {
            vector<Point> p, q;
            for (int j = 0; j < P1.size(); j++) {
                Point newp = P1[j]*i / (float)frames + P2[j]*(frames - i) / (float)frames;
                Point newq = Q1[j]*i / (float)frames + Q2[j]*(frames - i) / (float)frames;
                p.push_back(newp);
                q.push_back(newq);
            }
            //CImgList<unsigned char> tlist;
            CImg<unsigned char> M1 = BNAlgorithm(_img1, P1, Q1, p, q, 1, 1, 3);
            CImg<unsigned char> M2 = BNAlgorithm(_img2, P2, Q2, p, q, 1, 1, 3);
            //tlist.push_back(M1); tlist.push_back(M2);
            //tlist.display();
            CImg<unsigned char> Mix(M1);
            cimg_forXYC(Mix, x, y, c) {
                Mix(x, y, c) = M1(x, y, c)*i / (float)frames + M2(x, y, c) * (frames - i) / (float)frames;
            }
            list.push_back(Mix);
        }
        list.push_back(_img1);
        CImgDisplay dis;
        int index = 0;
        bool flag = true;
        dis.display(list[index]);
        while (!dis.is_closed()) {
            dis.wait(1000/frames);
            if (flag) index++;
            else index--;
            dis.display(list[index]);
            if (index >= list.size() - 1 || index <= 0) {
                dis.wait(1000);
                flag = !flag;
            }
            
        }
        for (int i = 0; i < list.size(); i++) {
            string ffp = "Data/dest/" + to_string(i + 1) + ".bmp";
            list[i].save_bmp(ffp.c_str());
        }
    }

    /*
    Unify the size of two images
    */
    void unifySize() {
        if (_img1.size() > _img2.size()) _img1.resize(_img2._width, _img2._height, _img1._depth, _img1._spectrum);
        else _img2.resize(_img1._width, _img1._height, _img2._depth, _img2._spectrum);
    }
    /*
    Press F to go on
    */
    void FtoSure(CImgDisplay &dis) {
        dis.wait();
        while (dis.is_moved()) {
            dis.wait();
            if (dis.is_keyF()) break;
        }
    }
    /*
    Bilinear interpolation for color image
    return a float array whose length is 3, the value in r, g, b channels
    */
    template<typename NT>
    float *Bilinear(const CImg<NT> &img, float x, float y) {
        float* ans = new float[img._spectrum];
        for (int i = 0; i < img._spectrum; i++) {
            float xu = ceil(x), xd = floor(x),
                yu = ceil(y), yd = floor(y);
            //检验范围
            if (xu < 0) xu = 0;
            if (xd < 0) xd = 0;
            if (yu < 0) yu = 0;
            if (yd < 0) yd = 0;
            if (xu >= img._width) xu = img._width - 1;
            if (xd >= img._width) xd = img._width - 1;
            if (yu >= img._height) yu = img._height - 1;
            if (yd >= img._height) yd = img._height - 1;

            //双线性插值
            float r1, r2;
            if (xu == xd) {
                r1 = img(xu, yd, i);
                r2 = img(xu, yu, i);
            }
            else {
                r1 = (xu - x) / (xu - xd)*img(xd, yd, i) + (x - xd) / (xu - xd)*img(xu, yd, i);
                r2 = (xu - x) / (xu - xd)*img(xd, yu, i) + (x - xd) / (xu - xd)*img(xu, yu, i);
            }
            if (yu == yd) {
                ans[i] = r1;
            }
            else {
                ans[i] = (yu - y) / (yu - yd)*r1 + (y - yd) / (yu - yd)*r2;
            }
        }
        return ans;
    }

    /*
    save P1 to "p1.txt"
    save Q1 to "q1.txt"
    save P2 to "p2.txt"
    save Q2 to "q2.txt"
    the first line of the .txt is a number which means the number of points
    then it has n line with two number which means x and y
    */
    void SavePointsToFile(string fp, vector<Point> data) {
        ofstream f;
        f.open(fp);
        if (f.is_open()) {
            f << data.size() << endl;
            for (int i = 0; i < data.size(); i++) {
                f << data[i].x << " " << data[i].y << endl;
            }
            f.close();
        }
    }

    //Load P1,Q1,P2,Q2 from file
    void LoadPQ(string fp) {
        LoadP1(fp);
        LoadQ1(fp);
        LoadP2(fp);
        LoadQ2(fp);
    }
    
private:
    CImg<T> _img1;
    CImg<T> _img2;
    vector<Point> P1, Q1, P2, Q2;
    void LoadP1(string fp) {
        ifstream f;
        f.open(fp + "p1.txt");
        if (f.is_open()) {
            int n;
            f >> n;
            for (int i = 0; i < n; i++) {
                double x, y;
                f >> x >> y;
                Point newp(x, y);
                P1.push_back(newp);
            }
            f.close();
        }
    }
    void LoadQ1(string fp) {
        ifstream f;
        f.open(fp + "q1.txt");
        if (f.is_open()) {
            int n;
            f >> n;
            for (int i = 0; i < n; i++) {
                double x, y;
                f >> x >> y;
                Point newp(x, y);
                Q1.push_back(newp);
            }
            f.close();
        }
    }
    void LoadP2(string fp) {
        ifstream f;
        f.open(fp + "p2.txt");
        if (f.is_open()) {
            int n;
            f >> n;
            for (int i = 0; i < n; i++) {
                double x, y;
                f >> x >> y;
                Point newp(x, y);
                P2.push_back(newp);
            }
            f.close();
        }
    }
    void LoadQ2(string fp) {
        ifstream f;
        f.open(fp + "q2.txt");
        if (f.is_open()) {
            int n;
            f >> n;
            for (int i = 0; i < n; i++) {
                double x, y;
                f >> x >> y;
                Point newp(x, y);
                Q2.push_back(newp);
            }
            f.close();
        }
    }
};


#endif