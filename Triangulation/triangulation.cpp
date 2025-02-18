#include <conio.h>
#include <iostream>
#include <vector>
#include <windows.h>
#include <windowsx.h>

int STARTX = 100;
int STARTY = 50;
int STOPX = 100 + 8 * 90;
int STOPY = 50 + 8 * 90;
int STEP = 8;

struct Point {
    int x;
    int y;
};

/*
* Она сравнивает координаты x и y обеих точек и возвращает true, если они равны, и false в противном случае
*/

bool BelongsPoints(const Point& point1, const Point& point2) {
    if (point1.x == point2.x && point1.y == point2.y) 
        return true;

    return false;
}

/*
* Эта функция проверяет, принадлежат ли два ребра, которые заданны точками point1, point2 и point3, point4, одному
* и тому же ребру. Она использует функцию BelongsPoints для проверки принадлежности точек и возвращает
* true, если ребра совпадают, и false в противном случае
*/

bool BelongsEdges(const Point& point1, const Point& point2, const Point& point3, const Point& point4) {
    if (BelongsPoints(point1, point3) && BelongsPoints(point2, point4)) 
        return true;
    if (BelongsPoints(point1, point4) && BelongsPoints(point2, point3)) 
        return true;

    return false;
}

/*
* Она использует функции BelongsPoints и BelongsEdges для проверки принадлежности точек и ребер и возвращает true,
если треугольники совпадают, и false в противном случае
*/

bool BelongsTriangles(const Point& point1_1, const Point& point1_2, const Point& point1_3, const Point& point2_1, const Point& point2_2, const Point& point2_3) {
    if (BelongsPoints(point1_1, point2_1)) 
        return BelongsEdges(point1_2, point1_3, point2_2, point2_3);
    if (BelongsPoints(point1_1, point2_2)) 
        return BelongsEdges(point1_2, point1_3, point2_1, point2_3);
    if (BelongsPoints(point1_1, point2_3)) 
        return BelongsEdges(point1_2, point1_3, point2_1, point2_2);

    return false;
}

/*
* Эта функция проверяет, принадлежит ли ребро, заданное точками point1 и point2, одному из треугольников в списке triangles
* Она использует функции BelongsPoints и BelongsEdge, которые описаны выше
*/

bool EdgeInTriangle(const std::vector<Point>& triangles, const Point& point1, const Point& point2) {

    int k = 0;
    for (size_t i = 0; i < triangles.size(); i += 3)
        for (size_t j = 0; j < 3; ++j, ++k) {
            if (BelongsPoints(triangles[k], point1)) {
                if (j == 0) {
                    if (BelongsPoints(triangles[k + 1], point2) || BelongsPoints(triangles[k + 2], point2)) 
                        return true;
                }
                if (j == 1) {
                    if (BelongsPoints(triangles[k - 1], point2) || BelongsPoints(triangles[k + 1], point2)) 
                        return true;
                }
                if (j == 2) {
                    if (BelongsPoints(triangles[k - 2], point2) || BelongsPoints(triangles[k - 1], point2)) 
                        return true;
                }
            }
        }
    return false;
}

/*
* Эта функция вычисляет значение правила Делоне для трех точек point1, point2 и point3
*/

double DelaunaysRule(const Point& point1, const Point& point2, const Point& point3) {

    double d = double(point3.y - point2.y) * (point2.x - point1.x) - double(point3.x - point2.x) * (point2.y - point1.y);
    if (d == 0) 
        return 1e99;

    return (double(point3.x - point2.x) * (point1.x - point3.x) + double(point3.y - point2.y) * (point1.y - point3.y)) / d;
}

/*
* Эта функция вычисляет косинус угла между векторами, образованными точками point1, point2 и point3
*/

double Cos(Point& point1, Point& point2, Point& point3) {

    int x0 = point1.x - point2.x;
    int y0 = point1.y - point2.y;
    int x1 = point3.x - point2.x;
    int y1 = point3.y - point2.y;

    double L = (x0 * x0 + y0 * y0) * (x1 * x1 + y1 * y1);
    if (L == 0) 
        return 0;

    double cos = (x0 * x1 + y0 * y1) / sqrt(L);

    if (cos > 1) 
        cos = 1;
    if (cos < -1) 
        cos = -1;

    return cos;
}

/*
* Эта функция рисует линию между двумя точками point1 и point2
*/

void Line(HDC hdc, const Point& point1, const Point& point2) {

    double steps = 0;
    int x = 0;
    int y = 0;
    int left_x = 0;
    int top_y = 0;
    Point p1 = point1;
    Point p2 = point2;

    if (p1.x > p2.x) {
        Point temp = p1;
        p1 = p2;
        p2 = temp;
    }
    steps = (double)(p2.y - p1.y) / (p2.x - p1.x);
    x = p1.x;
    y = p1.y;
    double double_y = p1.y;
    while (x <= p2.x) {
        left_x = STARTX + x * STEP;
        top_y = STOPY - (y + 1) * STEP;
        Rectangle(hdc, left_x + 1, top_y + 1, left_x + STEP, top_y + STEP);
        double_y += steps;
        y = round(double_y);
        x++;
    }

    p1 = point1;
    p2 = point2;
    if (p1.y > p2.y) {
        Point temp = p1;
        p1 = p2;
        p2 = temp;
    }
    steps = (double)(p2.x - p1.x) / (p2.y - p1.y);
    x = p1.x;
    y = p1.y;
    double double_x = p1.x;
    while (y <= p2.y) {
        left_x = STARTX + x * STEP;
        top_y = STOPY - (y + 1) * STEP;
        Rectangle(hdc, left_x + 1, top_y + 1, left_x + STEP, top_y + STEP);
        double_x += steps;
        x = round(double_x);
        y++;
    }
}

/*
* Эта функция находит первое ребро в списке точек points. Она находит самую левую точку a0 и затем
* выбирает точку с наименьшим угловым коэффициентом
*/

void FindFirstEdge(const std::vector<Point>& points, int& i0, int& i1) {
    // Сначала найдем a0 - самую левую точку
    int min_x = points[0].x;
    int min_i = 0;
    for (size_t i = 1; i < points.size(); ++i) {
        int x = points[i].x;
        if (x < min_x) {
            min_x = x;
            min_i = i;
        }
    }
    i0 = min_i;
    // Теперь среди остальных вершин выбираем идущую круче всех вверх -
    // то есть с наименьшим угловым коэффициентом (dy/dx - минимум)
    double k_min = 2e99;
    int min_diff = 999999999;
    for (size_t i = 0; i < points.size(); ++i) {

        if (i == i0) 
            continue;

        int dx = points[i].x - points[i0].x;			// Так как (index0) самая левая точка, то dx>=0
        int dy = points[i].y - points[i0].y;            // dy - любое
        double k = 0;
        if (dx != 0) {
            k = dy / (double)dx;
        }
        else 
            k = -2e99;
        
        if (k < k_min) {
            k_min = k;
            min_diff = abs(dx) + abs(dy);
            min_i = i;
            continue;
        }
        // Если же коэффициенты равны, выберем более короткий вектор:
        if (k == k_min) {
            int diff = abs(dx) + abs(dy);
            if (diff < min_diff) {
                min_diff = diff;
                min_i = i;
            }
        }
    }
    i1 = min_i;
}

/*
* Эта функция находит самое короткое ребро в списке ребер edges.
  Она вычисляет длину каждого ребра и выбирает ребро с наименьшей длиной
*/

void FindShortEdge(const std::vector<Point>& edges, int& index) {

    int len = 2147483647; // INT_MAX;
    for (size_t i = 0; i < edges.size(); i += 2) {

        int a = edges[i].x - edges[i + 1].x;
        int b = edges[i].y - edges[i + 1].y;
        int sum_sq = a * a + b * b;
        if (sum_sq < len) {
            index = i;
            len = sum_sq;
        }
    }
}

/*
* Эта функция находит вершины, присоединенные к ребру, заданному точками point1 и point2, в списке точек points
* Она проверяет каждую точку из списка и выбирает те, которые удовлетворяют условию Делоне
*/

void PutDot(const Point& point1, const Point& point2, const std::vector<Point>& points, std::vector<Point>& points_for_edge) {

    double min_Two_T = 1e99;
    points_for_edge.clear();
    for (size_t i = 0; i < points.size(); ++i) {

        Point dot3;
        dot3.x = points[i].x;
        dot3.y = points[i].y;

        int triangle_sq = (point2.x - point1.x) * (dot3.y - point1.y) - (point2.y - point1.y) * (dot3.x - point1.x);
        if (triangle_sq <= 0) {
            continue;
        }

        double Two_T = -DelaunaysRule(point1, point2, dot3);
        if (Two_T == min_Two_T) {
            points_for_edge.push_back(dot3); // Еще одна точка на одной и той же окружности
        }

        if (Two_T < min_Two_T) {
            min_Two_T = Two_T;
            points_for_edge.clear();
            points_for_edge.push_back(dot3);
        } // Нашли новую окружность, меньше прежней
    }
}

/*
* Эта функция добавляет ребро, заданное точками point1 и point2, в список ребер edges
* Кроме того она проверяет, не является ли это ребро ребром треугольника в списке triangles
  Если ребро уже существует или является ребром треугольника, функция ничего не делает
*/

void PutEdge(std::vector<Point>& edges, const std::vector<Point>& triangles, const Point& point1, const Point& point2) {
    for (size_t i = 0; i < edges.size(); i += 2) {
        if (BelongsEdges(point1, point2, edges[i], edges[i + 1]))
        {
            edges.erase(edges.begin() + i, edges.begin() + i + 2);
            return;
        }
    }
    if (EdgeInTriangle(triangles, point1, point2)) 
        return;	                        // Если есть такой отрезок, то ничего не делаем.
    edges.push_back(point1);
    edges.push_back(point2);
}

/*
* Эта функция добавляет треугольник, заданный точками point1, point2 и point3, в список треугольников triangles
* И еще проверяет, не является ли этот треугольник уже существующим треугольником в списке, если треугольник существует, то ничего не далаем
*/

void PutTriangle(std::vector<Point>& triangles, const Point& point1, const Point& point2, const Point& point3) {
    for (size_t i = 0; i < triangles.size(); i += 3) {
        if (BelongsTriangles(triangles[i], triangles[i + 1], triangles[i + 2], point1, point2, point3))
            return;
    }
    // Такого треугольника нет в списке, добавим его.			
    triangles.push_back(point1);
    triangles.push_back(point2);
    triangles.push_back(point3);
}

/*
* Эта функция выполняет триангуляцию Делоне для списка точек points и отображает результат
*/

void DelaynauTriangulation(HDC hdc, const std::vector<Point> points) {
    std::vector<Point> triangles;           // Список полученных треугольников по Треангуляции Делоне
    std::vector<Point> edges;               // Список ребер.
    std::vector<Point> points_for_edge;     // Список вершин, присоединенных к ребру.
    Point point1, point2, point3;

    if (points.size() < 3) {
        return;
    }

    // Вставляем начальное ребро.
    int i0 = 0, i1 = 0;
    FindFirstEdge(points, i0, i1);
    edges.push_back(points[i0]);
    edges.push_back(points[i1]);

    while (edges.size() > 0) { // Пока есть ребра, которые мы не обработали

        // Находим самое короткое ребро в списке ребер edges
        int index = -1;
        FindShortEdge(edges, index);
        point1 = edges[index];
        point2 = edges[index + 1];
        edges.erase(edges.begin() + index, edges.begin() + index + 2);
        PutDot(point1, point2, points, points_for_edge);

        // Если мы не нашли вершину, удовлетворяющую условию Делоне, то просто идем дальше.
        if (points_for_edge.size() == 0) continue;

        // Если мы нашли только одну вершину, удовлетворяющую правилу условию Делоне.
        if (points_for_edge.size() == 1) {
            point3 = points_for_edge[0];
            int triangle_area = (point2.x - point1.x) * (point3.y - point1.y) - (point2.y - point1.y) * (point3.x - point1.x);
            if (triangle_area <= 0) return;

            // Добавим ребра в список, если они не являются ребрами треугольника.
            PutEdge(edges, triangles, point1, point3);
            PutEdge(edges, triangles, point3, point2);
            PutTriangle(triangles, point1, point2, point3);
            continue;
        }

        // Остался случай, когда присоедиенных вершин несколько. Сортировка происходит по возрастанию косинуса угла.
        // Значит, изначально занесется самое большое ребро, лежащее напротив изначального.
        for (size_t i = 1; i < points_for_edge.size(); ++i) {
            double phi1 = Cos(point2, point1, points_for_edge[i - 1]);
            double phi2 = Cos(point2, point1, points_for_edge[i]);
            if (phi1 < phi2) {
                std::swap(points_for_edge[i - 1].x, points_for_edge[i].x);
                std::swap(points_for_edge[i - 1].y, points_for_edge[i].y);
                if (i > 1) i -= 2;
            }

        }

        // Добавляем полученные треугольники по очереди.
        PutEdge(edges, triangles, points_for_edge[0], point2);
        PutTriangle(triangles, point1, point2, points_for_edge[0]);

        for (size_t i = 1; i < points_for_edge.size(); ++i) {
            PutEdge(edges, triangles, points_for_edge[i], points_for_edge[i - 1]);
            PutTriangle(triangles, point1, points_for_edge[i], points_for_edge[i - 1]);
        }

        PutEdge(edges, triangles, point1, points_for_edge[points_for_edge.size() - 1]);
    }

    // Теперь переставим все вершины в треугольнике так, чтобы пройдя по ним мы смогли построить треугольник.
    for (size_t i = 0; i < triangles.size(); i += 3) {
        int triangle_area = (triangles[i + 1].x - triangles[i].x) * (triangles[i + 2].y - triangles[i].y) -
            (triangles[i + 1].y - triangles[i].y) * (triangles[i + 2].x - triangles[i].x);
        if (triangle_area < 0) {
            std::swap(triangles[i].x, triangles[i + 1].x);
            std::swap(triangles[i].y, triangles[i + 1].y);
        }
    }

    // Выведем Триангуляцию на экран.
    for (size_t i = 0; i < triangles.size(); i += 3) {

        Line(hdc, triangles[i], triangles[i + 1]);
        _getch();
        Line(hdc, triangles[i], triangles[i + 2]);
        Line(hdc, triangles[i + 1], triangles[i + 2]);
    }
}

int main() {
    HWND hwnd = GetConsoleWindow();
    HDC hdc = GetDC(hwnd);
    HBRUSH purpleBrush = CreateSolidBrush(RGB(0, 255, 255));
    SelectBrush(hdc, purpleBrush);
    HPEN dcPen = CreatePen(PS_SOLID, 1, RGB(10, 10, 10));
    SelectPen(hdc, dcPen);
    int n;
    std::cin >> n;
    
    // for example
    std::vector<Point> points = { {78, 65}, {65, 30}, {32, 54}, {48, 24}, {43, 48}, {93, 57}, {22, 20}, {120, 9}, {60,38 } };

    for (size_t i = 0; i < points.size(); ++i)
    {
        int x = STARTX + points[i].x * STEP;
        int y = STOPY - (points[i].y + 1) * STEP;
        Rectangle(hdc, x + 1, y + 1, x + STEP, y + STEP);
    }

    DelaynauTriangulation(hdc, points);

    ReleaseDC(hwnd, hdc);
    return 0;
 
}