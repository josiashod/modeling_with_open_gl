/**
 * @file main.cpp
 * @brief 3D Point Cloud Visualization
 * @version 1.0 (version ECM)
 * 
 * This program visualizes a 3D point cloud using OpenGL and GLUT. It reads point cloud data from a file,
 * calculates the bounding box and barycenter, and displays the points along with a voxel grid representation.
 * 
 * The program includes the following functionalities:
 * - Reading point cloud data from a file
 * - Calculating the minimum, maximum, and barycenter of the point cloud
 * - Initializing a 2D voxel grid and distributing points into the grid
 * - Drawing the bounding box and voxel grid
 * - Handling user input for rotation, perspective, and other display settings
 * 
 * @dependencies
 * - OpenGL
 * - GLUT
 * 
 * @usage
 * - Use the arrow keys to rotate the view
 * - Use the '+' and '-' keys to adjust the rotation speed
 * - Use the '0' and '1' keys to adjust the perspective
 * - Press 'ESC' to exit the program
 * 
 * @note
 * - Ensure that the file "Facade.xyz" is present in the same directory as the executable
 * - The program uses a fixed grid size of 20x20 for the voxel grid
 * 
 * @functions
 * - void DrawBoundingBox(): Draws the bounding box of the point cloud
 * - void DrawVoxels(): Draws the voxel grid
 * - void InitVoxelGrid(): Initializes the voxel grid and distributes points into it
 * - void Display(): Handles the display of the point cloud and voxel grid
 * - void Init(): Initializes OpenGL settings
 * - void Reshape(int width, int height): Handles window resizing
 * - void Idle(): Handles idle state updates
 * - void Keyboard(unsigned char key, int x, int y): Handles keyboard input
 * - void ShowInfos(): Displays information about the current settings
 * - void CalculerMinMaxBarycentre(): Calculates the minimum, maximum, and barycenter of the point cloud
 * 
 * @mainpage
 * The main function initializes GLUT, loads the point cloud data, calculates the necessary parameters,
 * and starts the GLUT main loop.
 */

//3D Point Cloud Visualization
//Filename: main.cpp
//Release: 1.0 (version ECM)
//*********************************************************
//Librairie STL
#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <math.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <string.h>
#include <stdio.h>
#include <iomanip>
#include <time.h>
#include <stdlib.h>
#include <map>
#include <numeric> // Pour std::accumulate

//********************************************************
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
//********************************************************
using namespace std;

template<class T>
class point3D
{
public:
    T x;
    T y;
    T z;
    T r, g, b; // Ajout des composantes RGB
    point3D();
    point3D(T& a, T& b, T& c, T& d, T& e, T& f);
};

template<class T>
point3D<T>::point3D() :
    x(0),
    y(0),
    z(0),
    r(0),
    g(0),
    b(0)
{
}

template<class T>
point3D<T>::point3D(T& a, T& b, T& c, T& d, T& e, T& f) :
    x(a),
    y(b),
    z(c),
    r(d),
    g(e),
    b(f)
{
}

//********************************************************
template<class T>
class point2D
{
public:
    T x;
    T y;
    point2D();
    point2D(T& a, T& b);
};

template<class T>
point2D<T>::point2D() :
    x(0),
    y(0)
{
}

template<class T>
point2D<T>::point2D(T& a, T& b) :
    x(a),
    y(b)
{
}

point3D<double> pointTemp, barycentre;
vector<point3D<double>> nuage;
double Xmin, Ymin, Zmin, Xmax, Ymax, Zmax;

int mouseX, mouseY; // Position de la souris
bool leftButtonPressed = false; // Clic gauche enfoncé
bool middleButtonPressed = false; // Clic milieu enfoncé
float zoomFactor = 1.0; // Facteur de zoom
float translateX = 0.0, translateY = 0.0; // Déplacement de la vue
float rotateX = 0.0, rotateY = 0.0; // Rotation de la vue

int display_mode = 0;

std::map<std::pair<double, double>, std::vector<double>> sommet_z;
std::map<std::pair<double, double>, double> sommet_z_moyen;


struct Color {
    float r, g, b;
    Color(float r_ = 0, float g_ = 0, float b_ = 0) : r(r_), g(g_), b(b_) {}
};

struct Voxel {
    vector<point3D<double>> points;
    double avgZ;
    Color avgColor;
};

vector<vector<Voxel>> voxelGrid2D;
const int GRID_SIZE = 150;

// Add these function declarations
void DrawBoundingBox();
void DrawVoxels();
void InitVoxelGrid();
void DrawHeatMap();
void DrawContourLines();
void ControlSharedPoints();
void Mouse(int button, int state, int x, int y);
void Motion(int x, int y);

// Add these function implementations
Color getRandomColor() {
    return Color(
        (float)rand() / RAND_MAX,
        (float)rand() / RAND_MAX,
        (float)rand() / RAND_MAX
    );
}

/**
 * @brief aroundie les valeurs pour qu'elle soit identique.
 * car les float sont souvent sensiblement différent
 * 
 * @param value 
 * @param precision 
 * @return double 
 */
double roundPrecision(double value, int precision = 5) {
    double scale = pow(10, precision);
    return round(value * scale) / scale;
}

/**
 * @brief divise notre nuage de point matrixe selon le GRID_SIZE demandé
 * 
 */
void InitVoxelGrid() {
    voxelGrid2D.resize(GRID_SIZE, vector<Voxel>(GRID_SIZE));
    double dx = (Xmax - Xmin) / GRID_SIZE;
    double dy = (Ymax - Ymin) / GRID_SIZE;

    // Distribution des points dans la grille
    for (const auto& point : nuage) {
        int i = min(GRID_SIZE - 1, (int)((point.x - Xmin) / dx));
        int j = min(GRID_SIZE - 1, (int)((point.y - Ymin) / dy));
        voxelGrid2D[i][j].points.push_back(point);
    }

    // Calcul des Z moyens et des couleurs moyennes
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (!voxelGrid2D[i][j].points.empty()) {
                double sumZ = 0;
                float sumR = 0, sumG = 0, sumB = 0;
                for (const auto& point : voxelGrid2D[i][j].points) {
                    sumZ += point.z;
                    sumR += point.r;
                    sumG += point.g;
                    sumB += point.b;

                }
                size_t size = voxelGrid2D[i][j].points.size();
                voxelGrid2D[i][j].avgZ = sumZ / size;
                voxelGrid2D[i][j].avgColor = Color(sumR / size, sumG / size, sumB / size);
            }
        }
    }
}

/**
 * @brief Dessine le box contenant de l'objet 3d à representé
 * 
 */
void DrawBoundingBox() {
    double scale = 50;
    glLineWidth(2.0);
    glBegin(GL_LINES);
    glColor3f(1.0, 1.0, 1.0);
    // Bottom face
    glVertex3f((Xmin - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    // Top face
    glVertex3f((Xmin - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    // Vertical edges
    glVertex3f((Xmin - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymin - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmax - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmin - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glVertex3f((Xmin - barycentre.x) / scale, (Zmax - barycentre.z) / scale, (Ymax - barycentre.y) / scale);
    glEnd();
}

/**
 * @brief Dessine chaque pavet de la grid de contenant notre nuage de points
 * 
 */
void DrawVoxels() {
    double scale = 50;
    double dx = (Xmax - Xmin) / GRID_SIZE;
    double dy = (Ymax - Ymin) / GRID_SIZE;
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (!voxelGrid2D[i][j].points.empty()) {
                double x = Xmin + i * dx;
                double y = Ymin + j * dy;

                // Générer une couleur aléatoire pour ce voxel
                Color color = voxelGrid2D[i][j].avgColor;

                glColor3f(color.r, color.g, color.b);

                double x1 = roundPrecision((x - barycentre.x) / scale);
                double z1 = roundPrecision((y - barycentre.y) / scale);
                double x2 = roundPrecision((x + dx - barycentre.x) / scale);
                double z2 = roundPrecision((y - barycentre.y) / scale);
                double x3 = roundPrecision((x + dx - barycentre.x) / scale);
                double z3 = roundPrecision((y + dy - barycentre.y) / scale);
                double x4 = roundPrecision((x - barycentre.x) / scale);
                double z4 = roundPrecision((y + dy - barycentre.y) / scale);

                // Récupération des valeurs Z synchronisées
                double y1 = sommet_z_moyen[{x1, z1}];
                double y2 = sommet_z_moyen[{x2, z2}];
                double y3 = sommet_z_moyen[{x3, z3}];
                double y4 = sommet_z_moyen[{x4, z4}];

                // Dessiner les deux triangles
                glBegin(GL_TRIANGLES);
                // Triangle 1
                glVertex3f(x1, y1, z1);
                glVertex3f(x2, y2, z2);
                glVertex3f(x3, y3, z3);

                // Triangle 2
                glVertex3f(x1, y1, z1);
                glVertex3f(x3, y3, z3);
                glVertex3f(x4, y4, z4);
                glEnd();
            }
        }
    }
    glDisable(GL_BLEND);
}

/**
 * @brief Dessine la carte de chaleur
 * 
 */
void DrawHeatMap() {
    double scale = 50;
    double dx = (Xmax - Xmin) / GRID_SIZE;
    double dy = (Ymax - Ymin) / GRID_SIZE;

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (!voxelGrid2D[i][j].points.empty()) {
                double z = voxelGrid2D[i][j].avgZ;

                // Couleur basée sur l'altitude
                double heightRatio = (z - Zmin) / (Zmax - Zmin);

                Color color;

                if (heightRatio < 0.5) {
                    // Transition du vert foncé (0,0.5,0) au jaune (1,0.8,0)
                    double factor = heightRatio * 2;  // Passe de 0 à 1
                    color.r = factor;  // Augmente R progressivement
                    color.g = 0.5 + factor * 0.3;  // Passe de 0.5 à 0.8
                    color.b = 0.0;
                } else {
                    // Transition du jaune (1,0.8,0) au rouge (1,0,0)
                    double factor = (heightRatio - 0.5) * 2;  // Passe de 0 à 1
                    color.r = 1.0;
                    color.g = 0.8 - factor * 0.8;  // Passe de 0.8 à 0.0
                    color.b = 0.0;
                }

                glColor3f(color.r, color.g, color.b);

                double x = Xmin + i * dx;
                double y = Ymin + j * dy;

                double x1 = roundPrecision((x - barycentre.x) / scale);
                double z1 = roundPrecision((y - barycentre.y) / scale);
                double x2 = roundPrecision((x + dx - barycentre.x) / scale);
                double z2 = roundPrecision((y - barycentre.y) / scale);
                double x3 = roundPrecision((x + dx - barycentre.x) / scale);
                double z3 = roundPrecision((y + dy - barycentre.y) / scale);
                double x4 = roundPrecision((x - barycentre.x) / scale);
                double z4 = roundPrecision((y + dy - barycentre.y) / scale);

                // Récupération des valeurs Z synchronisées
                double y1 = sommet_z_moyen[{x1, z1}];
                double y2 = sommet_z_moyen[{x2, z2}];
                double y3 = sommet_z_moyen[{x3, z3}];
                double y4 = sommet_z_moyen[{x4, z4}];

                glBegin(GL_TRIANGLES);
                // Triangle ABC
                glVertex3f(x1, y1, z1);
                glVertex3f(x2, y2, z2);
                glVertex3f(x3, y3, z3);

                // Triangle ACD
                glVertex3f(x1, y1, z1);
                glVertex3f(x3, y3, z3);
                glVertex3f(x4, y4, z4);
                glEnd();
            }
        }
    }
}

void DrawContourLines() {
    double scale = 50;
    double dx = (Xmax - Xmin) / GRID_SIZE;
    double dy = (Ymax - Ymin) / GRID_SIZE;

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (!voxelGrid2D[i][j].points.empty()) {
                double x = Xmin + i * dx;
                double y = Ymin + j * dy;
                double z = voxelGrid2D[i][j].avgZ;

                // Tracer des lignes de contour
                glColor3f(1, 1, 1);
                glBegin(GL_LINE_LOOP);
                glVertex3f((x - barycentre.x) / scale, (z - barycentre.z) / scale, (y - barycentre.y) / scale);
                glVertex3f((x + dx - barycentre.x) / scale, (z - barycentre.z) / scale, (y - barycentre.y) / scale);
                glVertex3f((x + dx - barycentre.x) / scale, (z - barycentre.z) / scale, (y + dy - barycentre.y) / scale);
                glVertex3f((x - barycentre.x) / scale, (z - barycentre.z) / scale, (y + dy - barycentre.y) / scale);
                glEnd();
            }
        }
    }
}

/**
 * @brief Controle souris
 * 
 * @param button 
 * @param state 
 * @param x 
 * @param y 
 */
void Mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            leftButtonPressed = true; // Clic gauche enfoncé
            mouseX = x;
            mouseY = y;
        } else if (state == GLUT_UP) {
            leftButtonPressed = false; // Clic gauche relâché
        }
    } else if (button == GLUT_MIDDLE_BUTTON) {
        if (state == GLUT_DOWN) {
            middleButtonPressed = true; // Clic milieu enfoncé
            mouseX = x;
            mouseY = y;
        } else if (state == GLUT_UP) {
            middleButtonPressed = false; // Clic milieu relâché
        }
    } else if (button == 3) { // Molette de la souris vers le haut (zoom in)
        zoomFactor *= 1.1; // Augmenter le zoom
    } else if (button == 4) { // Molette de la souris vers le bas (zoom out)
        zoomFactor /= 1.1; // Diminuer le zoom
    }
    glutPostRedisplay(); // Redessiner la scène
}

/**
 * @brief controle mouvement de la souris
 * 
 * @param x 
 * @param y 
 */
void Motion(int x, int y) {
    if (leftButtonPressed) { // Rotation
        rotateX += (y - mouseY) * 0.1; // Rotation autour de l'axe X
        rotateY += (x - mouseX) * 0.1; // Rotation autour de l'axe Y
        mouseX = x;
        mouseY = y;
    } else if (middleButtonPressed) { // Déplacement
        translateX += (x - mouseX) * 0.01; // Déplacement horizontal
        translateY -= (y - mouseY) * 0.01; // Déplacement vertical
        mouseX = x;
        mouseY = y;
    }
    glutPostRedisplay(); // Redessiner la scène
}

//*********************************************
// fonctions d'affichage des données avec GLUT
//*********************************************
void Display();
void Init();
void Reshape(int width, int height);
void Idle();
void Keyboard(unsigned char, int, int);
void ShowInfos();
void CalculerMinMaxBarycentre();
void CollecterZParSommet();
void CalculerZMoyen();

// vue de face
float angle = 92;
float vitesse = 0;
char axe = 'y';
int perspective = 20;

void Display()
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(4, 3, 3, 0, 0, 0, 0, 1, 0);

    // Appliquer le zoom
    glScalef(zoomFactor, zoomFactor, zoomFactor);

    // Appliquer le déplacement
    glTranslatef(translateX, translateY, 0.0);

    // Appliquer la rotation
    glRotatef(rotateX, 1.0, 0.0, 0.0); // Rotation autour de l'axe X
    glRotatef(rotateY, 0.0, 1.0, 0.0); // Rotation autour de l'axe Y


    switch (axe)
    {
        case 'x': glRotated(angle, 1, 0, 0); break;
        case 'y': glRotated(angle, 0, 1, 0); break;
        case 'z': glRotated(angle, 0, 0, 1); break;
    }

    // Affichage du nuage de points
    double scale = 50;
    glEnable(GL_COLOR_MATERIAL);
    glPointSize(1);


    // Dessiner le centre
    glPointSize(5);
    glBegin(GL_POINTS);
    glColor3d(1, 0, 0);
    glVertex3f(0, 0, 0);
    glEnd();

    // Dessiner la boîte englobante, les voxels, la carte de chaleur et les courbes de niveau
    DrawBoundingBox();
    switch(display_mode)
    {
        case 0:
            DrawVoxels();
        break;
        case 1:
            DrawHeatMap();
        break;
        case 2:
            // Afficher le nuage de point
            glBegin(GL_POINTS);
            glColor3d(0, 0, 1);
            for (vector<point3D<double>>::iterator j = nuage.begin(); j != nuage.end(); ++j) {
                pointTemp = *j;
                glVertex3f((pointTemp.x - barycentre.x) / scale, (pointTemp.z - barycentre.z) / scale, (pointTemp.y - barycentre.y) / scale);
            }
            glEnd();
        break;
    }
    // DrawContourLines();

    glutSwapBuffers();
}

void Init()
{
    glEnable(GL_DEPTH_TEST); // activation du test de Z-Buffering
    glutSetCursor(GLUT_CURSOR_NONE); // curseur invisible
}

void Reshape(int width, int height)
{
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(perspective, (float)width / height, 1.0 /*0.001*/, 10 /*100000*/);
}

void Idle()
{
    angle = angle + vitesse;
    if (angle > 360) angle = 0;
    glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case 27: // 'ESC'
        printf("\nFermeture en cours...\n");
        exit(0);
        break;
    case 120: // 'x'
        axe = 'x';
        ShowInfos();
        break;
    case 121: // 'y'
        axe = 'y';
        ShowInfos();
        break;
    case 122: // 'z'
        axe = 'z';
        ShowInfos();
        break;
    case 42: // '*'
        vitesse = -vitesse;
        ShowInfos();
        break;
    case 43: // '+'
        if (vitesse < 359.9) vitesse += 0.1;
        ShowInfos();
        break;
    case 45: // '-'
        if (vitesse > 0.1) vitesse -= 0.1;
        ShowInfos();
        break;
    case 48: // '0'
        perspective--;
        Reshape(640, 480);
        ShowInfos();
        break;
    case 49: // '1'
        perspective++;
        Reshape(640, 480);
        ShowInfos();
        break;
    case 110: // 'n'
        display_mode = (++display_mode) % 3;
        ShowInfos();
        break;
    case 78: // 'N'
        display_mode = (++display_mode) % 3;
        ShowInfos();
        break;
    case 127: // 'DEL'
        axe = 'y';
        vitesse = 0.1;
        perspective = 40; /*45*/
        Reshape(640, 480);
        ShowInfos();
        break;
    }
}

void ShowInfos()
{
    system("cls");
    printf("+ augmente la vitesse\n");
    printf("- diminue la vitesse\n");
    printf("n afficher la carte de maniere differente \n");
    printf("1 augmente la perspective\n");
    printf("0 diminue la perspective\n");
    printf("* inverse le sens de rotation\n");
    printf("x, y et z pour changer l'axe de rotation\n");
    printf("DEL pour restaurer les parametres par defaut\n");
    printf("ESC pour quitter\n\n");
    printf("Vitesse de rotation : %.1f\n", vitesse);
    printf("Axe de rotation : %c\n", axe);
    printf("Perspective : %d\n", perspective);
}

void CalculerMinMaxBarycentre()
{
    if (nuage.empty()) return; // Si le nuage est vide, on quitte la fonction.
    Xmin = Ymin = Zmin = numeric_limits<double>::max();
    Xmax = Ymax = Zmax = numeric_limits<double>::lowest();
    double sommeX = 0, sommeY = 0, sommeZ = 0;
    // Parcours du nuage pour calculer les min, max et barycentre
    for (const auto& point : nuage) {
        if (point.x < Xmin) Xmin = point.x;
        if (point.x > Xmax) Xmax = point.x;
        if (point.y < Ymin) Ymin = point.y;
        if (point.y > Ymax) Ymax = point.y;
        if (point.z < Zmin) Zmin = point.z;
        if (point.z > Zmax) Zmax = point.z;
        sommeX += point.x;
        sommeY += point.y;
        sommeZ += point.z;
    }
    // Calcul du barycentre
    int nbPoints = nuage.size();
    barycentre.x = sommeX / nbPoints;
    barycentre.y = sommeY / nbPoints;
    barycentre.z = sommeZ / nbPoints;
    cout << "Coordonnées du barycentre : (" << barycentre.x << ", " << barycentre.y << ", " << barycentre.z << ")" << endl;
    cout << "Xmin : " << Xmin << ", Ymin : " << Ymin << ", Zmin : " << Zmin << endl;
    cout << "Xmax : " << Xmax << ", Ymax : " << Ymax << ", Zmax : " << Zmax << endl;
}

void CollecterZParSommet() {
    double scale = 50;
    double dx = (Xmax - Xmin) / GRID_SIZE;
    double dy = (Ymax - Ymin) / GRID_SIZE;

    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (!voxelGrid2D[i][j].points.empty()) {
                double x = Xmin + i * dx;
                double y = Ymin + j * dy;
                double z = voxelGrid2D[i][j].avgZ;

                double x1 = roundPrecision((x - barycentre.x) / scale);
                double z1 = roundPrecision((y - barycentre.y) / scale);
                double x2 = roundPrecision((x + dx - barycentre.x) / scale);
                double z2 = roundPrecision((y - barycentre.y) / scale);
                double x3 = roundPrecision((x + dx - barycentre.x) / scale);
                double z3 = roundPrecision((y + dy - barycentre.y) / scale);
                double x4 = roundPrecision((x - barycentre.x) / scale);
                double z4 = roundPrecision((y + dy - barycentre.y) / scale);

                sommet_z[{x1, z1}].push_back((z - barycentre.z) / scale);
                sommet_z[{x2, z2}].push_back((z - barycentre.z) / scale);
                sommet_z[{x3, z3}].push_back((z - barycentre.z) / scale);
                sommet_z[{x4, z4}].push_back((z - barycentre.z) / scale);
            }
        }
    }
}

void CalculerZMoyen() {
    
    for (auto& [coord, z_values] : sommet_z) {
        double sommeZ = std::accumulate(z_values.begin(), z_values.end(), 0.0);
        double z_moyen = sommeZ / z_values.size();
        sommet_z_moyen[coord] = z_moyen;
    }
}


int main(int argc, char **argv)
{
    // Déclare variables
    glutInit(&argc, argv);
    point3D<double> pointTemp;
    int compteur = 0;
    string ligne;
    double a, b, c, d, e, f, g, h;
    nuage.reserve(400000);
    
    ifstream fichier;
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <nom_du_fichier.xyz>" << endl;
        return 1;
    }
    fichier.open(argv[1]);

    if (!fichier) {
        cout << "Fichier inexistant" << endl;
    }
    else {
        while (!fichier.eof()) {
            getline(fichier, ligne);
            sscanf(ligne.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf", &a, &b, &c, &d, &e, &f, &g, &h);
            d /= 255.0;
            e /= 255.0;
            f /= 255.0;
            pointTemp = point3D<double>(a, b, c, d, e, f);
            nuage.push_back(pointTemp);
            ++compteur;
        }
    }
    fichier.close();
    cout << "Fichier " << argv[1] <<" chargé" << endl;
    cout << "Nombre de points dans le nuage : " << compteur << endl;
    ShowInfos();
    // Calcul des min, max et barycentre
    CalculerMinMaxBarycentre();
    // ControlSharedPoints();
    InitVoxelGrid();
    CollecterZParSommet();
    CalculerZMoyen();
    // Initialisation de GLUT
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(640, 480);
    glutInitWindowPosition(40, 40);
    glutCreateWindow("3D Point Cloud Visualization (ECM)");
    Init();
    glutDisplayFunc(Display);
    glutReshapeFunc(Reshape);
    glutIdleFunc(Idle);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);
    glutMotionFunc(Motion);
    glutMainLoop();
    return 0;
}
