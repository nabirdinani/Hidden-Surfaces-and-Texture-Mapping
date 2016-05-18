//Nabir Dinani
#include <stdlib.h>
#include <gl\glut.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include <list>
#include <vector>/

using namespace std;

/******************************************************************
	Notes:	Image size is 400 by 400.
	This is a LEFT handed coordinate system.  That is, x is in the
 horizontal direction starting from the left, y is the
 vertical direction starting from the bottom, and z is pointing
 INTO the screen (so bigger z values mean farther away).  Think
 of it as inverting the z-axis if that helps.
	Your view vector is ALWAYS [0 0 1] (i.e. you are infinitely far away,
 looking in the POSITIVE Z direction).  So, from any point on
 a triangle, the direction vector to your eye is [0 0 -1]
	This file already contains code to load in data (triangles,
 lights, textures).  You just need to access the data stored
 in the trianglelist, lightlist, and texturelist.  Some
 other helpful routines are also included.
	Call setFramebuffer to set a pixel.  This should be the only
 routine you use to set the color.  Use the getTextureRGB to
 get a texture value.  drawit() will cause the current
 framebuffer to be displayed.
	You can create separate routines, global variables, etc. as
 necessary.  You'll probably want to define a global Z buffer.
 *****************************************************************/

#define ImageW 400
#define ImageH 400

float framebuffer[ImageH][ImageW][3];
float ZMAX = 10000.0;	// NOTE: Assume no point has a Z value greater than 10000.0

char* sourcefile="//Users//NabirDinani//Desktop//csce 441//assignment5//assignment5//triangle.dat";

struct color {
    float r, g, b;		// Color (R,G,B values)
};

float zdepthbuffer[ImageH][ImageW];

struct vertex {
    float x,y,z;		// x, y, z coordinates
    float nx,ny,nz;		// Normal at the vertex
    float u,v;			// Texture coordinates
};

struct triangle {
    int whichtexture;	// The index number of the corresponding texture to apply
    // Note: Use the color returned by the texture for the
    // ambient, diffuse, and specular color, scaled by the
    // coefficients of ambient, diffuse, and specular reflection
    vertex v[3];		// The three vertices
    float kamb;			// The coefficient of ambient reflection
    float kdiff;		// The coefficient of diffuse reflection
    float kspec;		// The coefficient of specular reflection
    int shininess;		// The exponent to use for Specular Phong Illumination
};

struct light {
    // Note: assume all lights are white
    float x,y,z;		// x, y, z coordinates of light
    color brightness;	// Level of brightness of light (0.0 - 1.0)
};

struct texture {
    // Note access using getTextureRGB provided below
    int xsize, ysize;	// The size of the texture in x and y
    float* elements;	// RGB values
};


int numtriangles;		// The number of triangles in the scene
int numlights;			// The number of lights (not including ambient) in the scene
int numtextures;		// The number of textures used in the scene

color ambientlight;		// The coefficient of ambient light

triangle* trianglelist;	// Array of triangles
light* lightlist;		// Array of lights
texture* texturelist;	// Array of textures

/* Pass in a pointer to the texture, t, and the texture coordinates, u and v
 Returns (in R,G,B) the color of the texture at those coordinates */
void getTextureRGB(texture* t, float u, float v, float& R, float& G, float& B) {
    int xval,yval;
    if (u<1.0)
        if (u>=0.0) xval = (int)(u*t->xsize);
        else xval = 0;
        else xval = t->xsize-1;
    if (v<1.0)
        if (v>=0.0) yval = (int)(v*t->ysize);
        else yval = 0;
        else yval = t->ysize-1;
    
    R = t->elements[3*(xval*t->ysize+yval)];
    G = t->elements[(3*(xval*t->ysize+yval))+1];
    B = t->elements[(3*(xval*t->ysize+yval))+2];
}


void swap_values(int i, int j) {
    int temp = i;
    i = j;
    j = temp;
}


// Draws the scene
void drawit(void)
{
    glDrawPixels(ImageW,ImageH,GL_RGB,GL_FLOAT,framebuffer);
    glFlush();
}

// Sets pixel x,y to the color RGB
void setFramebuffer(int x, int y, float R, float G, float B)
{
    if (R<=1.0)
        if (R>=0.0)
            framebuffer[x][y][0]=R;
        else
            framebuffer[x][y][0]=0.0;
        else
            framebuffer[x][y][0]=1.0;
    if (G<=1.0)
        if (G>=0.0)
            framebuffer[x][y][1]=G;
        else
            framebuffer[x][y][1]=0.0;
        else
            framebuffer[x][y][1]=1.0;
    if (B<=1.0)
        if (B>=0.0)
            framebuffer[x][y][2]=B;
        else
            framebuffer[x][y][2]=0.0;
        else
            framebuffer[x][y][2]=1.0;
}

// Normalizes the vector passed in
void normalize(float& x, float& y, float& z) {
    float temp = sqrt(x*x+y*y+z*z);
    if (temp > 0.0) {
        x /= temp;
        y /= temp;
        z /= temp;
    } else {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
}

// Returns dot product of two vectors
float dot(float x1, float y1, float z1, float x2, float y2, float z2) {
    return (x1*x2+y1*y2+z1*z2);
}

// Returns angle between two vectors (in radians)
float angle(float x1, float y1, float z1, float x2, float y2, float z2) {
    normalize(x1,y1,z1);
    normalize(x2,y2,z2);
    return  acos(dot(x1,y1,z1,x2,y2,z2));
}

void drawtriangle(triangle tri) {
    
    int num[3];
    color C;
    
    if (tri.v[0].y <= tri.v[1].y && tri.v[0].y <= tri.v[2].y){
        num[0] = 0; num[1] = 1; num[2] = 2;
        
        if(tri.v[0].y > tri.v[2].y){
            num[0] = 2; num[1] = 0; num[2] = 1;
        }
    }
    else if (tri.v[1].y <= tri.v[2].y){
        num[0] = 1; num[1] = 0; num[2] = 2;
    }
    else{
        num[0] = 2; num[1] = 1; num[2] = 0;
    }
    
    if (tri.v[num[0]].y == tri.v[num[1]].y){
        
        swap_values(num[1], num[2]);
    }
    
    if (tri.v[num[0]].y == tri.v[num[2]].y && tri.v[num[0]].x > tri.v[num[2]].x){
        
        swap_values(num[2], num[0]);
    }
    if (tri.v[num[2]].y - tri.v[num[0]].y == 0){
        
    }
    else if ((tri.v[num[1]].y - tri.v[num[0]].y) == 0){
        
        int temp = num[2];
        num[2] = num[1];
        num[1] = temp;
    }
    else if ((tri.v[num[1]].x - tri.v[num[0]].x) / (tri.v[num[1]].y - tri.v[num[0]].y) > (tri.v[num[2]].x - tri.v[num[0]].x) / (tri.v[num[2]].y - tri.v[num[0]].y)){
        
        int temp = num[2];
        num[2] = num[1];
        num[1] = temp;
    }
    
    GLfloat i, j, k, curr_x, curr_y, curr_z, max_value, min_value, jz, kz;
    
    GLfloat jnx, knx, jny, kny, jnz, knz, ju, ku, jv, kv;
    
    max_value = max(tri.v[num[1]].y, tri.v[num[2]].y);
    min_value = min(tri.v[num[1]].y, tri.v[num[2]].y);
    
    for (i = tri.v[num[0]].y;  i < max_value; i++){
        
        GLfloat ab = (i - tri.v[num[0]].y) / (tri.v[num[1]].y - tri.v[num[0]].y);
        GLfloat ac = (i - tri.v[num[0]].y) / (tri.v[num[2]].y - tri.v[num[0]].y);
        GLfloat bc = (i - tri.v[num[1]].y) / (tri.v[num[2]].y - tri.v[num[1]].y);
        
        GLfloat valx1 = (tri.v[num[1]].x - tri.v[num[0]].x) * ab;
        GLfloat valx2 = (tri.v[num[2]].x - tri.v[num[0]].x) * ac;
        GLfloat valx3 = (tri.v[num[2]].x - tri.v[num[1]].x) * bc;
        
        GLfloat valz1 = (tri.v[num[1]].z - tri.v[num[0]].z) * ab;
        GLfloat valz2 = (tri.v[num[2]].z - tri.v[num[0]].z) * ac;
        GLfloat valz3 = (tri.v[num[2]].z - tri.v[num[1]].z) * bc;
        
        GLfloat valjnx1 = (tri.v[num[1]].nx - tri.v[num[0]].nx) * ab;
        GLfloat valjny2 = (tri.v[num[1]].ny - tri.v[num[0]].ny) * ab;
        GLfloat valjnz3 = (tri.v[num[1]].nz - tri.v[num[0]].nz) * ab;
        GLfloat valjnx4 = (tri.v[num[2]].nx - tri.v[num[1]].nx) * bc;
        GLfloat valjny5 = (tri.v[num[2]].ny - tri.v[num[1]].ny) * bc;
        GLfloat valjnz6 = (tri.v[num[2]].nz - tri.v[num[1]].nz) * bc;
        
        GLfloat valknx1 = (tri.v[num[2]].nx - tri.v[num[0]].nx) * ac;
        GLfloat valkny2 = (tri.v[num[2]].ny - tri.v[num[0]].ny) * ac;
        GLfloat valknz3 = (tri.v[num[2]].nz - tri.v[num[0]].nz) * ac;
        GLfloat valknx4 = (tri.v[num[2]].nx - tri.v[num[1]].nx) * bc;
        GLfloat valkny5 = (tri.v[num[2]].ny - tri.v[num[1]].ny) * bc;
        GLfloat valknz6 = (tri.v[num[2]].nz - tri.v[num[1]].nz) * bc;
        
        GLfloat valju1 = (tri.v[num[1]].u - tri.v[num[0]].u) * ab;
        GLfloat valjv2 = (tri.v[num[1]].v - tri.v[num[0]].v) * ab;
        GLfloat valju3 = (tri.v[num[2]].u - tri.v[num[1]].u) * bc;
        GLfloat valjv4 = (tri.v[num[2]].v - tri.v[num[1]].v) * bc;
        
        GLfloat valku1 = (tri.v[num[2]].u - tri.v[num[0]].u) * ac;
        GLfloat valkv2 = (tri.v[num[2]].v - tri.v[num[0]].v) * ac;
        GLfloat valku3 = (tri.v[num[2]].u - tri.v[num[1]].u) * bc;
        GLfloat valkv4 = (tri.v[num[2]].v - tri.v[num[1]].v) * bc;
        
        if(i <  min_value){
            j = tri.v[num[0]].x + (valx1);
            k = tri.v[num[0]].x + (valx2);
            
            jz = tri.v[num[0]].z + valz1;
            kz = tri.v[num[0]].z + valz2;
            
            jnx = tri.v[num[0]].nx + valjnx1;
            knx = tri.v[num[0]].nx + valknx1;
            jny = tri.v[num[0]].ny + valjny2;
            kny = tri.v[num[0]].ny + valkny2;
            jnz = tri.v[num[0]].nz + valjnz3;
            knz = tri.v[num[0]].nz + valknz3;
            
            ju = tri.v[num[0]].u + valju1;
            ku = tri.v[num[0]].u + valku1;
            jv = tri.v[num[0]].v + valjv2;
            kv = tri.v[num[0]].v + valkv2;
            
        }
        else if(tri.v[num[1]].y < tri.v[num[2]].y){
            j = tri.v[num[1]].x + (valx3);
            k = tri.v[num[0]].x + (valx2);
            
            jz = tri.v[num[1]].z + valz3;
            kz = tri.v[num[0]].z + valz2;
            
            jnx = tri.v[num[1]].nx + valjnx4;
            knx = tri.v[num[0]].nx + valknx1;
            jny = tri.v[num[1]].ny + valjny5;
            kny = tri.v[num[0]].ny + valkny2;
            jnz = tri.v[num[1]].nz + valjnz6;
            knz = tri.v[num[0]].nz + valknz3;
            
            ju = tri.v[num[1]].u + valju3;
            ku = tri.v[num[0]].u + valku1;
            jv = tri.v[num[1]].v + valjv4;
            kv = tri.v[num[0]].v + valkv2;

        }
        else{
            j = tri.v[num[0]].x + (valx1);
            k = tri.v[num[1]].x + (valx3);
            
            jz = tri.v[num[0]].z + valz1;
            kz = tri.v[num[1]].z + valz3;
            
            jnx = tri.v[num[0]].nx + valjnx1;
            knx = tri.v[num[1]].nx + valknx4;
            jny = tri.v[num[0]].ny + valjny2;
            kny = tri.v[num[1]].ny + valkny5;
            jnz = tri.v[num[0]].nz + valjnz3;
            knz = tri.v[num[1]].nz + valknz6;
            
            ju = tri.v[num[0]].u + valju1;
            ku = tri.v[num[1]].u + valku3;
            jv = tri.v[num[0]].v + valjv2;
            kv = tri.v[num[1]].v + valkv4;
            
        }
        
        curr_y = i;
        curr_x = j;
    
        
        GLfloat curr_nx, curr_ny, curr_nz, curr_u, curr_v;
        
        while (curr_x < k){
            
            GLfloat val1 = (k - j);
            GLfloat val2 = (curr_x - j);
            
            curr_z = jz + ((kz - jz) / val1) * val2;
            curr_nx = jnx + ((knx - jnx) / val1) * val2;
            curr_ny = jny + ((kny - jny) / val1) * val2;
            curr_nz = jnz + ((knz - jnz) / val1) * val2;
            curr_u = ju + ((ku - ju) / val1) * val2;
            curr_v = jv + ((kv - jv) / val1) * val2;
            
            int index1 = (int)curr_x;
            int index2 = (int)curr_y;
            
            if (curr_z < zdepthbuffer[index1][index2]){
            
                zdepthbuffer[index1][index2] = curr_z;
                
                getTextureRGB(&texturelist[tri.whichtexture], curr_u, curr_v, C.r, C.g, C.b);
                
                color newcolor;
                int n = 0;
                
                // Ambient lighting
                newcolor.r = tri.kamb * ambientlight.r * C.r;
                newcolor.g = tri.kamb * ambientlight.g * C.g;
                newcolor.b = tri.kamb * ambientlight.b * C.b;
                
                GLfloat left_X, left_Y, left_Z, right_X, right_Y, right_Z;
                
                while(n < numlights){
                    
                    left_X = lightlist[n].x - curr_x;
                    left_Y = lightlist[n].y - curr_y;
                    left_Z = lightlist[n].z - curr_z;
                    
                    normalize(left_X, left_Y, left_Z);
                    normalize(curr_nx, curr_ny, curr_nz);
                    
                    GLfloat temp = dot(left_X, left_Y, left_Z, curr_nx, curr_ny, curr_nz);
                    
                    GLfloat dot_product1 = tri.kdiff * temp;
                    
                    right_X = 2 * dot_product1 * curr_nx - left_X;
                    right_Y = 2 * dot_product1 * curr_ny - left_Y;
                    right_Z = 2 * dot_product1 * curr_nz - left_Z;
                    
                    GLfloat dot_product2 = dot(right_X, right_Y, right_Z, 0, 0, -1);
                    
                    if (!(dot_product1 > 0)){
                        dot_product1 = 0;
                    }
                    if (!(dot_product2 > 0)){
                        dot_product2 = 0;
                    }
                    
                    GLfloat lightingFactor = (tri.kdiff * dot_product1) + (tri.kspec * pow(dot_product2, tri.shininess));
                    
                    newcolor.r = newcolor.r + C.r * lightlist[n].brightness.r * lightingFactor;
                    newcolor.g = newcolor.g + C.g * lightlist[n].brightness.g * lightingFactor;
                    newcolor.b = newcolor.b + C.b * lightlist[n].brightness.b * lightingFactor;
                    
                    n++;
                    
                }
                setFramebuffer(curr_y, curr_x, newcolor.r, newcolor.g, newcolor.b);
            }
            
            curr_x++;
        }
    }
}
void display(void)
{
    /* Your routine here */
    for(int i = 0; i < numtriangles; i++){
        drawtriangle(trianglelist[i]);
    }
    
    drawit();

}

void init(void)
{
    int i,j,k;
    
    // Initialize framebuffer to clear
    for(i=0;i<ImageH;i++) {
        for (j=0;j<ImageW;j++) {
            framebuffer[i][j][0] = 0.0;
            framebuffer[i][j][1] = 0.0;
            framebuffer[i][j][2] = 0.0;
            zdepthbuffer[i][j] = ZMAX;
        }
    }
    
    // Load in data
    ifstream infile(sourcefile);
    if (!infile) {
        cout << "Error! Input file " << sourcefile << " does not exist!" << endl;
        exit(-1);
    }
    infile >> numtriangles >> numlights >> numtextures;
    
    // First read triangles
    trianglelist = new triangle[numtriangles];
    for(i=0;i<numtriangles;i++) {
        infile >> trianglelist[i].whichtexture;
        infile >> trianglelist[i].kamb >> trianglelist[i].kdiff >> trianglelist[i].kspec;
        infile >> trianglelist[i].shininess;
        for(j=0;j<3;j++) {
            infile >> trianglelist[i].v[j].x >> trianglelist[i].v[j].y >> trianglelist[i].v[j].z;
            infile >> trianglelist[i].v[j].nx >> trianglelist[i].v[j].ny >> trianglelist[i].v[j].nz;
            infile >> trianglelist[i].v[j].u >> trianglelist[i].v[j].v;
        }
    }
    
    // Now read lights
    lightlist = new light[numlights];
    infile >> ambientlight.r >> ambientlight.g >> ambientlight.b;
    for(i=0;i<numlights;i++) {
        infile >> lightlist[i].x >> lightlist[i].y >> lightlist[i].z;
        infile >> lightlist[i].brightness.r >> lightlist[i].brightness.g >> lightlist[i].brightness.b;
    }
    
    // Now read textures
    texturelist = new texture[numtextures];
    for(i=0;i<numtextures;i++) {
        infile >> texturelist[i].xsize >> texturelist[i].ysize;
        texturelist[i].elements = new float[texturelist[i].xsize*texturelist[i].ysize*3];
        for(j=0;j<texturelist[i].xsize;j++) {
            for (k=0;k<texturelist[i].ysize;k++) {
                infile >> texturelist[i].elements[3*(j*texturelist[i].ysize+k)];
                infile >> texturelist[i].elements[3*(j*texturelist[i].ysize+k)+1];
                infile >> texturelist[i].elements[3*(j*texturelist[i].ysize+k)+2];
            }
        }
    }
    
    infile.close();
}

int main(int argc, char** argv)
{
	std::cout << "argc = " << argc << std::endl;
	for(int i = 0 ; i < argc ; i ++)
	{
		std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
	}
	std::cout << std::endl;
	int a = 0;

    glutInit(&argc,argv);

	if(argc == 2)
		sourcefile=argv[1];

    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
    glutInitWindowSize(ImageW,ImageH);
    glutInitWindowPosition(100,100);
    glutCreateWindow("Nabir Dinani - Assignment 5");
    init();	
    glutDisplayFunc(display);
    glutMainLoop();
    return 0;
}

