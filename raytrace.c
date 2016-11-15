#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "3dmath.h"


//#define DEBUG

typedef struct //Struct Redefined for lights
{
  int kind; // 0 = plane, 1 = sphere, 2 = camera, 3 = light
  union {
    struct {
      double diffuse_color[3];
      double specular_color[3];
      double position[3];
      double normal[3];
    } plane;
    struct {
      double diffuse_color[3];
      double specular_color[3];
      double position[3];
      int radius;
    } sphere;
    struct {
      double width;
      double height;
    } camera;
    struct {
      double color[3];
      double direction[3];
      double position[3];
      double radial_a2;
      double radial_a1;
      double radial_a0;
      double angular_a0;
      double theta;
    } light;
  };
} Object;

typedef struct Pixel
  {
    unsigned char r, g, b;
  } Pixel;

Object** parseScene(char* input);
int nextChar(FILE* json);
int getC(FILE* json);
int checkNextChar(FILE* json, int val);
char* nextString(FILE* json);
char* checkNextString(FILE* json, char* value);
double* nextVector(FILE* json);
double nextNumber(FILE* json);
Pixel* raycast(Object** objects, int pxW, int pxH);
double planeIntersect(Object* object, double* rO, double* rD);
int imageWriter(Pixel* image, char* input, int pxW, int pxH);


int line = 1;

static inline double sqr(double v)
{
  return v*v;
}

static inline void normalize(double* v)
{
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

static inline double distance(double* v1,double* v2)
{
  double len = sqrt(sqr(v2[0]-v1[0]) + sqr(v2[1]-v1[1]) + sqr(v2[2]-v1[2]));
  if(len <0)
  {
    len = -len;
  }
  return len;
}

static inline double v3_length(double* v)
{
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  return len;
}

static inline double clamp(double in,double l,double u)
{
  if(in < l)
  in = l;
  if(in > u)
  in = u;
  return in;
}
static inline void printObjects(Object** objects)
{
  for(int i = 0;objects[i]!= NULL;i+=1)
  {
    switch(objects[i]->kind) // 0 = plane, 1 = sphere, 2 = camera, 3 = light
    {
       case 0:
       printf("%i: plane\n", i);
       printf("diffuse color: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->plane.diffuse_color[j]);
       }
       printf("\n");
       printf("specular color: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->plane.specular_color[j]);
       }
       printf("\n");
       printf("position: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->plane.position[j]);
       }
       printf("\n");
       printf("normal: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->plane.normal[j]);
       }
       printf("\n");
       break;
       case 1:
       printf("%i: sphere\n", i);
       printf("diffuse color: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->sphere.diffuse_color[j]);
       }
       printf("\n");
       printf("specular color: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->sphere.specular_color[j]);
       }
       printf("\n");
       printf("position: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->sphere.position[j]);
       }
       printf("\n");
       printf("radius: ");
       printf("%i\n", objects[i]->sphere.radius);
       break;
       case 2:
       printf("%i: camera\n", i);
       printf("width: ");
       printf("%lf\n", objects[i]->camera.width);
       printf("height: ");
       printf("%lf\n", objects[i]->camera.height);
       break;
       case 3:
       printf("%i: light\n", i);
       printf("color: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->light.color[j]);
       }
       printf("\n");
       printf("direction: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->light.direction[j]);
       }
       printf("\n");
       printf("position: ");
       for(int j = 0;j<3;j+=1)
       {
         printf("%lf ", objects[i]->light.position[j]);
       }
       printf("\n");
       printf("radial a2: ");
       printf("%lf\n", objects[i]->light.radial_a2);
       printf("radial a1: ");
       printf("%lf\n", objects[i]->light.radial_a1);
       printf("radial a0: ");
       printf("%lf\n", objects[i]->light.radial_a0);
       printf("angular a0: ");
       printf("%lf\n", objects[i]->light.angular_a0);
       break;
       default:
       printf("%i: unknown: %i \n", i, objects[i]->kind);
    }
  }
}
int main (int c, char** argv)
{
  if(c != 5)
  {
    fprintf(stderr, "Error: Invalid number of arguments\n");
    exit(1);
  }
  Object** r = parseScene(argv[3]);
  //printObjects(r);
  int i = 0;
  int pxW = atoi(argv[1]);
  int pxH = atoi(argv[2]);
  Pixel* p = raycast(r,pxW,pxH);
  int q = imageWriter(p, argv[4],pxW, pxH);
  return q;
}

int imageWriter(Pixel* image, char* input, int pxW, int pxH)
{
  FILE* fw = fopen(input, "w"); // File Write as P3

  fprintf(fw, "P3\n");
  fprintf(fw, "%i ", pxH);
  fprintf(fw, "%i\n", pxW);
  fprintf(fw,"%i\n",255);

  int row, col;
  for (row = 0; row < pxH; row += 1) //itterate through image array
    {
    for (col = 0; col < pxW; col += 1)
      {
      fprintf(fw,"%i ", image[pxW*row + col].r);
      fprintf(fw, "%i ", image[pxW*row + col].g);
      fprintf(fw, "%i\n", image[pxW*row + col].b);
      }
    }
  fclose(fw);
}

double planeIntersect(Object* object, double* rO, double* rD)
{
  double* norm = object->plane.normal;
  double* pnt = object->plane.position;
  //a(x-x0) + b(y-y0)+c(z-z0) = 0
  //a((rOx + t*rDx - x0)) + b((rOy + t*rDy - y0)) + c((rOz + t*rDz - z0))
  //a*rOx + t*a*rDx - a*x0 + b*rOy + t*b*rDy - b*y0 + c*rOz + t*c*rDz - c*z0
  //t(a*rDx+ b*rDy+ c*rDz) + (a*rOx + b*roy + c*rOz - a*xO - b*yO - c*zO) = 0
  //t = -(a*rOx + b*roy + c*rOz - a*xO - b*yO - c*zO) / (a*rDx+ b*rDy+ c*rDz)
  //form: mt+b = 0
  double m = norm[0]*rD[0] + norm[1]*rD[1] + norm[2]*rD[2];
  double b = norm[0]*rO[0] + norm[1]*rO[1] + norm[2]*rO[2] - norm[0]*pnt[0] - norm[1]*pnt[1] - norm[2]*pnt[2];
  double t = (-1.0*b)/m;
  if(t >= 0)
  {
    return t;
  }
  else
  {
    return -1;
  }
}

double sphereIntersect(Object* object, double* rO, double* rD)
{
  double r = object->sphere.radius;
  double* pnt = object->sphere.position;
  // (x-h)^2 + (y-j)^2 + (z-k)^2 = r^2
  // (rOx - t*rDx-x0)^2 + (rOy - t*rDy-y0)^2 + (rOz - t*rDz-z0)^2 = r^2

  //rDx^2t^2 - 2rDxrOxt + 2rDxtx0+ rOx^2 - 2rOx x0 + x0^2
  // + rDy^2t^2 - 2rDyrOyt + 2rDyty0+ rOy^2 - 2rOyy0 + y0^2
  // + rDZ^2t^2 - 2rDZrOZt + 2rDZtZ0+ rOZ^2 - 2rOzz0 + z0^2
  // - r^2 = 0

  //quadratic function
  double a = sqr(rD[0] - rO[0]) + sqr(rD[1]- rO[1]) + sqr(rD[2]- rO[2]);

  double b = 2*((rD[0]* (rO[0] - pnt[0])) + (rD[1]* (rO[1]- pnt[1])) + (rD[2]*(rO[2]- pnt[2])));

  double c = sqr(rO[0]- pnt[0]) + sqr(rO[1]-pnt[1]) + sqr(rO[2] -pnt[2])- sqr(r);

  double det = sqr(b) - 4.0 * a * c;

  if (det < 0)
    return -1;

  det = sqrt(det);

  double t0 = (-b - det) / (2.0*a);
  if (t0 >= 0)
    return t0;

  double t1 = (-b + det) / (2.0*a);
  if (t1 >= 0)
    return t1;

  return -1;
}

Pixel* raycast(Object** objects, int pxW, int pxH)
{
  double cx = 0; //camera location
  double cy = 0;
  double ch = 0; //initialize camera frame size
  double cw = 0;

  int i = 0;
  while (objects[i] != NULL) {
    if(objects[i]->kind == 2)
    {
      cw = objects[i]->camera.width;
      ch = objects[i]->camera.height;
      i++;
      break;
    }
  }
  if(cw == 0 || ch == 0)
  {
    fprintf(stderr, "Error: No camera defined ");
    exit(1);
  }
  double pixHeight = ch / (double) pxH; //size of the pixels
  double pixWidth = cw / (double) pxW;

  double rO[3] = {cx, cy, 0}; //origin of ray

  printf("prepare memory for image data\n");
  Pixel* image;
  image = malloc(sizeof(Pixel) * pxW * pxH); //Prepare memory for image data
  for (int y = pxH; y >= 0; y -= 1) {
    for (int x = 0; x < pxW; x += 1) {
      double rD[3] = {cx - (cw/2.0) + pixWidth * (x + 0.5),cy - (ch/2.0) + pixHeight * (y + 0.5),1.0}; //location of current pixel

      normalize(rD);

      double bestT = INFINITY; //initialize the best intersection
      int bestO = -1;
      int i = 0;
      while(objects[i] != NULL) // check all objects for intersection
      {
	       double t = 0;
	        switch(objects[i]->kind) // Added support for Lights
          {
	           case 0:
	            t = planeIntersect(objects[i],rO, rD);
	           break;
             case 1:
              t = sphereIntersect(objects[i],rO, rD);
             break;
             case 2:
             break;
             case 3:
             break;
             default:
             // Horrible error
             fprintf(stderr, "Error: invalid object %i\n", objects[i]->kind);
              exit(1);
	        }
          if (t > 0 && t < bestT) // if the object is closer, replace as new best
          {
            bestT = t;
            bestO = i;
          }
          i++;
        }
        if (bestT > 0 && bestT != INFINITY) // Collect color data // ADD LIGHTS HERE
        {
          //printf("hit\n");

          double* color = malloc(sizeof(double)*3);
          color[0] = 0; // ambient_color[0];
          color[1] = 0; // ambient_color[1];
          color[2] = 0; // ambient_color[2];
          int lightI = 0;
          while(objects[lightI] != NULL)
          {
            //printf("shadow test\n");
             //Shadow test
            if(objects[lightI]->kind != 3) // if Not a light, move on
            {
              lightI+=1;
              continue;
            }

              // IF LIGHT DO WORK!
              double* rOn = malloc(sizeof(double)*3); //ro +(t) rd = intersection
              v3_scale(rD,bestT,rOn);
              v3_add(rOn,rO,rOn);

              double* rDn = malloc(sizeof(double)*3); // dl = |light position - intersection|
              v3_subtract(objects[lightI]->light.position,rOn,rDn);

              double dl = v3_length(rDn);
              normalize(rDn);

              int closest_shadow_object = -1;
              double closest_shadow_object_distance = INFINITY;
              double* closest_position = malloc(sizeof(double)*3);

              int k = 0;
              int shadow = 0;
              while(objects[k] != NULL)
              { //printf("check for shadow %i\n", k);

                if (k == bestO) //skip self
                  {
                    k+=1;
                  continue;
                }

                switch(objects[k]->kind) // Added support for Lights
                {
      	           case 0:
      	            closest_shadow_object_distance = planeIntersect(objects[k],rOn, rDn);
      	           break;
                   case 1:
                    closest_shadow_object_distance = sphereIntersect(objects[k],rOn, rDn);
                   break;
                   case 2:
                   k+=1; //skip cameras
                 continue;
                   break;
                   case 3:
                   k+=1; // skip lights
                 continue;
                   break;
                   default:
                   // Horrible error
                   fprintf(stderr, "Error: invalid object %i\n", objects[k]->kind);
                    exit(1);
      	        }
                if (closest_shadow_object_distance > dl) {
                  k+=1;
                  //printf("beyond light\n");
                  continue;
              	}
                if (closest_shadow_object_distance < dl && closest_shadow_object_distance > 0)
                {
                  closest_shadow_object = k;
                  //printf("found shadow: %i\n", k);
                  shadow = 1;
                }
                k+=1;
              }
              if (closest_shadow_object == -1 && shadow == 0)
              {
                //printf("calc light\n");
              	// N, L, R, V
                double* N = malloc(sizeof(double)*3);
                double* L = malloc(sizeof(double)*3);
                double* R = malloc(sizeof(double)*3);
                double* V = malloc(sizeof(double)*3);
                double* diffuse = malloc(sizeof(double)*3);
                double* specular = malloc(sizeof(double)*3);
                double* position = malloc(sizeof(double)*3);
                switch(objects[bestO]->kind) // Added support for Lights
                {
      	           case 0:
      	            N = objects[bestO]->plane.normal; // plane
                    diffuse = objects[bestO]->plane.diffuse_color;
                  	specular = objects[bestO]->plane.specular_color;
                    position = objects[bestO]->plane.position;
      	           break;
                   case 1:
                    v3_subtract(rOn,objects[bestO]->sphere.position, N); // sphere
                    normalize(N);
                    diffuse = objects[bestO]->sphere.diffuse_color;
                  	specular = objects[bestO]->sphere.specular_color;
                    position = objects[bestO]->sphere.position;
                   break;
                   case 2:
                   break;
                   case 3:
                   break;
                   default:
                   // Horrible error
                   fprintf(stderr, "Error: invalid object %i\n", objects[i]->kind);
                    exit(1);
      	        }

              	L = rDn;

                //v3_scale(rDn,-1,L);

                double dot = v3_dot(N, L);
                v3_scale(N,2.0*(dot),R);
                v3_subtract(R,L,R);//reflection of L;

                v3_scale(rD,-1,V);

                double frad = 1.0/(objects[lightI]->light.radial_a2*sqr(dl)
                + objects[lightI]->light.radial_a1*dl
                + objects[lightI]->light.radial_a0);

                double fang;
                double* vlight = malloc(sizeof(double)*3);
                vlight = objects[lightI]->light.direction;
                if(objects[lightI]->light.theta == 0 || objects[lightI]->light.angular_a0 == 0) //if point light
                {
                  fang = 1.0;
                }
                else
                {
                  double targetVal = cos(objects[lightI]->light.theta * 3.14159265 / 180.0); // if spotlight
                  if(targetVal > v3_dot(rDn,vlight))
                  {
                    fang = 0.0;
                  }
                  else{
                    fang = pow(v3_dot(rDn,vlight),objects[lightI]->light.angular_a0);
                  }
                }


                double* diffusecalc = malloc(sizeof(double)*3);
                double* specularcalc = malloc(sizeof(double)*3);
                if(v3_dot(N,L) > 0)
                {
                  v3_mult(objects[lightI]->light.color,diffuse,diffusecalc); // add diffuse component
                  v3_scale(diffusecalc,v3_dot(N,L),diffusecalc);


                  if(v3_dot(V,R) > 0)
                  {
                  v3_mult(objects[lightI]->light.color,specular,specularcalc); // add specular component
                  v3_scale(specularcalc,pow(v3_dot(V,R),20),specularcalc);
                  }
                  else{
                    for(int q = 0;q < 3; q +=1)
                    {
                      specularcalc[q] = 0.0;
                    }
                  }
                }
              	color[0] += frad * fang * clamp(diffusecalc[0] + specularcalc[0],0,1); // add components
              	color[1] += frad * fang * clamp(diffusecalc[1] + specularcalc[1],0,1);
              	color[2] += frad * fang * clamp(diffusecalc[2] + specularcalc[2],0,1);
              }
              lightI+=1;
          }
          color[0] = clamp(color[0],0,255);
          color[1] = clamp(color[1],0,255);
          color[2] = clamp(color[2],0,255);

          color[0] = color[0]*255.0;
          color[1] = color[1]*255.0;
          color[2] = color[2]*255.0;


          image[pxH*(pxH - y-1) + x].r = (unsigned char)(color[0]); // store color data
          image[pxH*(pxH - y-1) + x].g = (unsigned char)(color[1]);
          image[pxH*(pxH - y-1) + x].b = (unsigned char)(color[2]);
        }
        else
        {

        }
      }
  }

  return image;
}
//Parsing JSON
Object** parseScene(char* input)
{

  int c;
  int objectI = 0;
  Object** objects;
  objects = malloc(sizeof(Object*)*129); //create object array

  FILE* json = fopen(input,"r"); // read file
  if (json == NULL)
  {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", input);
    exit(1);
  }

  checkNextChar(json, '['); //first char should be [

  while(1)
  {
    c = getC(json);
    if (c == ']') // if end of file - done
    {
      fprintf(stderr, "Error: Empty scene file.\n");
      fclose(json);
      exit(1);
    }
    if (c == '{') // if new object
    {
      checkNextString(json,"type");
      checkNextChar(json,':');
      char* value = nextString(json);
      objects[objectI] = malloc(sizeof(Object));
      if (strcmp(value, "camera") == 0) //Get type of object
      {
        objects[objectI]->kind = 2;
      }
      else if (strcmp(value, "sphere") == 0)
      {
        objects[objectI]->kind = 1;
      }
      else if (strcmp(value, "plane") == 0)
      {
        objects[objectI]->kind = 0;
      }
      else if (strcmp(value, "light") == 0)
      {
        objects[objectI]->kind = 3;
      }
      else
      {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        fclose(json);
        exit(1);
      }

      while (1)
      {
        c = nextChar(json);
        if (c == '}') // if end of object
        {
      	  // stop parsing this object
          objectI++; // move to next object index
      	  break;
      	}
        else if (c == ',') // if there is another parameter
        {
      	  // read another field
      	  char* key = nextString(json);
      	  checkNextChar(json, ':');
      	  if ((strcmp(key, "width") == 0)) //save values
              {
      	    double value = nextNumber(json);
            if (objects[objectI]->kind == 2) {
              objects[objectI]->camera.width = value;
            }
            else
            {
              fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
              fclose(json);
              exit(1);
            }
      	      }
          else if ((strcmp(key, "height") == 0))
              {
      	    double value = nextNumber(json);
              if (objects[objectI]->kind == 2) {
                objects[objectI]->camera.height = value;
              }
              else
              {
                fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                fclose(json);
                exit(1);
              }
      	      }
          else if ((strcmp(key, "radius") == 0))
              {
          	    double value = nextNumber(json);
                if (objects[objectI]->kind == 1)
                {
                  objects[objectI]->sphere.radius = value;
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
      	      }
          else if ((strcmp(key, "radial-a2") == 0))
              {
          	    double value = nextNumber(json);
                if (objects[objectI]->kind == 3)
                {
                  objects[objectI]->light.radial_a2 = value;
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
      	      }
          else if ((strcmp(key, "radial-a1") == 0))
              {
          	    double value = nextNumber(json);
                if (objects[objectI]->kind == 3)
                {
                  objects[objectI]->light.radial_a1 = value;
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
      	      }
          else if ((strcmp(key, "radial-a0") == 0))
              {
          	    double value = nextNumber(json);
                if (objects[objectI]->kind == 3)
                {
                  objects[objectI]->light.radial_a0 = value;
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
      	      }
          else if ((strcmp(key, "angular-a0") == 0))
              {
          	    double value = nextNumber(json);
                if (objects[objectI]->kind == 3)
                {
                  objects[objectI]->light.angular_a0 = value;
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
      	      }
          else if ((strcmp(key, "theta") == 0))
              {
          	    double value = nextNumber(json);
                if (objects[objectI]->kind == 3)
                {
                  objects[objectI]->light.theta = value;
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
      	      }
          //Vectors
          else if ((strcmp(key, "color") == 0))
               {
      	          double* value = nextVector(json);
                  if (objects[objectI]->kind == 3)
                  {
                    for(int i = 0;i<3;i++)
                    {
                    objects[objectI]->light.color[i] = value[i];
                    }
                  }
                  else
                  {
                    fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                    fclose(json);
                    exit(1);
                  }
      	       }
           else if ((strcmp(key, "specular_color") == 0))
                {
       	          double* value = nextVector(json);
                   if (objects[objectI]->kind == 1)
                   {
                     for(int i = 0;i<3;i++)
                     {
                     objects[objectI]->sphere.specular_color[i] = value[i];
                     }
                   }
                   else if (objects[objectI]->kind == 0)
                   {
                     for(int i = 0;i<3;i++)
                     {
                     objects[objectI]->plane.specular_color[i] = value[i];
                     }
                   }
                   else
                   {
                     fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                     fclose(json);
                     exit(1);
                   }
       	       }
           else if ((strcmp(key, "diffuse_color") == 0))
                {
       	          double* value = nextVector(json);
                   if (objects[objectI]->kind == 1)
                   {
                     for(int i = 0;i<3;i++)
                     {
                     objects[objectI]->sphere.diffuse_color[i] = value[i];
                     }
                   }
                   else if (objects[objectI]->kind == 0)
                   {
                     for(int i = 0;i<3;i++)
                     {
                     objects[objectI]->plane.diffuse_color[i] = value[i];
                     }
                   }
                   else
                   {
                     fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                     fclose(json);
                     exit(1);
                   }
       	       }
          else if ((strcmp(key, "position") == 0))
              {
                double* value = nextVector(json);
                if (objects[objectI]->kind == 1)
                {
                  for(int i = 0;i<3;i++)
                  {
                  objects[objectI]->sphere.position[i] = value[i];
                  }
                }
                else if (objects[objectI]->kind == 0)
                {
                  for(int i = 0;i<3;i++)
                  {
                  objects[objectI]->plane.position[i] = value[i];
                  }
                }
                else if (objects[objectI]->kind == 3)
                {
                  for(int i = 0;i<3;i++)
                  {
                  objects[objectI]->light.position[i] = value[i];
                  }
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
              }
         else if ((strcmp(key, "normal") == 0))
              {
                double* value = nextVector(json);
                if (objects[objectI]->kind == 0)
                {
                  for(int i = 0;i<3;i++)
                  {
                  objects[objectI]->plane.normal[i] = value[i];
                  }
                }
                else
                {
                  fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                  fclose(json);
                  exit(1);
                }
              }
          else if ((strcmp(key, "direction") == 0))
               {
                 double* value = nextVector(json);
                 if (objects[objectI]->kind == 3)
                 {
                   for(int i = 0;i<3;i++)
                   {
                   objects[objectI]->light.direction[i] = value[i];
                   }
                 }
                 else
                 {
                   fprintf(stderr, "Error: Invalid property, \"%s\", on line number %d.\n", key, line);
                   fclose(json);
                   exit(1);
                 }
              }

          else {
      	    fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n", key, line);
            fclose(json);
              exit(1);
      	    //char* value = next_string(json);
      	       }

      	}
        else
        {
      	  fprintf(stderr, "Error: Unexpected value on line %d\n", line);
          fclose(json);
          exit(1);
      	}
      }
      c = nextChar(json);
      if (c == ',')
      {
        // noop
      }
      else if (c == ']') //if end of file
      {
        objects[objectI] = NULL;
        fclose(json);
        return objects;
      }
      else
      {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        fclose(json);
        exit(1);
      }
    }
  }
}

int getC(FILE* json)
{
  int c = fgetc(json); // get the character

  if (c == '\n') // count the number of lines
  {
    line += 1;
  }

  if (c == EOF)
  {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}
// Grabs the next non whitespace character from the file and returns it.
int nextChar(FILE* json)
{
  int c = getC(json);
  while (isspace(c))
  {
    c = getC(json);
  }
  return c;
}
//check the character that is next, throw error if not correct
int checkNextChar(FILE* json, int val)
{
  int c = nextChar(json);
  if(c==val)
  {
    return c;
  }
  else
  {
    fprintf(stderr, "Error: Expected '%c' on line %d.\n", val, line);
    exit(1);
  }
}
//collect the characters until the next "
char* nextString(FILE* json)
{
  char buffer[129];
  int c = checkNextChar(json,'"');
  c = nextChar(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = nextChar(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}
//check if the string is the one defined.
char* checkNextString(FILE* json, char* value)
{
  char* key = nextString(json);
  if (strcmp(key, value) != 0)
  {
    fprintf(stderr, "Error: Expected %s key on line number %d.\n", value, line);
    exit(1);
  }
  else
  {
    return key;
  }
}
//parse the values in a vector
double* nextVector(FILE* json)
{
  double* v = malloc(3*sizeof(double));
  checkNextChar(json, '[');
  v[0] = nextNumber(json);
  checkNextChar(json, ',');
  v[1] = nextNumber(json);
  checkNextChar(json, ',');
  v[2] = nextNumber(json);
  checkNextChar(json, ']');
  return v;
}
//collect next numeric value
double nextNumber(FILE* json)
{
  float value;
  fscanf(json, "%f", &value);
  return value;
}
