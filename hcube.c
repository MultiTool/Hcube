/* **********************************************************************
 Hypercube map maker by John Atkeson, February 2002.

 compile with
 gcc -o hcube2 hcube.c -lm
 unless using the allegro or gd graphics libraries.

 When using the gd image file IO library, compile with
 gcc -o hcube hcube.c -lgd -lpng -lm

********************************************************************** */
#define allegro_here 0
#define gd_here 0
#define svg_here 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#if allegro_here
  #define grx 1 /* grx=1 means turn graphics output on. */
#else
  #define grx 0
#endif
#if allegro_here
  #include "allegro.h"
#endif

#if gd_here
  #include "gd.h"
  #include "gdfontg.h"
  #include "gdfonts.h"
#endif

#define picwdt 600
#define pichgt picwdt

#define uchar unsigned char
#define uint unsigned int

#define sysint int
#define sysuint unsigned sysint

#ifndef bool
  #define bool unsigned char
  #define true (1==1)
  #define false (1==0)
#endif
#define infinity (1.0e200)
#define fudge (1.0/infinity) /* fudge is to avoid divide-by-zero */
#define rattle (0.5)

#ifndef PI
  #define PI 3.141592653589
#endif
#define twopi (PI*2.0)
#define full_circle (twopi)
#define half_circle (PI)
#define quarter_circle (half_circle/2.0)

#define maxn(a,b) ((a>b)?a:b)
#define minn(a,b) ((a<b)?a:b)

/* C++ compatible mem allocation */
#define allocsafe(num_elements,dattype)\
  ((dattype*)calloc(num_elements,sizeof(dattype)) )
#define freesafe(dat) {if (dat!=NULL){free(dat);dat=NULL;} }
#define reallocsafe(ptr,num_elements,dattype)\
  ((dattype*)realloc((ptr),(num_elements)*sizeof(dattype)))

#define array_chunksize 10
#define expand_array_1(array,oldsize,vartype)\
  if (((oldsize) % array_chunksize)==0){array=reallocsafe(array,(oldsize)+array_chunksize,vartype*);}
#define shrink_array_1(array,newsize,vartype)\
  if ((((newsize) % array_chunksize)==0) && ((newsize)>0)){array=reallocsafe(array,(newsize),vartype*);}


typedef struct pointT{/* for drawing */
  double x,y;
}pointT;

typedef pointT lineT[100];


typedef struct nodeT{
  double x,y;/* location */
  uint address; /* binary designation of which vertex of the cube I am. */
}nodeT;

#if gd_here
gdImagePtr im_out;

void image_init(gdImagePtr *im_out_ptr)
/* **********************************************************************
********************************************************************** */
{
  /* Input and output images */
  gdImagePtr im_out = 0;
  /* Color indexes */
  int white,blue,red,green;

  im_out = *im_out_ptr;

  /* Create output image, 256 by 256 pixels, true color. */
  im_out = gdImageCreateTrueColor(picwdt, pichgt);
  /* First color allocated is background. */
  white = gdImageColorAllocate(im_out, 255, 255, 255);
  /* Set transparent color. */
  gdImageColorTransparent(im_out, white);

  (*im_out_ptr) = im_out;
}/* image_init */
void image_close(gdImagePtr im_out,char *fname)
/* **********************************************************************
********************************************************************** */
{
  FILE *out;
  /* Make output image interlaced (progressive, in the case of JPEG) */
  gdImageInterlace(im_out, 1);
  out = fopen(fname, "wb");
  /* Write PNG */
  gdImagePng(im_out, out);
  fclose(out);
  gdImageDestroy(im_out);
}/* image_close */

void color_line(gdImagePtr im_out,int x0,int y0,int x1,int y1,uchar r,uchar g,uchar b)
/* **********************************************************************
********************************************************************** */
{
  int color;
  color = gdImageColorAllocate(im_out, r, g, b);
  gdImageLine(im_out, x0, y0, x1, y1, color);
}/* color_line */
#endif

void initgraph()
/* **********************************************************************
********************************************************************** */
{
#if allegro_here
  PALLETE pal;
  int error,ccnt;
#if grx
  /* you should always do this at the start of Allegro programs */
  allegro_init();

  /* set VGA graphics mode 13h (sized 320x200) */
  // error=set_gfx_mode(GFX_AUTODETECT, 800, 600, 0, 0);
  error=set_gfx_mode(GFX_VGA, 320, 200, 0, 0);
  if (error<0){
    exit(1);
  }
  /* set the color pallete */
  /* color 0 = black */
  // pal[0].r = pal[0].g = pal[0].b = 0;
  /* set up a greyscale pallete */
  pal[0].r = pal[0].g = pal[0].b = 0;

  for (ccnt=1; ccnt<256; ccnt++){
    pal[ccnt].r = (ccnt/4);
    pal[ccnt].g = (63);
    pal[ccnt].b = (ccnt/4);
  }
  set_pallete(pal);
#endif
#endif
}/* initgraph */

void poly_line(pointT *polyline,int num_points,int color)
/* **********************************************************************
*  Output a poly-line to the screen, or to a file.
********************************************************************** */
{
  int pcnt;
  int xorg,yorg;
  double oldx,oldy,x,y,wdt,hgt;

//  num_points--;
  xorg=100; yorg=100;
//  wdt=100.0; hgt=100.0;
  wdt=60.0; hgt=60.0;

  xorg=picwdt/2; yorg=pichgt/2;
  // wdt=600.0; hgt=600.0;
  wdt=picwdt*(3.0/10.0); hgt=pichgt*(3.0/10.0);
  x=polyline[0].x; y=polyline[0].y;
  for (pcnt=1;pcnt<num_points;pcnt++){
    oldx=x; oldy=y;
    x=polyline[pcnt].x;
    y=polyline[pcnt].y;
#if grx
    /* line(screen,
      (1+xorg+oldx*wdt), (1+yorg+oldy*hgt),
      (1+xorg+x*wdt), (1+yorg+y*hgt),
      255);*/
    line(screen,
      (1+yorg+oldy*hgt),(1+xorg+oldx*wdt),
      (1+yorg+y*hgt),(1+xorg+x*wdt),
      255-(1+color*50));//255);
#endif
#if gd_here
    color_line(im_out,
      (1+yorg+oldy*hgt),(1+xorg+oldx*wdt), (1+yorg+y*hgt),(1+xorg+x*wdt),
      255,255,255);
#endif
  }
  /* This would be the node. */
//    putpixel(screen, 1+xorg+x*wdt, 1+yorg+y*hgt, 0);
}/* poly_line */


double asin(double sine)
/* **********************************************************************
*  Returns the arc sine.
********************************************************************** */
{
  double cos_sq,cosine,tangeant;

  cos_sq=1.0-(sine*sine);
  cos_sq=fabs(cos_sq);
  cosine = sqrt(cos_sq);
  if (cosine==0.0){cosine=fudge;}
  tangeant = sine/cosine;
  return(atan(tangeant));
}/* asin */

void svg_start(FILE *outfile)
/* **********************************************************************
*  Output header for svg (Scalable Vector Graphics) file.
********************************************************************** */
{
  char outstr[256];
  sprintf(outstr,"<!DOCTYPE svg ><svg width=\"595\" height=\"841\" >\n");
  fprintf(outfile,outstr);
}/* svg_start */

void svg_line(FILE *outfile,pointT *polyline,int num_points)
/* **********************************************************************
*  Output an svg polyline structure.
********************************************************************** */
{
  int pcnt;

  fprintf(outfile," <polyline points=\"");
  for (pcnt=0;pcnt<num_points;pcnt++){
    // fprintf(outfile,"%li,%li ",(int)polyline[pcnt].x,(int)polyline[pcnt].y);
    fprintf(outfile,"%lf,%lf ",polyline[pcnt].x,polyline[pcnt].y);
  }
  // 104,96 171,59 230,128 279,64 343,169 156,193 100,183
  fprintf(outfile,"\"");
  fprintf(outfile," style=\"fill:none ; stroke:rgb(0,0,0) ; stroke-width:1\" />\n");
}/* svg_line */

void svg_end(FILE *outfile)
/* **********************************************************************
*  Output the closer to the svg file.
********************************************************************** */
{
  char outstr[256];
  sprintf(outstr,"</svg>\n");
  fprintf(outfile,outstr);
}/* svg_end */

void project(double x,double y,double cx,double cy,double sradius,double *newx,double *newy)
/* **********************************************************************
*  Project a point onto the surface of a ball, then flatten out the
* ball w/point to relocate that point in a proper warp.
********************************************************************** */
{
  double radius0,radius1,ratio;

  x-=cx; y-=cy;
  radius0=hypot(x,y);
//  frac_around=asin(radius0)/twopi;
//  radius1=asin(radius0)*twopi*sradius;
  radius1=asin(radius0)*sradius;
//  radius1=twopi*sradius;
  ratio=radius1/radius0;
//  ratio=1.0;
  (*newx)=cx+x*ratio; (*newy)=cy+y*ratio;
  /* Here we could normalize to a standard max radius if we want */
}/* project */

void plotline(FILE *outfile,double x0,double y0,double x1,double y1,int num_steps,double cx,double cy,double sradius,int color)
/* **********************************************************************
*  Plot a line segment as a poly-line, warped to be perpendicular to
* the circle rim at both ends.
********************************************************************** */
{
  double newx,newy,xstep,ystep,straightx,straighty,xlen,ylen;
  int cnt;
  pointT *polyline;

  xstep=(x1-x0)/((double)num_steps);  ystep=(y1-y0)/((double)num_steps);

  num_steps+=1;/* add one to include both endpoints */
  polyline=allocsafe(num_steps,pointT);
  for (cnt=0;cnt<num_steps;cnt++){
    straightx = x0+(((double)cnt)*xstep); straighty = y0+(((double)cnt)*ystep);
    project(straightx,straighty,cx,cy,sradius,&newx,&newy);
    /* now output the line. */
    polyline[cnt].x=newx; polyline[cnt].y=newy;
    // polyline[cnt].x=straightx; polyline[cnt].y=straighty;
  }
  /* output the polyline here. */

  /* plot to allegro or gd png file. */
  poly_line(polyline,num_steps,color);
#if svg_here
  svg_line(outfile,polyline,num_steps);
#endif
  freesafe(polyline);
}/* plotline */


int sgn(int num)
/* **********************************************************************
********************************************************************** */
{
  if (num<0){
    return (-1);
  }else{
    return (1);
  }
}/* sgn */
#define flipped 1

void shape_create(FILE *outfile,nodeT *array,double minangle,double maxangle,int mindex,int maxdex,int depth)
/* **********************************************************************
*  This creates the shape so that the array is ordered around the circle.
* IE counting from array index 0 to the end of the array goes
* counterclockwise around the circle, pure and simple.
*   NOTE: this way the array index does not correspond to the
* hypercube vertex binary coordinates (index 101 is NOT vertex 101).
********************************************************************** */
{
  double medangle,border;
  int meddex,cnt,num_points;
  nodeT *node0,*node1;

  border=((maxangle-minangle)/2.0)*0.03;
//  border=0.0;
  depth-=1;/* recursion depth */
  num_points=maxdex-mindex;
  meddex=(mindex+maxdex)/2;
  medangle=(minangle+maxangle)/2.0;
  if (depth>=0){
    shape_create(outfile,array,minangle+border,medangle-border,mindex,meddex,depth);
    shape_create(outfile,array,medangle+border,maxangle-border,meddex+1,maxdex,depth);
  }else{
    node0=&(array[meddex]);
    /* assumed to be centered on 0,0 with a radius of 1.0 */
    node0->x=cos(medangle); node0->y=sin(medangle);
    node0->address=0; /* not worked out yet. */
  }
  /* connect all nodes from min to med with max downto (med+1) */
  num_points/=2;
  for (cnt=0;cnt<=num_points;cnt++){
    node0=&(array[mindex+cnt]);
    node1=&(array[maxdex-cnt]);
    plotline(outfile,node0->x,node0->y,node1->x,node1->y,30, 0.0,0.0, 1.0, depth);
  }
}/* shape_create */


void shape_create2(FILE *outfile,nodeT *array,double minangle,double maxangle,int mindex,int maxdex,int depth)
/* **********************************************************************
*  This creates the shape so that array index 000 is
*  vertex coordinate 000, and array index 111 is vertex 111, etc.
********************************************************************** */
{
  double medangle,border;
  int meddex,cnt,num_points;
  nodeT *node0,*node1;

  border=((maxangle-minangle)/2.0)*0.03;
  // border=0.0;
  depth-=1;/* recursion depth */
  meddex=(mindex+maxdex+0)/2;
  medangle=(minangle+maxangle)/2.0;
  if (depth>=0){
    shape_create2(outfile,array,minangle+border,medangle-border,mindex,meddex,depth);
    shape_create2(outfile,array,maxangle-border,medangle+border,meddex+1,maxdex,depth);
  }else{
    node0=&(array[meddex]);
    /* assumed to be centered on 0,0 with a radius of 1.0 */
    node0->x=cos(medangle); node0->y=sin(medangle);
    node0->address=0; /* not worked out yet. */
  }
  /* connect all nodes from min to med with max downto (med+1) */
  num_points=1+abs(maxdex-mindex);
  num_points/=2;
  for (cnt=0;cnt<num_points;cnt++){
    node0=&(array[mindex+cnt]);
    node1=&(array[1+meddex+cnt]);
    plotline(outfile,node0->x,node0->y,node1->x,node1->y,30, 0.0,0.0, 1.0, depth);
  }
}/* shape_create2 */


int main()
/* **********************************************************************
********************************************************************** */
{
  FILE *outfile;
  nodeT *array;
  int num_dims,num_elements,cnt;

//  printf("%lf\n",asin(PI)); getchar();

  initgraph();
#if gd_here
  image_init(&im_out);
#endif

#if svg_here
  outfile = fopen("hcube.svg","w"); svg_start(outfile);
#endif

  num_dims=2;
  num_dims=7;
  // num_dims=10;
  num_elements=(1<<num_dims);
  array=allocsafe(num_elements,nodeT);
  for (cnt=0;cnt<num_elements;cnt++){
    array[cnt].x=-123.0; array[cnt].y=-123.4;
  }
  shape_create2(outfile,array,0.0,twopi,0,num_elements-1,num_dims);
  freesafe(array);

#if svg_here
  svg_end(outfile); fclose(outfile);
#endif

#if gd_here
  image_close(im_out,"hcube.png");
#endif

  return(0);
}/* main */

#if allegro_here
END_OF_MAIN();
#endif

#if 0

  Junkyard

<!DOCTYPE svg ><svg width="595" height="841" >
 <polyline points="104,96 171,59 230,128 279,64 343,169 156,193 100,183 " style="fill:none ; stroke:rgb(0,0,0) ; stroke-width:1" />
</svg>

<!DOCTYPE svg ><svg width="595" height="841" >
 <polyline points="104,96 171,59 230,128 279,64 343,169 156,193 100,183 " style="fill:none ; stroke:rgb(0,0,0) ; stroke-width:1" />
 <polyline points="104,96 171,59 230,128 279,64 343,169 156,193 100,183 " style="fill:none ; stroke:rgb(0,0,0) ; stroke-width:1" />
</svg>


#include <stdio.h>
#include "gd.h"
#include "gdfontg.h"
#include "gdfonts.h"


/* Try to load demoin.png and paste part of it into the
output image. */
in = fopen("scale.png", "rb");
if (!in) {
  fprintf(stderr, "Can't load source image; this demo\n");
  fprintf(stderr, "is much more impressive if demoin.png\n");
  fprintf(stderr, "is available.\n");
  im_in = 0;
} else {
  im_in = gdImageCreateFromPng(in);
  fclose(in);
  /* Now copy, and magnify as we do so */
  gdImageCopyResized(im_out, im_in,32, 32, 0, 0, 700, 700, 1000, 1000);
}
red = gdImageColorAllocate(im_out, 255, 0, 0);
green = gdImageColorAllocate(im_out, 0, 255, 0);
blue = gdImageColorAllocate(im_out, 0, 0, 255);
/* Rectangle */
gdImageLine(im_out, 16, 16, 240, 16, green);
gdImageLine(im_out, 240, 16, 240, 240, green);
gdImageLine(im_out, 240, 240, 16, 240, green);
gdImageLine(im_out, 16, 240, 16, 16, green);
/* Circle */
gdImageArc(im_out, 128, 128, 60, 20, 0, 720, blue);
/* Arc */
gdImageArc(im_out, 128, 128, 40, 40, 90, 270, blue);
/* Flood fill: doesn't do much on a continuously
variable tone jpeg original. */
gdImageFill(im_out, 8, 8, blue);
/* Polygon */
points[0].x = 64;
points[0].y = 0;
points[1].x = 0;
points[1].y = 128;
points[2].x = 128;
points[2].y = 128;
gdImageFilledPolygon(im_out, points, 3, green);

/* Text */
gdImageString(im_out, gdFontGiant, 32, 32,(unsigned char *) "hi", red);
gdImageStringUp(im_out, gdFontSmall, 64, 64,(unsigned char *) "hi", red);




#endif
