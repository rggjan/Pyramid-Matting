#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <math.h>

using namespace std;

unsigned char*
load_image (const char* filename, int dimx, int dimy, int num_colors) {
  unsigned char* data = new unsigned char[dimx*dimy*num_colors];

  FILE *fp = fopen (filename, "rb");

  fread(data, 1, dimx*dimy*num_colors, fp);
  
  fclose (fp);

  return data;
}


void
save_image (const char* filename, int dimx, int dimy, int num_colors, unsigned char* data) {
  FILE *fp = fopen (filename, "wb");
  
  if (num_colors == 3)
    {
      fprintf (fp, "P6\n%d %d\n255\n", dimx, dimy);
    }
  else if (num_colors == 1)
    {
      fprintf (fp, "P5\n%d %d\n255\n", dimx, dimy);
    }
  else
    {
      printf ("Problem!\n");
      exit (1);
    }

  fwrite(data, 1, dimx*dimy*num_colors, fp);
  
  fclose (fp);
}

// a1 and f1 known
double solve_equations(unsigned char* f0, unsigned char* f1, unsigned char* f2,
    unsigned char* b0, unsigned char* b1, unsigned char* b2,
    unsigned char* a0, unsigned char* a1, unsigned char* a2,
    unsigned char* c1, unsigned char* c2) {
  int nb1[3];
  int nf2[3];
  int nb2[3];
  a2[0] = 2*a0[0]-a1[0];
  for (int c=0; c<3; c++) {
    nb1[c] = (c1[c] - f1[c]*a1[0]/255.)/(1-a1[0]/255.);
    nf2[c] = (2*f0[c]*a0[0]/255. - f1[c]*a1[0]/255.)/(2*a0[0]/255.-a1[0]/255.);
    nb2[c] = (c2[c] - (2*f0[c]*a0[0]-f1[c]*a1[0]))/(1-2*a0[0]+a1[0]);
  }
  
  double quality = 0;
  for (int c=0; c<3; c++) {

    double diff = nb1[c]-255;
    if (diff > 0) {
      quality -= diff*diff;
      nb1[c] = 255;
    } else if (nb1[c] < 0) {
      quality -= nb1[c]*nb1[c];
      nb1[c] = 0;
    }

    diff = nf2[c]-255;
    if (diff > 0) {
      quality -= diff*diff;
      nf2[c] = 255;
    } else if (nf2[c] < 0) {
      quality -= nf2[c]*nf2[c];
      nf2[c] = 0;
    }
    
    diff = nb2[c]-255;
    if (diff > 0) {
      quality -= diff*diff;
      nb2[c] = 255;
    } else if (nb2[c] < 0) {
      quality -= nb2[c]*nb2[c];
      nb2[c] = 0;
    }

    b1[c] = nb1[c];
    f2[c] = nf2[c];
    b2[c] = nb2[c];
  }

  if (quality == 0) {
    for (int c=0; c<3; c++) {
      int diff = f1[c]-f0[c];
      quality += diff*diff;
      
      diff = f2[c]-f0[c];
      quality += diff*diff;
      
      diff = b1[c]-b0[c];
      quality += diff*diff;
      
      diff = b2[c]-b0[c];
      quality += diff*diff;
    }
      
    int diff = a1[0]-a0[0];
    quality += diff*diff;

    diff = a2[0]-a0[0];
    quality += diff*diff;
  }
  return quality;
}

// a1 and f1 unknown
void optimize(unsigned char* f0, unsigned char* f1, unsigned char* f2,
    unsigned char* b0, unsigned char* b1, unsigned char* b2,
    unsigned char* a0, unsigned char* a1, unsigned char* a2,
    unsigned char* c1, unsigned char* c2) {

  f1[0] = f0[0];
  f1[1] = f0[1];
  f1[2] = f0[2];
  a1[0] = 128;
//  int alpha = 128;

  for (int i=0; i<100; i++) {
    // alpha
    int qplus, qminus;
    if (a1[0] != 255)
      a1[0]++;
    qplus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
    if (a1[0] != 255)
      a1[0]--;

    if (a1[0] != 0)
      a1[0]--;
    qminus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
    if (a1[0] != 0)
      a1[0]++;

    //cout << (int)f1[0] << "/" << (int)f1[1] << "/" << (int)f1[2] << "/" << (int)a1[0]+1 << ": " << qplus << endl;

    if (qplus > 0)
      if (qminus > 0)
        if (qplus < qminus)
          a1[0]++;
        else
          a1[0]--;
      else
        a1[0]++;
    else
      if (qminus > 0)
        a1[0]--;
      else
        if (qminus > qplus)
          a1[0]--;
        else
          a1[0]++;

    for (int c=0; c<3; c++) {
      int qplus, qminus;
      if (f1[c] != 255)
        f1[c]++;
      qplus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
      if (f1[c] != 255)
        f1[c]--;

      if (f1[c] != 0)
        f1[c]--;
      qminus = solve_equations(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);
      if (f1[c] != 0)
        f1[c]++;

      if (qplus > 0)
        if (qminus > 0)
          if (qplus < qminus)
            f1[c]++;
          else
            f1[c]--;
        else
          f1[c]++;
      else
        if (qminus > 0)
          f1[c]--;
        else
          if (qminus > qplus)
            f1[c]--;
          else
            f1[c]++;
      }
  }
  
  //exit(1);

  //cout << "alpha 0/1/2: " << (int)a0[0] << "/" << (int)a1[0] << "/" << (int)a2[0] << "\n";
}

// Project the Point P onto the line from A-B.
// alpha_pointer receives the calculated alpha value between 0 and 1
// where (0=A, 1=B)
// returns the squared distance of the best projected point to P
float projection (unsigned char B[3], unsigned char A[3], unsigned char P[3], unsigned char* alpha_pointer) {
  int ABx = B[0] - A[0];
  int ABy = B[1] - A[1];
  int ABz = B[2] - A[2];

  int APx = P[0] - A[0];
  int APy = P[1] - A[1];
  int APz = P[2] - A[2];

  float len_AB = sqrtf(ABx * ABx + ABy * ABy + ABz * ABz);
  float len_AP = sqrtf(APx * APx + APy * APy + APz * APz);
  float alpha = (ABx * APx + ABy * APy + ABz * APz) / (len_AB * len_AP);

  float PPx, PPy, PPz;

  PPx = P[0] - (A[0] + alpha * ABx);
  PPy = P[1] - (A[1] + alpha * ABy);
  PPz = P[2] - (A[2] + alpha * ABz);

  alpha = alpha > 1 ? 1 : (alpha < 0 ? 0 : alpha);

  if (alpha_pointer) {
    *alpha_pointer = (unsigned char)(alpha*255);
  }

  // Normalize, so that it is in the unit cube of colors
  return (PPx * PPx + PPy * PPy + PPz * PPz) / (255. * 255.);
}

#define RESULTS "results/"

int main() {
  // ensure the "result" directory exists
  struct stat result_dir;
  if (stat(RESULTS, &result_dir) < 0) {
    int err = errno;
    if (err == ENOENT) {
      // create directory
      mkdir(RESULTS, 0777);
    } else {
      cerr << "Error while testing existence of 'results' directory: " << strerror(err) << '\n';
      return 1;
    }
  }
  if (!S_ISDIR(result_dir.st_mode)) {
    cerr << "The file 'results' is not a directory.\n";
    return 1;
  }

  int width = 512;
  int height = 512;
  int raise = 9;

  unsigned char* mask = load_image("trimap.pnm", width, height, 1);

  unsigned char* originals[raise+1][2];
  unsigned char* new_alphas[raise+1][2];
  unsigned char* foregrounds[raise+1][2];
  unsigned char* new_foregrounds[raise+1][2];
  double* foreground_ps[raise+1][2];
  unsigned char* backgrounds[raise+1][2];
  unsigned char* new_backgrounds[raise+1][2];
  double* background_ps[raise+1][2];

  originals[9][0] = load_image("test.ppm", width, height, 3);
  foregrounds[9][0] = new unsigned char[width*height*3];
  foreground_ps[9][0] = new double[width*height];
  backgrounds[9][0] = new unsigned char[width*height*3];
  background_ps[9][0] = new double[width*height];

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      if (mask[width*y+x] == 255) {
        foreground_ps[9][0][width*y+x] = 1;
        background_ps[9][0][width*y+x] = 0;
        for (int c=0; c<3; c++) {
          backgrounds[9][0][(width*y+x)*3+c] = 0;
          foregrounds[9][0][(width*y+x)*3+c] =
            originals[9][0][(width*y+x)*3+c];
        }
      } else if (mask[width*y+x] == 0) {
        foreground_ps[9][0][width*y+x] = 0;
        background_ps[9][0][width*y+x] = 1;
        for (int c=0; c<3; c++) {
          foregrounds[9][0][(width*y+x)*3+c] = 0;
          backgrounds[9][0][(width*y+x)*3+c] =
            originals[9][0][(width*y+x)*3+c];
        }
      } else {
        foreground_ps[9][0][width*y+x] = 0;
        background_ps[9][0][width*y+x] = 0;
        for (int c=0; c<3; c++) {
          foregrounds[9][0][(width*y+x)*3+c] = 0;
          backgrounds[9][0][(width*y+x)*3+c] = 0;
        }
      }
    }
  }

  save_image(RESULTS "originals_9.ppm", width, height, 3, originals[9][0]);
  save_image(RESULTS "foregrounds_9.ppm", width, height, 3, foregrounds[9][0]);
  save_image(RESULTS "backgrounds_9.ppm", width, height, 3, backgrounds[9][0]);

  while (raise > 0) {
    raise--;

    // Height half
    height = height/2;
    originals[raise][1] = new unsigned char[width*height*3];
    foregrounds[raise][1] = new unsigned char[width*height*3];
    foreground_ps[raise][1] = new double[width*height];
    backgrounds[raise][1] = new unsigned char[width*height*3];
    background_ps[raise][1] = new double[width*height];
  
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        foreground_ps[raise][1][y*width+x] =
          (foreground_ps[raise+1][0][2*y*width+x] +
           foreground_ps[raise+1][0][(2*y+1)*width+x])/2;
        background_ps[raise][1][y*width+x] =
          (background_ps[raise+1][0][2*y*width+x] +
           background_ps[raise+1][0][(2*y+1)*width+x])/2;

        for (int c=0; c<3; c++) {
          originals[raise][1][(y*width+x)*3+c] =
            (originals[raise+1][0][(2*y*width+x)*3+c] +
             originals[raise+1][0][((2*y+1)*width+x)*3+c])/2;

          foregrounds[raise][1][(y*width+x)*3+c] =
            (foregrounds[raise+1][0][(2*y*width+x)*3+c]*
             foreground_ps[raise+1][0][2*y*width+x] +
             foregrounds[raise+1][0][((2*y+1)*width+x)*3+c]*
             foreground_ps[raise+1][0][(2*y+1)*width+x])/
             foreground_ps[raise][1][y*width+x]/2;
          
          backgrounds[raise][1][(y*width+x)*3+c] =
            (backgrounds[raise+1][0][(2*y*width+x)*3+c]*
             background_ps[raise+1][0][2*y*width+x] +
             backgrounds[raise+1][0][((2*y+1)*width+x)*3+c]*
             background_ps[raise+1][0][(2*y+1)*width+x])/
             background_ps[raise][1][y*width+x]/2;
        }
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "originals_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise][1]);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, foregrounds[raise][1]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, backgrounds[raise][1]);

    // Width half
    width = width/2;
    originals[raise][0] = new unsigned char[width*height*3];
    foregrounds[raise][0] = new unsigned char[width*height*3];
    foreground_ps[raise][0] = new double[width*height];
    
    backgrounds[raise][0] = new unsigned char[width*height*3];
    background_ps[raise][0] = new double[width*height];

    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        foreground_ps[raise][0][y*width+x] =
          (foreground_ps[raise][1][y*width*2+x*2] +
           foreground_ps[raise][1][y*width*2+x*2+1])/2;

        background_ps[raise][0][y*width+x] =
          (background_ps[raise][1][y*width*2+x*2] +
           background_ps[raise][1][y*width*2+x*2+1])/2;

        for (int c=0; c<3; c++) {
          originals[raise][0][(y*width+x)*3+c] =
            (originals[raise][1][(y*width*2+2*x)*3+c] +
             originals[raise][1][(y*width*2+2*x+1)*3+c])/2;

          foregrounds[raise][0][(y*width+x)*3+c] =
            (foregrounds[raise][1][(y*width*2+x*2)*3+c]*
             foreground_ps[raise][1][y*width*2+x*2] +
             foregrounds[raise][1][(y*width*2+x*2+1)*3+c]*
             foreground_ps[raise][1][y*width*2+x*2+1])/
             foreground_ps[raise][0][y*width+x]/2;
          
          backgrounds[raise][0][(y*width+x)*3+c] =
            (backgrounds[raise][1][(y*width*2+x*2)*3+c]*
             background_ps[raise][1][y*width*2+x*2] +
             backgrounds[raise][1][(y*width*2+x*2+1)*3+c]*
             background_ps[raise][1][y*width*2+x*2+1])/
             background_ps[raise][0][y*width+x]/2;
        }
      }
    }

    snprintf(buffer, 100, RESULTS "originals_%i.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise][0]);
    
    snprintf(buffer, 100, RESULTS "foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, foregrounds[raise][0]);
    
    snprintf(buffer, 100, RESULTS "backgrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, backgrounds[raise][0]);
  }

  // Upscaling
  raise = 0;
  width = 1;
  height = 1;

  new_foregrounds[0][0] = new unsigned char[3];
  new_backgrounds[0][0] = new unsigned char[3];
  new_alphas[0][0] = new unsigned char[1];

  for (int c=0; c<3; c++) {
    new_foregrounds[0][0][c] = foregrounds[0][0][c];
    new_backgrounds[0][0][c] = backgrounds[0][0][c];
  }

  // Calculate alpha
  unsigned char merged_point[3];
  for (int c=0; c<3; c++) {
    merged_point[c] = (originals[0][0][c]
      -foreground_ps[0][0][0]*foregrounds[0][0][c]
      -background_ps[0][0][0]*backgrounds[0][0][c])
      /(1-foreground_ps[0][0][0]-background_ps[0][0][0]);
  }
  projection(new_foregrounds[0][0], new_backgrounds[0][0], merged_point, &(new_alphas[0][0][0]));

  save_image(RESULTS "new_foregrounds_0.ppm", 1, 1, 3, new_foregrounds[0][0]);

  

  while (raise <= 9) {
    cout << "first " << raise << endl;
    width *= 2;
    new_foregrounds[raise][1] = new unsigned char[width*height*3];
    new_backgrounds[raise][1] = new unsigned char[width*height*3];
    new_alphas[raise][1] = new unsigned char[width*height];

    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x+=2) {
        unsigned char* f0 = &(new_foregrounds[raise][0][(y*width/2+x/2)*3]);
        unsigned char* b0 = &(new_backgrounds[raise][0][(y*width/2+x/2)*3]);
        unsigned char* a0 = &(new_alphas[raise][0][(y*width/2+x/2)]);

        unsigned char* f1 = &(new_foregrounds[raise][1][(y*width+x)*3]);
        unsigned char* f2 = &(new_foregrounds[raise][1][(y*width+x+1)*3]);
        
        unsigned char* b1 = &(new_backgrounds[raise][1][(y*width+x)*3]);
        unsigned char* b2 = &(new_backgrounds[raise][1][(y*width+x+1)*3]);
        
        unsigned char* a1 = &(new_alphas[raise][1][y*width+x]);
        unsigned char* a2 = &(new_alphas[raise][1][y*width+x+1]);
  
        // Subtract known foreground / background 
        unsigned char* original1 = &(originals[raise][1][(y*width+x)*3]);
        unsigned char* original2 = &(originals[raise][1][(y*width+x+1)*3]);
       

        unsigned char* foreground1 = &(foregrounds[raise][1][(y*width+x)*3]);
        unsigned char* foreground2 = &(foregrounds[raise][1][(y*width+x+1)*3]);
        
        unsigned char* background1 = &(backgrounds[raise][1][(y*width+x)*3]);
        unsigned char* background2 = &(backgrounds[raise][1][(y*width+x+1)*3]);

        double* foreground_ps1 = &(foreground_ps[raise][1][(y*width+x)]);
        double* foreground_ps2 = &(foreground_ps[raise][1][(y*width+x+1)]);
        
        double* background_ps1 = &(background_ps[raise][1][(y*width+x)]);
        double* background_ps2 = &(background_ps[raise][1][(y*width+x+1)]);

        unsigned char c1[3];
        unsigned char c2[3];

        for (int c=0; c<3; c++) {
          c1[c] = (original1[c] - foreground_ps1[0]*foreground1[c]
              -background_ps1[0]*background1[c])
              /(1-foreground_ps1[0]-background_ps1[0]);
          c2[c] = (original2[c] - foreground_ps2[0]*foreground2[c]
              -background_ps2[0]*background2[c])
              /(1-foreground_ps2[0]-background_ps2[0]);
        }
        
//        cout << foreground_ps1[0] << endl << endl;
//        cout << foreground_ps2[0] << endl << endl;
/*        for (int c=0; c<3; c++) {
          cout << (int)c1[c] <<  "," << (int)c2[c] << endl;
        }*/

        optimize(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);

        /*cout << "f0\tf1\tf2\tb0\tb1\tb2\tc1\tc2\n";
        for (int c=0; c<3; c++) {
          cout << (int)f0[c] << "\t" << (int)f1[c] << "\t" << (int)f2[c] << "\t" <<
            (int)b0[c] << "\t" << (int)b1[c] << "\t" << (int)b2[c] << "\t" <<
            (int)c1[c] << "\t" << (int)c2[c] << "\t" << "\n";
        }
        cout << (int)a1[0] << "," << (int)a2[0] << endl;*/
      }
    }

    static char buffer[100];
    snprintf(buffer, 100, RESULTS "final_%i_h.ppm", raise);
    unsigned char *tmp = new unsigned char[width*height*3];
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int c=0; c<3; c++) {
          tmp[(y*width+x)*3+c] =
            foregrounds[raise][1][(y*width+x)*3+c]*
            foreground_ps[raise][1][(y*width+x)]+

            (new_foregrounds[raise][1][(y*width+x)*3+c]
            *new_alphas[raise][1][y*width+x]/255.)*

            (1-foreground_ps[raise][1][(y*width+x)]-
            background_ps[raise][1][(y*width+x)]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, new_foregrounds[raise][1]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i_h.ppm", raise);
    save_image(buffer, width, height, 1, new_alphas[raise][1]);
    
    cout << "second " << raise << endl;
    // Stretch height
    height *= 2;
    new_foregrounds[raise+1][0] = new unsigned char[width*height*3];
    new_backgrounds[raise+1][0] = new unsigned char[width*height*3];
    new_alphas[raise+1][0] = new unsigned char[width*height];

    for (int y=0; y<height; y+=2) {
      for (int x=0; x<width; x++) {
        unsigned char* f0 = &(new_foregrounds[raise][1][((y/2)*width+x)*3]);
        unsigned char* b0 = &(new_backgrounds[raise][1][((y/2)*width+x)*3]);
        unsigned char* a0 = &(new_alphas[raise][1][((y/2)*width+x)]);

        unsigned char* f1 = &(new_foregrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* f2 = &(new_foregrounds[raise+1][0][((y+1)*width+x)*3]);
        
        unsigned char* b1 = &(new_backgrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* b2 = &(new_backgrounds[raise+1][0][((y+1)*width+x)*3]);
        
        unsigned char* a1 = &(new_alphas[raise+1][0][y*width+x]);
        unsigned char* a2 = &(new_alphas[raise+1][0][(y+1)*width+x]);

        // Subtract known foreground / background 
        unsigned char* original1 = &(originals[raise+1][0][(y*width+x)*3]);
        unsigned char* original2 = &(originals[raise+1][0][((y+1)*width+x)*3]);
       
        unsigned char* foreground1 = &(foregrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* foreground2 = &(foregrounds[raise+1][0][((y+1)*width+x)*3]);
        
        unsigned char* background1 = &(backgrounds[raise+1][0][(y*width+x)*3]);
        unsigned char* background2 = &(backgrounds[raise+1][0][((y+1)*width+x)*3]);

        double* foreground_ps1 = &(foreground_ps[raise+1][0][((y+1)*width+x)]);
        double* foreground_ps2 = &(foreground_ps[raise+1][0][((y+1)*width+x)]);
        
        double* background_ps1 = &(background_ps[raise+1][0][((y+1)*width+x)]);
        double* background_ps2 = &(background_ps[raise+1][0][((y+1)*width+x)]);

        unsigned char c1[3];
        unsigned char c2[3];

        for (int c=0; c<3; c++) {
          c1[c] = (original1[c] - foreground_ps1[0]*foreground1[c]
              -background_ps1[0]*background1[c])
              /(1-foreground_ps1[0]-background_ps1[0]);
          c2[c] = (original2[c] - foreground_ps2[0]*foreground2[c]
              -background_ps2[0]*background2[c])
              /(1-foreground_ps2[0]-background_ps2[0]);
        }
        
        optimize(f0, f1, f2, b0, b1, b2, a0, a1, a2, c1, c2);

/*        cout << "f0\tf1\tf2\tb0\tb1\tb2\tc1\tc2\n";
        for (int c=0; c<3; c++) {
          cout << (int)f0[c] << "\t" << (int)f1[c] << "\t" << (int)f2[c] << "\t" <<
            (int)b0[c] << "\t" << (int)b1[c] << "\t" << (int)b2[c] << "\t" <<
            (int)c1[c] << "\t" << (int)c2[c] << "\t" << "\n";
        }*/
      }
    }

    snprintf(buffer, 100, RESULTS "final_%i.ppm", raise+1);
    tmp = new unsigned char[width*height*3];
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
        for (int c=0; c<3; c++) {
          tmp[(y*width+x)*3+c] =
            foregrounds[raise+1][0][(y*width+x)*3+c]*
            foreground_ps[raise+1][0][(y*width+x)]+

            (new_foregrounds[raise+1][0][(y*width+x)*3+c]
            *new_alphas[raise+1][0][y*width+x]/255.)*

            (1-foreground_ps[raise+1][0][(y*width+x)]-
            background_ps[raise+1][0][(y*width+x)]);
        }
      }
    }
    save_image(buffer, width, height, 3, tmp);
    delete[] tmp;

    snprintf(buffer, 100, RESULTS "new_foregrounds_%i.ppm", raise+1);
    save_image(buffer, width, height, 3, new_foregrounds[raise+1][0]);
    
    snprintf(buffer, 100, RESULTS "new_alphas_%i.ppm", raise+1);
    save_image(buffer, width, height, 1, new_alphas[raise+1][0]);

    raise++;
  }
}
