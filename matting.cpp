#include <stdio.h>
#include <stdlib.h>


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

  int dot_AB = ABx * ABx + ABy * ABy + ABz * ABz;
  float alpha = (float)(ABx * APx + ABy * APy + ABz * APz) / dot_AB;

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

int main() {
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

  save_image("results/final_9.ppm", width, height, 3, originals[9][0]);
  save_image("results/foregrounds_9.ppm", width, height, 3, foregrounds[9][0]);
  save_image("results/backgrounds_9.ppm", width, height, 3, backgrounds[9][0]);

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
    snprintf(buffer, 100, "results/final_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise][1]);
    
    snprintf(buffer, 100, "results/foregrounds_%i_h.ppm", raise);
    save_image(buffer, width, height, 3, foregrounds[raise][1]);
    
    snprintf(buffer, 100, "results/backgrounds_%i_h.ppm", raise);
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

    snprintf(buffer, 100, "results/final_%i.ppm", raise);
    save_image(buffer, width, height, 3, originals[raise][0]);
    
    snprintf(buffer, 100, "results/foregrounds_%i.ppm", raise);
    save_image(buffer, width, height, 3, foregrounds[raise][0]);
    
    snprintf(buffer, 100, "results/backgrounds_%i.ppm", raise);
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

/*  while (raise <= 9) {
    for (int y=0; y<height; y++) {
      for (int x=0; x<width; x++) {
      }
    }*/


  printf("Alpha: %i\n", new_alphas[0][0][0]);


/*  unsigned char* final = new unsigned char[width*height*3];*/

  /*

  int fg[3];
  fg[0] = 0;
  fg[1] = 0;
  fg[2] = 0;
  int fg_num = 0;

  int bg[3];
  bg[0] = 0;
  bg[1] = 0;
  bg[2] = 0;
  int bg_num = 0;

  for (int y=0; y<height; y++) {
    for (int x=0; x<width; x++) {
      int index = y*width+x;
      if (mask[index] == 0) {
        bg_num++;
        for (int color=0; color<3; color++) {
          bg[color] += original[index*3+color];
        }
      } else if (mask[index] == 255) {
        fg_num++;
        for (int color=0; color<3; color++) {
          fg[color] += original[index*3+color];
        }
      }
    }
  }

  printf("fg: %i/%i/%i\n", fg[0]/fg_num, fg[1]/fg_num, fg[2]/fg_num);
  printf("%%fg: %f\n", ((double)fg_num)/((double)width*height)*100);

  printf("bg: %i/%i/%i\n", bg[0]/bg_num, bg[1]/bg_num, bg[2]/bg_num);
  printf("%%bg: %f\n", ((double)bg_num)/((double)width*height)*100);

  save_image("final.ppm", width, height, 3, final);*/

  /*delete[] original;
  delete[] mask;
  delete[] final;*/
}
