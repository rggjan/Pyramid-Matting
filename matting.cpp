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

int main() {
  int width = 512;
  int height = 512;
  int raise = 9;

  unsigned char* mask = load_image("trimap.pnm", width, height, 1);

  unsigned char* originals[raise][2];
  unsigned char* foregrounds[raise][2];
  double* foreground_ps[raise][2];
  unsigned char* backgrounds[raise][2];
  double* background_ps[raise][2];

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
