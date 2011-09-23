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

  unsigned char* original = load_image("test.ppm", width, height, 3);
  unsigned char* mask = load_image("trimap.pnm", width, height, 1);
  unsigned char* final = new unsigned char[width*height*3];

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

  save_image("final.ppm", width, height, 3, final);

  delete[] original;
  delete[] mask;
  delete[] final;
}
