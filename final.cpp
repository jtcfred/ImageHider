#include "matrix.h"
#include <OpenImageIO/imageio.h>
#include <iostream>
#include <GL/glut.h>
#include <math.h>
#include <vector>
 
using namespace std;
OIIO_NAMESPACE_USING
 
 
struct Pixel{ // defines a pixel structure
	unsigned char r,g,b,a;
}; 
 
//
// Global variables and constants
//
const int DEFAULTWIDTH = 600;	// default window dimensions if no image
const int DEFAULTHEIGHT = 600;
 
int WinWidth, WinHeight;	// window width and height
int ImWidth, ImHeight, ImWidth2, ImHeight2;		// image width and height for both
int scaledWidth, scaledHeight;
int newImWidth, newImHeight;   //warped image width and height
int ImChannels;           // number of channels per image pixel
 
int VpWidth, VpHeight;		// viewport width and height
int Xoffset, Yoffset;     // viewport offset from lower left corner of window
 
Pixel **pixmap = NULL;  // the image pixmap used for OpenGL display
Pixel **pixmap2 = NULL;  //the pixmap used for the hidden image
Pixel **scaledHidden = NULL; // pixmap used for the hidden image post scaling
Pixel **warpedPixmap = NULL;  //pixmap used for warping
Pixel **separated = NULL; //pixmap for the extracted image
int pixformat; 			// the pixel format used to correctly  draw the image
 
int shiftedX = 0;
int shiftedY = 0; //these hold information on how much the origin has shifted
 
int scaleShiftedX = 0;
int scaleShiftedY = 0;
 
//
//  Routine to cleanup the memory.   
//
void destroy(Pixel **pixels){
 if (pixels){
     delete pixels[0];
	 delete pixels;  
  }
}
 
 
//
//  Routine to read an image file and store in a pixmap
//  returns the size of the image in pixels if correctly read, or 0 if failure
//
int readImage(string infilename){
  // Create the oiio file handler for the image, and open the file for reading the image.
  // Once open, the file spec will indicate the width, height and number of channels.
  std::unique_ptr<ImageInput> infile = ImageInput::open(infilename);
  if(!infile){
    cerr << "Could not input image file " << infilename << ", error = " << geterror() << endl;
    return 0;
  }
 
  // Record image width, height and number of channels in global variables
  if(pixmap == NULL) {
    ImWidth = infile->spec().width;
    ImHeight = infile->spec().height;
    ImChannels = infile->spec().nchannels;
    unsigned char tmp_pixels[ImWidth * ImHeight * ImChannels];
    int scanlinesize = ImWidth * ImChannels * sizeof(unsigned char);
    if(!infile->read_image(TypeDesc::UINT8, &tmp_pixels[0] + (ImHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
      cerr << "Could not read image from " << infilename << ", error = " << geterror() << endl;
      return 0;
    }

    //This is if there is no image already, this becomes the first image
    pixmap = new Pixel*[ImHeight];
    if(pixmap != NULL)
    pixmap[0] = new Pixel[ImWidth * ImHeight];
    for(int i = 1; i < ImHeight; i++)
    pixmap[i] = pixmap[i - 1] + ImWidth;
  
    //  assign the read pixels to the the data structure
    int index;
    for(int row = 0; row < ImHeight; ++row) {
      for(int col = 0; col < ImWidth; ++col) {
        index = (row*ImWidth+col)*ImChannels;
        
        if (ImChannels==1){ 
          pixmap[row][col].r = tmp_pixels[index];
          pixmap[row][col].g = tmp_pixels[index];
          pixmap[row][col].b = tmp_pixels[index];
          pixmap[row][col].a = 255;
        }
        else{
          pixmap[row][col].r = tmp_pixels[index];
          pixmap[row][col].g = tmp_pixels[index+1];
          pixmap[row][col].b = tmp_pixels[index+2];			
          if (ImChannels <4) // no alpha value is present so set it to 255
            pixmap[row][col].a = 255; 
          else // read the alpha value
            pixmap[row][col].a = tmp_pixels[index+3];			
        }
      }
    }
  } else {
    ImWidth2 = infile->spec().width;
    ImHeight2 = infile->spec().height;
    ImChannels = infile->spec().nchannels;
    unsigned char tmp_pixels[ImWidth2 * ImHeight2 * ImChannels];
    int scanlinesize = ImWidth2 * ImChannels * sizeof(unsigned char);
    if(!infile->read_image(TypeDesc::UINT8, &tmp_pixels[0] + (ImHeight2 - 1) * scanlinesize, AutoStride, -scanlinesize)){
      cerr << "Could not read image from " << infilename << ", error = " << geterror() << endl;
      return 0;
    }

    pixmap2 = new Pixel*[ImHeight2];
    if(pixmap2 != NULL)
    pixmap2[0] = new Pixel[ImWidth2 * ImHeight2];
    for(int i = 1; i < ImHeight2; i++)
    pixmap2[i] = pixmap2[i - 1] + ImWidth2;
 
    //  assign the read pixels to the the data structure
    int index;
    for(int row = 0; row < ImHeight2; ++row) {
      for(int col = 0; col < ImWidth2; ++col) {
        index = (row*ImWidth2+col)*ImChannels;
        
        if (ImChannels==1){ 
          pixmap2[row][col].r = tmp_pixels[index];
          pixmap2[row][col].g = tmp_pixels[index];
          pixmap2[row][col].b = tmp_pixels[index];
          pixmap2[row][col].a = 255;
        }
        else{
          pixmap2[row][col].r = tmp_pixels[index];
          pixmap2[row][col].g = tmp_pixels[index+1];
          pixmap2[row][col].b = tmp_pixels[index+2];			
          if (ImChannels <4) // no alpha value is present so set it to 255
            pixmap2[row][col].a = 255; 
          else // read the alpha value
            pixmap2[row][col].a = tmp_pixels[index+3];			
        }
      }
    }
  }
 
  // close the image file after reading, and free up space for the oiio file handler
  infile->close();
  
  // set the pixel format to GL_RGBA and fix the # channels to 4  
  pixformat = GL_RGBA;  
  ImChannels = 4;
 
  // return image size in pixels
  return ImWidth * ImHeight;
}


//copies pixels from im1 onto im2
void copyPixels(Pixel**& im1, Pixel**& im2) {
  for(int row = 0; row < newImHeight; row++) {
    for(int col = 0; col < newImWidth; col++) {
      im2[row][col].r = im1[row][col].r;
      im2[row][col].g = im1[row][col].g;
      im2[row][col].b = im1[row][col].b;
      im2[row][col].a = 255;
    }
  }
}

//
// Routine to display a pixmap in the current window
//
void displayImage(){
  // if the window is smaller than the image, scale it down, otherwise do not scale
  if(WinWidth < newImWidth  || WinHeight < newImHeight)
    glPixelZoom(float(VpWidth) / newImWidth, float(VpHeight) / newImHeight);
  else
    glPixelZoom(1.0, 1.0);
  
  // display starting at the lower lefthand corner of the viewport
  glRasterPos2i(0, 0);
 
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glDrawPixels(newImWidth, newImHeight, pixformat, GL_UNSIGNED_BYTE, warpedPixmap[0]);
}
 
 
//
// Routine to write the current framebuffer to an image file
//
void writeImage(string outfilename, int flag){
  // make a pixmap that is the size of the window and grab OpenGL framebuffer into it
  // alternatively, you can read the pixmap into a 1d array and export this 
   //unsigned char local_pixmap[WinWidth * WinHeight * ImChannels];
   //glReadPixels(0, 0, WinWidth, WinHeight, pixformat, GL_UNSIGNED_BYTE, local_pixmap);
  
  // create the oiio file handler for the image
  std::unique_ptr<ImageOutput> outfile = ImageOutput::create(outfilename);
  if(!outfile){
    cerr << "Could not create output image for " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // Open a file for writing the image. The file header will indicate an image of
  // width WinWidth, height WinHeight, and ImChannels channels per pixel.
  // All channels will be of type unsigned char
  ImageSpec spec(WinWidth, WinHeight, ImChannels, TypeDesc::UINT8);
  if(!outfile->open(outfilename, spec)){
    cerr << "Could not open " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  Pixel **temp = new Pixel*[newImHeight];
  if(temp != NULL){
    temp[0] = new Pixel[newImWidth * newImHeight];
  }
  for(int i = 1; i < newImHeight; i++){
    temp[i] = temp[i - 1] + newImWidth;
  }
  //making a temp pixmap to flip and write since openimageio is weird
  if(flag == 0) {
    copyPixels(warpedPixmap,temp);
  } else {
    copyPixels(separated,temp);
  }
  for(int i = 0; i < newImHeight / 2; i++) {
    for(int j = 0; j < newImWidth; j++) {
        swap(temp[i][j], temp[newImHeight - 1 - i][j]);
    }
  }




  
  // Write the image to the file. All channel values in the pixmap are taken to be
  // unsigned chars. While writing, flip the image upside down by using negative y stride, 
  // since OpenGL pixmaps have the bottom scanline first, and oiio writes the top scanline first in the image file.
  //int scanlinesize = WinWidth * ImChannels * sizeof(unsigned char);
  //if(!outfile->write_image(TypeDesc::UINT8, warpedPixmap[0] + (WinHeight - 1) * scanlinesize, AutoStride, -scanlinesize)){
      if(!outfile->write_image(TypeDesc::UINT8, temp[0])) {  
    cerr << "Could not write image to " << outfilename << ", error = " << geterror() << endl;
    return;
  }
  
  // close the image file after the image is written and free up space for the
  // ooio file handler
  outfile->close();
  destroy(temp);
}
 
//handles the display for the repaired window
void handleNewDisplay() {
  glClearColor(0,0,0,1);
  glClear(GL_COLOR_BUFFER_BIT);

  if(newImWidth > 0 && newImHeight > 0) {
    if(WinWidth < newImWidth  || WinHeight < newImHeight)
      glPixelZoom(float(VpWidth) / newImWidth, float(VpHeight) / newImHeight);
    else
      glPixelZoom(1.0, 1.0);
    
    // display starting at the lower lefthand corner of the viewport
    glRasterPos2i(0, 0);
  
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glDrawPixels(newImWidth, newImHeight, pixformat, GL_UNSIGNED_BYTE, separated[0]);
  }

  glFlush();
} 
 
//
//   Display Callback Routine: clear the screen and draw the current image
//
void handleDisplay(){
  
  // specify window clear (background) color to be opaque black
  glClearColor(0.7, 0, 0, 1);
  // clear window to background color
  glClear(GL_COLOR_BUFFER_BIT);  
  
  // only draw the image if it is of a valid size
  if(newImWidth > 0 && newImHeight > 0)
    displayImage();
  
  // flush the OpenGL pipeline to the viewport
  glFlush();
}
 
//determines the size of the warped pixmap and allocates memory
//(steps 1 and 2)
void allocateBounds(Matrix3D &M, int ImHeight, int ImWidth, int scaleFlag) {
  double minX, minY, maxX, maxY;
  //botleft
  Vector3D bottomLeft(0,0,1);
  Vector3D newbotLeft = M * bottomLeft;
 
  //normalize
  newbotLeft.x /= newbotLeft.z;
  newbotLeft.y /= newbotLeft.z;
 
  //set mins and maxes to first values
  minX = maxX = newbotLeft.x;
  minY = maxY = newbotLeft.y;
 
 
  //topleft
  Vector3D topLeft(0, ImHeight, 1);
  Vector3D newtopLeft = M * topLeft;
 
  //normalize
  newtopLeft.x /= newtopLeft.z;
  newtopLeft.y /= newtopLeft.z;
 
  if(newtopLeft.x < minX) {
    minX = newtopLeft.x;
  }if(newtopLeft.x > maxX) {
    maxX = newtopLeft.x;
  }if(newtopLeft.y < minY) {
    minY = newtopLeft.y;
  }if(newtopLeft.y > maxY) {
    maxY = newtopLeft.y;
  }
  //botright
  Vector3D bottomRight(ImWidth, 0 , 1);
  Vector3D newbotRight = M * bottomRight;
 
  //normalize
  newbotRight.x /= newbotRight.z;
  newbotRight.y /= newbotRight.z;
 
  if(newbotRight.x < minX) {
    minX = newbotRight.x;
  }if(newbotRight.x > maxX) {
    maxX = newbotRight.x;
  }if(newbotRight.y < minY) {
    minY = newbotRight.y;
  }if(newbotRight.y > maxY) {
    maxY = newbotRight.y;
  }
  //topright
  Vector3D topRight(ImWidth, ImHeight, 1);
  Vector3D newtopRight = M * topRight;
 
  //normalize
  newtopRight.x /= newtopRight.z;
  newtopRight.y /= newtopRight.z;
 
  if(newtopRight.x < minX) {
    minX = newtopRight.x;
  }if(newtopRight.x > maxX) {
    maxX = newtopRight.x;
  }if(newtopRight.y < minY) {
    minY = newtopRight.y;
  }if(newtopRight.y > maxY) {
    maxY = newtopRight.y;
  }
  //make extrema ints
  minX = ceil(minX);
  minY = ceil(minY);
  maxX = ceil(maxX);
  maxY = ceil(maxY);
 
  if(scaleFlag == 0) {
    //determine how much to move the origin
    shiftedX = minX;
    shiftedY = minY;
  
    //determine new image dimensions
    newImWidth = maxX - minX;
    newImHeight = maxY - minY;
  
    //allocated memory for new image dimensions
    warpedPixmap = new Pixel*[newImHeight];
    if(warpedPixmap != NULL) {
      warpedPixmap[0] = new Pixel[newImWidth * newImHeight];
    }
  
    for(int i = 1; i < newImHeight; i++){
      warpedPixmap[i] = warpedPixmap[i - 1] + newImWidth;
    }
  } else if(scaleFlag == 1){
    //determine how much to move the origin
    scaleShiftedX = minX;
    scaleShiftedY = minY;

    //determine new image dimensions
    /*
    scaledWidth = maxX - minX;
    scaledHeight = maxY - minY;
    */
    scaledWidth = ::ImWidth;
    scaledHeight = ::ImHeight;
    cout << "We are now in the Allocate Bounds Method\n";
    cout << "The maxX is " << maxX << " while the minX is " << minX << endl;
    cout << "The maxY is " << maxY << " while the minY is " << minY << endl;

    //allocated memory for new image dimensions
    scaledHidden = new Pixel*[scaledHeight];
    if(scaledHidden != NULL) {
      scaledHidden[0] = new Pixel[scaledWidth * scaledHeight];
    }
  
    for(int i = 1; i < scaledHeight; i++){
      scaledHidden[i] = scaledHidden[i - 1] + scaledWidth;
    }
  }
}

void applyMatrix(Matrix3D &IM, Pixel** source, Pixel** warped, int shiftedX, int shiftedY, int sourceWidth, int sourceHeight, int newWidth, int newHeight) {
  for(int y = 0; y < newHeight; y++) {
    for(int x = 0; x < newWidth; x++) {
      Vector3D outpixel(x + shiftedX, y + shiftedY, 1);
      Vector3D inpixel = IM * outpixel;
 
      //normalize
      double u = inpixel.x / inpixel.z;
      double v = inpixel.y / inpixel.z;
      
     int V = (int)round(v);
     int U = (int)round(u);

     //prevent invalid input pixel locations
      if((V >= 0 && U >= 0 && V < sourceHeight && U < sourceWidth )) {
        //inverse map
        warped[y][x] = source[((int)V)][((int)U)];
        //cout << "Pixel at " << y << "," << x << " comes from input " << (int)round(V) << "," << (int)round(U) << endl;
      }
      
    }
  }
}

/*
Multiply M by a rotation matrix of angle theta
*/
 
void Rotate(Matrix3D &M, float theta) {
 
	Matrix3D R;  // this initializes R to identity  
	double rad = PI * theta / 180.0; // convert degrees to radians
 
 
	R[0][0] = cos(rad);
	R[0][1] = -sin(rad);
	R[1][0] = sin(rad);
	R[1][1] = cos(rad); 
 
	M = R * M; //append the rotation to your transformation matrix
 
}
 
void Scale(Matrix3D &M, float scaleX, float scaleY) {
  Matrix3D R;
 
  R[0][0] = scaleX;
  R[1][1] = scaleY;
 
  M = R * M;
}
 
void Translate(Matrix3D &M, int moveX, int moveY) {
  Matrix3D R;
 
  R[0][2] = moveX;
  R[1][2] = moveY;
 
  M = R * M;
}
 
void Shear(Matrix3D &M, float shearX, float shearY) {
  Matrix3D R;
 
  R[0][1] = shearX;
  R[1][0] = shearY;
 
  M = R * M;
}
 
void Flip(Matrix3D &M, int flipX, int flipY) {
  Matrix3D R;
 
  if(flipY == 1) {
    R[0][0] = -1;
  }
  if(flipX == 1) {
    R[1][1] = -1;
  }
 
  M = R * M;
}
 
void Perspective(Matrix3D &M, float a31, float a32) {
  Matrix3D R;

  R[2][0] = a31;
  R[2][1] = a32;
 
  M = R * M;
}
 
/*
Build a transformation matrix from input text
*/
void read_input(Matrix3D &M) {
	string cmd;
	
	/* prompt for user input */
	do
	{
    cout << "Would you like to perform warps on the merged image?\n";
    cout << "Enter a command: r , s , t , h, f, p, d\n";
		cout << "> ";
		cin >> cmd;
		if (cmd.length() != 1){
			cout << "invalid command, enter r, s, t, h, f, p, d\n";
		}
		else {
			switch (cmd[0]) {
				case 'r':		/* Rotation, accept angle in degrees */
					float theta;
          cout << "Enter a theta in degrees (counterclockwise):";
					cin >> theta;
					if (cin) {
						cout << "calling rotate\n";
						Rotate(M, theta);
					}
					else {
						cerr << "invalid rotation angle\n";
						cin.clear();
					}						
					break;
				case 's':		/* scale, accept scale factors */
          float scaleX, scaleY;
          cout << "Enter X axis scale factor: ";
          cin >> scaleX;
          cout << "Enter Y axis scale factor: ";
          cin >> scaleY;
          if(scaleX == 0 || scaleY == 0) {
            cout << "You tried to enter a 0 scale factor :(";
          } else {
            cout << "calling scale\n";
            Scale(M, scaleX, scaleY);
          }
					break;
				case 't':		/* Translation, accept translations */
          int moveX, moveY;
          cout << "Enter the X axis translate amount: ";
          cin >> moveX;
          cout << "Enter the Y axis translate amount: ";
          cin >> moveY;
          cout << "calling translate\n";
          Translate(M, moveX, moveY);
					break;
				case 'h':		/* Shear, accept shear factors */
          float shearX, shearY;
          cout << "Enter the X axis shear factor: ";
          cin >> shearX;
          cout << "Enter the Y axis shear factor: ";
          cin >> shearY;
          cout << "calling shear\n";
          Shear(M, shearX, shearY);
					break;
				case 'f':		/* Flip, accept flip factors */
          int flipX, flipY;
          cout << "Flip the x axis? 1 for yes: ";
          cin >> flipX;
          cout <<"Flip the y axis? 1 for yes: ";
          cin >> flipY;
          cout << "calling flip\n";
          Flip(M, flipX, flipY);
					break;
				case 'p':		/* Perspective, accept perspective factors */
          float a31, a32;
          cout << "Enter the first factor: ";
          cin >> a31;
          cout << "Enter the second factor: ";
          cin >> a32;
          cout << "calling perspective\n";
          Perspective(M, a31, a32);
					break;		
				case 'd':		/* Done, that's all for now */
					break;
				default:
					cout << "invalid command, enter r, s, t, h, f, p, d\n";
			}
		}
	} while (cmd.compare("d")!=0);
 
}
 

void handleKey(unsigned char key, int x, int y){
  
  switch(key){
    //cout << "Would you like to extract the image?";
    //cin >> 
	case 'q':		// q or ESC - quit
    case 'Q':
    case 27:
      destroy(pixmap);
      exit(0);
      
    default:		// not a valid key -- just ignore it
      return;
  }
}
 
 
//
//  Reshape Callback Routine: If the window is too small to fit the image,
//  make a viewport of the maximum size that maintains the image proportions.
//  Otherwise, size the viewport to match the image size. In either case, the
//  viewport is centered in the window.
//
void handleReshape(int w, int h){
  float imageaspect = (float)ImWidth / (float)ImHeight;	// aspect ratio of image
  float newaspect = (float)w / (float)h; // new aspect ratio of window
  
  // record the new window size in global variables for easy access
  WinWidth = w;
  WinHeight = h;
  
  // if the image fits in the window, viewport is the same size as the image
  if(w >= newImWidth && h >= newImHeight){
    Xoffset = (w - newImWidth) / 2;
    Yoffset = (h - newImHeight) / 2;
    VpWidth = newImWidth;
    VpHeight = newImHeight;
  }
  // if the window is wider than the image, use the full window height
  // and size the width to match the image aspect ratio
  else if(newaspect > imageaspect){
    VpHeight = h;
    VpWidth = int(imageaspect * VpHeight);
    Xoffset = int((w - VpWidth) / 2);
    Yoffset = 0;
  }
  // if the window is narrower than the image, use the full window width
  // and size the height to match the image aspect ratio
  else{
    VpWidth = w;
    VpHeight = int(VpWidth / imageaspect);
    Yoffset = int((h - VpHeight) / 2);
    Xoffset = 0;
  }
  
  // center the viewport in the window
  glViewport(Xoffset, Yoffset, VpWidth, VpHeight);
  
  // viewport coordinates are simply pixel coordinates
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(0, VpWidth, 0, VpHeight);
  glMatrixMode(GL_MODELVIEW);
}
 
 void swap(Pixel *a, Pixel *b)
{
    Pixel temp = *a;
    *a = *b;
    *b = temp;
}

//this function clears int x number of least significant bits from pix
void clearBits(Pixel &pix, int x) {
  pix.r >>= x;
  pix.r <<= x;
  pix.g >>= x;
  pix.g <<= x;
  pix.b >>= x;
  pix.b <<= x;
  if(ImChannels == 4) {
    pix.a >>= x;
    pix.a <<= x;
  }
}

//this function gets the x significant bits for the hidden image
void hiddenBits(Pixel &pix, int x) {
  pix.r >>= (8 - x);
  pix.g >>= (8 - x);
  pix.b >>= (8 - x);
  if(ImChannels == 4) {
    pix.a >>= (8 - x);
  }
}

//this function loads most significant bits into the least significant position
//of the cleared large pixel
void mergeBits(Pixel &big, Pixel &hidden) {
  big.r = big.r | hidden.r;
  big.g = big.g | hidden.g;
  big.b = big.b | hidden.b;
  if(ImChannels == 4) {
    big.a = big.a | hidden.a;
  }
}

//loads the least significant positions into a vector to be 
//unloaded as significant bits in a pixmap for the extracted image
void clearSigBits(Pixel &pix, int x, vector<unsigned char> &rgb) {

  //load vector with the cleared bits
  rgb.push_back(pix.r << (8 - x));
  rgb.push_back(pix.g << (8 - x));
  rgb.push_back(pix.b << (8 - x));
  if(ImChannels == 4) {
    rgb.push_back(pix.a << (8 - x));
  } 
}

//This function goes through both images and merges them into the large
void mergeImages(Pixel** big, Pixel** hidden, int x) {

  //Testing check: checks to see if the width and height are
  //now the same after scaling
  if(ImHeight == scaledHeight && ImWidth == scaledWidth) {
    cout << "Dimensions now match :)" << endl;
  } else {
    cout << "The two image dimensions do not match\n" <<endl;
    cout << "Big Width is : " << ImWidth << endl;
    cout << "Scaled width is : " << scaledWidth << endl;
    cout << "Big Height is : " << ImHeight << endl;
    cout << "Scaled Height is : " << scaledHeight << endl << endl;
  }


  //loop through bigger image to clear
  for(int r = 0; r < ImHeight; r++) {
    for(int c = 0; c < ImWidth; c++) {
      clearBits(big[r][c], x);
    }
  }
  //loop through hidden image to clear
  for(int r = 0; r < ImHeight; r++) {
    for(int c = 0; c < ImWidth; c++) {
      hiddenBits(hidden[r][c], x);
    }
  }
  //merge the bits
  for(int r = 0; r < ImHeight; r++) {
    for(int c = 0; c < ImWidth; c++) {
      mergeBits(big[r][c], hidden[r][c]);
    }
  }
}

void bilinearInterpolation(const float &u, const float &v,
  unsigned char &red, unsigned char &green, unsigned char &blue) {
  float j = floor(u);
  float k = floor(v);

  float a = u - j;
  float b = v - k;

  Matrix3D bVect;
  bVect[0][0] = 1-b;
  bVect[0][1] = b;
  Vector3D aVect(1-a, a, 1);
  Matrix3D midR, midG, midB;
  if(k + 1 < ImHeight2 && j + 1 < ImHeight2) {
  //red
  midR[0][0] = pixmap2[(int)k][(int)j].r;
  midR[0][1] = pixmap2[(int)k][(int)j + 1].r;
  midR[1][0] = pixmap2[(int)k + 1][(int)j].r;
  midR[1][1] = pixmap2[(int)k + 1][(int)j + 1].r;
  //green
  midG[0][0] = pixmap2[(int)k][(int)j].g;
  midG[0][1] = pixmap2[(int)k][(int)j + 1].g;
  midG[1][0] = pixmap2[(int)k + 1][(int)j].g;
  midG[1][1] = pixmap2[(int)k + 1][(int)j + 1].g;
  //blue
  midB[0][0] = pixmap2[(int)k][(int)j].b;
  midB[0][1] = pixmap2[(int)k][(int)j + 1].b;
  midB[1][0] = pixmap2[(int)k + 1][(int)j].b;
  midB[1][1] = pixmap2[(int)k + 1][(int)j + 1].b;

  //calculate
  Matrix3D temp = bVect * midR;
  red = (temp * aVect).x;

  temp = bVect * midG;
  green = (temp * aVect).x;

  temp = bVect * midB;
  blue = (temp * aVect).x;
  }
}

void applyBI() {
  for(int y = 0; y < ImHeight; y++) {
    for(int x = 0; x < ImHeight; x++) {
      float u = (float)x;
      float v = (float)y;
      unsigned char r, g, b;
      bilinearInterpolation(u, v, r, g, b);
      scaledHidden[y][x].r = r;
      scaledHidden[y][x].g = g;
      scaledHidden[y][x].b = b;
    }
  }
}

//this function scales the size of the hidden image so that
//it can fit inside of the larger image
void scaleToFit() {
  Matrix3D S;

  float scaleX = ((float)ImWidth) / ImWidth2;
  float scaleY = ((float)ImHeight) / ImHeight2;
  cout << "\nThe hidden width is " << ImWidth2 <<endl;
  cout << "The hidden height is " << ImHeight2 << endl;
  cout << "\nThe scalex is " << scaleX <<endl;
  cout <<"The scaley is " << scaleY <<endl <<endl;
  Scale(S, scaleX, scaleY);

  allocateBounds(S, ImHeight2, ImWidth2, 1);
  Matrix3D inverseScaled = S.inverse();
  applyMatrix(inverseScaled, pixmap2, scaledHidden, scaleShiftedX, scaleShiftedY, ImWidth2, ImHeight2, scaledWidth, scaledHeight);
  //applyBI();
}

void extractImage(int bits) {
  vector<unsigned char> rgb;

  //allocate memory for the separated image
  separated = new Pixel*[newImHeight];
  if(separated != NULL) {
    separated[0] = new Pixel[newImWidth * newImHeight];
  }

  for(int i = 1; i < newImHeight; i++){
    separated[i] = separated[i - 1] + newImWidth;
  }

  //loop image
  for(int r = 0; r < newImHeight; r++) {
    for(int c = 0; c < newImWidth; c++) {
      //extract lsb from image
      clearSigBits(warpedPixmap[r][c], bits, rgb);
      //apply lsb as msb in separated pixmap for r,g,b,/a
      separated[r][c].r = rgb[0];
      rgb.erase(rgb.begin());
      separated[r][c].g = rgb[0];
      rgb.erase(rgb.begin());
      separated[r][c].b = rgb[0];
      rgb.erase(rgb.begin());
      if(ImChannels == 4) {
        separated[r][c].a = rgb[0];
        rgb.erase(rgb.begin());
      }
      //just incase vector isnt empty for some reason
      rgb.clear();
    }
  }
}

//
// Main program to scan the commandline, set up GLUT and OpenGL, and start Main Loop
//
int main(int argc, char* argv[]){
  // scan command line and process
  // only one parameter allowed, an optional image filename and extension
  if(argc > 3 || argc < 2){
    cout << "usage: ./final [BiggerImage.ext] [HiddenImage.ext]" << endl;
    cout << "extraction usage: ./final [MergedImage.ext]" <<endl;
    exit(1);
  }
 
  // initialize transformation matrix to identity
	Matrix3D M;
  // set up the default window and empty pixmap if no image or image fails to load
  WinWidth = DEFAULTWIDTH;
  WinHeight = DEFAULTHEIGHT;
  ImWidth = 0;
  ImHeight = 0;
  
  // load the image if present, and size the window to match
  if(argc == 3){
    if(readImage(argv[1])){
      // :)
    } else {
      cout << "Failure reading image 1" <<endl;
      exit(1);
    }
    if(readImage(argv[2])) {
      // :)
    } else {
      cout << "Failure reading image 2" <<endl;
      exit(1);
    }
    //scale the hidden image to fit inside the larger image
    scaleToFit();
    //merge the scaled hidden image into the larger image
    int bits = 4;
    cout << "How many bits would you like to hide in the image?";
    cin >> bits;
    mergeImages(pixmap, scaledHidden, bits);

    //build the transformation matrix based on user input
    read_input(M);
 
    // your code to perform inverse mapping (4 steps)
    allocateBounds(M, ImHeight, ImWidth, 0);
    WinWidth = newImWidth;
    WinHeight = newImHeight;
    Matrix3D inversed = M.inverse();

    applyMatrix(inversed, pixmap, warpedPixmap, shiftedX, shiftedY, ImWidth, ImHeight, newImWidth, newImHeight);
    writeImage("MergedImage.png", 0);
    extractImage(bits);
    writeImage("ExtractedImage1.png", 1);

    // start up GLUT
    glutInit(&argc, argv);
    
    // create the graphics window, giving width, height, and title text
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(WinWidth, WinHeight);
    glutCreateWindow("Merged Image");
    
    // set up the callback routines
    glutDisplayFunc(handleDisplay); // display update callback
    glutKeyboardFunc(handleKey);	  // keyboard key press callback
    glutReshapeFunc(handleReshape); // window resize callback
    

    //routines for extracted image
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(WinWidth, WinHeight);
    glutCreateWindow("Extracted Image");

    glutDisplayFunc(handleNewDisplay);
    glutKeyboardFunc(handleKey);
    glutReshapeFunc(handleReshape);
    // Enter GLUT's event loop
    glutMainLoop();

  } else if(argc == 2) {
    if(readImage(argv[1])) {
      int bits = 4;
      cout << "How many bits were hidden in the merged image?";
      cin >> bits;
      newImWidth = ImWidth;
      newImHeight = ImHeight;
      WinWidth = ImWidth;
      WinHeight = ImHeight;
      warpedPixmap = new Pixel*[newImHeight];
      if(warpedPixmap != NULL) {
        warpedPixmap[0] = new Pixel[newImWidth * newImHeight];
      }
    
      for(int i = 1; i < newImHeight; i++){
        warpedPixmap[i] = warpedPixmap[i - 1] + newImWidth;
      }
      copyPixels(pixmap, warpedPixmap);

      extractImage(bits);
      writeImage("ExtractedImage.png", 1);
      // start up GLUT
      glutInit(&argc, argv);
      
      // create the graphics window, giving width, height, and title text
      glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
      glutInitWindowSize(WinWidth, WinHeight);
      glutCreateWindow("Extracted Image");
      
      // set up the callback routines
      glutDisplayFunc(handleNewDisplay); // display update callback
      glutKeyboardFunc(handleKey);	  // keyboard key press callback
      glutReshapeFunc(handleReshape); // window resize callback

      // Enter GLUT's event loop
      glutMainLoop();
    }

  }
  
  
  return 0;
}
 
