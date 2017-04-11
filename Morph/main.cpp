#include "CImg.h"
#include "Morph.h"


using namespace cimg_library;

int main() {
    CImg<unsigned char> img1("Data/1.bmp"), img2("Data/2.bmp");
    Morph<unsigned char> m(img1, img2);
    //m.LoadPQ("Data/temp/");
    m.drawPQ();
    m.ShowMorph(11);
}