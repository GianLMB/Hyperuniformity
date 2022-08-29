//  The code has the purpose of computing the scaled variance of the number of particles inside circles of radii R in the interval [0.1, L/2], where L is the size of the system. \


#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <omp.h>
using namespace std;


// The structure matpt contains two vectors (x and y), which are respectively the x and y positions of the "particles" in the hyperuniform system.\
Each vector has thus dimension (L^2)x1
struct matpt {
    vector<float> x;
    vector<float> y;
};


// The function initializes the particles' position vectors x and y. If the displacement delta is zero a square grid is obtained (with lattice parameter 1),\
while if delta is different from zero the square grid is randomly shuffled, displacing each particle by a uniform random number in the range [-delta/2,delta/2] in both x\
and y directions. In order to take into account periodc boundary conditions, we take particles back when they exit from the lattice.\
Inputs: \
    - grid : sctructure cointaining the x and y positions' vectors. \
    - L : lattice size \
    - delta : particle random displacement \
Output: \
    - grid : is the sctructure cointaining the x and y positions' vectors.

void grid_make(matpt* grid, int L, float delta) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {  // j is the y index (position).
            float xpos = float(i); // i is the x index (position).
            float ypos = float(j); // j is the y index (position).
            if (delta != 0) {
                xpos += ((float) rand()/RAND_MAX-0.5)*delta; // when delta is different from zero each particle is randomly displaced
                if (xpos < 0 || xpos >= L) {xpos = fmod(xpos, L);}  // the modulo operation is needed to tak into account periodic boundary conditions. 
                ypos += ((float) rand()/RAND_MAX-0.5)*delta;
                if (ypos < 0 || xpos >= L) {ypos = fmod(ypos, L);}
            }
            grid->x.push_back(xpos);  // each particle position is inserted in the matpt structure
            grid->y.push_back(ypos);
        }
    }
}

// The function computes the distance between two points P1 = (x1, y1) and P2 = (x2, y2), taking into account periodic boundarcy conditions on a square lattice of size L. \
Inputs: \
    - x1, y1, x2, y2 : x and y positions of the points P1 and P2. \
    - L : lattice size. \
Output: \
    - d : distance between P1 and P2.

float distance(float x1, float y1, float x2, float y2, int L) {
        float dx = abs(x2-x1); // the usual distance along x and y is computed as |x1-x2| and |y1-y2|
        float dy = abs(y2-y1);
        if (dx > (float) L/2) {dx = L - fmod(dx,L);} // due to periodic boundary conditions, whenever dx (or dy) are larger than L/2 the distance is taken as L-dx (or L-dy).\
        The modulo operation is needed to take into account periodic boundary conditions.
        if (dy > (float) L/2) {dy = L - fmod(dy,L);}
        float d = sqrt(pow(dx,2)+pow(dy,2)); // the distance between P1 and P2 is, using pitagora's theorem, d=sqrt(dx^2+dy^2)
    return d;
}

// The function computes the number of particles inside a circle of radius R centered in C = (x0, y0). A particle in position P = (x, y) is inside the circle if its distance d \
from C is less or equal to R. Then cycling over all the particles and computing the distance from C we obtain the number of particles inside the circle. \
Inputs: \
    - x0, y0 : x and y positions of the center C. \
    - grid : sctructure cointaining the x and y positions' vectors. \
    - R : radius of the circle. \
    - L : lattice size. \
Output: \
    - count : number of particles inside the circle.

int number_points(float x0, float y0, matpt* grid, float R, float L) {
    int count = 0;
    for (int i = 0; i < pow(L,2); i++) {  // for every particle in the lattice we compute its distance d from C and we increment count by 1 if d<=R.
        if (distance(grid->x[i], grid->y[i], x0, y0, L) <= R) {count++;}
    }
    return count;
}

// The function compute the scaled variance scvar of the number of particles inside a circle of radius R, defined as the variance divided by R^2. In order to do it the function \
randomly extracts npoints centers C, and compute the number of particles inside the circles of radius R centered in them, computing then the variance of that sample. \
Inputs: \
    - grid : sctructure cointaining the x and y positions' vectors. \
    - R : radius of the circle. \
    - L : lattice size. \
    - npoints : number of centers extracted. \
Output: \
    - scvar : scaled variance.

float scaled_variance(matpt* grid, float R, int L, int npoints) {
    float x0;  // x0 and y0 are the x and y position of the centers C.
    float y0;
    int num = 0;  // num is the running sum of the number of points inside each circle, summed over all the random centers.
    int sqnum = 0;  // sqnum is the running sum of the square of the number of points inside each circle, summed over all the random centers.
    x0 =  (float) rand()/RAND_MAX*(L-1);  // a random center C is extracted
    y0 = (float) rand()/RAND_MAX*(L-1);
    int N0 = number_points(x0, y0, grid, R, L);  // in order to avoid roundoff and overflow errors, the value of the number of particles N for each extraction is shifted back \
    by the value of the first extraction N0. In this way all the values are small, and when summing them we avoid overflow errors.

    // the code is parallelized over the random extractions
    #pragma omp parallel for reduction(+:num) reduction(+:sqnum)
    for (int i = 1; i < npoints; i++) {
        x0 =  (float) rand()/RAND_MAX*(L-1);
        y0 = (float) rand()/RAND_MAX*(L-1);
        int N = number_points(x0, y0, grid, R, L);
        N -= N0;  // the value of N is shifted back by the value of the first extraction N0.
        num += N;  // num and sqnum are increased by summing the values obtained in each extraction.
        sqnum += pow(N,2);
    }
    float mean = (float) num/npoints;  // the average of the number of particle inside the circle of radius R is computed.
    float sqmean = (float) sqnum/npoints;  // the second moment is computed.
    float scvar = (sqmean - pow(mean,2)) / pow(R,2);  // the scaled variance is computed.
    return scvar;
}



int main() {
    srand(time(NULL));  // the seed of the pseudo random number generator is initialized.
    int L;
    float delta;
    cin >> L >> delta;  // L and delta are given from input.
    matpt grid;

    grid_make(&grid, L, delta);  // the grid structure is initialized.

    // the scaled variance is computed for each radius R in the interval [0.1, L/2], using a geometric sequence with common ratio r such that 200 values of R are obtained.
    float R = 0.1;
    float r = pow((float) L/(2*0.1), (float) 1/(200-1));  // the common ratio r is computed in order to have 200 values of R in the interval [0.1, L/2].
    cout << "L = " << L << "         delta = " << delta << endl;
    while (R <= (float) L/2) {
        cout << R << "    " << scaled_variance(&grid, R, L, 10000) << endl;  // the scaled variance is computed and outputted together with the corresponding value of R.
        R *= r;
    }

    return 0;
}
