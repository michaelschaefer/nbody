#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <boost/numeric/ublas/vector.hpp>


using namespace std;


typedef unsigned int uint;
typedef boost::numeric::ublas::vector<double> vec;


void rhs(vec x, vec& res, vec& gm, double eps) {
  uint N = x.size() / 4;
  assert(x.size() == res.size() && gm.size() == N);

  double ax, ay, denominator, dx, dy, xi, yi;
  for (uint i = 0; i < N; ++i) {
    ax = 0;
    ay = 0;
    xi = x[2*i];
    yi = x[2*i+1];

    for (uint j = 0; j < N; ++j) {
      if (i == j)
	continue;
      dx = x[2*j] - xi;
      dy = x[2*j+1] - yi;
      denominator = pow(dx*dx + dy*dy + eps, 1.5);
      ax += gm[j] * dx / denominator;
      ay += gm[j] * dy / denominator;
    }

    res[2*i] = x[2*(N+i)];
    res[2*i+1] = x[2*(N+i)+1];
    res[2*(N+i)] = ax;
    res[2*(N+i)+1] = ay;
  }
}


void euler(vec& x, double dt, vec& gm, double eps) {
  vec f(x.size());
  rhs(x, f, gm, eps);
  x += dt * f;
}


void runge_kutta(vec& x, double dt, vec& gm, double eps) {
  uint N = x.size();
  vec f1(N), f2(N), f3(N), f4(N);

  rhs(x, f1, gm, eps);
  rhs(x + 0.5 * dt * f1, f2, gm, eps);
  rhs(x + 0.5 * dt * f2, f3, gm, eps);
  rhs(x + dt * f3, f4, gm, eps);
  
  x += dt / 6 * (f1 + 2*f2 + 2*f3 + f4);
}


void gravity_center(vec& x, vec& gm) {
  uint N = x.size() / 4;
  assert(gm.size() == N);

  double M = 0, sum_x = 0, sum_y = 0;
  for (uint i = 0; i < N; ++i) {
    M += gm[i];
    sum_x += (gm[i] * x[2*i]);
    sum_y += (gm[i] * x[2*i+1]);
  }

  cout << "center: " << (sum_x / M) << ", " << (sum_y / M) << endl;
}


void make_periodic(vec& x, double dx, double dy) {
  uint N = x.size() / 4;
  for (uint i = 0; i < N; ++i) {
    if (x[2*i] >= dx) x[2*i] -= 2*dx;
    if (x[2*i] < -dx) x[2*i] += 2*dx;
    if (x[2*i+1] >= dy) x[2*i+1] -= 2*dy;
    if (x[2*i+1] < -dy) x[2*i+1] += 2*dy;
  }
}


void write_to_gnuplot(vec& x, vec& gm, string filename, double dx, double dy) {
  uint N = x.size() / 4;
  assert(gm.size() == N);

  double gm_max = 0;
  ofstream data(filename + ".txt");
  for (uint i = 0; i < N; ++i) {
    if (gm[i] > gm_max) gm_max = gm[i];
    data << x[2*i] << " " << x[2*i+1] << " " << gm[i] << endl;
  }
  data.close();

  ofstream plot(filename + ".plot"); 
  plot << "unset key" << endl;
  /*plot << "unset colorbox" << endl;
    plot << "unset border" << endl;*/
  plot << "set palette defined (0 \"blue\", 1 \"red\")" << endl;
  plot << "set cbrange [0:" << gm_max << "]" << endl;
  plot << "set xrange [" << (-dx) << ":" << dx << "]" << endl;
  plot << "set yrange [" << (-dy) << ":" << dy << "]" << endl;
  plot << "plot '" << filename << ".txt' using 1:2:3 with points palette pt 7 ps 2" << endl;
  plot.close();
}


int main(int argc, char* argv[]) {
  uint N;
  if (argc == 1)
    N = 2;
  else
    N = atoi(argv[1]);

  double eps = 0.01;

  double dt = 0.001;
  double t = 0.0;
  double tf = 100.0;
  uint n_iter = 0;

  double dgm = 1, dx = 5, dy = 5;
  std::uniform_real_distribution<double> dist_x(-dx, dx);
  std::uniform_real_distribution<double> dist_y(-dy, dy);
  std::uniform_real_distribution<double> dist_gm(0, dgm);
  std::default_random_engine re;
  re.seed(123);
  
  vec gm(N);
  vec x(4*N);

  for (uint i = 0; i < N; ++i) {
    gm[i] = exp(dist_gm(re));
    x[2*i] = dist_x(re);
    x[2*i+1] = dist_y(re);
    x[2*(N+i)] = 0;
    x[2*(N+1)+1] = 0;
  }

  write_to_gnuplot(x, gm, "start", dx, dy);

  while (t < tf) {
    runge_kutta(x, dt, gm, eps);
    make_periodic(x, 5, 5);
    t = n_iter * dt;
    n_iter++;
    
    cout << "Progress: " << (int)(100*t/tf) << "%\r" << flush;
  }

  write_to_gnuplot(x, gm, "end", dx, dy);
  
  return 0;
}
