#include "LegendreFET.h"
#include "MonteCarloSimulation.h"
#include <cmath>
#include <gtest/gtest.h>

using namespace nonte_fonte;

#define N_SAMPLES 1000000

TEST(LegendreFETTest, Construction1D) {
  Domain domain(-1.0, 1.0);
  LegendreFET tally(domain, {5});
  EXPECT_EQ(tally.orders()[0], 5);
  EXPECT_EQ(tally.coefficients().size(), 6);
}

TEST(LegendreFETTest, Construction2D) {
  Domain domain(-1.0, 1.0, -1.0, 1.0);
  LegendreFET tally(domain, {3, 4});
  EXPECT_EQ(tally.orders()[0], 3);
  EXPECT_EQ(tally.orders()[1], 4);
  EXPECT_EQ(tally.coefficients().size(), 4 * 5);
}

TEST(LegendreFETTest, InvalidConstruction) {
  Domain domain(-1.0, 1.0);
  EXPECT_THROW(LegendreFET(domain, {5, 5}), std::invalid_argument);
  EXPECT_THROW(LegendreFET(domain, {-1}), std::invalid_argument);
}

TEST(LegendreFETTest, UniformPDF1D) {
  Domain domain(-1.0, 1.0);
  LegendreFET tally(domain, {5});
  auto pdf = [](const std::vector<double> &) { return 0.5; };
  MCSimulation mc(domain, pdf, N_SAMPLES, 42);
  mc.addTally(&tally);
  mc.run();
  const auto &coeffs = tally.coefficients();
  EXPECT_NEAR(coeffs[0], 0.5, 0.01);
  for (size_t i = 1; i < coeffs.size(); ++i)
    EXPECT_NEAR(coeffs[i], 0.0, 0.01);
}

TEST(LegendreFETTest, LinearPDF) {
  Domain domain(-1.0, 1.0);
  LegendreFET tally(domain, {2});
  auto pdf = [](const std::vector<double> &x) { return (x[0] + 1.0) / 2.0; };
  MCSimulation mc(domain, pdf, N_SAMPLES, 42);
  mc.addTally(&tally);
  mc.run();
  const auto &coeffs = tally.coefficients();
  EXPECT_NEAR(coeffs[0], 0.5, 0.02);

}


TEST(LegendreFETTest, DomainMapping) {
  Domain domain(0.0, 2.0);
  LegendreFET tally(domain, {5});
  auto pdf = [](const std::vector<double> &) { return 0.5; };
  MCSimulation mc(domain, pdf, N_SAMPLES, 42);
  mc.addTally(&tally);
  mc.run();
  EXPECT_NEAR(tally.reconstruct({0.5}), 0.5, 0.01);
  EXPECT_NEAR(tally.reconstruct({1.0}), 0.5, 0.01);
  EXPECT_NEAR(tally.reconstruct({1.5}), 0.5, 0.01);
}

TEST(LegendreFETTest, Normalization) {
  Domain domain(-1.0, 1.0);
  LegendreFET tally(domain, {5});
  auto pdf = [](const std::vector<double> &x) { return 0.5; };
  MCSimulation mc(domain, pdf, N_SAMPLES, 42);
  mc.addTally(&tally);
  mc.run();
  int n_samples = 1000;
  double dx = 2.0 / n_samples;
  double integral = 0.0;
  for (int i = 0; i < n_samples; ++i) {
    double x = -1.0 + i * dx;
    integral += tally.reconstruct({x}) * dx;
  }
  EXPECT_NEAR(integral, 1.0, 0.01);
}