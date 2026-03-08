/**
 * @file test_domain.cpp
 * @brief Unit tests for Domain class
 */

#include "Domain.h"
#include <gtest/gtest.h>
#include <random>

using namespace nonte_fonte;

// Test 1D domain construction
TEST(DomainTest, Construction1D) {
  Domain domain(0.0, 1.0);

  EXPECT_EQ(domain.dim(), 1);
  EXPECT_DOUBLE_EQ(domain.measure(), 1.0);

  auto [min, max] = domain.bounds(0);
  EXPECT_DOUBLE_EQ(min, 0.0);
  EXPECT_DOUBLE_EQ(max, 1.0);
}

// Test 2D domain construction
TEST(DomainTest, Construction2D) {
  Domain domain(-1.0, 1.0, -2.0, 2.0);

  EXPECT_EQ(domain.dim(), 2);
  EXPECT_DOUBLE_EQ(domain.measure(), 8.0); // 2 * 4

  auto [x_min, x_max] = domain.bounds(0);
  EXPECT_DOUBLE_EQ(x_min, -1.0);
  EXPECT_DOUBLE_EQ(x_max, 1.0);

  auto [y_min, y_max] = domain.bounds(1);
  EXPECT_DOUBLE_EQ(y_min, -2.0);
  EXPECT_DOUBLE_EQ(y_max, 2.0);
}

// Test 3D domain construction
TEST(DomainTest, Construction3D) {
  Domain domain(0.0, 1.0, 0.0, 2.0, 0.0, 3.0);

  EXPECT_EQ(domain.dim(), 3);
  EXPECT_DOUBLE_EQ(domain.measure(), 6.0); // 1 * 2 * 3
}

// Test invalid domain (min >= max)
TEST(DomainTest, InvalidDomain1D) {
  EXPECT_THROW(Domain(1.0, 0.0), std::invalid_argument);
  EXPECT_THROW(Domain(1.0, 1.0), std::invalid_argument);
}

TEST(DomainTest, InvalidDomain2D) {
  EXPECT_THROW(Domain(1.0, 0.0, 0.0, 1.0), std::invalid_argument);
  EXPECT_THROW(Domain(0.0, 1.0, 1.0, 0.0), std::invalid_argument);
}

// Test contains method
TEST(DomainTest, Contains1D) {
  Domain domain(-1.0, 1.0);

  EXPECT_TRUE(domain.contains({0.0}));
  EXPECT_TRUE(domain.contains({-1.0}));
  EXPECT_TRUE(domain.contains({1.0}));
  EXPECT_FALSE(domain.contains({-1.1}));
  EXPECT_FALSE(domain.contains({1.1}));
  EXPECT_FALSE(domain.contains({0.0, 0.0})); // Wrong dimension
}

TEST(DomainTest, Contains2D) {
  Domain domain(-1.0, 1.0, -2.0, 2.0);

  EXPECT_TRUE(domain.contains({0.0, 0.0}));
  EXPECT_TRUE(domain.contains({-1.0, -2.0}));
  EXPECT_TRUE(domain.contains({1.0, 2.0}));
  EXPECT_FALSE(domain.contains({1.1, 0.0}));
  EXPECT_FALSE(domain.contains({0.0, 2.1}));
  EXPECT_FALSE(domain.contains({0.0})); // Wrong dimension
}

// Test sampling
TEST(DomainTest, Sampling1D) {
  Domain domain(0.0, 1.0);
  std::mt19937 rng(42);

  for (int i = 0; i < 1000; ++i) {
    auto point = domain.sample(rng);
    ASSERT_EQ(point.size(), 1);
    EXPECT_GE(point[0], 0.0);
    EXPECT_LE(point[0], 1.0);
  }
}

TEST(DomainTest, Sampling2D) {
  Domain domain(-1.0, 1.0, -2.0, 2.0);
  std::mt19937 rng(42);

  for (int i = 0; i < 1000; ++i) {
    auto point = domain.sample(rng);
    ASSERT_EQ(point.size(), 2);
    EXPECT_GE(point[0], -1.0);
    EXPECT_LE(point[0], 1.0);
    EXPECT_GE(point[1], -2.0);
    EXPECT_LE(point[1], 2.0);
  }
}

// Test sampling distribution uniformity
TEST(DomainTest, SamplingUniformity) {
  Domain domain(0.0, 1.0);
  std::mt19937 rng(42);

  int n_samples = 10000;
  int n_bins = 10;
  std::vector<int> bins(n_bins, 0);

  for (int i = 0; i < n_samples; ++i) {
    auto point = domain.sample(rng);
    int bin = static_cast<int>(point[0] * n_bins);
    if (bin >= 0 && bin < n_bins) {
      bins[bin]++;
    }
  }

  // Each bin should have ~1000 samples (±10%)
  int expected = n_samples / n_bins;
  for (int count : bins) {
    EXPECT_NEAR(count, expected, expected * 0.1);
  }
}