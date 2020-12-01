// Copyright 2019 dimakirol <your_email>

#include "header.hpp"
#include <gtest/gtest.h>

TEST(Example, EmptyTest)
{
  std::ifstream fin("matrixEZ");
    if (!fin)
    {
        std::cout << "Could't read" << std::endl;
        return;
    }

    SqMatrix A(CURRENT_SIZE);
    SqMatrix B(CURRENT_SIZE);
    SqMatrix C(CURRENT_SIZE);
    A.fill(fin);
    B.fill(fin);
    std::cout << std::endl;
    std::cout << std::endl;
    SqMatrix D(CURRENT_SIZE);
    clock_t t_s = clock();
    SqMatrix::simple_multiply(D, A, B, 0, 0, 0, 0, CURRENT_SIZE);
    t_s = clock() - t_s;
    D.print();
    clock_t t_b = clock();
    SqMatrix::block_multiply(C, A, B, 0, 0, 0, 0, CURRENT_SIZE);
    t_b = clock() - t_b;
    std::cout << std::endl;
    C.print();
    Check(C, D);
    std::cout << "Time taken by simple = " << t_s << "    Time taken by block = " << t_b << std::endl;
    std::cout << "k = " << log(k)/log(7) << std::endl;
    std::cout << "Hello, World!" << std::endl;
  EXPECT_TRUE(true);
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
