// Copyright 2019 ScissorsS <file>

#ifndef INCLUDE_HEADER_HPP_
#define INCLUDE_HEADER_HPP_

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#define CURRENT_SIZE 17
#define MINIMUM_SIZE 2

int k = 0;

class SqMatrix {
public:
    explicit SqMatrix(uint size_in)
    {
        size = size_in;
        //std::cout << "Taking memory" << std::endl;
        content = new double* [size];
        for (uint i = 0; i < size; ++i)
        {
            content[i] = new double[size];
            for (uint j = 0; j < size; ++j)
            {
                content[i][j] = 0;
            }
        }
    }

    SqMatrix(const SqMatrix& old_matrix,
             uint size_in, uint start_pos_line, uint start_pos_column)
    {
        size = size_in;
        content = new double*[size];
        for (uint i = 0; i < size; ++i)
        {
            content[i] = new double[size];
            for (uint j = 0; j < size; ++j)
            {
                content[i][j] = 0;
            }
        }

        for (uint i = 0; i < size; ++i)
        {
            for (uint j = 0; j < size; ++j)
            {
                content[i][j] = old_matrix.content[start_pos_line+i][start_pos_column+j];
            }
        }
    }

    SqMatrix (const SqMatrix& matrix)
    {
        size = matrix.size;
        content = new double*[size];
        for (uint i = 0; i < size; ++i)
        {
            content[i] = new double[size];
            for (uint j = 0; j < size; ++j)
            {
                content[i][j] = 0;
            }
        }
        for (uint i = 0; i < size; ++i)
        {
            for (uint j = 0; j < size; ++j)
            {
                content[i][j] = matrix.content[i][j];
            }
        }
    }

    void fill(std::ifstream &fin) const
    {
        for (uint i = 0; i < size; ++i)
        {
            for (uint j = 0; j < size; ++j)
            {
                fin >> content[i][j];
            }
        }
    }

    static void/*SqMatrix*/ simple_multiply(SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                            uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size) {
        //SqMatrix C(cur_size);
        //std::cout << "Multiply: " << std::endl;
        for (uint i = 0; i < cur_size; ++i)
        {
            for (uint j = 0; j < cur_size; ++j)
            {
                uint b = 0;
                for (uint a = 0; (a < cur_size) && (b < cur_size); ++a, ++b)
                {
                    //std::cout << "A: " << A(AposL + i, AposC + a) << "  B: " << B(BposL + b, BposC + j) << std::endl;
                    C(i, j) += A(AposL + i, AposC + a) * B(BposL + b, BposC + j);
                }
            }
        }
        /*return C;*/
    }

    static SqMatrix PNadd(const SqMatrix& A, const SqMatrix& B,
                          uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        SqMatrix C(cur_size);
        //std::cout << "PNAdd: " << std::endl;
        for (uint i = 0; i < cur_size; ++i)
        {
            for (uint j = 0; j < cur_size; ++j)
            {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }

    static SqMatrix PNsubtract(const SqMatrix& A, const SqMatrix& B,
                               uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        SqMatrix C(cur_size);
        for (uint i = 0; i < cur_size; ++i)
        {
            for (uint j = 0; j < cur_size; ++j)
            {
                C(i, j) = A(i, j) - B(i, j);
            }
        }
        return C;
    }

    static SqMatrix add(const SqMatrix& A, const SqMatrix& B,
                        uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        SqMatrix C(cur_size);
        //std::cout << "Add: " << std::endl;
        for (uint i = 0; i < cur_size; ++i)
        {
            for (uint j = 0; j < cur_size; ++j)
            {
                //std::cout << "A: " << A(AposL+i,AposC+j) << "  B: " << B(BposL+i,BposC+j) << std::endl;
                C(i, j) = A(AposL+i, AposC+j) + B(BposL+i, BposC+j);
            }
        }
        return C;
    }

    static SqMatrix subtract(const SqMatrix& A, const SqMatrix& B,
                             uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        SqMatrix C(cur_size);
        for (uint i = 0; i < cur_size; ++i)
        {
            for (uint j = 0; j < cur_size; ++j)
            {
                C(i, j) = A(AposL+i, AposC+j) - B(BposL+i, BposC+j);
            }
        }
        return C;
    }

    static void/*SqMatrix*/ P1 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P1: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //(A11+A22)*(B11+B22)
            /*return*/ SqMatrix::simple_multiply(C, SqMatrix::add(A, A, AposL, AposC, AposL+cur_size, AposC+cur_size, cur_size),
                                                 SqMatrix::add(B, B, BposL, BposC, BposL+cur_size, BposC+cur_size, cur_size),
                                                 0, 0, 0, 0, cur_size);

        else
            /*return*/ SqMatrix::block_multiply(C, SqMatrix::add(A, A, AposL, AposC, AposL+cur_size, AposC+cur_size, cur_size),
                                                SqMatrix::add(B, B, BposL, BposC, BposL+cur_size, BposC+cur_size, cur_size),
                                                0, 0, 0, 0, cur_size);
    }

    static void/*SqMatrix*/ P2 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P2: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //(A21+A22)*B11
            /*return*/ SqMatrix::simple_multiply(C, SqMatrix::add(A, A, AposL+cur_size, AposC, AposL+cur_size, AposC+cur_size, cur_size),
                                                 B,
                                                 0, 0, BposL, BposC, cur_size);
        else
            /*return*/ SqMatrix::block_multiply(C, SqMatrix::add(A, A, AposL+cur_size, AposC, AposL+cur_size, AposC+cur_size, cur_size),
                                                B,
                                                0, 0, BposL, BposC, cur_size);
    }

    static void/*SqMatrix*/ P3 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P3: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //A11*(B12-B22)
            /*return*/ SqMatrix::simple_multiply(C, A,
                                                 SqMatrix::subtract(B, B, BposL, BposC+cur_size, BposL+cur_size, BposC+cur_size, cur_size),
                                                 AposL, AposC, 0, 0, cur_size);
        else
            /*return*/ SqMatrix::block_multiply(C, A,
                                                SqMatrix::subtract(B, B, BposL, BposC+cur_size, BposL+cur_size, BposC+cur_size, cur_size),
                                                AposL, AposC, 0, 0, cur_size);
    }

    static void/*SqMatrix*/ P4 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P4: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //A22*(B21-B11)
            /*return*/ SqMatrix::simple_multiply(C, A,
                                                 SqMatrix::subtract(B, B, BposL+cur_size, BposC, BposL, BposC, cur_size),
                                                 AposL+cur_size, AposC+cur_size, 0, 0, cur_size);
        else
            /*return*/ SqMatrix::block_multiply(C, A,
                                                SqMatrix::subtract(B, B, BposL+cur_size, BposC, BposL, BposC, cur_size),
                                                AposL+cur_size, AposC+cur_size, 0, 0, cur_size);
    }

    static void/*SqMatrix*/ P5 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P5: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //(A11+A12)*B22
            /*return*/ SqMatrix::simple_multiply(C, SqMatrix::add(A, A, AposL, AposC, AposL, AposC+cur_size, cur_size),
                                                 B,
                                                 0, 0, BposL+cur_size, BposC+cur_size, cur_size);
        else
            /*return*/ SqMatrix::block_multiply(C, SqMatrix::add(A, A, AposL, AposC, AposL, AposC+cur_size, cur_size),
                                                B,
                                                0, 0, BposL+cur_size, BposC+cur_size, cur_size);
    }

    static void/*SqMatrix*/ P6 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P6: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //(A21-A11)*(B11+B12)
            /*return*/ SqMatrix::simple_multiply(C, SqMatrix::subtract(A, A, AposL+cur_size, AposC, AposL, AposC, cur_size),
                                                 SqMatrix::add(B, B, BposL, BposC, BposL, BposC+cur_size, cur_size),
                                                 0, 0, 0, 0, cur_size);

        else
            /*return*/ SqMatrix::block_multiply(C, SqMatrix::subtract(A, A, AposL+cur_size, AposC, AposL, AposC, cur_size),
                                                SqMatrix::add(B, B, BposL, BposC, BposL, BposC+cur_size, cur_size),
                                                0, 0, 0, 0, cur_size);
    }

    static void/*SqMatrix*/ P7 (SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //std::cout << "P7: " << cur_size << std::endl;
        if (cur_size == MINIMUM_SIZE)              //(A12-A22)*(B21+B22)
            /*return*/ SqMatrix::simple_multiply(C, SqMatrix::subtract(A, A, AposL, AposC+cur_size, AposL+cur_size, AposC+cur_size, cur_size),
                                                 SqMatrix::add(B, B, BposL+cur_size, BposC, BposL+cur_size, BposC+cur_size, cur_size),
                                                 0, 0, 0, 0, cur_size);
        else
            /*return*/ SqMatrix::block_multiply(C, SqMatrix::subtract(A, A, AposL, AposC+cur_size, AposL+cur_size, AposC+cur_size, cur_size),
                                                SqMatrix::add(B, B, BposL+cur_size, BposC, BposL+cur_size, BposC+cur_size, cur_size),
                                                0, 0, 0, 0, cur_size);
    }

    static void/*SqMatrix*/ block_multiply(SqMatrix& C, const SqMatrix& A, const SqMatrix& B,
                                           uint AposL, uint AposC, uint BposL, uint BposC, uint cur_size)
    {
        //SqMatrix C(B.size);
        //std::cout << "BlocK: " << cur_size << std::endl;
        k++;
        if (cur_size == MINIMUM_SIZE)
            return SqMatrix::simple_multiply(C, A, B, AposL, AposC, BposL, BposC, cur_size);

        uint temp_size = 2;
        while (temp_size < cur_size)
        {
            temp_size *= 2;
        }

        if (temp_size > cur_size) {
            SqMatrix newA(temp_size);
            SqMatrix newB(temp_size);
            SqMatrix newC(temp_size);
            /*SqMatrix pn1(temp_size/2);
            SqMatrix pn2(temp_size/2);
            SqMatrix pn3(temp_size/2);
            SqMatrix pn4(temp_size/2);*/
            for (uint i = 0; i < cur_size; ++i) {
                for (uint j = 0; j < cur_size; ++j) {
                    newA(i, j) = A(i, j);
                    newB(i, j) = B(i, j);
                }
            }
            std::cout << "Inside: " << cur_size << std::endl;
            block_multiply(newC, newA, newB, AposL, AposC, BposL, BposC, temp_size);
            for (uint i = 0; i < cur_size; ++i) {
                for (uint j = 0; j < cur_size; ++j) {
                    C(i, j) = newC(i, j);
                }
            }
            /*P1(pn1, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            P4(pn2, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            P5(pn3, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            P7(pn4, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            SqMatrix C11 = SqMatrix::PNadd(SqMatrix::PNsubtract(SqMatrix::PNadd(pn1, pn2, AposL, AposC, BposL, BposC, temp_size/2),
                                                            pn3, AposL, AposC, BposL, BposC, temp_size/2),
                                         pn4, AposL, AposC, BposL, BposC, temp_size/2);
            pn4 = 0;
            P3(pn4, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            SqMatrix C12 = SqMatrix::PNadd(pn3, pn4, AposL, AposC, BposL, BposC, temp_size/2);
            pn3 = 0;
            P2(pn3, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            SqMatrix C21 = SqMatrix::PNadd(pn2, pn3, AposL, AposC, BposL, BposC, temp_size/2);
            pn2 = 0;
            P6(pn2, newA, newB, AposL, AposC, BposL, BposC, temp_size/2);
            SqMatrix C22 = SqMatrix::PNadd(SqMatrix::PNadd(SqMatrix::PNsubtract(pn1, pn3, AposL, AposC, BposL, BposC, temp_size/2),
                                                       pn2, AposL, AposC, BposL, BposC, temp_size/2),
                                         pn4, AposL, AposC, BposL, BposC, temp_size/2);
            uint point = temp_size/2;
            for (uint i = 0; i < cur_size; ++i)
            {
                for (uint j = 0; j < cur_size; ++j)
                {
                    C(i, j) = C11(i, j);
                    C(point + i, j) = C21(i, j);
                    C(i, point + j) = C12(i, j);
                    C(point + i, point + j) = C22(i, j);
                }
            }*/
            //C.print();
            //return C;
        } else {
            //std::cout << "Outside: " << cur_size << std::endl;
            SqMatrix pn1(temp_size/2);
            SqMatrix pn2(temp_size/2);
            SqMatrix pn3(temp_size/2);
            SqMatrix pn4(temp_size/2);
            P1(pn1, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P1: " << std::endl;
            //pn1.print();
            //std::cout << std::endl;


            P4(pn2, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P4: " << std::endl;
            //pn2.print();
            //std::cout << std::endl;


            P5(pn3, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P5: " << std::endl;
            //pn3.print();
            //std::cout << std::endl;


            P7(pn4, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P7: " << std::endl;
            //pn4.print();
            //std::cout << std::endl;

            SqMatrix C11 = SqMatrix::PNadd(SqMatrix::PNsubtract(SqMatrix::PNadd(pn1, pn2, AposL, AposC, BposL, BposC, temp_size/2),
                                                            pn3, AposL, AposC, BposL, BposC, temp_size/2),
                                         pn4, AposL, AposC, BposL, BposC, temp_size/2);

            pn4 = 0;
            P3(pn4, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P3: " << std::endl;
            //pn4.print();
            //std::cout << std::endl;

            SqMatrix C12 = SqMatrix::PNadd(pn3, pn4, AposL, AposC, BposL, BposC, temp_size/2); // P3 + P5

            pn3 = 0;
            P2(pn3, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P2: " << std::endl;
            //pn3.print();
            //std::cout << std::endl;

            SqMatrix C21 = SqMatrix::PNadd(pn2, pn3, AposL, AposC, BposL, BposC, temp_size/2); // P2 + P4

            pn2 = 0;
            P6(pn2, A, B, AposL, AposC, BposL, BposC, temp_size/2);

            //std::cout << "P6: " << std::endl;
            //pn2.print();
            //std::cout << std::endl;

            SqMatrix C22 = SqMatrix::PNadd(SqMatrix::PNadd(SqMatrix::PNsubtract(pn1, pn3, AposL, AposC, BposL, BposC, temp_size/2),
                                                       pn2, AposL, AposC, BposL, BposC, temp_size/2),
                                         pn4, AposL, AposC, BposL, BposC, temp_size/2); // P1 - P2 + P3 + P6

            //std::cout << "qq" << std::endl;
            uint point = cur_size/2;
            //std::cout << "qq" << std::endl;
            for (uint i = 0; i < cur_size/2; ++i)
            {
                for (uint j = 0; j < cur_size/2; ++j)
                {
                    C(i, j) = C11(i, j);
                    C(point + i, j) = C21(i, j);
                    C(i, point + j) = C12(i, j);
                    C(point + i, point + j) = C22(i, j);
                    //std::cout << "qq: " << point + i << "  " << point + j << "  " << cur_size << std::endl;
                }
            }
            //C.print();
            //return C;
        }
    }




    void print() const
    {
        std::cout.setf((std::ios::fixed));
        for (uint i = 0; i < size; ++i)
        {
            for (uint j = 0; j < size; ++j)
            {
                std::cout << content[i][j] << "\t\t";
            }
            std::cout << std::endl;
        }
    }

    const double& operator()(uint i, uint j) const;

    double& operator()(uint i, uint j);

    SqMatrix& operator=(uint i)
    {
        for (uint j = 0; j < size; ++j)
        {
            for (uint l = 0; l < size; ++l)
            {
                content[j][l] = i;
            }
        }
        return *this;
    }

    /*   friend SqMatrix operator+(const SqMatrix& A, const SqMatrix& B)
    {
        return SqMatrix::add(A, B);
    }
    friend SqMatrix operator-(const SqMatrix& A, const SqMatrix& B)
    {
        return SqMatrix::subtract(A, B);
    }*/

    ~SqMatrix()
    {
        for (uint i = 0; i < size; ++i)
        {
            if (content[i] != nullptr)
            {
                delete[] content[i];
            }
        }
        delete[] content;
        //std::cout << "Cleared" << std::endl;
    }
    uint size;
    double** content;
};

const double & SqMatrix::operator()(uint i, uint j) const {
    {
        return content[i][j];
    }
}

double & SqMatrix::operator()(uint i, uint j) {
    {
        return content[i][j];
    }
}

void Check(const SqMatrix &A, const SqMatrix &B)
{
    bool checker = true;
    double a;
    for (uint i = 0; i < CURRENT_SIZE; ++i)
    {
        for (uint j = 0; j < CURRENT_SIZE; ++j)
        {
            //Accuracy == eight digits after point
            /*a = ( ( (double)( (int)( ( A(i, j) - ( (int) A(i, j) ) ) * 10000) ) ) / 10000 ) + (int)A(i, j);
            b = ( ( (double)( (int)( ( B(i, j) - ( (int) B(i, j) ) ) * 10000) ) ) / 10000 ) + (int)B(i, j);*/
            a = std::abs(A(i, j) - B(i, j));
            //std::cout << "a = " << a << std::endl;
            if (a > 0.00000001)
            {
                checker = false;
            }
        }
    }

    if (checker)
        std::cout << "Equal!!!" << std::endl;
    else
        std::cout << "Not equal" << std::endl;
}

void rec ()
{
    SqMatrix A(4);
    rec();
}

int main() {
    std::ifstream fin("matrixEZ");
    if (!fin)
    {
        std::cout << "Could't read" << std::endl;
        return -1;
    }

    SqMatrix A(CURRENT_SIZE);
    SqMatrix B(CURRENT_SIZE);
    SqMatrix C(CURRENT_SIZE);
    /*for (uint i = 0; i < A.size; ++i)
    {
        for (uint j = 0; j < A.size; ++j)
        {
            A(i, j) = 4.24*i + 5.632*j;
            B(i, j) = 6.72*i + 3.876*j;
        }
    }*/
    A.fill(fin);
    B.fill(fin);
    //A.print();
    std::cout << std::endl;
    //B.print();
    std::cout << std::endl;
    //(SqMatrix::subtract(A, B, 0, 0, 0, 0, 4)).print();
    SqMatrix D(CURRENT_SIZE);
    clock_t t_s = clock();
    SqMatrix::simple_multiply(D, A, B, 0, 0, 0, 0, CURRENT_SIZE);
    t_s = clock() - t_s;
    D.print();
    clock_t t_b = clock();
    SqMatrix::block_multiply(C, A, B, 0, 0, 0, 0, CURRENT_SIZE); // SEG FAULT <<REDO WITHOUT USING NEW CLASS OBJECTS>>
    t_b = clock() - t_b;
    std::cout << std::endl;
    C.print();
    Check(C, D);
    std::cout << "Time taken by simple = " << t_s << "    Time taken by block = " << t_b << std::endl;
    std::cout << "k = " << log(k)/log(7) << std::endl;
    std::cout << "Hello, World!" << std::endl;
    return 0;
}
#endif // INCLUDE_HEADER_HPP_
