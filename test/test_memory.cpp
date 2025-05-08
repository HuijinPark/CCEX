#include "../include/memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// Include your memory management header here
// #include "your_memory_management_header.h"

void testDouble1d() {
    size_t m = 5, m_new = 10;
    double* array = allocDouble1d(m);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        array[i] = i * 1.0;
    }

    // Reallocate and test
    reallocDouble1d(&array, m_new);
    for (size_t i = m-1; i < m_new; i++)
    {
        array[i] = i * 1.0;
    }
    
    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        assert(array[i] == i * 1.0); 
    }

    freeDouble1d(array);
    printf("testDouble1d passed.\n");
}

void testDouble2d() {
    size_t m = 2, n = 3, m_new = 4;
    double** array = allocDouble2d(m, n);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            array[i][j] = i + j * 1.0;
        }
    }

    // Reallocate and test
    reallocDouble2d(&array, m, m_new, n);
    for (size_t i = m-1; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            array[i][j] = i + j * 1.0;
        }
    }

    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            assert(array[i][j] == i + j * 1.0); 
        }
    }

    freeDouble2d(array, m_new);
    printf("testDouble2d passed.\n");
}

void testDouble3d() {
    size_t m = 2, n = 2, l = 2, m_new = 3;
    double*** array = allocDouble3d(m, n, l);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                array[i][j][k] = i + j + k * 1.0;
            }
        }
    }

    // Reallocate and test
    reallocDouble3d(&array, m, m_new, n, l);
    for (size_t i = m-1; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                array[i][j][k] = i + j + k * 1.0;
            }
        }
    }

    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                assert(array[i][j][k] == i + j + k * 1.0); 
            }
        }
    }

    freeDouble3d(array, m_new, n);
    printf("testDouble3d passed.\n");
}

void testDoubleComplex1d(){
    size_t m = 5, m_new = 10;
    DoubleComplex* array = allocDoubleComplex1d(m);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        array[i].real = i * 1.0;
        array[i].imag = i * 1.0;
    }

    // Reallocate and test
    reallocDoubleComplex1d(&array, m_new);
    for (size_t i = m-1; i < m_new; i++)
    {
        array[i].real = i * 1.0;
        array[i].imag = i * 1.0;
    }
    
    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        assert(array[i].real == i * 1.0); 
        assert(array[i].imag == i * 1.0); 
    }

    freeDoubleComplex1d(array);
    printf("testDoubleComplex1d passed.\n");
}

void testDoubleComplex2d(){
    size_t m = 2, n = 3, m_new = 4;
    DoubleComplex** array = allocDoubleComplex2d(m, n);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            array[i][j].real = i + j * 1.0;
            array[i][j].imag = i + j * 1.0;
        }
    }

    // Reallocate and test
    reallocDoubleComplex2d(&array, m, m_new, n);
    for (size_t i = m-1; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            array[i][j].real = i + j * 1.0;
            array[i][j].imag = i + j * 1.0;
        }
    }

    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            assert(array[i][j].real == i + j * 1.0); 
            assert(array[i][j].imag == i + j * 1.0); 
        }
    }

    freeDoubleComplex2d(array, m_new);
    printf("testDoubleComplex2d passed.\n");

}

void testDoubleComplex3d(){
    size_t m = 2, n = 2, l = 2, m_new = 3;
    DoubleComplex*** array = allocDoubleComplex3d(m, n, l);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                array[i][j][k].real = i + j + k * 1.0;
                array[i][j][k].imag = i + j + k * 1.0;
            }
        }
    }

    // Reallocate and test
    reallocDoubleComplex3d(&array, m, m_new, n, l);
    for (size_t i = m-1; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                array[i][j][k].real = i + j + k * 1.0;
                array[i][j][k].imag = i + j + k * 1.0;
            }
        }
    }

    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                assert(array[i][j][k].real == i + j + k * 1.0); 
                assert(array[i][j][k].imag == i + j + k * 1.0); 
            }
        }
    }

    freeDoubleComplex3d(array, m_new, n);
    printf("testDoubleComplex3d passed.\n");

}

void testDoubleTensor1d(){
    size_t m = 5, m_new = 10;
    DoubleTensor* array = allocDoubleTensor1d(m);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        array[i].xx = i * 1.0;
        array[i].xy = i * 1.0;
        array[i].xz = i * 1.0;
        array[i].yx = i * 1.0;
        array[i].yy = i * 1.0;
        array[i].yz = i * 1.0;
        array[i].zx = i * 1.0;
        array[i].zy = i * 1.0;
        array[i].zz = i * 1.0;
    }

    // Reallocate and test
    reallocDoubleTensor1d(&array, m_new);
    for (size_t i = m-1; i < m_new; i++)
    {
        array[i].xx = i * 1.0;
        array[i].xy = i * 1.0;
        array[i].xz = i * 1.0;
        array[i].yx = i * 1.0;
        array[i].yy = i * 1.0;
        array[i].yz = i * 1.0;
        array[i].zx = i * 1.0;
        array[i].zy = i * 1.0;
        array[i].zz = i * 1.0;
    }
    
    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        assert(array[i].xx == i * 1.0); 
        assert(array[i].xy == i * 1.0); 
        assert(array[i].xz == i * 1.0); 
        assert(array[i].yx == i * 1.0); 
        assert(array[i].yy == i * 1.0); 
        assert(array[i].yz == i * 1.0); 
        assert(array[i].zx == i * 1.0); 
        assert(array[i].zy == i * 1.0); 
        assert(array[i].zz == i * 1.0); 
    }

    freeDoubleTensor1d(array);
    printf("testDoubleTensor1d passed.\n");
}

void testDoubleTensor2d(){
    size_t m = 2, n = 3, m_new = 4;
    DoubleTensor** array = allocDoubleTensor2d(m, n);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            array[i][j].xx = i + j * 1.0;
            array[i][j].xy = i + j * 1.0;
            array[i][j].xz = i + j * 1.0;
            array[i][j].yx = i + j * 1.0;
            array[i][j].yy = i + j * 1.0;
            array[i][j].yz = i + j * 1.0;
            array[i][j].zx = i + j * 1.0;
            array[i][j].zy = i + j * 1.0;
            array[i][j].zz = i + j * 1.0;
        }
    }

    // Reallocate and test
    reallocDoubleTensor2d(&array, m, m_new, n);
    for (size_t i = m-1; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            array[i][j].xx = i + j * 1.0;
            array[i][j].xy = i + j * 1.0;
            array[i][j].xz = i + j * 1.0;
            array[i][j].yx = i + j * 1.0;
            array[i][j].yy = i + j * 1.0;
            array[i][j].yz = i + j * 1.0;
            array[i][j].zx = i + j * 1.0;
            array[i][j].zy = i + j * 1.0;
            array[i][j].zz = i + j * 1.0;
        }
    }

    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            assert(array[i][j].xx == i + j * 1.0);
            assert(array[i][j].xy == i + j * 1.0);
            assert(array[i][j].xz == i + j * 1.0);
            assert(array[i][j].yx == i + j * 1.0);
            assert(array[i][j].yy == i + j * 1.0);
            assert(array[i][j].yz == i + j * 1.0);
            assert(array[i][j].zx == i + j * 1.0);
            assert(array[i][j].zy == i + j * 1.0);
            assert(array[i][j].zz == i + j * 1.0);
        }
    }
    freeDoubleTensor2d(array, m_new);
    printf("testDoubleTensor2d passed.\n");
}

void testDoubleTensor3d(){
    size_t m = 2, n = 2, l = 2, m_new = 3;
    DoubleTensor*** array = allocDoubleTensor3d(m, n, l);

    // Fill the array
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                array[i][j][k].xx = i + j + k * 1.0;
                array[i][j][k].xy = i + j + k * 1.0;
                array[i][j][k].xz = i + j + k * 1.0;
                array[i][j][k].yx = i + j + k * 1.0;
                array[i][j][k].yy = i + j + k * 1.0;
                array[i][j][k].yz = i + j + k * 1.0;
                array[i][j][k].zx = i + j + k * 1.0;
                array[i][j][k].zy = i + j + k * 1.0;
                array[i][j][k].zz = i + j + k * 1.0;
            }
        }
    }

    // Reallocate and test
    reallocDoubleTensor3d(&array, m, m_new, n, l);
    for (size_t i = m-1; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                array[i][j][k].xx = i + j + k * 1.0;
                array[i][j][k].xy = i + j + k * 1.0;
                array[i][j][k].xz = i + j + k * 1.0;
                array[i][j][k].yx = i + j + k * 1.0;
                array[i][j][k].yy = i + j + k * 1.0;
                array[i][j][k].yz = i + j + k * 1.0;
                array[i][j][k].zx = i + j + k * 1.0;
                array[i][j][k].zy = i + j + k * 1.0;
                array[i][j][k].zz = i + j + k * 1.0;
            }
        }
    }

    // Check old values
    for (size_t i = 0; i < m_new; i++) {
        for (size_t j = 0; j < n; j++) {
            for (size_t k = 0; k < l; k++) {
                assert(array[i][j][k].xx == i + j + k * 1.0);
                assert(array[i][j][k].xy == i + j + k * 1.0);
                assert(array[i][j][k].xz == i + j + k * 1.0);
                assert(array[i][j][k].yx == i + j + k * 1.0);
                assert(array[i][j][k].yy == i + j + k * 1.0);
                assert(array[i][j][k].yz == i + j + k * 1.0);
                assert(array[i][j][k].zx == i + j + k * 1.0);
                assert(array[i][j][k].zy == i + j + k * 1.0);
                assert(array[i][j][k].zz == i + j + k * 1.0);
            }
        }
    }

    freeDoubleTensor3d(array, m_new, n);
    printf("testDoubleTensor3d passed.\n");
}

int main() {

    printf("\n");
    testDouble1d();
    testDouble2d();
    testDouble3d();
    printf("\n");
    testDoubleComplex1d();
    testDoubleComplex2d();
    testDoubleComplex3d();
    printf("\n");
    testDoubleTensor1d();
    testDoubleTensor2d();
    testDoubleTensor3d();
    printf("\n");
    return 0;
}
